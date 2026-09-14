#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第一段: 全N=6..40 SSBバッチの忠実再走行＋終端Z保存（2026-09-14 木原承認）。

正本 run_ssb_batch_N6_N40_20260912.py の忠実コピー。変更は次の3点のみで、物理式・
シード生成(default_rng(1000*N+seed))・停止条件・診断計算は一切変更しない。
  (1) STATES を絶対パスで指定（配置場所が異なるための環境変更）
  (2) 出力先を本フォルダ（per_N/, ssb_batch_summary.csv）にする
  (3) 各Nで終端状態 Z_terminal（NSEED×M complex128）と診断量を terminal_Z/N{N}_terminal_Z.npz に保存する
元ランナー・正本フォルダには一切書き込まない。
"""
import csv
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# 変更(1): 正本と同一の床データを絶対パスで参照（読み取りのみ）
STATES = ('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/'
          'マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/'
          'seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/'
          'make_parent型初期値_自己無撞着構造_20260907/full_N3_N40_sweep/states')
NS = list(range(6, 41))
NSEED = 5
SEED_AMP = 1e-8
RATIO_TOL = 1.005      # 等振幅（max/min）判定
RATE_TOL = 3e-4        # 相対位相変化率（ゲージ不変）rad/step
MAXSTEP = 2600


def adjacency(N):
    ea, eb = np.triu_indices(N, 1); M = len(ea); A = np.zeros((M, M))
    for e in range(M):
        sh = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e]); sh[e] = False; A[e, sh] = 1.0
    return A


def one_step(z, A, den):
    u = np.exp(1j * np.angle(z)); H = A * (np.conj(u)[:, None] * u[None, :]); np.fill_diagonal(H, 0)
    H = 1j * np.imag(H); w, V = np.linalg.eigh(H); ph = np.exp(-1j * (2 * math.pi / den) * w)
    return V @ (ph * (V.conj().T @ z))


def evolve_to_lock(z0, A, den):
    z = z0.copy(); prev_rel = None
    for t in range(1, MAXSTEP + 1):
        z = one_step(z, A, den)
        if t % 20 == 0 or t == MAXSTEP:
            r = np.abs(z); ratio = r.max() / r.min()
            rel = np.angle(z * np.conj(z[0]))
            if prev_rel is not None:
                rate = np.median(np.abs(((rel - prev_rel + np.pi) % (2 * np.pi) - np.pi))) / 20.0
                if ratio < RATIO_TOL and rate < RATE_TOL:
                    return z, t
            prev_rel = rel
    return z, MAXSTEP


def main():
    os.makedirs(os.path.join(HERE, 'per_N'), exist_ok=True)
    os.makedirs(os.path.join(HERE, 'terminal_Z'), exist_ok=True)   # 変更(3)
    summ = []
    csv_path = os.path.join(HERE, 'ssb_batch_summary.csv')
    with open(csv_path, 'w', newline='') as cf:
        w = csv.writer(cf)
        w.writerow(['N', 'M', 'onset_ok', 'all_equal_amp', 'all_lock_2piN', 'vacuum_spread_med_deg',
                    'vacuum_spread_max_deg', 'SSB_confirmed'])
        for N in NS:
            A = adjacency(N); M = N * (N - 1) // 2; den = N; ideal = -2 * math.pi / N
            Z0 = np.asarray(np.load(os.path.join(STATES, f'hm_N{N}_den_{N}_states_500.npz'))['Z'][0], np.complex128)
            p = Z0.real / np.linalg.norm(Z0.real); q = Z0.imag - np.dot(Z0.imag, p) * p; q /= np.linalg.norm(q)
            perp = lambda z: (lambda zp: np.vdot(zp, zp).real / np.vdot(z, z).real)(
                z - p * np.dot(p, z) - q * np.dot(q, z))
            h0 = perp(Z0)   # 床の初期H⊥/H（≈0のはず）
            per = []
            rels = []
            zfs = []        # 変更(3): 終端状態の収集
            for seed in range(NSEED):
                rng = np.random.default_rng(1000 * N + seed)
                z = Z0 + SEED_AMP * (rng.standard_normal(M) + 1j * rng.standard_normal(M))
                z = z / np.linalg.norm(z) * np.linalg.norm(Z0)
                zf, ls = evolve_to_lock(z, A, den)
                r = np.abs(zf); z1 = one_step(zf, A, den); dth = float(np.angle(np.vdot(zf, z1)))
                per.append(dict(seed=seed, lock_step=ls, Hperp=float(perp(zf)), amp_ratio=float(r.max()/r.min()),
                                delta=dth, derr=abs(dth-ideal)))
                rels.append(np.angle(zf * np.conj(zf[0])))
                zfs.append(zf)   # 変更(3)
            # 変更(3): 終端Zと診断量の保存（診断は per と同一値）
            np.savez_compressed(
                os.path.join(HERE, 'terminal_Z', f'N{N}_terminal_Z.npz'),
                Z_terminal=np.asarray(zfs, np.complex128),
                seed_ids=np.arange(NSEED, dtype=np.int64),
                lock_step=np.asarray([pp['lock_step'] for pp in per], np.int64),
                Hperp=np.asarray([pp['Hperp'] for pp in per], np.float64),
                amp_ratio=np.asarray([pp['amp_ratio'] for pp in per], np.float64),
                delta=np.asarray([pp['delta'] for pp in per], np.float64),
                derr=np.asarray([pp['derr'] for pp in per], np.float64),
                N=np.int64(N), denominator=np.int64(den), M=np.int64(M),
                seed_amp=np.float64(SEED_AMP))
            # 個別分析（正本と同一）
            onset_ok = all(pp['Hperp'] > 10 * (h0 + 1e-12) and pp['Hperp'] > 1e-3 for pp in per)  # 床から明確に離脱
            all_eq = all(pp['amp_ratio'] < 1.01 for pp in per)
            all_lock = all(pp['derr'] < 1e-3 for pp in per)
            base = rels[0]; spreads = []
            for s in range(1, NSEED):
                d = np.degrees(np.abs(((rels[s] - base + np.pi) % (2 * np.pi)) - np.pi)); spreads.append(np.median(d))
            vac_med = float(np.median(spreads)); vac_max = float(max(np.degrees(np.abs(((rels[s]-base+np.pi)%(2*np.pi))-np.pi)).max() for s in range(1, NSEED)))
            ssb = onset_ok and all_eq and all_lock and vac_med > 1.0
            # 個別ファイル
            with open(os.path.join(HERE, 'per_N', f'N{N}.txt'), 'w', encoding='utf-8') as pf:
                pf.write(f'SSB検証 N={N} den={N} M={M}  床H⊥/H={h0:.2e}  −2π/N={ideal:+.6f}\n')
                pf.write(f'{"seed":>4}{"lock_step":>10}{"H⊥/H":>9}{"amp比":>9}{"Δ":>11}{"|Δ+2π/N|":>11}\n')
                for pp in per:
                    pf.write(f'{pp["seed"]:>4}{pp["lock_step"]:>10}{pp["Hperp"]:>9.4f}{pp["amp_ratio"]:>9.4f}'
                             f'{pp["delta"]:>11.6f}{pp["derr"]:>11.1e}\n')
                pf.write(f'真空の向き(相対位相配置)のシード間差[deg] 中央={vac_med:.1f} 最大={vac_max:.1f}\n')
                pf.write(f'判定: inflation全シード={onset_ok} 等振幅={all_eq} Δ=−2π/N={all_lock} 真空相異={vac_med>1.0}'
                         f'  → SSB={"確認" if ssb else "要確認"}\n')
            w.writerow([N, M, onset_ok, all_eq, all_lock, f'{vac_med:.1f}', f'{vac_max:.1f}', ssb]); cf.flush()
            summ.append((N, ssb, vac_med))
            print(f'N={N:>3} M={M:>4}: inflation={onset_ok} 等振幅={all_eq} Δ=−2π/N={all_lock} '
                  f'真空差中央={vac_med:.1f}° → SSB={"確認" if ssb else "要確認"}', flush=True)
    n_ok = sum(s for _, s, _ in summ)
    print(f'\n全{len(summ)}本中 SSB確認: {n_ok}/{len(summ)}')


if __name__ == '__main__':
    main()
