#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SENS-2 用終端窓状態の生成（stage1 ランナーの忠実コピー、2026-09-14）。

変更点は3つのみ（物理式・シード生成・停止条件は不変）:
  (1) NS を感度検査代表 N {6,10,16,20,30,40} に限定（パラメータ部分集合）
  (2) 20-step チェックポイントの状態を巡回バッファに保持し、ロック時に
      直近3点（lock-40, lock-20, lock）を保存
  (3) 出力を本フォルダ window_states/ に保存
内蔵対照: 各 (N,seed) の最終状態が stage1 正本 terminal_Z と一致することを検証し
window_crosscheck.json に記録する。
"""
import math
import os
import json

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STATES = ('/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/'
          'マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/電子の反跳実験/'
          'seed除去による準安定相追試/chatgpt追試/干渉保存力学_資格審査とシード無し系列_20260831/'
          'make_parent型初期値_自己無撞着構造_20260907/full_N3_N40_sweep/states')
STAGE1 = os.path.abspath(os.path.join(HERE, '..', '..', 'stage1_ssb_terminal_Z_v1', 'terminal_Z'))
NS = [6, 10, 16, 20, 30, 40]
NSEED = 5
SEED_AMP = 1e-8
RATIO_TOL = 1.005
RATE_TOL = 3e-4
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


def evolve_to_lock_with_window(z0, A, den):
    z = z0.copy(); prev_rel = None; buf = []
    for t in range(1, MAXSTEP + 1):
        z = one_step(z, A, den)
        if t % 20 == 0 or t == MAXSTEP:
            buf.append((t, z.copy()))
            if len(buf) > 3:
                buf.pop(0)
            r = np.abs(z); ratio = r.max() / r.min()
            rel = np.angle(z * np.conj(z[0]))
            if prev_rel is not None:
                rate = np.median(np.abs(((rel - prev_rel + np.pi) % (2 * np.pi) - np.pi))) / 20.0
                if ratio < RATIO_TOL and rate < RATE_TOL:
                    return z, t, buf
            prev_rel = rel
    return z, MAXSTEP, buf


def main():
    os.makedirs(os.path.join(HERE, 'window_states'), exist_ok=True)
    xc = {}
    for N in NS:
        A = adjacency(N); M = N * (N - 1) // 2; den = N
        Z0 = np.asarray(np.load(os.path.join(STATES, f'hm_N{N}_den_{N}_states_500.npz'))['Z'][0], np.complex128)
        ref = np.load(os.path.join(STAGE1, f'N{N}_terminal_Z.npz'))
        Zw, Tw = [], []
        maxdiff = 0.0
        for seed in range(NSEED):
            rng = np.random.default_rng(1000 * N + seed)
            z = Z0 + SEED_AMP * (rng.standard_normal(M) + 1j * rng.standard_normal(M))
            z = z / np.linalg.norm(z) * np.linalg.norm(Z0)
            zf, ls, buf = evolve_to_lock_with_window(z, A, den)
            maxdiff = max(maxdiff, float(np.abs(zf - ref['Z_terminal'][seed]).max()))
            assert ls == int(ref['lock_step'][seed]), (N, seed, ls)
            Zw.append(np.array([b[1] for b in buf]))
            Tw.append(np.array([b[0] for b in buf], np.int64))
        np.savez_compressed(os.path.join(HERE, 'window_states', f'N{N}_window.npz'),
                            Z_window=np.asarray(Zw, np.complex128),
                            steps=np.asarray(Tw, np.int64),
                            N=np.int64(N), M=np.int64(M))
        xc[str(N)] = {'max_abs_diff_final_vs_stage1': maxdiff,
                      'lock_steps_match': True}
        print(f'N={N}: window saved, final-vs-stage1 maxdiff={maxdiff:.2e}', flush=True)
    with open(os.path.join(HERE, 'window_crosscheck.json'), 'w') as f:
        json.dump(xc, f, indent=2)
    print('WINDOW DONE')


if __name__ == '__main__':
    main()
