#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全N≥6：写像の「面の乗り換え」を調べる（読出しのみ・物理無変更）。

木原指示（2026-09-12）: 終点を「円」と読むのはつまらない。元の課題＝この写像が「ある面から別の面へ
乗り換える」の観察。全N≥6でこれを定量化せよ。

測定（各Nはロック到達済みデータを使用: 小N・N23-29=500step, 未飽和N(22,30-40)=2000step, N64=2000step）:
  - H⊥/H_end: 初期の90度床平面（step0のp,q）に直交する成分の割合の終端値。
      初期面からの回転面の傾き＝tilt = arcsin(√(H⊥/H))[deg]。乗り換えの大きさ。
  - transient_s3/s1: 遷移(区間[0:min(500,T)])の状態空間SVD上位第3特異値/第1。
      ≳0.1 なら (s1,s2)+(s3,s4) の2回転面＝rank4＝S³候補（思考実験⑥）。
  - tail_s3/s1: ロック後(末尾200)のSVD第3/第1。≈0 なら単一回転面＝rank2＝S¹。
状態空間埋め込み X(t)=[Re z; Im z]∈R^{2M}。
出力: plane_switch_allN.csv, fig_plane_switch_allN.png
"""
import csv
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907'))
_spec = importlib.util.spec_from_file_location('orig', os.path.join(GEN, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(orig)
PL = os.path.abspath(os.path.join(HERE, '..'))  # 位相ロック全N確認フォルダ


def path(N):
    if N == 64:
        return os.path.join(GEN, 'N64_2000step走行_20260907', 'results', 'hm_N64_den_64_states_2000.npz')
    if N == 22:
        return os.path.join(PL, 'N22_2000step検証走行_20260912', 'results', 'hm_N22_den_22_states_2000.npz')
    if N == 40:
        return os.path.join(PL, 'N40_2000step検証走行_20260912', 'results', 'hm_N40_den_40_states_2000.npz')
    if 30 <= N <= 39:
        return os.path.join(PL, '残り未飽和N_2000step検証走行_20260912', 'results', f'hm_N{N}_den_{N}_states_2000.npz')
    return os.path.join(orig.STATES, f'hm_N{N}_den_{N}_states_500.npz')


def perp_frac(S):
    z0 = S[0]; p = z0.real / np.linalg.norm(z0.real)
    q = z0.imag - np.dot(z0.imag, p) * p; q /= np.linalg.norm(q)
    out = np.empty(len(S))
    for t in range(len(S)):
        z = S[t]; a = np.dot(p, z); b = np.dot(q, z); zp = z - p * a - q * b
        out[t] = np.vdot(zp, zp).real / np.vdot(z, z).real
    return out


def s3_over_s1(X):
    W = X - X.mean(0); s = np.linalg.svd(W, compute_uv=False); return s[2] / s[0]


def main():
    rows = []
    for N in list(range(6, 41)) + [64]:
        p = path(N)
        if not os.path.exists(p):
            continue
        S = np.asarray(np.load(p)['Z'], np.complex128); T = len(S)
        src = 2000 if T > 1000 else 500
        hp = perp_frac(S); hend = float(np.median(hp[-100:]))
        tilt = float(np.degrees(np.arcsin(np.sqrt(min(hend, 1.0)))))
        X = np.concatenate([S.real, S.imag], 1)
        rows.append(dict(N=N, src=src, Hperp_end=hend, tilt_deg=tilt,
                         transient_s3s1=float(s3_over_s1(X[:min(500, T)])),
                         tail_s3s1=float(s3_over_s1(X[-200:]))))
    with open(os.path.join(HERE, 'plane_switch_allN.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"{'N':>3}{'src':>6}{'Hperp_end':>10}{'tilt°':>7}{'transient_s3/s1':>16}{'tail_s3/s1':>11}")
    for r in rows:
        print(f"{r['N']:>3}{r['src']:>6}{r['Hperp_end']:>10.4f}{r['tilt_deg']:>7.1f}"
              f"{r['transient_s3s1']:>16.4f}{r['tail_s3s1']:>11.4f}")
    tilts = [r['tilt_deg'] for r in rows]
    print(f"\ntilt°: 範囲 {min(tilts):.1f}–{max(tilts):.1f}, 中央 {np.median(tilts):.1f}")
    print("普遍則: transient で第2面出現(rank4/S³候補), tail で単一面(rank2/S¹)へ再ロック, 面は初期床から傾く")

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        Ns = [r['N'] for r in rows]
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        ax[0].plot(Ns, [r['tilt_deg'] for r in rows], 'o-')
        ax[0].axhline(np.median(tilts), ls=':', c='gray', label=f'median {np.median(tilts):.1f} deg')
        ax[0].set_xlabel('N'); ax[0].set_ylabel('tilt of locked plane from initial floor [deg]')
        ax[0].set_title('plane switch: tilt of final rotation plane vs N'); ax[0].legend(fontsize=8)
        ax[1].plot(Ns, [r['transient_s3s1'] for r in rows], 's-', label='transient s3/s1 (rank4/S3)')
        ax[1].plot(Ns, [r['tail_s3s1'] for r in rows], 'o-', label='locked-tail s3/s1 (rank2/S1)')
        ax[1].axhline(0, ls=':', c='r'); ax[1].set_xlabel('N'); ax[1].set_ylabel('s3/s1')
        ax[1].set_title('transient two-plane (S3) -> locked single-plane (S1)'); ax[1].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_plane_switch_allN.png'), dpi=110)
        print('wrote fig_plane_switch_allN.png')
    except Exception as e:
        print('figure skipped:', e)


if __name__ == '__main__':
    main()
