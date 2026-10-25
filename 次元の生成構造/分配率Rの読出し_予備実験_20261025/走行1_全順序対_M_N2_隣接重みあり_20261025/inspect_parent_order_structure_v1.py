#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""走行1 の親（parents_all/、全順序対 M=N²）の値を配列の並び順（e = i·N + j）で N×N に置き、
並び順から隣接の構造が読めるかを測る（読み出しのみ。木原指示 2026-10-25）。

測る量（N ごと）
  振幅  : 対角 (i,i) と非対角 (i,j) の振幅の平均・幅、|V_ij| と |V_ji| の差の最大、行ごとの振幅の幅
  位相  : 相対位相（v/v[0]）の束の数（1e-6° で同一視）と束の角度、転置差 θ_ji − θ_ij の取り得る値の数
  不変性: 転置 V^T、共役転置 conj(V^T)、同時巡回シフト (i,j)→(i+s,j+s) に対する大域位相を揃えた相対距離
  自己相関: C(s) = |Σ_e conj(v_e) v_{e+s mod M}| / Σ|v_e|²。s = 1, N, N+1 と上位ピーク
  読めるか: 全対 (e≠f) の |sin(θ_f−θ_e)| を A_ef=1 と A_ef=0 に分け、分布（平均、|sin|<1e-9 の割合）と、
            |sin| の大小で A を言い当てる AUC（Mann–Whitney）。AUC=0.5 なら値から隣接は読めない
出力: results_parent_order_structure/{parent_order_structure_v1.json, parent_order_structure_v1.md, fig_parent_order_structure_N{N}.png}"""
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parents_dense_v1 as P

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, 'results_parent_order_structure')
NS = [int(x) for x in sys.argv[1:]] or [3, 4, 5, 6, 8, 10, 12, 15, 20, 25]
FIG_NS = {5, 10, 20}


def aligned_rel_dist(X, Y):
    phi = np.angle(np.vdot(Y, X))
    return float(np.linalg.norm(X - np.exp(1j * phi) * Y) / np.linalg.norm(Y))


def n_clusters_deg(angles_deg, tol=1e-6):
    a = np.sort(np.mod(angles_deg, 360.0))
    if a.size == 0:
        return 0, []
    reps = [a[0]]
    for x in a[1:]:
        if x - reps[-1] > tol:
            reps.append(x)
    if len(reps) > 1 and (360.0 - reps[-1] + reps[0]) <= tol:
        reps.pop()
    return len(reps), reps


def auc_mann_whitney(x_pos, x_neg):
    # P(x_pos > x_neg) + 0.5 P(=)
    x = np.concatenate([x_pos, x_neg]); r = np.argsort(np.argsort(x, kind='mergesort'), kind='mergesort').astype(float) + 1.0
    # 同順位の平均ランク
    order = np.argsort(x, kind='mergesort'); xs = x[order]; ranks = np.empty_like(xs); i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        ranks[i:j + 1] = (i + j) / 2.0 + 1.0; i = j + 1
    r = np.empty_like(ranks); r[order] = ranks
    n1, n0 = len(x_pos), len(x_neg)
    return float((r[:n1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


rows = []
for N in NS:
    d = np.load(os.path.join(BASE, 'parents_all', f'parent_dense_all_N{N:05d}_20261025.npz'))
    v = d['v'].astype(np.complex128); ea, eb = d['ea'], d['eb']; M = v.size
    assert np.array_equal(ea, np.repeat(np.arange(N), N)) and np.array_equal(eb, np.tile(np.arange(N), N))
    V = v.reshape(N, N)
    amp = np.abs(V); th = np.degrees(np.angle(v / v[0]))
    diag = np.diag(amp); off = amp[~np.eye(N, dtype=bool)]
    nb, reps = n_clusters_deg(th)
    T_diff = np.degrees(np.angle(V.T / V))
    nT, _ = n_clusters_deg(T_diff[~np.eye(N, dtype=bool)].ravel())
    row_spread = float(np.mean(amp.max(axis=1) - amp.min(axis=1)))
    # 不変性
    d_T = aligned_rel_dist(V.T, V); d_cT = aligned_rel_dist(np.conj(V.T), V)
    d_shift = [aligned_rel_dist(np.roll(np.roll(V, -s, axis=0), -s, axis=1), V) for s in range(1, N)]
    d_shift_amp = [float(np.linalg.norm(np.roll(np.roll(amp, -s, axis=0), -s, axis=1) - amp) / np.linalg.norm(amp)) for s in range(1, N)]
    # 自己相関
    C = np.array([abs(np.vdot(v, np.roll(v, -s))) for s in range(M)]) / float(np.vdot(v, v).real)
    top = np.argsort(-C[1:])[:5] + 1
    # 読めるか
    A = P.adjacency(ea, eb)
    th_r = np.angle(v)
    S = np.abs(np.sin(th_r[None, :] - th_r[:, None]))
    iu = np.triu_indices(M, k=1)
    s_all = S[iu]; a_all = A[iu]
    s_pos, s_neg = s_all[a_all == 1], s_all[a_all == 0]
    auc = auc_mann_whitney(s_pos, s_neg)
    # 面積 |Im(conj v_e v_f)| も同じく
    W = np.abs(np.imag(np.conj(v)[:, None] * v[None, :]))[iu]
    auc_area = auc_mann_whitney(W[a_all == 1], W[a_all == 0])
    rec = dict(N=N, M=M, sigma_planes=int(len(d['sigma'])),
               amp_diag_mean=float(diag.mean()), amp_diag_spread=float(diag.max() - diag.min()),
               amp_off_mean=float(off.mean()), amp_off_spread=float(off.max() - off.min()),
               amp_transpose_asym_max=float(np.max(np.abs(amp - amp.T))), amp_row_spread_mean=row_spread,
               phase_bundles=nb, phase_bundle_angles_deg=[round(float(x), 6) for x in reps][:12],
               transpose_phase_diff_values=nT,
               dist_transpose=d_T, dist_conj_transpose=d_cT,
               dist_cyclic_shift_min=float(min(d_shift)) if d_shift else None, dist_cyclic_shift_max=float(max(d_shift)) if d_shift else None,
               dist_cyclic_shift_amp_min=float(min(d_shift_amp)) if d_shift_amp else None,
               autocorr_C1=float(C[1]), autocorr_CN=float(C[N]), autocorr_CN1=float(C[N + 1]) if N + 1 < M else None,
               autocorr_top=[(int(s), round(float(C[s]), 4)) for s in top],
               n_pairs_adj=int(a_all.sum()), n_pairs_nonadj=int((a_all == 0).sum()),
               sin_mean_adj=float(s_pos.mean()), sin_mean_nonadj=float(s_neg.mean()),
               frac_sin_zero_adj=float(np.mean(s_pos < 1e-9)), frac_sin_zero_nonadj=float(np.mean(s_neg < 1e-9)),
               auc_sin_predicts_adjacency=auc, auc_area_predicts_adjacency=auc_area)
    rows.append(rec)
    print(json.dumps({k: (round(v_, 4) if isinstance(v_, float) else v_) for k, v_ in rec.items()}, ensure_ascii=False), flush=True)
    if N in FIG_NS:
        fig, axs = plt.subplots(1, 4, figsize=(20, 4.6))
        im = axs[0].imshow(amp, cmap='viridis'); axs[0].set_title(f'N={N}: |V_ij| (row i, col j)'); plt.colorbar(im, ax=axs[0], fraction=0.046)
        im = axs[1].imshow(np.mod(th.reshape(N, N), 360.0), cmap='hsv', vmin=0, vmax=360); axs[1].set_title('arg(V_ij / V_00) [deg]'); plt.colorbar(im, ax=axs[1], fraction=0.046)
        axs[2].bar(np.arange(1, M), C[1:], width=1.0); axs[2].set_title('autocorrelation C(s) over array index'); axs[2].set_xlabel('shift s'); axs[2].axvline(N, color='r', lw=0.8, ls='--', label='s=N'); axs[2].axvline(N + 1, color='g', lw=0.8, ls='--', label='s=N+1'); axs[2].legend(fontsize=7)
        axs[3].hist(s_neg, bins=50, range=(0, 1), alpha=0.6, label=f'A=0 ({len(s_neg)})', density=True); axs[3].hist(s_pos, bins=50, range=(0, 1), alpha=0.6, label=f'A=1 ({len(s_pos)})', density=True)
        axs[3].set_title(f'|sin(theta_f - theta_e)| by adjacency, AUC={auc:.3f}'); axs[3].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, f'fig_parent_order_structure_N{N}.png'), dpi=150); plt.close(fig)

json.dump({'rows': rows}, open(os.path.join(OUT, 'parent_order_structure_v1.json'), 'w'), indent=1, ensure_ascii=False)
md = ['| N | 平面数 | 対角振幅 平均(幅) | 非対角振幅 平均(幅) | max\\|\\|V_ij\\|−\\|V_ji\\|\\| | 位相の束 | 転置差の値数 | d(転置) | d(共役転置) | d(巡回) min/max | C(1) / C(N) / C(N+1) | 隣接対 / 非隣接対 | mean\\|sin\\| 隣接/非隣接 | \\|sin\\|<1e-9 の割合 隣接/非隣接 | AUC(sin) | AUC(面積) |',
      '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
for r in rows:
    md.append(f"| {r['N']} | {r['sigma_planes']} | {r['amp_diag_mean']:.4f} ({r['amp_diag_spread']:.1e}) | {r['amp_off_mean']:.4f} ({r['amp_off_spread']:.1e}) | {r['amp_transpose_asym_max']:.1e} | {r['phase_bundles']} | {r['transpose_phase_diff_values']} | {r['dist_transpose']:.2f} | {r['dist_conj_transpose']:.2f} | "
              f"{r['dist_cyclic_shift_min']:.2f}/{r['dist_cyclic_shift_max']:.2f} | {r['autocorr_C1']:.3f} / {r['autocorr_CN']:.3f} / {r['autocorr_CN1'] if r['autocorr_CN1'] is None else round(r['autocorr_CN1'],3)} | {r['n_pairs_adj']} / {r['n_pairs_nonadj']} | "
              f"{r['sin_mean_adj']:.3f}/{r['sin_mean_nonadj']:.3f} | {r['frac_sin_zero_adj']:.3f}/{r['frac_sin_zero_nonadj']:.3f} | {r['auc_sin_predicts_adjacency']:.3f} | {r['auc_area_predicts_adjacency']:.3f} |")
open(os.path.join(OUT, 'parent_order_structure_v1.md'), 'w').write('\n'.join(md) + '\n')
print('\n'.join(md)); print('wrote', OUT)
