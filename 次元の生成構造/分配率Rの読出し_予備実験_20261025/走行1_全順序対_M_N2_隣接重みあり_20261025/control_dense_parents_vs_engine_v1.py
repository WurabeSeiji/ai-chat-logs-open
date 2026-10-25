#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""親生成器の対照（読み出し＋親生成）: 密行列版 parents_dense_v1 を**現構造（上三角、M=N(N−1)/2）**で走らせ、
旧エンジン（run_n_scaling_lowrank_v1.py）が作った親（対照実験_stage123_忠実コピー_20261025/parents/）と比べる。

測る量（N ごと）
  K の同一性   : 旧エンジンの親 v_旧 を密行列 K で評価した残差（旧エンジンの記録値と同じ桁なら K_of は旧 K と同じ）
  収束         : 密行列版の残差、反復時間
  σ スペクトル : σ_max の差、本数、全体の max|Δσ|（同じ本数のとき）
  同一解か     : 大域位相を揃えた距離 d = min_φ ‖v_密 − e^{iφ} v_旧‖（φ = arg(v_旧† v_密)）、振幅分布（ソート）の差、
                 相対位相の 90° 格子からの最大ずれ（両者）
  種           : ‖g_密 − g_旧‖、旧の種が K の値域に持つ成分 ‖Q Qᵀ g_旧‖、閉塞 |gᵀg|
  Z0           : 大域位相を揃えた距離
出力: results_control_parents/control_dense_parents_vs_engine_v1.json と同名 .md。解釈は含まない。"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parents_dense_v1 as P

BASE = os.path.dirname(os.path.abspath(__file__))
ENGINE_PARENTS = os.path.join(BASE, '..', '対照実験_stage123_忠実コピー_20261025', 'parents')
OUT = os.path.join(BASE, 'results_control_parents')
SAVE = os.path.join(BASE, 'control_parents_triu')
NS = [int(x) for x in sys.argv[1:]] or [3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40]
os.makedirs(OUT, exist_ok=True); os.makedirs(SAVE, exist_ok=True)


def aligned_dist(a, b):
    phi = np.angle(np.vdot(b, a))
    return float(np.linalg.norm(a - np.exp(1j * phi) * b)), float(phi)


def dev90(v):
    rel = np.degrees(np.angle(v / v[0]))
    return float(np.max(np.abs(((rel + 45.0) % 90.0) - 45.0)))


rows = []
path_json = os.path.join(OUT, 'control_dense_parents_vs_engine_v1.json')
if os.path.exists(path_json):
    rows = json.load(open(path_json))['rows']
done = {r['N'] for r in rows}
for N in NS:
    if N in done:
        continue
    eng = np.load(os.path.join(ENGINE_PARENTS, f'parent_static_N{N:05d}_makeparent_20260905.npz'))
    ea, eb = P.edges_triu(N); A = P.adjacency(ea, eb)
    res_eng_under_dense_K = P.eigenmode_residual(eng['v'].astype(np.complex128), A)
    t0 = time.time()
    d = P.build_parent(N, P.edges_triu)
    np.savez_compressed(os.path.join(SAVE, f'parent_dense_triu_N{N:05d}_20261025.npz'), structure='triu', **d)
    v, g, Z0, sig = d['v'], d['g'], d['Z0'], d['sigma']
    ve, ge, Z0e, sige = eng['v'], eng['g'], eng['Z0'], eng['sigma']
    dist_v, phi = aligned_dist(v, ve)
    dist_Z0, _ = aligned_dist(Z0, Z0e)
    _, Q = P.zero_closure_seed_dense(np.angle(ve), A, np.random.default_rng(0))
    row = dict(N=N, M=int(d['M']),
               residual_engine_recorded=float(eng['residual']), residual_engine_under_dense_K=res_eng_under_dense_K,
               residual_dense=float(d['residual']), elapsed_s=float(d['elapsed_s']),
               sigma_max_engine=float(np.max(sige)), sigma_max_dense=float(np.max(sig)),
               planes_engine=int(len(sige)), planes_dense=int(len(sig)),
               max_abs_dsigma=(float(np.max(np.abs(np.sort(sig) - np.sort(sige)))) if len(sig) == len(sige) else None),
               aligned_dist_v=dist_v, align_phase_deg=float(np.degrees(phi)),
               max_abs_damp_sorted=float(np.max(np.abs(np.sort(np.abs(v)) - np.sort(np.abs(ve))))),
               dev90_engine_deg=dev90(ve), dev90_dense_deg=dev90(v),
               seed_dist=float(np.linalg.norm(g - ge)), seed_engine_in_rangeK=float(np.linalg.norm(Q @ (Q.T @ ge))),
               seed_closure_dense=float(abs(g @ g)), seed_closure_engine=float(abs(ge @ ge)),
               aligned_dist_Z0=dist_Z0, rank_K=int(d['rank_K']))
    rows.append(row)
    print(json.dumps({k: (float('%.4g' % v) if isinstance(v, float) else v) for k, v in row.items()}, ensure_ascii=False), flush=True)
    json.dump({'rows': rows}, open(path_json, 'w'), indent=1, ensure_ascii=False)

rows.sort(key=lambda r: r['N'])
md = ['| N | M | 旧残差 | 旧親を密K で評価 | 密残差 | σ_max 旧/密 | 本数 旧/密 | d(v) 位相揃え | Δ振幅分布 | 90°ずれ 旧/密 | ‖Δg‖ | 旧g の range(K) 成分 | d(Z0) | 秒 |',
      '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
for r in rows:
    md.append(f"| {r['N']} | {r['M']} | {r['residual_engine_recorded']:.1e} | {r['residual_engine_under_dense_K']:.1e} | {r['residual_dense']:.1e} | "
              f"{r['sigma_max_engine']:.6f}/{r['sigma_max_dense']:.6f} | {r['planes_engine']}/{r['planes_dense']} | {r['aligned_dist_v']:.1e} | "
              f"{r['max_abs_damp_sorted']:.1e} | {r['dev90_engine_deg']:.1e}/{r['dev90_dense_deg']:.1e} | {r['seed_dist']:.1e} | "
              f"{r['seed_engine_in_rangeK']:.1e} | {r['aligned_dist_Z0']:.1e} | {r['elapsed_s']:.1f} |")
open(os.path.join(OUT, 'control_dense_parents_vs_engine_v1.md'), 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))
print('wrote', path_json)
