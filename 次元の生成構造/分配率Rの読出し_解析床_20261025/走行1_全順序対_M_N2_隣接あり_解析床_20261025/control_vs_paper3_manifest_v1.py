#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照（読み出しのみ）：本フォルダで再生成した上三角の床（floors_triu/manifest_floor_triu.json）を、論文3 の manifest
（original_reference_manifest_90deg_floor_analytic.json、2026-09-09 生成）の検算値と比べる。
比べる量：σ、|Σz|、|Σz²|、相対平衡残差、象限数、振幅種数、反復回数。σ と |Σz| は解析量なので一致が期待される。"""
import json, os, sys
BASE = os.path.dirname(os.path.abspath(__file__))
ours = json.load(open(os.path.join(BASE, 'floors_triu', 'manifest_floor_triu.json')))
ref = json.load(open(os.path.join(BASE, 'original_reference_manifest_90deg_floor_analytic.json')))
rows = ['| N | σ 本/論文3 | Δσ | \\|Σz\\| 本/論文3 | \\|Σz²\\| 本/論文3 | 残差 本/論文3 | 象限 | 振幅種 | 反復 |', '|---|---|---|---|---|---|---|---|---|']
ok = True
for name, c in sorted(ours.items()):
    N = c['N']; r = ref[f'parent_90deg_floor_N{N:05d}_analytic.npz']
    ds = abs(c['sigma'] - r['sigma'])
    same_struct = (c['n_quadrant'] == r['n_quadrant'] and c['n_distinct_amp'] == r['n_distinct_amp'])
    ok &= (ds < 1e-9 and abs(c['centroid_abs'] - r['centroid_abs']) < 1e-9 and c['closure_abs'] < 1e-9 and c['releq_residual'] < 1e-9 and same_struct)
    rows.append(f"| {N} | {c['sigma']:.9f} / {r['sigma']:.9f} | {ds:.1e} | {c['centroid_abs']:.6f} / {r['centroid_abs']:.6f} | {c['closure_abs']:.1e} / {r['closure_abs']:.1e} | "
                f"{c['releq_residual']:.1e} / {r['releq_residual']:.1e} | {c['n_quadrant']}/{r['n_quadrant']} | {c['n_distinct_amp']}/{r['n_distinct_amp']} | {c['iters']}/{r['iters']} |")
out = '\n'.join(rows) + f"\n\nCONTROL vs paper3 manifest: {'PASS' if ok else 'FAIL'}\n"
open(os.path.join(BASE, 'results_control', 'control_vs_paper3_manifest.md'), 'w').write(out); print(out)
sys.exit(0 if ok else 1)
