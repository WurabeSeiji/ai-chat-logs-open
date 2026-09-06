#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""現行の親（make_parent 静的親 Z0, N=3..40）の位相族と振幅の構造（読み出しのみ）。
木原の問い（2026-09-06）: 極限対称親は「全て等振幅・相対位相差90°」。しかし M=N(N−1)/2 が
奇数のとき 0°族と 90°族の本数を等しくできない。現行の親はこの場合どうしているか。
測定:
 (1) 4回対称秩序変数 |⟨e^{4iu}⟩| (u=arg z mod π) と主軸 φ=¼arg Σe^{4iu}: 位相が2直交軸に乗るか
 (2) 各波を軸A(φ)／軸B(φ+90°)に分類: 本数 n_A,n_B、軸からの最大ずれ(deg)
 (3) 振幅: 全体CV・12桁で異なる|z|の個数・族別の平均振幅 r_A,r_B と族内CV
 (4) 閉塞の分担: 族別パワー P_A=Σ_A|z|², P_B=Σ_B|z|²（軸上なら Σz²=0 ⇔ P_A=P_B）、実 |Σz²|
入力: ../N3_N40_stage123_sweep_20260905/parents/parent_static_N{N:05d}_makeparent_20260905.npz（SHA台帳照合）
出力: parent_phase_families_v1.csv / check_parent_phase_families_v1.json"""
import csv, hashlib, json, os
import numpy as np
BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
led = {}
for line in open(os.path.join(PKG, 'SHA256SUMS.txt')):
    p = line.split()
    if len(p) == 2: led[p[1]] = p[0]
rows, out = [], {}
for N in range(3, 41):
    rel = f'parents/parent_static_N{N:05d}_makeparent_20260905.npz'
    path = os.path.join(PKG, rel)
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    assert led[rel] == h, f'INPUT GATE FAIL {rel}'
    d = np.load(path)
    z = np.asarray(d['Z0'], dtype=np.complex128); M = z.size
    r = np.abs(z); u = np.angle(z) % np.pi
    op4 = np.mean(np.exp(4j * u)); phi = (np.angle(op4) / 4.0) % (np.pi / 2)
    dev = ((u - phi + np.pi / 4) % (np.pi / 2)) - np.pi / 4       # 最近軸(φ or φ+90°)からのずれ
    fam = np.round((u - phi) / (np.pi / 2)).astype(int) % 2         # 0=軸A, 1=軸B
    A, B = fam == 0, fam == 1
    PA, PB = float(np.sum(r[A] ** 2)), float(np.sum(r[B] ** 2))
    rec = dict(N=N, M=M, M_odd=bool(M % 2), order4=float(abs(op4)), axis_deg=float(np.degrees(phi)),
               n_A=int(A.sum()), n_B=int(B.sum()), max_axis_dev_deg=float(np.degrees(np.abs(dev).max())),
               amp_cv_all=float(r.std() / r.mean()), n_distinct_amp_12dig=int(len(set(np.round(r, 12)))),
               r_A_mean=float(r[A].mean()) if A.any() else None, r_B_mean=float(r[B].mean()) if B.any() else None,
               amp_cv_A=float(r[A].std() / r[A].mean()) if A.sum() > 1 else 0.0,
               amp_cv_B=float(r[B].std() / r[B].mean()) if B.sum() > 1 else 0.0,
               P_A=PA, P_B=PB, P_A_over_P_B=(PA / PB if PB > 0 else None),
               abs_sum_z2=float(abs(np.sum(z ** 2))), sum_r2=float(np.sum(r ** 2)),
               seed_delta=float(d['delta']), parent_residual=float(d['residual']))
    out[N] = rec; rows.append(rec)
    print(f"N={N:2d} M={M:3d}{'奇' if M%2 else '偶'} 4回秩序={rec['order4']:.4f} n_A={rec['n_A']:3d} n_B={rec['n_B']:3d} "
          f"軸ずれmax={rec['max_axis_dev_deg']:6.2f}° ampCV={rec['amp_cv_all']:.3f} 異なる|z|={rec['n_distinct_amp_12dig']:3d} "
          f"r_A/r_B={ (rec['r_A_mean']/rec['r_B_mean']) if rec['r_B_mean'] else float('nan'):.4f} P_A/P_B={rec['P_A_over_P_B']:.6f} |Σz²|={rec['abs_sum_z2']:.1e}")
with open(os.path.join(BASE, 'parent_phase_families_v1.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
json.dump(out, open(os.path.join(BASE, 'check_parent_phase_families_v1.json'), 'w'), indent=2)
print('ALL DONE')
