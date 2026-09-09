#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=7/N=6 【90度解析床】親の正規化コピー作成（2026-09-09、木原指示・実験21を90度解析床で追試）。
make_normalized_parents_N7N6_v1.py の忠実コピーに最小変更のみ: 元データを staged理論床 → 90度解析床（実験24）へ差替。正規化ロジックは無変更。

元データ: 90度理論床_整数K解析解_全N_20260909/parents_90deg_floor_analytic/parent_90deg_floor_N0000{7,6}_analytic.npz（各 H=1・位相0/90/180/270）。
正規化: 密度一致条件 H_i/M_i = H9/M9 = 1/36 を満たすよう Z0 を振幅スケール
  c7 = sqrt((21/36)/H7_orig), c6 = sqrt((15/36)/H6_orig)
で一様スケール（位相・構造は無変更）。合計 H = 21/36 + 15/36 = 1 = H9。
Z0 以外のキーは元 npz からそのままコピーし、scale/H_target/source_sha256 を追記。
出力: parents_mixed_N7N6/parent_norm_N0000{7,6}_HtoM_over36_20260908.npz
"""
import hashlib
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PD = os.path.join(BASE, '..', '90度理論床_整数K解析解_全N_20260909', 'parents_90deg_floor_analytic')
OUT = os.path.join(BASE, 'parents_mixed_N7N6')
os.makedirs(OUT, exist_ok=True)

M9 = 36
for N in (7, 6):
    src_path = os.path.join(PD, f'parent_90deg_floor_N{N:05d}_analytic.npz')
    src_sha = hashlib.sha256(open(src_path, 'rb').read()).hexdigest()
    d = dict(np.load(src_path, allow_pickle=False))
    z0 = np.array(d['Z0'], dtype=np.complex128, copy=True)
    M = z0.size
    H_orig = float(np.sum(np.abs(z0) ** 2))
    H_target = M / M9
    c = np.sqrt(H_target / H_orig)
    d['Z0'] = (c * z0).astype(np.complex128)
    d['scale'] = np.float64(c)
    d['H_target'] = np.float64(H_target)
    d['H_orig'] = np.float64(H_orig)
    d['source_sha256'] = np.str_(src_sha)
    out_path = os.path.join(OUT, f'parent_norm_N{N:05d}_HtoM_over36_90degAnalytic_20260909.npz')
    np.savez_compressed(out_path, **d)
    H_new = float(np.sum(np.abs(d['Z0']) ** 2))
    print(f'N={N}: M={M} H_orig={H_orig:.15f} scale={c:.15f} H_new={H_new:.15f} '
          f'(目標 {H_target:.15f}) -> {os.path.basename(out_path)}', flush=True)
print('ALL DONE', flush=True)
