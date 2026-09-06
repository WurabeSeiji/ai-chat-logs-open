#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重心ゼロ系列（v3、仕様書 §2v3〜§2v3d・N=3,4,5 は B案）を本走行フォルダの初期値データへ配置する。

木原指示（2026-09-06）: 「インフレーションが失敗した実験系のフォルダをそのまま使い、
厳密に重心が0の系列のデータを初期値データにコピーする」。

- 状態の正本: ../最も対称性の高い初期値_20260906/plot_v3_centroid_series_v1.py の build_v3(N)
  （図 fig_v3_centroid_zero_series_N3_N40.png と同一の状態・同一の辺順ブロック割当）。
- 書き込み先: parents_symmetric_staged/parent_static_N{N:05d}_makeparent_20260905.npz
  （走行プログラムが読む旧名。v2 親のコピーを上書き。v2 原本は
   ../最も対称性の高い初期値_20260906/parents_symmetric/ に無傷で残る）。
- キーは v2 親と同じ（v, g, Z0, sigma, residual, n, seed, delta, tol, iters, family, theta, rule）。
  sigma/residual は v2 生成プログラムと同一式（LowRankSystem の K・Rayleigh 商）で実測して格納。
- 書き込み前に各 N で4条件（‖Z‖²=1・|Σz|・|Σa²−Σb²|・|Σab| ≤1e-12）をアサート。
- 対応表と SHA256 を parents_v3_centroid_zero_manifest.json に記録。"""
import hashlib
import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
SERIES = os.path.dirname(BASE)
sys.path.insert(0, os.path.join(SERIES, '最も対称性の高い初期値_20260906'))
ROOT6 = os.path.abspath(os.path.join(SERIES, '..', '..', '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT6, '時間軸Q軸とフェルミオンの生成構造', '検証_対照実験', '第5論文原本_自発的分裂予備実験_v1'))
from plot_v3_centroid_series_v1 import build_v3   # 状態の正本
from run_n_scaling_lowrank_v1 import LowRankSystem

DST = os.path.join(BASE, 'parents_symmetric_staged')
manifest = {}
for N in range(3, 41):
    Z, kind = build_v3(N)
    Z = np.asarray(Z, dtype=np.complex128)
    M = N * (N - 1) // 2
    assert Z.size == M
    a, b = Z.real, Z.imag
    assert abs(float(np.vdot(Z, Z).real) - 1) < 1e-12
    assert abs(Z.sum()) < 1e-12 and abs(float(a @ a - b @ b)) < 1e-12 and abs(float(a @ b)) < 1e-12
    sys_lr = LowRankSystem(N); sys_lr.set_theta(np.angle(Z))
    I = np.eye(M)
    K = np.column_stack([sys_lr.kmatvec(I[:, j]) for j in range(M)])
    lam = complex(np.vdot(Z, K @ Z) / np.vdot(Z, Z))
    resid = float(np.linalg.norm(K @ Z - lam * Z) / np.linalg.norm(Z))
    sigma = float(abs(lam))
    name = f'parent_static_N{N:05d}_makeparent_20260905.npz'
    path = os.path.join(DST, name)
    np.savez_compressed(path, v=Z, g=np.zeros(M, dtype=np.complex128), Z0=Z,
                        sigma=np.array([sigma]), residual=np.array(resid), n=np.array(N),
                        seed=np.array(-1), delta=np.array(0.0), tol=np.array(0.0), iters=np.array(0),
                        family=np.zeros(M, dtype=np.int64), theta=np.angle(Z),
                        rule=np.array('v3_centroid_zero'))
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    manifest[name] = {'source': 'plot_v3_centroid_series_v1.build_v3', 'kind': kind,
                      'sigma': sigma, 'residual': resid, 'sha256': h}
    print(f'N={N:2d} {kind:24s} σ={sigma:10.6f} 残差={resid:.1e} SHA={h[:12]}…')
with open(os.path.join(BASE, 'parents_v3_centroid_zero_manifest.json'), 'w') as f:
    json.dump(manifest, f, indent=2, ensure_ascii=False)
print('WROTE 38 parents + manifest')
