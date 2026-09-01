#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=40 の初期状態を生成
標準的なランダム複素ベクトルを生成し、正規化する。
N=40 対応の states_treatment.npz を作成"""
import numpy as np
import os

N = 40
M = N * (N - 1) // 2  # エッジ数 = 780

# シード固定：再現性確保
rng = np.random.default_rng(40)

# ランダム複素ベクトル生成
v = rng.standard_normal(M) + 1j * rng.standard_normal(M)

# 正規化
norm = np.linalg.norm(v)
v = v / norm

# 検証
print(f"N={N}, M={M}")
print(f"norm(v) = {np.linalg.norm(v):.10f}")
print(f"v[0:5] = {v[0:5]}")

# states_treatment.npz として保存（1ステップだけのダミー）
output_dir = os.path.dirname(__file__)
os.makedirs(os.path.join(output_dir, 'hm_N40'), exist_ok=True)

# 初期状態のみを保存（L=1）
Z = np.array([v], dtype=np.complex128)
npz_path = os.path.join(output_dir, 'hm_N40', 'states_treatment.npz')
np.savez_compressed(npz_path, Z=Z)

print(f"✓ 保存完了: {npz_path}")
print(f"  Z shape: {Z.shape}")
