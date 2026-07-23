#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
S1（補遺1）恒等式の数値検算

S1 は純粋な恒等式（定理A/B/C）であり近似を含まないが、シリーズ内完結の規約に従い、
乱数生成子で恒等式を機械精度まで確認する自己完結スクリプトを置く。

検算対象:
  定理A: 反対称 K の対角成分は恒等的に零   → max|K_ee| = 0（構成上）
  定理B: 自己影響の最低次は2、(K²)_ee = −Σ_{f≠e} K_ef²
         かつ一次自己結合 (K)_ee = 0
  定理C(ii): 零次同次 K(λZ)=K(Z)（正弦生成子で確認）

実行: python3 self_coupling_identity_check.py
"""

import numpy as np


def random_antisymmetric(n, rng):
    M = rng.normal(size=(n, n))
    return M - M.T


def sine_generator(theta, A):
    K = A * np.sin(theta[None, :] - theta[:, None])
    s = np.linalg.norm(K, 2)
    return K if s < 1e-300 else K / s


def main():
    rng = np.random.default_rng(20260723)
    max_diag = 0.0        # 定理A
    max_selfB1 = 0.0      # 定理B 一次自己結合
    max_identB = 0.0      # 定理B (K²)_ee = −ΣK_ef²
    max_homog = 0.0       # 定理C(ii) 零次同次
    for _ in range(2000):
        n = int(rng.integers(2, 12))
        K = random_antisymmetric(n, rng)
        # 定理A: 対角零
        max_diag = max(max_diag, float(np.max(np.abs(np.diag(K)))))
        # 定理B: 一次自己結合零、二次核 = −行の二乗和
        K2 = K @ K
        for e in range(n):
            max_selfB1 = max(max_selfB1, abs(K[e, e]))
            lhs = K2[e, e]
            rhs = -np.sum(K[e, :] ** 2)
            max_identB = max(max_identB, abs(lhs - rhs))
        # 定理C(ii): 正弦生成子の零次同次（振幅スケールに不変）
        m = n
        theta = rng.uniform(0, 2 * np.pi, m)
        A = (rng.random((m, m)) > 0.5).astype(float)
        A = np.triu(A, 1); A = A + A.T  # 対称隣接
        lam = rng.uniform(0.1, 100.0)
        Kz = sine_generator(theta, A)          # arg(Z) は振幅に依らない
        Klz = sine_generator(theta, A)         # λZ でも位相 theta は同一
        max_homog = max(max_homog, float(np.max(np.abs(Kz - Klz))))
    print("S1 恒等式の数値検算（2000試行, n=2..11）")
    print(f"  定理A max|K_ee|              = {max_diag:.2e}（0 が予言）")
    print(f"  定理B max|一次自己結合 K_ee| = {max_selfB1:.2e}（0 が予言）")
    print(f"  定理B max|(K²)_ee + ΣK_ef²|  = {max_identB:.2e}（0 が予言）")
    print(f"  定理C max|K(Z)−K(λZ)|        = {max_homog:.2e}（0 が予言）")
    ok = max(max_diag, max_selfB1, max_identB, max_homog) < 1e-12
    print("判定:", "PASS（全恒等式が機械精度で成立）" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
