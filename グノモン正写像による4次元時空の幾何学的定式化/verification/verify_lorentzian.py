#!/usr/bin/env python3
"""
グノモン正写像による誘導計量の記号的検証（Lorentzian版）

検証項目:
  P-8: 署名 (+,+,+,-) のLorentzian計量で
       G_μν + (3/R²)g_μν = 0 が成立するか

  l² = R² + x² + y² + z² - t²
  g_μν = R²/l² (η_μν - x̃_μ x̃_ν / l²)
  ただし η = diag(1,1,1,-1), x̃_μ = η_μν x^ν

依存: sympy (pip install sympy)
実行: python3 verify_lorentzian.py

著者: 木原 範昭 (Noriaki Kihara)
日付: 2026-04-06
"""

from sympy import *
import time
import sys


def gnomonic_verify_lorentzian(verbose=True):
    """
    4次元 Lorentzian 版グノモン計量の記号検証。
    署名 (+,+,+,-), l² = R² + x² + y² + z² - t²

    Returns:
        dict: 各検証項目の成否
    """
    t0 = time.time()
    n = 4

    if verbose:
        print(f"\n{'='*65}")
        print(f"  P-8: Lorentzian版グノモン計量の記号検証")
        print(f"  署名: (+,+,+,-), l² = R² + x² + y² + z² - t²")
        print(f"{'='*65}")

    # --- 座標と基本量 ---
    x, y, z, t = symbols('x y z t')
    coords = [x, y, z, t]
    R = Symbol('R', positive=True)

    # 署名 η_μν = diag(1, 1, 1, -1)
    eta = [1, 1, 1, -1]

    # l² = R² + x² + y² + z² - t²
    l2 = R**2 + x**2 + y**2 + z**2 - t**2

    # x̃_μ = η_μν x^ν（下付き添字）
    x_low = [eta[i] * coords[i] for i in range(n)]

    # =============================================
    # P-1: Lorentzian 計量テンソル
    # =============================================
    if verbose:
        print(f"\n  ── P-1: Lorentzian計量テンソル ──")

    # g_μν = R²/l² (η_μν - x̃_μ x̃_ν / l²)
    g = zeros(n, n)
    for i in range(n):
        for j in range(n):
            eta_ij = eta[i] if i == j else 0
            g[i, j] = R**2 / l2 * (eta_ij - x_low[i] * x_low[j] / l2)

    # 逆計量: g^μν = l²/R² (η^μν + x^μ x^ν / R²)
    ginv = zeros(n, n)
    for i in range(n):
        for j in range(n):
            eta_ij = eta[i] if i == j else 0
            ginv[i, j] = l2 / R**2 * (eta_ij + coords[i] * coords[j] / R**2)

    # 検証: g · g⁻¹ = I
    product = simplify(g * ginv)
    p1_ok = (product == eye(n))

    if verbose:
        print(f"  g · g⁻¹ = I → {'✅' if p1_ok else '❌'}")
        print(f"  ({time.time()-t0:.1f}s)")

    # =============================================
    # P-2: Christoffel 記号
    # =============================================
    if verbose:
        print(f"\n  ── P-2: Christoffel 記号 ──")

    dg = [[[diff(g[i, j], coords[k])
            for k in range(n)] for j in range(n)] for i in range(n)]

    Gamma = [[[None] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for m in range(n):
            for nu in range(m, n):
                val = sum(
                    ginv[a, s] * (dg[s][m][nu] + dg[s][nu][m] - dg[m][nu][s])
                    for s in range(n)
                )
                result = simplify(val / 2)
                Gamma[a][m][nu] = result
                Gamma[a][nu][m] = result

    if verbose:
        print(f"  計算完了 ({time.time()-t0:.1f}s)")

    # =============================================
    # P-3: Ricci テンソル
    # =============================================
    if verbose:
        print(f"\n  ── P-3: Ricci テンソル ──")

    Ric = zeros(n, n)
    for m in range(n):
        for nu in range(m, n):
            val = 0
            for a in range(n):
                val += diff(Gamma[a][m][nu], coords[a])
                val -= diff(Gamma[a][m][a], coords[nu])
                for b in range(n):
                    val += Gamma[a][a][b] * Gamma[b][m][nu]
                    val -= Gamma[a][nu][b] * Gamma[b][m][a]
            s = simplify(val)
            Ric[m, nu] = s
            Ric[nu, m] = s

    if verbose:
        print(f"  計算完了 ({time.time()-t0:.1f}s)")

    # 検証: R_μν = 3/R² · g_μν
    p3r_ok = True
    for i in range(n):
        for j in range(i, n):
            d = simplify(Ric[i, j] - Rational(3, 1) / R**2 * g[i, j])
            if d != 0:
                if verbose:
                    print(f"    R_{i+1}{j+1} - 3/R²·g ≠ 0")
                p3r_ok = False

    if verbose:
        print(f"  R_μν = 3/R² · g_μν → {'✅' if p3r_ok else '❌'}")

    # Ricci スカラー
    R_scalar = simplify(
        sum(ginv[i, j] * Ric[i, j] for i in range(n) for j in range(n))
    )
    p3s_ok = simplify(R_scalar - Rational(12, 1) / R**2) == 0

    if verbose:
        print(f"  R = {R_scalar}, 期待: 12/R² → {'✅' if p3s_ok else '❌'}")

    # =============================================
    # P-4: Einstein 方程式
    # =============================================
    if verbose:
        print(f"\n  ── P-4: G_μν + (3/R²)g_μν = 0 ──")

    Lambda_val = Rational(3, 1) / R**2
    G_tensor = Ric - Rational(1, 2) * R_scalar * g
    check = G_tensor + Lambda_val * g

    p4_ok = True
    for i in range(n):
        for j in range(i, n):
            v = simplify(check[i, j])
            if v != 0:
                if verbose:
                    print(f"    [G+Λg]_{i+1}{j+1} = {v}")
                p4_ok = False

    if verbose:
        print(f"  G_μν + (3/R²)g_μν = 0 → {'✅' if p4_ok else '❌'}")

    # =============================================
    # 結果
    # =============================================
    elapsed = time.time() - t0

    if verbose:
        print(f"\n{'='*65}")
        print(f"  P-8 最終結果: Lorentzian版 (署名 +,+,+,-)")
        print(f"{'='*65}")
        print(f"  P-1 計量:     {'✅' if p1_ok else '❌'}")
        print(f"  P-3 Ricci:    {'✅' if p3r_ok else '❌'}")
        print(f"  P-3 Scalar:   {'✅' if p3s_ok else '❌'}")
        print(f"  P-4 Einstein: {'✅' if p4_ok else '❌'}")
        status = '✅ ALL PASS' if all([p1_ok, p3r_ok, p3s_ok, p4_ok]) else '❌ FAIL'
        print(f"\n  {status}")
        print(f"  総時間: {elapsed:.1f}秒")

    return {
        'P1_metric': p1_ok,
        'P3_ricci': p3r_ok,
        'P3_scalar': p3s_ok,
        'P4_einstein': p4_ok,
        'time': elapsed,
    }


def main():
    result = gnomonic_verify_lorentzian(verbose=True)
    ok = all([result['P1_metric'], result['P3_ricci'],
              result['P3_scalar'], result['P4_einstein']])
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
