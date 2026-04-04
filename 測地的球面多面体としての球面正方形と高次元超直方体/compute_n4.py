#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_n4.py  ─  n=4（5次元埋め込み）計算実験
===================================================

方針：
  ① d_{4,m} 係数を論文の一般項で厳密計算（主要結果）
  ② 収束半径を解析的に導出
  ③ L=1 のみ nquad で厳密値を計算（小領域なので高速）
  ④ L=2 は級数のみ（収束境界であることを注記）
  ⑤ L=3 は係数の交代符号で発散傾向を解析的に説明

本スクリプトの結果は最終的に update_tables.py に統合する。
"""

import math
from fractions import Fraction
from scipy.integrate import nquad

# ══════════════════════════════════════════════════════════
# ① 係数 d_{n,m} の一般項（論文 7.2 節）
# ══════════════════════════════════════════════════════════

def weak_compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for a in range(total + 1):
        for rest in weak_compositions(total - a, parts - 1):
            yield (a,) + rest

def coeff_tuples(dim, n_terms=8):
    """論文の一般項 d_{dim,m} を厳密分数で返す"""
    out = []
    for m in range(n_terms):
        total = Fraction(0, 1)
        for alpha in weak_compositions(m, dim):
            mn = math.factorial(m)
            for a in alpha: mn //= math.factorial(a)
            denom = 1
            for a in alpha: denom *= (2 * a + 1)
            total += Fraction(mn, denom)
        poch = Fraction(1, 1)
        for j in range(m):
            poch *= Fraction(dim + 1, 2) + j
        c = Fraction((-1)**m, 1) * poch * total / Fraction(math.factorial(m) * 4**m, 1)
        out.append((c.numerator, c.denominator))
    return out

def series_val(R, L, coeffs, dim):
    """d_{dim,m} 係数による級数近似値"""
    x2 = L**2 / R**2
    return L**dim * sum(num / den * x2**m for m, (num, den) in enumerate(coeffs))

# ══════════════════════════════════════════════════════════
# ② 収束半径の解析的導出
# ══════════════════════════════════════════════════════════

def convergence_bound(n, R):
    """超直方体の角で二項展開が収束する辺長の上限 L < 2R/√n"""
    return 2 * R / math.sqrt(n)

# ══════════════════════════════════════════════════════════
# メイン
# ══════════════════════════════════════════════════════════

R   = 2.0
DIM = 4
N_TERMS = 8

print("=" * 65)
print(f"n={DIM}（5次元埋め込み: R^{DIM} -> S^{DIM}(R) ⊂ R^{DIM+1}）")
print(f"球半径 R = {R}")
print("=" * 65)

# ─── ① 係数表 ───────────────────────────────────────────
print(f"\n--- 係数 d_{{4,m}}（{N_TERMS}項、厳密分数）---")
coeffs = coeff_tuples(DIM, N_TERMS)
for m, (num, den) in enumerate(coeffs):
    sign = "+" if num >= 0 else "-"
    print(f"  m={m}: {num}/{den}  (= {num/den:+.8f})")

# ─── ② 収束半径 ──────────────────────────────────────────
Lmax = convergence_bound(DIM, R)
print(f"\n--- 収束半径の解析的導出 ---")
print(f"  条件: √n × (L/2R) < 1  ⟺  L < 2R/√n")
print(f"  n=2: L < {convergence_bound(2, R):.4f}  （L=1,2,3 はすべて収束域内）")
print(f"  n=4: L < {Lmax:.4f}  （L=1 のみ収束域内, L=2 は境界, L=3 は域外）")

# ─── ③ L=1: nquad で厳密値と比較 ─────────────────────────
print(f"\n--- L=1: nquad 厳密値 vs d_{{4,m}} 8項級数 ---")
L = 1.0
h = L / 2.0

def integrand(y4, y3, y2, y1, _h=h, _R=R):
    return _R**5 / (_R**2 + y1**2 + y2**2 + y3**2 + y4**2)**2.5

print("  nquad 実行中（L=1, 小領域）...", flush=True)
exact_L1, quad_err_L1 = nquad(
    integrand,
    [[-h, h], [-h, h], [-h, h], [-h, h]],
    opts={'limit': 60, 'epsabs': 1e-8, 'epsrel': 1e-8}
)
series_L1 = series_val(R, L, coeffs, DIM)
flat_L1   = L**DIM
rho_L1    = exact_L1 / flat_L1
ae_L1     = abs(series_L1 - exact_L1)
re_L1     = ae_L1 / exact_L1

print(f"  exact  = {exact_L1:.15f}  （quad_err={quad_err_L1:.3e}）")
print(f"  series = {series_L1:.15f}")
print(f"  abs_err= {ae_L1:.3e},  rel_err= {re_L1:.3e}")
print(f"  L^4    = {flat_L1:.1f},  充填率 ρ_4 = {rho_L1:.8f}")

# ─── ④ L=2: 級数のみ（境界） ────────────────────────────
print(f"\n--- L=2: d_{{4,m}} 8項級数のみ（収束境界, 厳密値省略）---")
L = 2.0
series_L2 = series_val(R, L, coeffs, DIM)
flat_L2   = L**DIM
rho_L2_approx = series_L2 / flat_L2
print(f"  series = {series_L2:.12f}")
print(f"  L^4    = {flat_L2:.1f},  級数による充填率推定 ρ_4 ≈ {rho_L2_approx:.6f}")
print(f"  注記: L=2 は収束境界（√4×1/2=1）であり級数精度が低下する可能性あり")

# ─── ⑤ L=3: 発散の解析的説明 ────────────────────────────
print(f"\n--- L=3: 発散傾向の解析的説明 ---")
L = 3.0
corner_ratio = math.sqrt(DIM) * (L / (2 * R))
print(f"  超直方体の角での展開パラメータ: √n × L/(2R) = {corner_ratio:.4f} > 1")
print(f"  → 二項展開の収束条件を満たさず、8項級数は参考値に過ぎない")
series_L3 = series_val(R, L, coeffs, DIM)
print(f"  series（参考） = {series_L3:.12f}")
print(f"  ※ 負値または不安定な値が出る場合は発散の証拠")

# ─── サマリ ──────────────────────────────────────────────
print()
print("=" * 65)
print("SUMMARY（MD 差し込み用）")
print("=" * 65)
print(f"n=4 (5次元), R={R}")
print(f"係数: {coeffs}")
print(f"L=1: exact={exact_L1:.12f}, series={series_L1:.12f}, rel_err={re_L1:.3e}, rho={rho_L1:.6f}")
print(f"L=2: series_only={series_L2:.12f}, rho_approx={rho_L2_approx:.6f}")
print(f"収束上限: L < {Lmax:.4f} (= 2R/√n)")
