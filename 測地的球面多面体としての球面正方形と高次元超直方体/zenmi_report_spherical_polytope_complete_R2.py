#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
zenmi_report_spherical_polytope_complete_R2.py

再現性情報
==========
本スクリプトは、以下の環境で実行・検証した。

Python
------
3.13.5

Platform
--------
Linux-4.4.0-x86_64-with-glibc2.41

主要ライブラリ
--------------
sympy == 1.14.0
scipy == 1.17.0
numpy == 2.3.5

研究目的
========
本スクリプトは、以下の 2 系統の幾何学的対象を同時に検証する。

1. 2 次元の内在的球面正方形
   - 球面そのものの上で、等辺・等角に定義される測地的正方形
   - 厳密式とテイラー展開係数を計算する

2. 高次元のグノモン型測地的球面超直方体
   - R^n の直方体を S^n 上へ持ち上げた像として定義する
   - 厳密数値積分と、級数展開による係数を計算する

本スクリプトが行うこと
======================
1. 内在的球面正方形のテイラー係数 c_n を (分子, 分母) のタプル列で生成
2. 内在的球面正方形の厳密値と級数近似値を比較
3. 高次元グノモン型測地的球面超直方体の係数 d_{n,m} を (分子, 分母) のタプル列で生成
4. 高次元版の厳密数値積分と級数近似値を比較
5. R = 2, 辺長 L = 1, 2, 3 の比較表を出力

主要テスト条件
==============
[共通]
- 球半径 R = 2

[内在的球面正方形]
- 一辺の測地線長 a = 1, 2, 3
- 展開次数 = 20（x^20 まで）
- a=3 の場合は x=a/R=1.5 が収束半径 π/2 に近く、20 項近似の精度低下を確認する

[高次元グノモン型]
- 各方向の辺長 L = 1, 2, 3
- 半幅 h = L/2
- 次元 n = 2
- 展開項数 = 8

ファイルの役割
==============
【本ファイル】zenmi_report_spherical_polytope_complete_R2.py
    - 本研究の「計算コア」スクリプトである。
    - 数値計算関数（係数生成・厳密値・級数近似・誤差）をすべて定義・実装している。
    - 単体で実行すると、計算結果を標準出力に表示するのみで、
      論文ファイル（MD）・図（PNG）には一切手を加えない。

【論文更新スクリプト】update_tables.py
    - 本ファイルの計算関数を呼び出した上で、以下を一括自動更新する。
        1. 論文 MD の 10.2 節「主要結果の整理」数表を計算値で上書き
        2. 論文 MD の 13.5 節「実行スクリプト全文」コードブロックを
           本ファイルの最新内容で自動同期
        3. 論文 MD に掲載している 3 枚の PNG 図（辺長 1/2/3）を再生成

実行方法
========
【計算結果だけ確認したい場合（論文ファイルを変更しない）】
    python3 zenmi_report_spherical_polytope_complete_R2.py

    標準出力に以下を表示する:
        - 内在的球面正方形の係数タプル列
        - 高次元版 n=2 の係数タプル列
        - R=2、辺長 1/2/3 に対する厳密値 / 級数値 / 誤差の比較表

【論文 MD と図を最新状態に更新する場合】
    python3 update_tables.py

    これを実行すると次の変更がすべて自動で行われる:
        ・zenmi_report_spherical_polytope_reproducible_full_R2.md
            → 数表（10.2.1 / 10.2.2 節）を計算値で上書き
            → スクリプト全文（13.5 節）を本ファイルの現在の内容で同期
        ・three_surfaces_new_R2{,_side2,_side3}.png の再生成

    ⚠ 注意: update_tables.py は MD ファイルを直接書き換える。
            論文本体を変更したくない場合は、上の
            python3 zenmi_report_spherical_polytope_complete_R2.py
            を実行すること。

理論的背景
==========
[A] 内在的球面正方形
    面積公式:
        A/R^2 = 8*atan(1/sqrt(cos(x))) - 2*pi,   x = a/R

    テイラー展開:
        A/R^2 = sum_{n>=1} c_n x^(2n)

    係数 c_n は、tan(x/2) と sec(x)^(1/2) の級数係数の畳み込みとして得られる。
    実装では SymPy を使って厳密分数として自動生成する。

[B] 高次元グノモン型測地的球面超直方体（2 次元面積版）
    面積公式:
        V_2(h,h;R)
          = ∫∫ R^3/(R^2 + y_1^2 + y_2^2)^(3/2) dy_1 dy_2,  y_i in [-h,h]

    等辺長 L の場合:
        V_2 = L^2 * sum_{m>=0} d_{2,m} * (L^2/R^2)^m

    係数 d_{2,m} は、
    - 二項展開
    - 多項展開
    - 直方体上のモーメント積分
    を組み合わせて厳密分数として計算する。
"""

import math
from math import atan, sqrt, cos, pi
from fractions import Fraction
import sympy as sp
from scipy.integrate import nquad


def intrinsic_coeff_tuples(n_terms=10):
    x = sp.symbols("x")
    f = 8 * sp.atan(1 / sp.sqrt(sp.cos(x))) - 2 * sp.pi
    ser = sp.series(f, x, 0, 2 * n_terms + 1).removeO()
    out = []
    for n in range(1, n_terms + 1):
        c = sp.simplify(ser.coeff(x, 2 * n))
        num, den = sp.fraction(c)
        out.append((int(num), int(den)))
    return out


def intrinsic_area_exact(a: float, R: float) -> float:
    x = a / R
    return (R ** 2) * (8.0 * atan(1.0 / sqrt(cos(x))) - 2.0 * pi)


def intrinsic_area_series(a: float, R: float, coeffs) -> float:
    x = a / R
    s = 0.0
    for n, (num, den) in enumerate(coeffs, start=1):
        s += (num / den) * (x ** (2 * n))
    return (R ** 2) * s


def weak_compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for a in range(total + 1):
        for rest in weak_compositions(total - a, parts - 1):
            yield (a,) + rest


def projected_equal_edge_coeff_tuples(dim: int, n_terms=8):
    out = []
    for m in range(n_terms):
        total = Fraction(0, 1)
        for alpha in weak_compositions(m, dim):
            multinomial = math.factorial(m)
            for a in alpha:
                multinomial //= math.factorial(a)
            denom = 1
            for a in alpha:
                denom *= (2 * a + 1)
            total += Fraction(multinomial, denom)

        poch = Fraction(1, 1)
        for j in range(m):
            poch *= Fraction(dim + 1, 2) + j

        coeff = Fraction((-1) ** m, 1) * poch * total / Fraction(math.factorial(m) * (4 ** m), 1)
        out.append((coeff.numerator, coeff.denominator))
    return out


def projected_area_exact_equal(R: float, L: float):
    h = L / 2.0

    def integrand(y1, y2):
        return R ** 3 / (R ** 2 + y1 * y1 + y2 * y2) ** (3 / 2)

    val, err = nquad(integrand, [[-h, h], [-h, h]])
    return val, err


def projected_area_series_equal(R: float, L: float, coeffs):
    x2 = (L * L) / (R * R)
    s = 0.0
    for m, (num, den) in enumerate(coeffs):
        s += (num / den) * (x2 ** m)
    return (L ** 2) * s


def rel_err(approx: float, exact: float) -> float:
    if exact == 0.0:
        return abs(approx - exact)
    return abs((approx - exact) / exact)


def main():
    R = 2.0
    lengths = [1.0, 2.0, 3.0]

    coeffs_intr = intrinsic_coeff_tuples(10)
    coeffs_proj = projected_equal_edge_coeff_tuples(2, 8)

    print("=== Coefficient tuples: intrinsic spherical square ===")
    print(coeffs_intr)
    print()

    print("=== Coefficient tuples: projected geodesic spherical orthotope (n=2) ===")
    print(coeffs_proj)
    print()

    print("=== TEST TABLE: R = 2 ===")
    print()

    print("[Intrinsic spherical square]")
    for a in lengths:
        exact = intrinsic_area_exact(a, R)
        approx = intrinsic_area_series(a, R, coeffs_intr)
        print(f"a={a:.1f} | exact={exact:.15f} | series20={approx:.15f} | abs_err={abs(approx-exact):.3e} | rel_err={rel_err(approx, exact):.3e}")
    print()

    print("[Projected gnomonic image, n=2]")
    for L in lengths:
        exact, quad_err = projected_area_exact_equal(R, L)
        approx = projected_area_series_equal(R, L, coeffs_proj)
        print(f"L={L:.1f} | exact={exact:.15f} | series8={approx:.15f} | abs_err={abs(approx-exact):.3e} | rel_err={rel_err(approx, exact):.3e} | quad_err={quad_err:.3e}")
    print()

    print("[Reference flat Euclidean square areas]")
    for L in lengths:
        print(f"L={L:.1f} | flat_area={L**2:.15f}")


if __name__ == "__main__":
    main()
