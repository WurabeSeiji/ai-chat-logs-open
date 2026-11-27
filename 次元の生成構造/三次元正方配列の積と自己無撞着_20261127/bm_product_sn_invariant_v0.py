#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 2: 完全無名 (S_N 不変) な 3 次元正方配列の BM 積の縮約方程式 (3 変数・5 変数) の
厳密導出・検算、および固定点・線形安定性・準周期軌道の有無・Σz^2=0 の解析。

5 つの同値類 (アソシエーション・スキームの 1 + 3 + 1):
    R0: i = j = k           -> a
    R1: j = k != i          -> b1   (1 番目の軸だけ違う)
    R2: i = k != j          -> b2   (2 番目の軸だけ違う)
    R3: i = j != k          -> b3   (3 番目の軸だけ違う)
    R4: i, j, k すべて異なる -> c

BM 積:
    T'_{ijk} = sum_{t=0}^{N-1} T_{itk} * T_{ijt} * T_{tjk}

乱数なし。まず N = 4, 5, 6 で 5 変数の解析式と実際の N^3 の BM 積が全成分一致することを
厳密整数で検算し、そのうえで固定点・周期解・線形安定性・Σz^2 = 0 を解析する。
"""
import cmath
import math
import os


def bm_5var_analytic(a, b1, b2, b3, c, N):
    """
    S_N 不変な 5 変数 (a, b1, b2, b3, c) に対する BM 積の厳密縮約方程式。
    各成分 (i, j, k) で t を場合分けして和を取る:

    1. (i, i, i):
       - t = i (1個):        T_{iii} T_{iii} T_{iii} = a^3
       - t != i (N-1個):     T_{iti} T_{iit} T_{tii} = b2 * b3 * b1
       => a' = a^3 + (N - 1) * b1 * b2 * b3

    2. (i, j, j) [i != j] (R1 -> b1'):
       - t = i (1個):        T_{iij} T_{iji} T_{ijj} = b3 * b2 * b1
       - t = j (1個):        T_{ijj} T_{ijj} T_{jjj} = b1 * b1 * a = a * b1^2
       - t != i, j (N-2個):  T_{itj} T_{ijt} T_{tjj} = c * c * b1 = b1 * c^2
       => b1' = a * b1^2 + b1 * b2 * b3 + (N - 2) * b1 * c^2

    3. (i, j, i) [i != j] (R2 -> b2'):
       - t = i (1個):        T_{iii} T_{iji} T_{iji} = a * b2^2
       - t = j (1個):        T_{iji} T_{ijj} T_{jji} = b2 * b1 * b3
       - t != i, j (N-2個):  T_{iti} T_{ijt} T_{tji} = b2 * c * c = b2 * c^2
       => b2' = a * b2^2 + b1 * b2 * b3 + (N - 2) * b2 * c^2

    4. (i, i, k) [i != k] (R3 -> b3'):
       - t = i (1個):        T_{iik} T_{iii} T_{iik} = b3 * a * b3 = a * b3^2
       - t = k (1個):        T_{ikk} T_{iik} T_{kik} = b1 * b3 * b2
       - t != i, k (N-2個):  T_{itk} T_{iit} T_{tik} = c * b3 * c = b3 * c^2
       => b3' = a * b3^2 + b1 * b2 * b3 + (N - 2) * b3 * c^2

    5. (i, j, k) [すべて異なる] (R4 -> c'):
       - t = i (1個):        T_{iik} T_{iji} T_{ijk} = b3 * b2 * c
       - t = j (1個):        T_{ijk} T_{ijj} T_{jjk} = c * b1 * b3
       - t = k (1個):        T_{ikk} T_{ijk} T_{kjk} = b1 * c * b2
       - t != i, j, k (N-3個): T_{itk} T_{ijt} T_{tjk} = c^3
       => c' = c * (b1*b2 + b2*b3 + b3*b1 + (N - 3) * c^2)
    """
    a_n = a**3 + (N - 1) * b1 * b2 * b3
    b1_n = a * b1**2 + b1 * b2 * b3 + (N - 2) * b1 * c**2
    b2_n = a * b2**2 + b1 * b2 * b3 + (N - 2) * b2 * c**2
    b3_n = a * b3**2 + b1 * b2 * b3 + (N - 2) * b3 * c**2
    c_n = c * (b1 * b2 + b2 * b3 + b3 * b1 + (N - 3) * c**2)
    return a_n, b1_n, b2_n, b3_n, c_n


def build_sn_tensor(a, b1, b2, b3, c, N):
    T = [[[0] * N for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if i == j == k:
                    T[i][j][k] = a
                elif j == k:
                    T[i][j][k] = b1
                elif i == k:
                    T[i][j][k] = b2
                elif i == j:
                    T[i][j][k] = b3
                else:
                    T[i][j][k] = c
    return T


def bm_prod_full(T):
    N = len(T)
    R = [[[0] * N for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                s = 0
                for t in range(N):
                    s += T[i][t][k] * T[i][j][t] * T[t][j][k]
                R[i][j][k] = s
    return R


def verify_5var_formula():
    """N = 4, 5, 6 で 5 変数の解析式と全成分 BM 積が厳密一致することを検算"""
    test_vals = (3, -5, 7, -11, 13)
    for N in (4, 5, 6):
        T = build_sn_tensor(*test_vals, N)
        R_full = bm_prod_full(T)
        an = bm_5var_analytic(*test_vals, N)
        R_expected = build_sn_tensor(*an, N)
        assert R_full == R_expected, f"N={N} で不一致!"
    return True


def jacobian_3var_ratio(x, y, N):
    """
    3 変数 (b1=b2=b3=b != 0) の比座標 x = a/b, y = c/b における写像:
      f(x, y) = (x^3 + N - 1) / (x + 1 + (N - 2) * y^2)
      g(x, y) = y * (3 + (N - 3) * y^2) / (x + 1 + (N - 2) * y^2)
    の厳密ヤコビアン [[df/dx, df/dy], [dg/dx, dg/dy]] と固有値を返す。
    """
    D = x + 1.0 + (N - 2) * y * y
    num_f = x**3 + (N - 1)
    num_g = y * (3.0 + (N - 3) * y * y)

    df_dx = (3.0 * x * x * D - num_f * 1.0) / (D * D)
    df_dy = (-num_f * 2.0 * (N - 2) * y) / (D * D)
    dg_dx = (-num_g * 1.0) / (D * D)
    dg_dy = ((3.0 + 3.0 * (N - 3) * y * y) * D - num_g * 2.0 * (N - 2) * y) / (D * D)

    tr = df_dx + dg_dy
    det = df_dx * dg_dy - df_dy * dg_dx
    disc = tr * tr - 4.0 * det
    sq = cmath.sqrt(disc)
    l1 = 0.5 * (tr + sq)
    l2 = 0.5 * (tr - sq)
    return (df_dx, df_dy, dg_dx, dg_dy), (l1, l2)


def sum_T2_5var(a, b1, b2, b3, c, N):
    """配列全成分の二乗和 Σ T_{ijk}^2 (複素数の場合は二乗の和、ノルム二乗ではない)"""
    return (
        N * (a**2)
        + N * (N - 1) * (b1**2 + b2**2 + b3**2)
        + N * (N - 1) * (N - 2) * (c**2)
    )


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("Step 2: 完全無名 (S_N 不変) 5 変数・3 変数縮約方程式の厳密解析")
    out("=" * 72)

    assert verify_5var_formula()
    out("[検算合格] 5 変数縮約方程式は N = 4, 5, 6 の全成分 BM 積と厳密一致。\n")

    out("--- 1. 3 変数 (b1 = b2 = b3 = b != 0) のすべての固定点と線形安定性 ---")
    out("比座標 x = a/b, y = c/b での固定点:")
    out("  (A) y = 0 (c = 0): x^3 + N - 1 = x(x + 1) => x^3 - x^2 - x + (N - 1) = 0")
    out("  (B) y = ±1 (c = ±b): x = 1 (a = b), 倍率 λ = N")
    out("  (C) y = ±sqrt(N+1): x = -(N-1), 倍率 λ = N(N-2)")
    out("")

    for N in (4, 8, 24):
        out(f"[N = {N}]")
        # (B) y = 1, x = 1 (y = -1 も同じ固有値)
        _, (l1, l2) = jacobian_3var_ratio(1.0, 1.0, N)
        out(f"  固定点 (B) x=1, y=±1 (a=b=±c): 固有値 = {l1.real:+.4f}{l1.imag:+.4f}i, {l2.real:+.4f}{l2.imag:+.4f}i  (|λ|={abs(l1):.4f}, {abs(l2):.4f})")
        # (C) y = sqrt(N+1), x = -(N-1)
        xC = -(N - 1.0)
        yC = math.sqrt(N + 1.0)
        _, (l1, l2) = jacobian_3var_ratio(xC, yC, N)
        out(f"  固定点 (C) x={xC:.1f}, y=±sqrt({N+1}): 固有値 = {l1.real:+.4f}{l1.imag:+.4f}i, {l2.real:+.4f}{l2.imag:+.4f}i  (|λ|={abs(l1):.4f}, {abs(l2):.4f})")
        # (A) y = 0 の複素根を含む 3 根: x^3 - x^2 - x + (N-1) = 0 (numpy なしでコンパニオンやカルダノで求める)
        # 実根 1 つ + 複素共役 2 つ (N >= 4 のとき判別式 < 0)
        # ニュートン法で実根を求める
        xr = -float(N)**(1/3)
        for _ in range(50):
            fx = xr**3 - xr**2 - xr + (N - 1)
            dfx = 3*xr**2 - 2*xr - 1
            xr -= fx / dfx
        # 残りの 2 次式: x^2 + (xr - 1)x + (xr^2 - xr - 1) = 0
        p = xr - 1.0
        q = xr**2 - xr - 1.0
        disc = p*p - 4.0*q
        xc1 = 0.5 * (-p + cmath.sqrt(disc))
        xc2 = 0.5 * (-p - cmath.sqrt(disc))
        for label, xv in [("実根", complex(xr, 0)), ("複素根+", xc1), ("複素根-", xc2)]:
            # 複素ヤコビアン (y = 0 なので対角化される!)
            D = xv + 1.0
            df_dx = (3*xv*xv*D - (xv**3 + N - 1)) / (D*D)
            dg_dy = 3.0 / D
            out(f"  固定点 (A: c=0, {label}) x={xv.real:+.4f}{xv.imag:+.4f}i: 固有値(x方向)={df_dx.real:+.4f}{df_dx.imag:+.4f}i (|λ|={abs(df_dx):.4f}), (y方向)={dg_dy.real:+.4f}{dg_dy.imag:+.4f}i (|λ|={abs(dg_dy):.4f})")
            s2 = sum_T2_5var(xv, 1, 1, 1, 0, N)
            out(f"    ΣT^2 (b=1) = {s2.real:+.4f}{s2.imag:+.4f}i")

    out("\n--- 2. 3 変数 (複素数) で ΣT_{ijk}^2 = 0 を満たす固定点の有無 ---")
    out("各固定点での ΣT^2 / b^2 = N x^2 + 3 N(N-1) + N(N-1)(N-2) y^2:")
    for N in (4, 8, 24):
        sB = N * 1 + 3*N*(N-1) + N*(N-1)*(N-2)*1
        sC = N * (N-1)**2 + 3*N*(N-1) + N*(N-1)*(N-2)*(N+1)
        out(f"  N={N}: (B) ΣT^2/b^2 = {sB}, (C) ΣT^2/b^2 = {sC} (どちらも正で 0 にならない)")
    out("  → 3 変数実数固定点 (B)(C) はすべて ΣT^2 > 0。複素数固定点 (A: c=0) も ΣT^2 != 0。")

    out("\n--- 3. 5 変数 (b1, b2, b3 を区別) の対称性の破れと新たな解 ---")
    out("5 変数方程式から、b_m' = b_m * (a * b_m + (N-2)*c^2) + b1*b2*b3 (m=1,2,3)。")
    out("注目点:")
    out("  もし b_m のどれか 1 つが 0 (例: b1 = 0, b2 != 0, b3 != 0) なら:")
    out("    b1*b2*b3 = 0 となり、b1' = 0 が厳密に保たれる (不変部分空間)!")
    out("    このとき:")
    out("      a'  = a^3")
    out("      b2' = b2 * (a * b2 + (N - 2) * c^2)")
    out("      b3' = b3 * (a * b3 + (N - 2) * c^2)")
    out("      c'  = c * (b2 * b3 + (N - 3) * c^2)")
    out("  さらに b2 != 0, b3 != 0, c != 0 で固定点 (λ = a^2) を解くと:")
    out("    b2' = a^2 b2 => a b2 + (N-2)c^2 = a^2")
    out("    b3' = a^2 b3 => a b3 + (N-2)c^2 = a^2 => b2 = b3")
    out("    c'  = a^2 c  => b2^2 + (N-3)c^2 = a^2")
    out("  b2 = b3 = b, a = 1 とおくと: b + (N-2)c^2 = 1 かつ b^2 + (N-3)c^2 = 1。")
    out("  ここから c^2 = (1 - b)/(N - 2) を代入すると:")
    out("    (N-2) b^2 - (N-3) b - 1 = 0 => (b - 1)((N-2) b + 1) = 0!")
    out("  → 非自明解: b2 = b3 = -1/(N-2), c^2 = (N-1)/(N-2)^2 (実数解あり!)")
    out("  さらに b1 = b2 = 0, b3 != 0 のとき:")
    out("    a' = a^3, b3' = b3*(a*b3 + (N-2)*c^2), c' = (N-3)*c^3")
    out("    c' = a^2 c => (N-3)*c^2 = a^2 => c = ± a / sqrt(N-3)")
    out("    b3' = a^2 b3 => a*b3 + (N-2)/(N-3)*a^2 = a^2 => b3 = -a / (N-3)")
    out("  そして最も重要な点: a = 0 の場合 (対角ゼロ):")
    out("    a'  = (N-1) b1 b2 b3  (b1=0 なら a'=0 も厳密に保たれる!)")
    out("    b1 = 0 かつ a = 0 の不変部分空間では:")
    out("      b2' = (N-2) b2 c^2")
    out("      b3' = (N-2) b3 c^2")
    out("      c'  = c (b2 b3 + (N-3) c^2)")
    out("    ここで b2/b3 は厳密な保存量! さらに b2*b3 = c^2 なら任意の (b2, b3, c) が固定点 (連続族)!")
    out("    しかも複素数で b2 = r e^{iθ}, b3 = r e^{-iθ} (または b2 = i c, b3 = -i c) にすると:")
    out("    ΣT^2 = N(N-1)(b2^2 + b3^2 + (N-2)c^2) を 0 にできるか確認:")
    out("      b2*b3 = c^2 かつ b2^2 + b3^2 + (N-2)c^2 = 0")
    out("      => (b2/c)^2 + (c/b2)^2 = -(N-2)")
    out("      u = (b2/c)^2 とおくと u^2 + (N-2)u + 1 = 0 => u は負の実数 (N >= 4)!")
    out("      よって b2/c = i * sqrt(|u|), b3/c = -i / sqrt(|u|) は厳密に純虚数解が存在し、")
    out("      ★ BM 積の厳密固定点 かつ ΣT^2 = 0 を同時に満たす! ★")

    # 実際に N = 4, 8 で上の「a=0, b1=0, ΣT^2=0 の厳密固定点」を数値検算
    out("\n  [検算: a = 0, b1 = 0, b2*b3 = c^2, ΣT^2 = 0 の厳密解]")
    for N in (4, 8):
        # u^2 + (N-2)u + 1 = 0 の負実根
        u = 0.5 * (-(N - 2) + math.sqrt((N - 2)**2 - 4))  # u < 0
        # b2 = i * sqrt(-u), b3 = -i / sqrt(-u) = 1 / b2, c = 1
        c_val = 1.0 + 0j
        b2_val = 1j * math.sqrt(-u)
        b3_val = c_val**2 / b2_val
        an = bm_5var_analytic(0j, 0j, b2_val, b3_val, c_val, N)
        lam = N - 2
        err = max(abs(an[m] - lam * v) for m, v in enumerate((0j, 0j, b2_val, b3_val, c_val)))
        s2 = sum_T2_5var(0j, 0j, b2_val, b3_val, c_val, N)
        out(f"    N={N}: b2={b2_val.imag:+.6f}i, b3={b3_val.imag:+.6f}i, c=1")
        out(f"      BM積誤差 max|T' - {lam}T| = {err:.2e}, ΣT^2 = {s2.real:+.2e}{s2.imag:+.2e}i")

    # 最後に、複素数 5 変数で「b2*b3 と c^2 の比」が回る準周期軌道があるかチェック
    out("\n--- 4. a = 0, b1 = 0 の不変部分空間における比 z = (b2*b3)/c^2 の 1 変数力学 ---")
    out("b2' = (N-2) b2 c^2, b3' = (N-2) b3 c^2, c' = c (b2 b3 + (N-3) c^2) より、")
    out("z = (b2*b3)/c^2 とおくと:")
    out("  z' = (b2' b3') / c'^2 = (N-2)^2 * z / (z + N - 3)^2")
    out("固定点 z' = z (z != 0) は (z + N - 3)^2 = (N - 2)^2 => z = 1 または z = -(2N - 5)。")
    for N in (4, 8):
        # dz'/dz at z = 1:
        # z' = (N-2)^2 * z * (z + N - 3)^(-2)
        # dz'/dz = (N-2)^2 * [(z+N-3) - 2z] / (z+N-3)^3 = (N-2)^2 * (N - 3 - z) / (z+N-3)^3
        for z0 in (1.0, -(2*N - 5.0)):
            dz = ((N - 2)**2 * (N - 3 - z0)) / ((z0 + N - 3)**3)
            out(f"  N={N}, z0={z0:+.1f}: dz'/dz = {dz:+.4f} (|dz'/dz| = {abs(dz):.4f})")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "bm_product_sn_invariant_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
