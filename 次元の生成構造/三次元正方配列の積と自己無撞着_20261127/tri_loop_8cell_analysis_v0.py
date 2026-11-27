#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
§16 で見つかった「8 マス溶け残り (2 枚の 2x2 複素回転行列)」の解析的解明と厳密検算。

解析定理:
  2 つの 2 要素集合 I = {i0, i1} と J = {j0, j1} (I と J は同じでも異なってもよい) の上で、
  j in J の各スライス M_j = T[I, j, I] が 2x2 複素回転行列
      M_j = [[a_j, b_j], [-b_j, a_j]]  <->  複素数 z_j = a_j + i b_j
  であり、それ以外の全成分が 0 のとき、三角形ループ積
      Tri(T)_{ijk} = sum_{p,q,r} T_{ipq} T_{qjr} T_{rpk}
  は厳密にこの 8 マス部分空間を保ち、
      M_j' = (M_{j0}^2 + M_{j1}^2) M_j   <=>   z_j' = (z_0^2 + z_1^2) z_j
  に従う。

  特に z_0, z_1 が同位相 (arg z_0 = arg z_1 = θ) のとき:
    1. ノルム比 ||Tri(U)|| = 1/2 = 0.5 が厳密に恒等的に保たれる (||U||=1)。
    2. 位相角は 1 ステップごとに 3 倍角写像 θ_{n+1} = 3 θ_n (mod 2π) に従う。
       (リャプノフ指数 = ln 3 ≈ 1.0986 の決定論的カオス回転。ただし θ が有理数倍
        2πk/(3^m - 1) なら厳密に周期 m の U^m = I になる。)
    3. 8 マスの外 (残り 504 マス) への微小摂動に対する安定性を、全 512 方向の
       線形化あるいは代表位置の摂動で検証する。
"""
import math
import os
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table


def tri_prod_np(T):
    return np.einsum("ipq,qjr,rpk->ijk", T, T, T, optimize=True)


def proj_dist(A, B):
    na = np.linalg.norm(A)
    nb = np.linalg.norm(B)
    if na == 0 or nb == 0:
        return 1.0
    ov = abs(np.vdot(A / na, B / nb))
    return float(np.sqrt(max(0.0, 1.0 - ov * ov)))


def extract_z(U, I, J):
    """I = (i0, i1), J = (j0, j1) から 2 枚の 2x2 回転行列の複素数 (z0, z1) と対称誤差を抽出"""
    i0, i1 = I
    j0, j1 = J
    zs = []
    max_asym = 0.0
    for j in (j0, j1):
        m00, m01 = U[i0, j, i0], U[i0, j, i1]
        m10, m11 = U[i1, j, i0], U[i1, j, i1]
        max_asym = max(max_asym, abs(m00 - m11), abs(m01 + m10))
        zs.append(complex(m00.real, m01.real))
    return zs[0], zs[1], max_asym


def run():
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out("8 マス溶け残り (2 枚の 2x2 複素回転行列) の解析定理と厳密検算")
    out("=" * 72)

    # 1. 解析定理 z_j' = (z_0^2 + z_1^2) z_j の直接検算 (N = 8, I={6,7}, J={0,1})
    z0_t, z1_t = 0.3 + 0.4j, -0.2 + 0.5j
    T_chk = np.zeros((8, 8, 8), dtype=np.float64)
    for idx_j, (j, z) in enumerate([(0, z0_t), (1, z1_t)]):
        a, b = z.real, z.imag
        T_chk[6, j, 6] = a
        T_chk[6, j, 7] = b
        T_chk[7, j, 6] = -b
        T_chk[7, j, 7] = a
    P_chk = tri_prod_np(T_chk)
    z0_p, z1_p, asym_p = extract_z(P_chk, (6, 7), (0, 1))
    z0_exp = (z0_t**2 + z1_t**2) * z0_t
    z1_exp = (z0_t**2 + z1_t**2) * z1_t
    err_thm = max(abs(z0_p - z0_exp), abs(z1_p - z1_exp), asym_p)
    # 8 マスの外が厳密に 0 であることも確認
    mask8 = np.zeros((8, 8, 8), dtype=bool)
    for i in (6, 7):
        for j in (0, 1):
            for k in (6, 7):
                mask8[i, j, k] = True
    err_out = float(np.max(np.abs(P_chk[~mask8])))
    assert err_thm < 1e-15 and err_out < 1e-15
    out(f"[検算合格] 解析定理 z_j' = (z_0^2 + z_1^2) z_j: 誤差 = {err_thm:.2e}, 8マス外最大値 = {err_out:.2e}\n")

    # 2. 四元数 (N=4) と八元数 (N=8) の (0,1,0)+ε 溶け残り (step 500〜510) で 3 倍角写像を直接確認
    for name, N, I, J in [("四元数 (H)", 4, (0, 1), (0, 1)), ("八元数 (O)", 8, (6, 7), (0, 1))]:
        out(f"--- {name}: N = {N}, I = {I}, J = {J} の step 500〜508 追跡 ---")
        T0 = np.array(make_mult_table(N), dtype=np.float64)
        T0[0, 1, 0] += 1e-4
        U = T0 / np.linalg.norm(T0)
        for _ in range(500):
            P = tri_prod_np(U)
            U = P / np.linalg.norm(P)

        prev_th = None
        for st in range(500, 509):
            z0, z1, asym = extract_z(U, I, J)
            th0 = math.atan2(z0.imag, z0.real)
            th1 = math.atan2(z1.imag, z1.real)
            # 位相差 (mod pi)
            dth_01 = ((th1 - th0 + math.pi / 2) % math.pi) - math.pi / 2
            P = tri_prod_np(U)
            nrm = float(np.linalg.norm(P))
            if prev_th is not None:
                # 3 * prev_th (mod 2pi) と th0 の差 (符号反転 ±pi を考慮して mod pi)
                pred = (3.0 * prev_th) % (2 * math.pi)
                actual = th0 % (2 * math.pi)
                diff3 = min(abs((actual - pred + math.pi) % (2 * math.pi) - math.pi),
                            abs((actual - pred) % (2 * math.pi) - math.pi))
                out(f"  step {st}: ||P||={nrm:.8f}, θ0={math.degrees(th0):+8.3f}°, θ1-θ0={math.degrees(dth_01):+.2e}°, 3倍角誤差={math.degrees(diff3):.2e}°, 回転行列誤差={asym:.2e}")
            else:
                out(f"  step {st}: ||P||={nrm:.8f}, θ0={math.degrees(th0):+8.3f}°, θ1-θ0={math.degrees(dth_01):+.2e}°, 回転行列誤差={asym:.2e}")
            prev_th = th0
            U = P / nrm

        # 3. 内部摂動 vs 外部摂動のリャプノフ指数・安定性比較 (step 500 の状態 U_500 に対して)
        U500 = T0 / np.linalg.norm(T0)
        for _ in range(500):
            P = tri_prod_np(U500)
            U500 = P / np.linalg.norm(P)

        out(f"  [step 500 溶け残りへの微小摂動 (1e-8) の 1〜10 ステップ後の伸び率]")
        # (a) 8 マス内部: (i0, j0, i1) つまり角度方向への摂動
        V_in = U500.copy()
        V_in[I[0], J[0], I[1]] += 1e-8
        V_in /= np.linalg.norm(V_in)
        U_ref = U500.copy()
        ratios_in = []
        d_prev = proj_dist(U_ref, V_in)
        for k_st in range(1, 11):
            U_ref = tri_prod_np(U_ref); U_ref /= np.linalg.norm(U_ref)
            V_in = tri_prod_np(V_in); V_in /= np.linalg.norm(V_in)
            d_cur = proj_dist(U_ref, V_in)
            ratios_in.append(d_cur / d_prev if d_prev > 1e-16 else 0.0)
            d_prev = d_cur
        out(f"    8マス内部 (角度方向): 1ステップ伸び率 (step 2..8) = {[round(x, 4) for x in ratios_in[1:8]]}  (理論値 = 3.0000)")

        # (b) 8 マス外部の代表位置 (全マスを 3 種類に分けて確認)
        #   - 完全外部 (i,j,k すべて I,J の外; N=8 なら (2,3,4))
        #   - 1 軸だけ共有 (例: N=8 なら (6,3,4), (2,0,4), (2,3,6))
        #   - 2 軸共有 (例: N=8 なら (6,0,4), (6,3,6), (2,0,6))
        if N == 8:
            ext_tests = [
                ("完全外部 (2,3,4)", (2, 3, 4)),
                ("1軸共有 (6,3,4)", (6, 3, 4)),
                ("2軸共有 (6,0,4)", (6, 0, 4)),
                ("2軸共有 (6,3,6)", (6, 3, 6)),
                ("2軸共有 (2,0,6)", (2, 0, 6)),
            ]
        else:
            ext_tests = [
                ("完全外部 (2,3,2)", (2, 3, 2)),
                ("1軸共有 (0,3,2)", (0, 3, 2)),
                ("2軸共有 (0,0,2)", (0, 0, 2)),
                ("2軸共有 (0,3,0)", (0, 3, 0)),
            ]
        for elabel, epos in ext_tests:
            V_ex = U500.copy()
            V_ex[epos] += 1e-6
            V_ex /= np.linalg.norm(V_ex)
            U_ref = U500.copy()
            ds = []
            for k_st in range(1, 11):
                U_ref = tri_prod_np(U_ref); U_ref /= np.linalg.norm(U_ref)
                V_ex = tri_prod_np(V_ex); V_ex /= np.linalg.norm(V_ex)
                ds.append(proj_dist(U_ref, V_ex))
            out(f"    外部摂動 {elabel:18s}: step 1,2,5,10 後距離 = {[f'{ds[m]:.2e}' for m in (0, 1, 4, 9)]}")
        out("")

    os.makedirs("results", exist_ok=True)
    out_path = os.path.join("results", "tri_loop_8cell_analysis_v0.txt")
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"結果を保存しました: {out_path}")


if __name__ == "__main__":
    run()
