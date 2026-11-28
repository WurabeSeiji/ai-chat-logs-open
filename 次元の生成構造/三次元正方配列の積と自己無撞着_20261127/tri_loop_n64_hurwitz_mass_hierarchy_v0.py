#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 6 速報（質量階層の代数的起源）用：
  Part 1: Hurwitz 4 段階局在系列 (O_8 [d=8], H_4 [d=4], C_2 [d=2], R_1 [d=1]) の
          同一条件 (eps=0.05, 0.10, τ=30) における
          - 衣/芯エネルギー比 E_cloud / E_core および 芯占有率 E_core / E_tot
          - 空間倍音スペクトル (k=0,1,2,3) と不確定性超過 (H_x + H_k - ln 8)
          - 粒子時計速度 phi_P / phi_0 とステップ距離 d1_P
  Part 2: 各 Hurwitz 局在系列 d=8, 4, 2, 1 の臨界崩壊振幅 ε_c(d) の比較
  Part 3: 直積準同型定理 Tri(A ⊗ B) == Tri(A) ⊗ Tri(B) の機械精度検証と
          K 段入れ子 (O_8^{⊗K}) による指数関数的質量階層 r^K の定量表
乱数なし。
"""
import os
import time
import numpy as np
from bm_product_algebra_tables_v0 import make_mult_table
from tri_loop_large_crystal_v0 import tri_prod_fast, tensor_prod_3d


def step_tri(U):
    P = tri_prod_fast(U)
    nrm = float(np.linalg.norm(P))
    if nrm < 1e-15:
        return np.zeros_like(U), 0.0
    return P / nrm, nrm


def walsh_hadamard_matrix(n=8):
    H = np.array([[1.0]])
    while H.shape[0] < n:
        H = np.block([[H, H], [H, -H]])
    H /= np.sqrt(n)
    orders = np.array([bin(i).count("1") for i in range(n)], dtype=int)
    return H, orders


def shannon_entropy(p):
    p = np.asarray(p, dtype=np.float64)
    p = p / np.sum(p)
    nz = p[p > 1e-15]
    return float(-np.sum(nz * np.log(nz)))


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    H4 = np.array(make_mult_table(4), dtype=np.float64)
    C2 = np.array(make_mult_table(2), dtype=np.float64)
    R1 = np.array(make_mult_table(1), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)
    H8, orders = walsh_hadamard_matrix(8)
    ln8 = float(np.log(8.0))

    log("=" * 80)
    log("第 6 速報：Hurwitz 4 段階局在系列 (O8, H4, C2, R1) 質量階層 & 直積指数則検証")
    log("=" * 80)

    # 第 1 層 H = U^(200)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, _ = step_tri(U)
    H = U.copy()
    H_mat = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    _, _, vt_svd = np.linalg.svd(H_mat, full_matrices=False)
    M8 = vt_svd[0].reshape(8, 8, 8)

    # 基準時計 phi_0 (τ=20..35 平均)
    Uh = H.copy()
    phi0_l = []
    for st in range(1, 36):
        Un, _ = step_tri(Uh)
        if st >= 20:
            c = float(np.clip(np.dot(Uh.ravel(), Un.ravel()), -1.0, 1.0))
            phi0_l.append(float(np.arccos(abs(c))))
        Uh = Un
    phi0 = float(np.mean(phi0_l))
    log(f"\n[Ref] 背景 H 単体 基準時計 phi_0 (τ=20..35) = {phi0:.6f} rad/step")

    def build_hurwitz_probe(sub_table, d):
        O_sub = np.zeros((8, 8, 8), dtype=np.float64)
        O_sub[:d, :d, :d] = sub_table
        P = np.zeros((64, 64, 64), dtype=np.float64)
        for A in range(8):
            for B in range(8):
                for C in range(8):
                    val = O_sub[A, B, C]
                    if val == 0.0:
                        continue
                    for a in [4, 5]:
                        for b in [2, 3]:
                            for c in [4, 5]:
                                P[A * 8 + a, B * 8 + b, C * 8 + c] += val * M8[a, b - 2, c]
        return P / np.linalg.norm(P)

    hurwitz_list = [
        ("O8 (d=8)", O8, 8),
        ("H4 (d=4)", H4, 4),
        ("C2 (d=2)", C2, 2),
        ("R1 (d=1)", R1, 1),
    ]
    probes = {name: (build_hurwitz_probe(tbl, d), d) for name, tbl, d in hurwitz_list}

    # ------------------------------------------------------------------
    # Part 1: Hurwitz 4 段階 (d=8, 4, 2, 1) の衣比・倍音比・不確定性・時計速度
    # ------------------------------------------------------------------
    for eps in [0.05, 0.10]:
        log("\n" + "-" * 80)
        log(f"[Part 1] Hurwitz 4 段階局在系列 (eps = {eps:.2f}, τ = 30 安定状態)")
        log("-" * 80)
        log(
            f"{'Algebra':>9s} | {'nz':>4s} {'E_tot':>9s} {'E_core':>9s} {'E_cloud':>9s} | "
            f"{'cloud/core':>10s} {'core/tot':>8s} | {'harm_E%':>7s} {'k0:k1:k2:k3':>23s} | "
            f"{'Hx+Hk-ln8':>9s} {'phi_P/phi0':>10s}"
        )
        for name, _, d in hurwitz_list:
            P_u, _ = probes[name]
            Ut = (H + eps * P_u)
            Ut /= np.linalg.norm(Ut)
            phi_l = []
            for st in range(1, 31):
                Un, _ = step_tri(Ut)
                if st >= 20:
                    Tp = Ut.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, 2:4, :].ravel()
                    Tc = Un.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, 2:4, :].ravel()
                    np_v, nc_v = float(np.linalg.norm(Tp)), float(np.linalg.norm(Tc))
                    if np_v > 1e-14 and nc_v > 1e-14:
                        c = float(np.clip(np.dot(Tp, Tc) / (np_v * nc_v), -1.0, 1.0))
                        phi_l.append(float(np.arccos(abs(c))))
                Ut = Un
            nz = int(np.sum(np.abs(Ut) > 1e-12))
            T_pr = Ut.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, 2:4, :]
            e_tot = float(np.sum(T_pr ** 2))
            # コア: 外側 (A,B,C) すべて < d の部分代数内部
            e_core = float(np.sum(T_pr[:d, :d, :d, :, :, :] ** 2))
            e_cloud = e_tot - e_core
            ratio_cc = e_cloud / e_core if e_core > 1e-18 else 0.0
            frac_core = e_core / e_tot if e_tot > 1e-18 else 0.0

            # 1軸射影エネルギー E(A) の WHT 倍音 (k=0..3) とエントロピー超過
            p_x = (
                np.sum(T_pr ** 2, axis=(1, 2, 3, 4, 5))
                + np.sum(T_pr ** 2, axis=(0, 2, 3, 4, 5))
                + np.sum(T_pr ** 2, axis=(0, 1, 3, 4, 5))
            ) / (3.0 * e_tot)
            wE = H8 @ p_x
            pE = (wE ** 2) / np.sum(wE ** 2)
            by_k = np.zeros(4)
            for idx in range(8):
                by_k[orders[idx]] += pE[idx]
            harm_pct = float(np.sum(by_k[1:]) * 100.0)

            amp_k = H8 @ np.sqrt(np.maximum(p_x, 0.0))
            exc = shannon_entropy(p_x) + shannon_entropy(amp_k ** 2) - ln8
            r_phi = float(np.mean(phi_l) / phi0)
            k_str = f"{by_k[0]:.3f}:{by_k[1]:.3f}:{by_k[2]:.3f}:{by_k[3]:.3f}"
            log(
                f"{name:>9s} | {nz:4d} {e_tot:9.6f} {e_core:9.6f} {e_cloud:9.6f} | "
                f"{ratio_cc:10.5f} {frac_core:8.5f} | {harm_pct:6.2f}% {k_str:>23s} | "
                f"{exc:+9.6f} {r_phi:10.5f}"
            )

    # ------------------------------------------------------------------
    # Part 2: 各 Hurwitz 局在系列 d=8, 4, 2, 1 の臨界崩壊振幅 ε_c(d) 比較
    # ------------------------------------------------------------------
    log("\n" + "-" * 80)
    log("[Part 2] 各 Hurwitz 局在系列 d=8, 4, 2, 1 の臨界転移・崩壊振幅 ε_c(d) (τ=50)")
    log("-" * 80)
    for name, _, d in hurwitz_list:
        P_u, _ = probes[name]
        # 二分探索的に E_bg(far) が消える（または ||P|| がジャンプする）臨界点を探す
        eps_scan = [0.10, 0.15, 0.20, 0.25, 0.275, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80, 1.00]
        first_jump = None
        for ep in eps_scan:
            Ut = (H + ep * P_u)
            Ut /= np.linalg.norm(Ut)
            nP = 0.0
            for _ in range(50):
                Ut, nP = step_tri(Ut)
            nz = int(np.sum(np.abs(Ut) > 1e-12))
            if nP > 0.06 and first_jump is None:
                first_jump = (ep, nz, nP)
                break
        if first_jump:
            log(f"  {name:>9s}: 臨界ジャンプ検出 at eps <= {first_jump[0]:.3f} -> nz={first_jump[1]:4d}, ||P||={first_jump[2]:.5f}")
        else:
            log(f"  {name:>9s}: eps <= 1.00 まで臨界崩壊なし (安定共存)")

    # ------------------------------------------------------------------
    # Part 3: 直積準同型定理 Tri(A ⊗ B) == Tri(A) ⊗ Tri(B) の機械精度検証
    # ------------------------------------------------------------------
    log("\n" + "-" * 80)
    log("[Part 3] 直積準同型定理 Tri(A ⊗ B) == Tri(A) ⊗ Tri(B) と K 段指数階層 r^K")
    log("-" * 80)
    # 8x8x8 と 8x8x8 で検証 (A = O8 + 0.1*pert, B = M8)
    A_test = O8.copy()
    A_test[0, 0, 0] += 0.13
    A_test /= np.linalg.norm(A_test)
    B_test = M8.copy()
    B_test /= np.linalg.norm(B_test)
    AB = tensor_prod_3d(A_test, B_test)
    lhs = tri_prod_fast(AB)
    rhs = tensor_prod_3d(tri_prod_fast(A_test), tri_prod_fast(B_test))
    homo_err = float(np.linalg.norm(lhs - rhs) / np.linalg.norm(lhs))
    log(f"  ||Tri(A ⊗ B) - Tri(A) ⊗ Tri(B)|| / ||Tri(A ⊗ B)|| = {homo_err:.3e}")

    # 指数関数的階層スケーリング表 (1 層比 r に対する K 段の r^{-K})
    log("  [K 段入れ子 O_8^{⊗K} による指数関数的階層増幅 r^{-K}]")
    for label, r_inv in [
        ("r = 1/3 (SO(3) スピン0/2 質量ギャップ)", 3.0),
        ("r = 1/7 (R1 局在コア占有率 ~1/7)", 7.0),
        ("r = 1/8 (1サイト空間サポート比 1/8)", 8.0),
    ]:
        k_ew = np.log(3.4e5) / np.log(r_inv)
        k_pl = np.log(1.0e16) / np.log(r_inv)
        log(f"    {label:40s}: 10^5.5 (e/top) -> K={k_ew:5.1f}段, 10^16 (EW/Pl) -> K={k_pl:5.1f}段")

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.2f}s.")

    with open("results/tri_loop_n64_hurwitz_mass_hierarchy_v0.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    run()
