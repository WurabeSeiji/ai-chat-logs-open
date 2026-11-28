#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 5 速報（本丸：重力曲率・局所時間遅延・信号伝播上限・ローレンツ因子）数値検証スクリプト
  Part 1: 局在質量 P4 (1サイト) / P5 (4サイト) の振幅 ε スイープによる
          - 局所固有時間進行率 (背景サイト A=0 vs A=1..7, および粒子セクター位相速度)
          - 空間エネルギー密度・背景計量変形 Δg_00(A) とドレスド雲の成長
  Part 2: スピン 2 (グラビトン候補: 内側 σ_x および外側 SO(3) 5重項 g_11-g_22) の
          背景時計・空間異方性への結合、および局在粒子 P4 との相互作用（吸収・変調）
  Part 3: 信号伝播の上限速度 c（ハミング距離 d_H=0,1,2,3 への因果伝播）と
          空間倍音エネルギー比 (運動量/空間変位) vs 内部時計速度 (ローレンツ因子検証)
  Part 4: 質量直上 (A=0) vs 遠方 (A=1..7) における三角形視野角 ρ = L3/L1（離散曲率）の歪み比較
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


def extract_sector_phase(T6, A_idx, b_slice):
    """
    指定外側サイト A_idx (i=A, j=0..7, k=A の対角和または特定成分) の
    内側 2x2 ブロック (a,c in {4,5}) の複素位相 θ = atan2(c_J, c_0) を返す。
    """
    # A=A_idx, B=0, C=A_idx かつ b=b_slice[0]
    B = T6[A_idx, 0, A_idx, 4:6, b_slice[0], 4:6]
    c0 = 0.5 * (B[0, 0] + B[1, 1])
    cJ = 0.5 * (B[1, 0] - B[0, 1])
    return float(np.arctan2(cJ, c0)), float(np.hypot(c0, cJ))


def unwrap_step_dtheta(th_prev, th_curr):
    """位相角のステップ差分 (-pi, pi]"""
    d = (th_curr - th_prev + np.pi) % (2 * np.pi) - np.pi
    return float(d)


def site_d1_and_energy(U_prev, U_curr, b_slice):
    """
    指定内側 b_slice において、外側サイト A=0..7 ごとの
    (1) エネルギー E(A)
    (2) 正規化 1 ステップ変化距離 d1(A) = ||v_curr - v_prev|| / (0.5*(||v_curr||+||v_prev||))
        および射影内積 cos(phi_A) = <v_prev, v_curr> / (||v_prev|| ||v_curr||) -> 固有回転角 phi_A = arccos(|cos|)
    を返す。
    """
    T_p = U_prev.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, b_slice, :]
    T_c = U_curr.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, b_slice, :]
    E = np.zeros(8)
    d1 = np.zeros(8)
    phi = np.zeros(8)
    for A in range(8):
        vp = T_p[A].ravel()
        vc = T_c[A].ravel()
        np_v = float(np.linalg.norm(vp))
        nc_v = float(np.linalg.norm(vc))
        E[A] = nc_v ** 2
        if np_v > 1e-14 and nc_v > 1e-14:
            # 符号反転は射影同一視せず生距離と角速度両方見る
            d1[A] = float(np.linalg.norm(vc / nc_v - vp / np_v))
            cval = float(np.clip(np.dot(vp, vc) / (np_v * nc_v), -1.0, 1.0))
            phi[A] = float(np.arccos(abs(cval)))
    return E, d1, phi


def measure_view_angle_stats(U, i_idx, j_idx, k_idx, thresh=1e-12):
    """指定成分 (i,j,k) に寄与する全三角形ループ (p,q,r) の視野角 ρ = L3/L1 統計を計測"""
    # Tri(U)_{ijk} = sum_{p,q,r} U_{ipq} U_{qjr} U_{rpk}
    a = np.abs(U[i_idx, :, :])   # (p, q)
    b = np.abs(U[:, j_idx, :])   # (q, r)
    c = np.abs(U[:, :, k_idx])   # (r, p)
    # 寄与の大きい (p,q,r) だけ高速抽出
    p_nz = np.where(np.max(a, axis=1) > thresh)[0]
    q_nz = np.where(np.max(a, axis=0) > thresh)[0]
    r_nz = np.where(np.max(b, axis=1) > thresh)[0]
    rhos = []
    weights = []
    for p in p_nz:
        for q in q_nz:
            apq = a[p, q]
            if apq <= thresh:
                continue
            for r in r_nz:
                bqr = b[q, r]
                crp = c[r, p]
                if bqr <= thresh or crp <= thresh:
                    continue
                prod = apq * bqr * crp
                if prod <= thresh ** 3:
                    continue
                L1 = max(apq, bqr, crp)
                L3 = min(apq, bqr, crp)
                rhos.append(L3 / L1)
                weights.append(prod)
    if not rhos:
        return 0, 0.0, 0.0, 0.0
    rhos = np.array(rhos)
    weights = np.array(weights)
    w_mean = float(np.sum(rhos * weights) / np.sum(weights))
    return len(rhos), float(np.mean(rhos)), w_mean, float(np.min(rhos))


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    H4 = np.array(make_mult_table(4), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)
    H8, orders = walsh_hadamard_matrix(8)

    log("=" * 76)
    log("第 5 速報（本丸）：重力曲率・局所時間遅延・因果律・ローレンツ因子 数値検証")
    log("=" * 76)

    # 第 1 層 H = U^(200)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, normP_H = step_tri(U)
    H = U.copy()
    H_mat = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    _, _, vt_svd = np.linalg.svd(H_mat, full_matrices=False)
    M8 = vt_svd[0].reshape(8, 8, 8)

    # 背景 H 単体の基準時計速度を計測 (τ=20..40 平均)
    Uh = H.copy()
    bg_phi_ref = []
    bg_d1_ref = []
    for st in range(1, 41):
        Uh_next, _ = step_tri(Uh)
        if st >= 20:
            _, d1_h, phi_h = site_d1_and_energy(Uh, Uh_next, slice(0, 2))
            bg_phi_ref.append(float(np.mean(phi_h)))
            bg_d1_ref.append(float(np.mean(d1_h)))
        Uh = Uh_next
    phi0_ref = float(np.mean(bg_phi_ref))
    d1_0_ref = float(np.mean(bg_d1_ref))
    log(f"\n[Ref] 背景 H 単体の基準時計 (τ=20..40 平均): phi_0 = {phi0_ref:.6f} rad/step, d1_0 = {d1_0_ref:.6f}")

    def build_b23_probe(O_part, inner_mode="M8"):
        P = np.zeros((64, 64, 64), dtype=np.float64)
        for A in range(8):
            for B in range(8):
                for C in range(8):
                    val = O_part[A, B, C]
                    if val == 0.0:
                        continue
                    for a in [4, 5]:
                        for b in [2, 3]:
                            for c in [4, 5]:
                                if inner_mode == "M8":
                                    mval = M8[a, b - 2, c]
                                elif inner_mode == "sigma_x":
                                    mval = 1.0 if (a != c) else 0.0
                                elif inner_mode == "sigma_z":
                                    mval = 1.0 if (a == 4 and c == 4) else (-1.0 if (a == 5 and c == 5) else 0.0)
                                P[A * 8 + a, B * 8 + b, C * 8 + c] += val * mval
        nrm = np.linalg.norm(P)
        return P / nrm if nrm > 0 else P

    O_d0 = np.zeros((8, 8, 8))
    O_d0[0, 0, 0] = 1.0
    P4_unit = build_b23_probe(O_d0, "M8")

    O_H4 = np.zeros((8, 8, 8))
    O_H4[:4, :4, :4] = H4
    P5_unit = build_b23_probe(O_H4, "M8")

    P3c_sx_unit = build_b23_probe(O8, "sigma_x")

    # ------------------------------------------------------------------
    # Part 1: 局在質量 P4 の振幅 ε スイープ -> 局所時間遅延 & 背景計量変形
    # ------------------------------------------------------------------
    log("\n" + "-" * 76)
    log("[Part 1] 局在質量 P4 (A=0) 振幅 ε スイープ：局所固有時間遅延 & 背景計量変形")
    log("-" * 76)
    log(
        f"{'eps':>6s} | {'E_P_tot':>9s} {'E_core':>9s} {'E_cloud':>9s} | "
        f"{'bg_E(0)/E(far)':>14s} | {'phi_bg(0)/phi0':>14s} {'phi_bg(far)/phi0':>16s} | "
        f"{'phi_P(0)/phi0':>13s} {'d1_P(0)':>8s}"
    )

    eps_list = [0.001, 0.01, 0.03, 0.05, 0.08, 0.12, 0.16, 0.20, 0.25, 0.30, 0.40, 0.50]
    p1_records = []
    U_p4_ref = None  # eps=0.20 の τ=30 状態を Part 4 曲率測定用に保持

    for eps in eps_list:
        Ut = H + eps * P4_unit
        Ut /= np.linalg.norm(Ut)
        bg_E0_l, bg_Ef_l = [], []
        bg_phi0_l, bg_phif_l = [], []
        p_phi0_l, p_d10_l = [], []
        E_core_l, E_cloud_l = [], []
        for st in range(1, 41):
            Un, _ = step_tri(Ut)
            if st == 30 and abs(eps - 0.20) < 1e-9:
                U_p4_ref = Un.copy()
            if st >= 20:
                E_bg, d1_bg, phi_bg = site_d1_and_energy(Ut, Un, slice(0, 2))
                E_pr, d1_pr, phi_pr = site_d1_and_energy(Ut, Un, slice(2, 4))
                bg_E0_l.append(E_bg[0])
                bg_Ef_l.append(float(np.mean(E_bg[1:])))
                bg_phi0_l.append(phi_bg[0])
                bg_phif_l.append(float(np.mean(phi_bg[1:])))
                p_phi0_l.append(phi_pr[0])
                p_d10_l.append(d1_pr[0])
                E_core_l.append(E_pr[0])
                E_cloud_l.append(float(np.sum(E_pr[1:])))
            Ut = Un
        e_core = float(np.mean(E_core_l))
        e_cloud = float(np.mean(E_cloud_l))
        e_tot = e_core + e_cloud
        g_ratio = float(np.mean(bg_E0_l) / np.mean(bg_Ef_l))
        r_bg0 = float(np.mean(bg_phi0_l) / phi0_ref)
        r_bgf = float(np.mean(bg_phif_l) / phi0_ref)
        r_p0 = float(np.mean(p_phi0_l) / phi0_ref)
        d1_p0 = float(np.mean(p_d10_l))
        p1_records.append((eps, e_tot, e_core, e_cloud, g_ratio, r_bg0, r_bgf, r_p0, d1_p0))
        log(
            f"{eps:6.3f} | {e_tot:9.6f} {e_core:9.6f} {e_cloud:9.6f} | "
            f"{g_ratio:14.6f} | {r_bg0:14.6f} {r_bgf:16.6f} | "
            f"{r_p0:13.6f} {d1_p0:8.4f}"
        )

    # ------------------------------------------------------------------
    # Part 2: スピン 2 (グラビトン候補 σ_x) と局在粒子 P4 の結合・時計変調
    # ------------------------------------------------------------------
    log("\n" + "-" * 76)
    log("[Part 2] スピン 2 グラビトン候補 (σ_x) の背景時計変調と P4 粒子との結合")
    log("-" * 76)
    for label, U_init in [
        ("σ_x only (eps=0.05)", H + 0.05 * P3c_sx_unit),
        ("σ_x only (eps=0.10)", H + 0.10 * P3c_sx_unit),
        ("σ_x only (eps=0.20)", H + 0.20 * P3c_sx_unit),
        ("P4(0.10) only", H + 0.10 * P4_unit),
        ("P4(0.10) + σ_x(0.05)", H + 0.10 * P4_unit + 0.05 * P3c_sx_unit),
        ("P4(0.10) + σ_x(0.10)", H + 0.10 * P4_unit + 0.10 * P3c_sx_unit),
    ]:
        Ut = U_init / np.linalg.norm(U_init)
        bg_phi_l, pr_phi_l, E_pr_l, E_bg_l, d1_pr_l = [], [], [], [], []
        for st in range(1, 41):
            Un, _ = step_tri(Ut)
            if st >= 20:
                E_bg, d1_bg, phi_bg = site_d1_and_energy(Ut, Un, slice(0, 2))
                E_pr, d1_pr, phi_pr = site_d1_and_energy(Ut, Un, slice(2, 4))
                bg_phi_l.append(float(np.mean(phi_bg)))
                pr_phi_l.append(float(np.mean(phi_pr)))
                E_pr_l.append(float(np.sum(E_pr)))
                E_bg_l.append(float(np.sum(E_bg)))
                d1_pr_l.append(float(np.mean(d1_pr)))
            Ut = Un
        log(
            f"  {label:22s}: E_bg={np.mean(E_bg_l):.5f}, E_b23={np.mean(E_pr_l):.5f} | "
            f"phi_bg/phi0={np.mean(bg_phi_l)/phi0_ref:.5f}, phi_b23/phi0={np.mean(pr_phi_l)/phi0_ref:.5f}, d1_b23={np.mean(d1_pr_l):.4f}"
        )

    # ------------------------------------------------------------------
    # Part 3: 因果伝播速度上限 c (ハミング距離 d_H=0,1,2,3) & 空間倍音 vs 固有時計速度
    # ------------------------------------------------------------------
    log("\n" + "-" * 76)
    log("[Part 3A] 信号伝播の因果構造：A=0 局所パルス (P4, eps=0.05) のハミング距離別エネルギー伝播 (τ=0..5)")
    log("-" * 76)
    Ut = H + 0.05 * P4_unit
    Ut /= np.linalg.norm(Ut)
    ham_dist = np.array([bin(A).count("1") for A in range(8)], dtype=int)
    p3a_rows = []
    for st in range(0, 6):
        T6 = Ut.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, 2:4, :]
        eA = np.sum(T6 ** 2, axis=(1, 2, 3, 4, 5))
        tot = float(np.sum(eA))
        by_d = [float(np.sum(eA[ham_dist == d]) / tot) if tot > 0 else 0.0 for d in range(4)]
        p3a_rows.append((st, tot, by_d))
        log(
            f"  tau={st}: E_b23={tot:.6e} | frac(d_H=0)={by_d[0]:.5f}, "
            f"d_H=1(3sites)={by_d[1]:.5f}, d_H=2(3sites)={by_d[2]:.5f}, d_H=3(1site)={by_d[3]:.5f}"
        )
        Ut, _ = step_tri(Ut)

    log("\n" + "-" * 76)
    log("[Part 3B] 空間倍音占有率 f_harm (空間運動量) vs 粒子固有時計速度 phi_P/phi_0（ローレンツ関係検証）")
    log("-" * 76)
    # 局在度を連続的に変えたファミリー: P(alpha) = (1-alpha)*O8 + alpha*delta_000 (alpha in [0, 1])
    log(
        f"{'alpha':>6s} | {'f_harm(k>=1)':>12s} {'sqrt(1-f_harm)':>14s} | "
        f"{'phi_P/phi0':>10s} {'cos_rel':>9s} {'v_space_est':>11s} | {'(phi/phi0)^2+f_h':>16s}"
    )
    p3b_records = []
    alphas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]
    for alpha in alphas:
        O_mix = (1.0 - alpha) * (O8 / np.linalg.norm(O8)) + alpha * O_d0
        P_mix = build_b23_probe(O_mix, "M8")
        Ut = H + 0.10 * P_mix
        Ut /= np.linalg.norm(Ut)
        phi_list = []
        fh_list = []
        for st in range(1, 36):
            Un, _ = step_tri(Ut)
            if st >= 20:
                E_pr, _, phi_pr = site_d1_and_energy(Ut, Un, slice(2, 4))
                # エネルギー重み付き固有時計速度
                w_phi = float(np.sum(E_pr * phi_pr) / np.sum(E_pr))
                phi_list.append(w_phi)
                # 空間倍音占有率 (振幅 sqrt(E_pr) の WHT)
                p_x = E_pr / np.sum(E_pr)
                amp_k = H8 @ np.sqrt(np.maximum(p_x, 0.0))
                p_k = amp_k ** 2
                fh_list.append(float(1.0 - p_k[0]))
            Ut = Un
        fh = float(np.mean(fh_list))
        r_phi = float(np.mean(phi_list) / phi0_ref)
        sq = float(np.sqrt(max(0.0, 1.0 - fh)))
        sum_sq = r_phi ** 2 + fh
        p3b_records.append((alpha, fh, sq, r_phi, sum_sq))
        log(
            f"{alpha:6.2f} | {fh:12.6f} {sq:14.6f} | "
            f"{r_phi:10.6f} {r_phi/sq if sq>0 else 0:9.5f} {np.sqrt(max(0, 1 - r_phi**2)):11.6f} | {sum_sq:16.6f}"
        )

    # ------------------------------------------------------------------
    # Part 4: 離散曲率（三角形視野角 ρ = L3/L1）の質量直上 (A=0) vs 遠方 (A=1..7) 比較
    # ------------------------------------------------------------------
    log("\n" + "-" * 76)
    log("[Part 4] 離散曲率：質量 P4 (eps=0.20, τ=30) 直上 A=0 vs 遠方 A=1,7 の視野角 ρ = L3/L1 比較")
    log("-" * 76)
    # 背景 H の代表セル (A, B, C) = (0,0,0), (1,0,1), (7,0,7) × 内側 (4,0,4)
    for tag, U_tgt in [("Vacuum H (no mass)", H), ("With Mass P4 (eps=0.20, tau=30)", U_p4_ref)]:
        log(f"  [{tag}]")
        for A_site in [0, 1, 7]:
            i_idx = A_site * 8 + 4
            j_idx = 0 * 8 + 0
            k_idx = A_site * 8 + 4
            cnt, r_mean, r_wmean, r_min = measure_view_angle_stats(U_tgt, i_idx, j_idx, k_idx)
            log(f"    Site A={A_site} ({i_idx:2d},{j_idx:2d},{k_idx:2d}): n_tri={cnt:4d}, mean_rho={r_mean:.6f}, w_mean_rho={r_wmean:.6f}, min_rho={r_min:.6f}")

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.2f}s.")

    with open("results/tri_loop_n64_gravity_time_lorentz_v0.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    np.savez(
        "results/tri_loop_n64_gravity_time_lorentz_v0.npz",
        phi0_ref=phi0_ref,
        d1_0_ref=d1_0_ref,
        p1_records=np.array(p1_records),
        p3b_records=np.array(p3b_records),
    )


if __name__ == "__main__":
    run()
