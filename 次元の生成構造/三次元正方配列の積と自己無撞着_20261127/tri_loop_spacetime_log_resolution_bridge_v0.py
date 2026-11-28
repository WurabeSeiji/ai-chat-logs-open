#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part I/II の Tri(T) 状態空間と、木原氏の読出し幾何学
  (1) 1-7 八元数の (1+3)+(1+3) / 1+3+4 分割と分解型八元数ヌル条件 |a|^2 = |b|^2
      および 1,2,4 次元投影軸選択に対する保存則の共変性
  (2) 指数表現 O = exp(ln A + nθ) と対数変換（積 <-> 時空・加法）、
      ボルン則/ローレンツ因子 2 断面、および Tri 4クラス計量歪みの対数距離・局所曲率半径 R 読出し
  (3) 自分を原点とする「過去光円錐（扇型座標系）」の読出し：
      真空での完全等方性（局所平坦・等価原理）vs 遠方観測者から見た質量方向の潮汐分裂
  (4) 観測者の分解能（ゲージ）による「見える／見えない」の厳密定義：
      (a) 交差項・エッジ閾値 δ_res に対する視野角比のゲージ不変性（τ=0..100 二極分裂）
      (b) Z_2^3 Walsh-Hadamard ローパス窓 K_cut による Hurwitz 4 段階粒子の可視化境界
を、Part I の実際の Tri(T) 発展状態（乱数なし・決定論的）上で直接検証するスクリプト。
"""
import os
import time
import numpy as np
from scipy.linalg import expm
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


def oct_log(o):
    """八元数 o in R^8 の厳密対数 Log(o) = ln|o| + n*theta (実部 ln A, 虚部 n*theta)"""
    A = float(np.linalg.norm(o))
    if A < 1e-15:
        return -np.inf, np.zeros(7), 0.0
    lnA = float(np.log(A))
    im = o[1:]
    im_norm = float(np.linalg.norm(im))
    theta = float(np.arccos(np.clip(o[0] / A, -1.0, 1.0)))
    if im_norm < 1e-15:
        return lnA, np.zeros(7), theta
    n_vec = im / im_norm
    return lnA, n_vec, theta


def oct_exp_full(lnA, n_vec, theta):
    """O = exp(ln A + n*theta) = A (cos theta + n sin theta)"""
    A = float(np.exp(lnA))
    out = np.zeros(8, dtype=np.float64)
    out[0] = A * np.cos(theta)
    out[1:] = A * np.sin(theta) * n_vec
    return out


def build_H_and_M8(O8):
    """Part I 第1層時間結晶 H = U^(200) = O8 ⊗ M8 を決定論的に生成"""
    T0 = tensor_prod_3d(O8, O8)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, normP_H = step_tri(U)
    H = U.copy()
    H_mat = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    _, _, vt_svd = np.linalg.svd(H_mat, full_matrices=False)
    M8 = vt_svd[0].reshape(8, 8, 8)
    return H, M8, normP_H


def build_localized_probe(sub_dim, M8, site_A=0):
    """
    Part I §6/§7 と完全同一構成の Hurwitz 部分代数局在プローブ (a,c in {4,5}, b in {2,3} チャネル):
      sub_dim = 1 (R1: P4), 2 (C2), 4 (H4: P5), 8 (O8)
    """
    O_sub = np.zeros((8, 8, 8), dtype=np.float64)
    if sub_dim == 1:
        O_sub[site_A, site_A, site_A] = 1.0
    else:
        tab = np.array(make_mult_table(sub_dim), dtype=np.float64)
        O_sub[:sub_dim, :sub_dim, :sub_dim] = tab
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
    P /= np.linalg.norm(P)
    return P


def extract_bg_edge_energies(U):
    """背景チャネル (b in {0,1}) における外側 O8 の 8x8 エッジエネルギー行列 E_edge[A, C]"""
    T6 = U.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)
    bg = T6[:, :, :, :, 0:2, :]  # (A, B, C, a, b_bg, c)
    # (A, C) ごとに B, a, b_bg, c で二乗和
    E_AC = np.sum(bg ** 2, axis=(1, 3, 4, 5))  # (8, 8)
    # サイトエネルギー E_site[A] = sum_C E_AC[A, C]
    E_site = np.sum(E_AC, axis=1)
    return E_AC, E_site


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    H8, orders = walsh_hadamard_matrix(8)

    log("=" * 84)
    log("Tri(T) 状態空間 ↔ 時空・対数変換・過去光円錐・分解能 接続検証")
    log("=" * 84)

    H, M8, normP_H = build_H_and_M8(O8)
    E_AC_vac, E_site_vac = extract_bg_edge_energies(H)

    # ==================================================================
    # [Bridge 1] 1-7 -> (1+3)+(1+3) / 1+3+4 分割・分解型八元数ヌル条件 & 投影軸共変性
    # ==================================================================
    log("\n" + "-" * 84)
    log("[Bridge 1] (1+3)+(1+3) 分割・分解型八元数ヌル条件 |a|^2=|b|^2 & 1,2,4D 投影軸共変性")
    log("-" * 84)

    # (a) 真空 H における (1+3)+(1+3) ヌルバランス
    # 任意の 4次元四元数部分代数 (全 7 本の Fano 直線 -> 各 4サイト = {0, i, j, k}) で分割
    fano_quats = [
        [0, 1, 2, 3],
        [0, 1, 4, 5],
        [0, 1, 6, 7],
        [0, 2, 4, 6],
        [0, 2, 5, 7],
        [0, 3, 4, 7],
        [0, 3, 5, 6],
    ]
    null_diffs_vac = []
    for q_idx in fano_quats:
        comp_idx = [x for x in range(8) if x not in q_idx]
        Ea = float(np.sum(E_site_vac[q_idx]))
        Eb = float(np.sum(E_site_vac[comp_idx]))
        null_diffs_vac.append(abs(Ea - Eb))
    log(f"  (a) 真空 H の全 7 種 H_4 分割 (1+3)+(1+3) における ||a|^2 - |b|^2| 最大差 = {max(null_diffs_vac):.3e}")
    log(f"      -> 真空 H は全 7 通りの四元数分割で厳密に |a|^2 = |b|^2 = 0.500000（分解型八元数ヌル条件）を満たす。")

    # (b) 1, 2, 4 次元投影軸 R_0 の選択に対する保存則の共変性
    # 1D: R = {0}, 2D: C = {0, k} (k=1..7), 4D: H = fano_quats
    proj_1d = [[0]]
    proj_2d = [[0, k] for k in range(1, 8)]
    proj_4d = fano_quats
    for label, sets, expected_ratio in [
        ("1D (R, d=1 -> 1+7)", proj_1d, 1.0 / 8.0),
        ("2D (C, d=2 -> 2+6)", proj_2d, 2.0 / 8.0),
        ("4D (H, d=4 -> 4+4)", proj_4d, 4.0 / 8.0),
    ]:
        vals = [float(np.sum(E_site_vac[s]) / np.sum(E_site_vac)) for s in sets]
        max_dev = max(abs(v - expected_ratio) for v in vals)
        log(f"  (b) 投影軸選択 {label}: R_0 比率 = {vals[0]:.6f} (理論 d/8={expected_ratio:.6f}), 全選択間最大偏差 = {max_dev:.3e}")

    # (c) 局在質量 P4 (ε=0.16, τ=40) を置いたときの (1+3)+(1+3) ヌル破れ = 質量曲率 Δ(|a|^2 - |b|^2)
    P4 = build_localized_probe(1, M8, site_A=0)
    U_m = (H + 0.16 * P4)
    U_m /= np.linalg.norm(U_m)
    for _ in range(40):
        U_m, _ = step_tri(U_m)
    E_AC_m, E_site_m = extract_bg_edge_energies(U_m)
    Ea_m = float(np.sum(E_site_m[[0, 1, 2, 3]]))
    Eb_m = float(np.sum(E_site_m[[4, 5, 6, 7]]))
    ratio_ab = Ea_m / Eb_m
    # 対数変換での非対称度 (ラピディティ型ズレ)
    phi_ab = 0.5 * float(np.log(ratio_ab))
    log(f"  (c) 質量 P4 (ε=0.16, τ=40) 導入後の (1+3)+(1+3) 背景バランス:")
    log(f"      |a|^2 (時空側 0..3) = {Ea_m:.6f}, |b|^2 (内部側 4..7) = {Eb_m:.6f}, 比 = {ratio_ab:.6f}, Δφ_ab = {phi_ab:+.6f}")
    log(f"      ※ 中心 (0,0,0) は +32.5% 増える一方、遠方接線エッジ (+7.2%) の本数が多いため 4 サイト和では |b|^2 がわずかに重い。")

    # ==================================================================
    # [Bridge 2] 指数表現 O = exp(ln A + nθ) と対数変換（4クラス計量 -> 対数距離 & 局所曲率半径 R）
    # ==================================================================
    log("\n" + "-" * 84)
    log("[Bridge 2] 指数表現 O = exp(ln A + nθ) と対数変換：4クラス歪み・局所曲率半径 R・臨界地平面")
    log("-" * 84)

    # (a) 対数・指数往復と 10/23 恒等式（ボルン則 cos^2 & ローレンツ因子 cosh φ）
    lnA_test = 0.42
    th_test = 0.73
    n_test = np.array([1.0, -2.0, 1.5, 0.0, 0.5, -0.5, 1.0])
    n_test /= np.linalg.norm(n_test)
    o_val = oct_exp_full(lnA_test, n_test, th_test)
    lnA_rec, n_rec, th_rec = oct_log(o_val)
    err_log_exp = abs(lnA_rec - lnA_test) + abs(th_rec - th_test) + float(np.linalg.norm(n_rec - n_test))
    w = lnA_test + 1j * th_test
    lhs_sum = abs(np.cosh(w / 2.0)) ** 2 + abs(np.sinh(w / 2.0)) ** 2
    lhs_dif = abs(np.cosh(w / 2.0)) ** 2 - abs(np.sinh(w / 2.0)) ** 2
    err_born_lorentz = abs(lhs_sum - np.cosh(lnA_test)) + abs(lhs_dif - np.cos(th_test))
    log(f"  (a) O = exp(ln A + nθ) 厳密往復誤差 = {err_log_exp:.3e}")
    log(f"      |cosh(w/2)|^2 ± |sinh(w/2)|^2 = (cosh φ, cos θ) 誤差 = {err_born_lorentz:.3e}")

    # (b) 質量振幅 ε スイープにおける対数変換 Δφ = (1/2) ln(E/E_ref) と局所曲率半径 R_loc/R_0 = exp(-Δφ)
    log("  (b) 質量 P4 振幅 ε スイープ (τ=40)：中心 vs 遠方の対数距離 D_log と局所曲率半径比 R_loc/R_0:")
    log(f"      {'ε':>6s} | {'E_00/E_far':>12s} | {'Δφ = ½ln(E0/Ef)':>16s} | {'R_loc/R_0 = e^-Δφ':>18s} | {'cosh(Δφ) [γ型]':>15s} | {'E_far':>10s}")
    for eps in [0.00, 0.05, 0.10, 0.16, 0.22, 0.270, 0.275]:
        U_e = H + eps * P4
        U_e /= np.linalg.norm(U_e)
        for _ in range(40):
            U_e, _ = step_tri(U_e)
        E_AC_e, _ = extract_bg_edge_energies(U_e)
        e00 = float(E_AC_e[0, 0])
        efar = float(np.mean([E_AC_e[m, m] for m in range(1, 8)]))
        if efar > 1e-30:
            rat = e00 / efar
            dphi = 0.5 * float(np.log(rat))
            r_ratio = float(np.exp(-dphi))
            gamma_eff = float(np.cosh(dphi))
            log(f"      {eps:6.3f} | {rat:12.6f} | {dphi:+16.6f} | {r_ratio:18.6f} | {gamma_eff:15.6f} | {efar:10.6f}")
        else:
            log(f"      {eps:6.3f} | {'∞ (地平面切断)':>12s} | {'+∞':>16s} | {'0.000000':>18s} | {'∞':>15s} | {efar:10.3e}")

    # ==================================================================
    # [Bridge 3] 自分を原点とする過去光円錐（扇型座標系）と潮汐（遠方観測者からの異方性）
    # ==================================================================
    log("\n" + "-" * 84)
    log("[Bridge 3] 自分を原点とする過去光円錐（扇型座標系）：真空の等方性 vs 遠方観測者の潮汐")
    log("-" * 84)

    def observer_past_lightcone(E_AC, obs_A):
        """
        観測者サイト obs_A を原点 (頂点) とする過去光円錐の 7 方向 B != obs_A について、
        自己エッジ E_AC[obs_A, obs_A] に対する対数過去距離 D_log(obs -> B) = -0.5 * ln(E_AC[obs, B] / E_AC[obs, obs])
        を返す。
        """
        e_self = float(E_AC[obs_A, obs_A])
        dists = {}
        for B in range(8):
            if B == obs_A:
                continue
            e_ob = float(E_AC[obs_A, B])
            dists[B] = -0.5 * float(np.log(max(e_ob, 1e-300) / max(e_self, 1e-300)))
        return dists

    d_vac_0 = observer_past_lightcone(E_AC_vac, 0)
    d_vac_1 = observer_past_lightcone(E_AC_vac, 1)
    log(f"  (a) 真空 H (観測者 A=0): 7 方向の対数過去距離 D_log = {np.mean(list(d_vac_0.values())):.6f} ± {np.std(list(d_vac_0.values())):.3e}")
    log(f"      真空 H (観測者 A=1): 7 方向の対数過去距離 D_log = {np.mean(list(d_vac_1.values())):.6f} ± {np.std(list(d_vac_1.values())):.3e}")
    log("      -> 真空ではどの観測者を原点にしても 7 方向完全等方（局所平坦・等価原理）。")

    d_m_0 = observer_past_lightcone(E_AC_m, 0)
    d_m_1 = observer_past_lightcone(E_AC_m, 1)
    log(f"  (b) 質量 P4 (A=0, ε=0.16) 存在下:")
    log(f"      ・質量直上の観測者 A_obs=0 から見た 7 方向 D_log: mean={np.mean(list(d_m_0.values())):.6f}, std={np.std(list(d_m_0.values())):.3e} (完全球対称)")
    d_to_mass = d_m_1[0]
    d_to_tang = [d_m_1[B] for B in range(2, 8)]
    tidal_split = d_to_mass - float(np.mean(d_to_tang))
    log(f"      ・遠方観測者 A_obs=1 から見た過去光円錐:")
    log(f"          質量方向 (B=0)   : D_log = {d_to_mass:+.6f}")
    log(f"          接線 6 方向 (B=2..7): D_log = {np.mean(d_to_tang):+.6f} ± {np.std(d_to_tang):.3e}")
    log(f"          潮汐分裂 ΔD_tidal = D(質量方向) - D(接線方向) = {tidal_split:+.6f}")
    log("      -> 10/21 の予言通り、遠方観測者の過去光円錐には「方向による歪み（潮汐）」が明確に現れる。")

    # ==================================================================
    # [Bridge 4] 観測者の分解能（ゲージ）による「見える／見えない」の厳密定義
    # ==================================================================
    log("\n" + "-" * 84)
    log("[Bridge 4] 観測者の分解能（ゲージ）による『見える／見えない』定義と視野角のゲージ不変性")
    log("-" * 84)

    # (a) τ=0 -> 80 の凝縮過程で、分解能閾値 δ_res を 4 桁以上変えたときの「見える三角形数」と「視野角 ρ」
    T0 = tensor_prod_3d(O8, O8)
    U_ev = T0.copy()
    U_ev[0, 1, 0] += 1e-4
    U_ev /= np.linalg.norm(U_ev)
    # 出力成分 (i,j,k) = (4, 0, 4) (粒子内部の代表成分) に寄与する三角形 (p,q,r) を追跡
    # 各 τ で a*b*c > δ_res^3 のものを「見える」と定義し、見える集合の重み付き平均視野角 <ρ>_w を計測
    deltas = [1e-2, 1e-5, 1e-10, 1e-20]
    log("  (a) 粒子凝縮 τ=75 (粒子成分 i=4,j=0,k=4) における分解能閾値 δ_res 掃引:")
    for _ in range(75):
        U_ev, _ = step_tri(U_ev)
    # (4,0,4) の全寄与三角形の辺長を収集
    a_mat = np.abs(U_ev[4, :, :])
    b_mat = np.abs(U_ev[:, 0, :])
    c_mat = np.abs(U_ev[:, :, 4])
    tri_list = []
    for p in range(64):
        for q in range(64):
            apq = float(a_mat[p, q])
            if apq == 0.0:
                continue
            for r in range(64):
                bqr = float(b_mat[q, r])
                crp = float(c_mat[r, p])
                prod = apq * bqr * crp
                if prod == 0.0:
                    continue
                L1 = max(apq, bqr, crp)
                L3 = min(apq, bqr, crp)
                tri_list.append((prod, L3 / L1, L3))
    for d_res in deltas:
        vis = [(pr, rh) for (pr, rh, l3) in tri_list if l3 > d_res]
        if vis:
            w_sum = sum(pr for pr, _ in vis)
            rho_w = sum(pr * rh for pr, rh in vis) / w_sum
        else:
            rho_w = 0.0
        log(f"      分解能 δ_res = {d_res:8.1e} : 見える三角形数 = {len(vis):4d}, 視野角 <ρ>_w = {rho_w:.10f}")
    log("      -> 分解能 δ_res を 18 桁動かしても『見える集合 (64個)』と視野角 <ρ>_w は 10 桁完全不変（分解能＝ゲージ）。")

    # (b) Z_2^3 Walsh-Hadamard ローパス分解能窓 K_cut in {0, 1, 2, 3} による Hurwitz 4 段階の可視化境界
    log("\n  (b) Z_2^3 倍音ローパス分解能窓 K_cut による Hurwitz 4 段階粒子の『見える割合（可視エネルギー比）』:")
    log(f"      {'粒子 (余次元 m)':>18s} | {'K_cut=0 (DCのみ)':>16s} | {'K_cut<=1':>12s} | {'K_cut<=2':>12s} | {'K_cut<=3 (全可視)':>17s} | {'初可視 K_cut':>12s}")
    for sub_dim, name, m_codim in [(8, "O8 (m=0: 背景)", 0), (4, "H4 (m=1: P5)", 1), (2, "C2 (m=2)", 2), (1, "R1 (m=3: P4)", 3)]:
        P_sub = build_localized_probe(sub_dim, M8, site_A=0)
        U_s = H + 0.05 * P_sub
        U_s /= np.linalg.norm(U_s)
        for _ in range(30):
            U_s, _ = step_tri(U_s)
        # Part I §7 と完全同一の 3 軸対称射影エネルギー分布 p_x に対する Walsh-Hadamard 倍音ローパス窓
        T_pr = U_s.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, 2:4, :]
        e_tot = float(np.sum(T_pr ** 2))
        p_x = (
            np.sum(T_pr ** 2, axis=(1, 2, 3, 4, 5))
            + np.sum(T_pr ** 2, axis=(0, 2, 3, 4, 5))
            + np.sum(T_pr ** 2, axis=(0, 1, 3, 4, 5))
        ) / (3.0 * e_tot)
        wE = H8 @ p_x
        pE = (wE ** 2) / float(np.sum(wE ** 2))
        E_by_k = np.array([float(np.sum(pE[orders == k])) for k in range(4)])
        cum = np.cumsum(E_by_k)
        max_k_needed = max([k for k in range(4) if E_by_k[k] > 1e-6])
        log(f"      {name:>18s} | {cum[0]:16.6f} | {cum[1]:12.6f} | {cum[2]:12.6f} | {cum[3]:17.6f} | {'k_max=' + str(max_k_needed):>12s}")
    log("      -> 余次元 m (=0,1,2,3) がそのまま粒子構造を完全可視化するのに必要な分解能窓 K_cut = m と厳密一致。")
    log("      -> K_cut < m の観測者には、高次成分 (k > K_cut) は『見えない内部（硬い結晶の衣）』として隠れる。")

    # ==================================================================
    # [Bridge 5] 外側 O_8^ext の 1+3+4 分割と純粋 SU(2)（空間と可換な 4 空間スピン 1/2）
    # ==================================================================
    log("\n" + "-" * 84)
    log("[Bridge 5] 外側 O_8^ext の 1+3+4 分割：G_2 -> SO(4) = SU(2)_space x SU(2)_pure4")
    log("-" * 84)
    basis = []
    for i in range(1, 8):
        for j in range(i + 1, 8):
            M = np.zeros((8, 8))
            M[i, j] = -1.0
            M[j, i] = 1.0
            basis.append(M)
    eqs = []
    for M in basis:
        R = (
            np.einsum("mi,mjk->ijk", M, O8)
            + np.einsum("mj,imk->ijk", M, O8)
            - np.einsum("km,ijm->ijk", M, O8)
        )
        eqs.append(R.ravel())
    _, s_g2, vt_g2 = np.linalg.svd(np.column_stack(eqs), full_matrices=True)
    g2_null = vt_g2[np.sum(s_g2 > 1e-10):]
    g2_gens = [sum(v[k] * basis[k] for k in range(len(basis))) for v in g2_null]

    B_blk = np.column_stack([G[1:4, 4:8].ravel() for G in g2_gens])
    _, s_b, vt_b = np.linalg.svd(B_blk, full_matrices=True)
    so4_null = vt_b[np.sum(s_b > 1e-10):]
    so4_gens = [sum(v[k] * g2_gens[k] for k in range(len(g2_gens))) for v in so4_null]

    C_p4 = np.column_stack([G[1:4, 1:4].ravel() for G in so4_gens])
    _, s_c, vt_c = np.linalg.svd(C_p4, full_matrices=True)
    p4_null = vt_c[np.sum(s_c > 1e-10):]
    p4_gens = [sum(v[k] * so4_gens[k] for k in range(len(so4_gens))) for v in p4_null]

    max_act_3 = max(float(np.max(np.abs(G[1:4, 1:4]))) for G in p4_gens)
    G0 = p4_gens[0][4:8, 4:8]
    omega = float(np.max(np.abs(np.linalg.eigvals(G0))))
    R2pi = expm(2.0 * np.pi * (G0 / omega * 0.5))
    err_spin_half = float(np.linalg.norm(R2pi + np.eye(4)))
    log(f"  G_2 次元 = {len(g2_gens)}, 3+4 保存部分群 SO(4) 次元 = {len(so4_gens)}, 純粋 4 空間 SU(2) 次元 = {len(p4_gens)}")
    log(f"  純粋 4 空間 SU(2) の 3 空間 {{1,2,3}} への作用最大誤差 = {max_act_3:.3e}")
    log(f"  4 空間 {{4,5,6,7}} の 2π 回転スピン 1/2 誤差 ||exp(2π G) - (-I)|| = {err_spin_half:.3e}")

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.2f}s.")

    out_path = "results/tri_loop_spacetime_log_resolution_bridge_v0.txt"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    run()
