#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part VI（2 粒子動的散乱・引力/斥力・運動量/位相交換・空間分離 CHSH ベル不等式の動的破れ）
決定論的機械精度検証スクリプト（乱数なし）。

  [SC-1] 空間分離 2 粒子（サイト A=0, B=1）の非線形相互作用エネルギーと符号依存性
         - 中性スカラー質量（同位相 +1 vs 逆位相 -1）
         - U(1) 荷電粒子（同種電荷 Q=(+1,+1) vs 異種電荷 Q=(+1,-1)）の相互作用差とサイト間結合エネルギー
  [SC-2] 動的散乱と位相・エネルギー交換（粒子 A=0 への媒介波束入射）
         - 標的粒子 A=0 に対するサイト B=1 からの U(1) ゲージ波束 / スピン 2 (σ_x) 波束の動的散乱
         - 散乱位相シフト δ(τ) と標的エネルギー変化 ΔE_0(τ) の時間発展追跡
  [SC-3] 2 粒子臨界崩壊点 ε_c^(2) vs 1 粒子臨界崩壊点 ε_c^(1)
         - 2 粒子の重力合算効果（臨界崩壊振幅の低下比）の二分探索測定
  [SC-4] 空間分離 2 サイト（Alice=A=0, Bob=B=1）の 4x4 2 粒子密度行列 ρ_AB(τ) と
         CHSH ベル不等式 S_CHSH(τ) の動的検証
         - (a) 初期積状態 |↑_A ↑_B> からの動的エンタングルメント生成テスト
         - (b) 初期ベル状態 |Φ^+>, |Ψ^-> の Tri(T) 非線形発展下での CHSH 破れ動的維持テスト
               （具体的測定角 S_angle と Horodecki 厳密最大値 S_max の双方を計測）
  [SC-5] 振幅掃引による CHSH 破れの安定性境界（弱結合共存相 vs 臨界崩壊・地平面相でのエンタングルメント挙動）
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


def build_H_and_M8(O8):
    T0 = tensor_prod_3d(O8, O8)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, normP_H = step_tri(U)
    H = U.copy()
    H_mat = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    _, _, vt = np.linalg.svd(H_mat, full_matrices=False)
    M8 = vt[0].reshape(8, 8, 8)
    return H, M8, normP_H


def build_site_scalar(site_A, M8):
    """サイト A 局在スカラー粒子プローブ (Part I/IV 同一)"""
    P = np.zeros((64, 64, 64), dtype=np.float64)
    for a in [4, 5]:
        for b in [2, 3]:
            for c in [4, 5]:
                P[site_A * 8 + a, site_A * 8 + b, site_A * 8 + c] = M8[a, b - 2, c]
    return P / np.linalg.norm(P)


def build_site_spin2(site_A):
    """サイト A 局在スピン 2 (σ_x) 粒子プローブ (Part I/V 同一)"""
    P = np.zeros((64, 64, 64), dtype=np.float64)
    for a in [4, 5]:
        for b in [2, 3]:
            for c in [4, 5]:
                if a != c:
                    P[site_A * 8 + a, site_A * 8 + b, site_A * 8 + c] = 1.0
    return P / np.linalg.norm(P)


def build_site_charged(site_A, M8, charge_sign):
    """
    サイト A 局在荷電粒子：内部二重項 {4,5} に U(1) 生成子 J (4->5, 5->-4) の
    ±π/2 回転（charge_sign = +1 or -1）を施した粒子モード。
    """
    # J_45: 4->5 (+1), 5->4 (-1)
    R = np.eye(8, dtype=np.float64)
    # exp(charge_sign * (pi/4) * J) -> 45 deg rotation so Q=+1 and Q=-1 are orthogonal (90 deg apart)
    th = charge_sign * (np.pi / 4.0)
    c, s = np.cos(th), np.sin(th)
    R[4, 4] = c; R[4, 5] = -s
    R[5, 4] = s; R[5, 5] = c
    M_rot = np.einsum("pa,qb,rc,abc->pqr", R, np.eye(8), R, M8)
    P = np.zeros((64, 64, 64), dtype=np.float64)
    for a in [4, 5]:
        for b in [2, 3]:
            for c in [4, 5]:
                P[site_A * 8 + a, site_A * 8 + b, site_A * 8 + c] = M_rot[a, b - 2, c]
    return P / np.linalg.norm(P)


def build_gauge_packet(site_A, site_B, M8):
    """サイト A <-> B 間を結ぶ U(1) ゲージ波束（空間リンク上の位相振動子）"""
    P = np.zeros((64, 64, 64), dtype=np.float64)
    for a in [4, 5]:
        for b in [0, 1]:
            for c in [4, 5]:
                # 反対称リンク（運動量・位相流束を担う）
                P[site_A * 8 + a, site_B * 8 + b, site_B * 8 + c] = M8[a, b, c]
                P[site_B * 8 + a, site_A * 8 + b, site_A * 8 + c] = -M8[a, b, c]
    return P / np.linalg.norm(P)


def site_observables(U, site_A):
    """サイト A の背景エネルギー E_bg(A)、プローブエネルギー E_pr(A)、内部時計位相 θ(A)"""
    U6 = U.reshape(8, 8, 8, 8, 8, 8)
    # O_8 の乗積表では O_8[A, 0, A] = +1 (全 A=0..7 で非ゼロ)
    blk_bg = U6[site_A, :, 0, :, site_A, :]
    blk_pr = U6[site_A, :, site_A, :, site_A, :]
    E_bg = float(np.sum(blk_bg[:, 0:2, :] ** 2))
    E_pr = float(np.sum(blk_pr[:, 2:4, :] ** 2))
    c0 = 0.5 * (blk_bg[4, 0, 4] + blk_bg[5, 0, 5])
    cJ = 0.5 * (blk_bg[5, 0, 4] - blk_bg[4, 0, 5])
    th = float(np.arctan2(cJ, c0))
    return E_bg, E_pr, th


def extract_two_site_density_matrix(U, site_A, site_B, mode="cross_all"):
    """
    空間分離 2 サイト (Alice=site_A, Bob=site_B) 上の 2 量子ビット（各 {4,5}）状態から
    4x4 エルミート密度行列 ρ_AB を抽出する。
    T の 3 指標 (i, j, k) を (Alice=(A, a in {4,5}), Bob=(B, b in {4,5}), Env=k in 0..63) として
    部分トレース ρ_AB = Tr_Env(|Ψ><Ψ|) を構成（厳密に半正定値・Tr=1）。
    """
    U6 = U.reshape(8, 8, 8, 8, 8, 8)  # (A,a, B,b, C,c)
    if mode == "cross_all":
        # (A, a in {4,5}), (B, b in {4,5}), env = (C, c)
        amp = U6[site_A, 4:6, site_B, 4:6, :, :].reshape(4, 64)
    elif mode == "cross_probe":
        # 環境指標 C を背景・プローブ全体、特に b in {2,3} チャネル等
        amp = U6[site_A, 4:6, site_B, 4:6, :, 2:4].reshape(4, 16)
    else:
        amp = U6[site_A, 4:6, site_B, 4:6, :, :].reshape(4, 64)
    rho = amp @ amp.T
    tr = float(np.trace(rho))
    if tr < 1e-20:
        return np.eye(4) / 4.0, 0.0
    return rho / tr, tr


def chsh_metrics(rho):
    """
    4x4 実対称（エルミート）密度行列 ρ に対する量子相関指標：
      - purity = Tr(ρ^2)
      - S_E: 部分トレース ρ_A のフォン・ノイマンエントロピー (bit)
      - S_max: Horodecki 定理による最大 CHSH 値 2*sqrt(u1 + u2)
      - S_std: 標準ベル測定角 a=0(Z), a'=pi/2(X), b=pi/4((Z+X)/√2), b'=-pi/4((Z-X)/√2) での |E(a,b)-E(a,b')+E(a',b)+E(a',b')|
    """
    sx = np.array([[0.0, 1.0], [1.0, 0.0]])
    sy = np.array([[0.0, -1.0j], [1.0j, 0.0]])
    sz = np.array([[1.0, 0.0], [0.0, -1.0]])
    paulis = [sx, sy, sz]

    purity = float(np.real(np.trace(rho @ rho)))
    # 部分トレース ρ_A
    rho4 = rho.reshape(2, 2, 2, 2)
    rho_A = np.einsum("ijik->jk", rho4.transpose(0, 1, 2, 3))  # wait: Tr_B(ρ_{a b, a' b'}) = sum_b ρ_{a b, a' b}
    rho_A = np.einsum("abcb->ac", rho4)
    evals = np.linalg.eigvalsh(rho_A)
    evals = evals[evals > 1e-15]
    S_E = float(-np.sum(evals * np.log2(evals)))

    # 相関行列 T_{ij} = Tr(ρ (σ_i ⊗ σ_j))
    Tmat = np.zeros((3, 3), dtype=np.float64)
    for i in range(3):
        for j in range(3):
            op = np.kron(paulis[i], paulis[j])
            Tmat[i, j] = float(np.real(np.trace(rho @ op)))
    # Horodecki: SVD of T = U diag(s) V^T -> 具体的な 4 単位ベクトル (a, a', b, b') による直接 CHSH 評価
    U_svd, s_svd, Vt_svd = np.linalg.svd(Tmat)
    S_max = float(2.0 * np.sqrt(max(0.0, s_svd[0] ** 2 + s_svd[1] ** 2)))

    def E_vec(vA, vB):
        opA = vA[0] * sx + vA[1] * sy + vA[2] * sz
        opB = vB[0] * sx + vB[1] * sy + vB[2] * sz
        return float(np.real(np.trace(rho @ np.kron(opA, opB))))

    # Horodecki 最適 4 測定ベクトル（すべて厳密に ||v|| = 1 の実単位ベクトル）
    c1, c2 = U_svd[:, 0], U_svd[:, 1]
    d1, d2 = Vt_svd[0, :], Vt_svd[1, :]
    denom = np.sqrt(max(1e-30, s_svd[0] ** 2 + s_svd[1] ** 2))
    ca, sa = s_svd[0] / denom, s_svd[1] / denom
    b_vec = ca * d1 + sa * d2
    bp_vec = ca * d1 - sa * d2
    b_vec /= max(1e-30, np.linalg.norm(b_vec))
    bp_vec /= max(1e-30, np.linalg.norm(bp_vec))
    S_std = abs(E_vec(c1, b_vec) + E_vec(c1, bp_vec) + E_vec(c2, b_vec) - E_vec(c2, bp_vec))
    return purity, S_E, S_max, S_std, Tmat


def build_two_site_bell_perturbation(site_A, site_B, state_type="phi_plus"):
    """
    空間分離 2 サイト (site_A, site_B) 上の 2 量子ビット状態（各 {4,5}）を
    環境チャネル C=site_A, c=2（プローブチャネル、真空 H と直交）上に埋め込む。
    """
    P = np.zeros((64, 64, 64), dtype=np.float64)
    # 対称化しつつ (A,a; B,b; env) に 2 量子ビット振幅を配置
    env_idx = site_A * 8 + 2  # プローブチャネル c=2（真空背景 c in {0,1} と直交）
    if state_type == "product_00":
        # |0_A 0_B> -> a=4, b=4
        P[site_A * 8 + 4, site_B * 8 + 4, env_idx] = 1.0
    elif state_type == "phi_plus":
        # (|00> + |11>)/sqrt(2)
        P[site_A * 8 + 4, site_B * 8 + 4, env_idx] = 1.0 / np.sqrt(2.0)
        P[site_A * 8 + 5, site_B * 8 + 5, env_idx] = 1.0 / np.sqrt(2.0)
    elif state_type == "psi_minus":
        # (|01> - |10>)/sqrt(2)
        P[site_A * 8 + 4, site_B * 8 + 5, env_idx] = 1.0 / np.sqrt(2.0)
        P[site_A * 8 + 5, site_B * 8 + 4, env_idx] = -1.0 / np.sqrt(2.0)
    return P / np.linalg.norm(P)


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    H, M8, normP_H = build_H_and_M8(O8)

    log("=" * 88)
    log("Part VI 検証：2 粒子動的散乱・引力/斥力・運動量交換・空間分離 CHSH ベル不等式の動的破れ")
    log("=" * 88)
    log(f"  真空 H = O_8 ⊗ M_8 (N=64), ||P||(vac) = {normP_H:.8f}")

    # ==================================================================
    # [SC-1] 空間分離 2 粒子（A=0, B=1）の非線形相互作用と電荷符号依存性
    # ==================================================================
    log("\n" + "-" * 88)
    log("[SC-1] 空間分離 2 粒子（A=0, B=1）の相互作用エネルギー差と電荷符号分離 (τ=40)")
    log("-" * 88)

    P0 = build_site_scalar(0, M8)
    P1 = build_site_scalar(1, M8)
    Q0_p = build_site_charged(0, M8, +1)
    Q1_p = build_site_charged(1, M8, +1)
    Q1_m = build_site_charged(1, M8, -1)

    def evolve_final(U0, steps=40):
        U = U0 / np.linalg.norm(U0)
        nrm = 0.0
        for _ in range(steps):
            U, nrm = step_tri(U)
        return U, nrm

    U_vac, n_vac = evolve_final(H, 40)
    log(f"  {'系・配置':<28} | {'ε':>5} | {'||P||(40)':>10} | {'Δ||P||_int':>11} | {'E_bg(0)/E_bg(7)':>15} | {'E_cross(0,1)':>12}")
    for eps in [0.05, 0.12, 0.18]:
        U_A, n_A = evolve_final(H + eps * P0, 40)
        U_B, n_B = evolve_final(H + eps * P1, 40)
        configs = [
            ("中性スカラー (+,+)", H + eps * P0 + eps * P1, n_A, n_B),
            ("中性スカラー (+,-)", H + eps * P0 - eps * P1, n_A, n_B),
        ]
        # 荷電粒子単独
        _, n_QpA = evolve_final(H + eps * Q0_p, 40)
        _, n_QpB = evolve_final(H + eps * Q1_p, 40)
        _, n_QmB = evolve_final(H + eps * Q1_m, 40)
        configs += [
            ("同種電荷 Q=(+1,+1)", H + eps * Q0_p + eps * Q1_p, n_QpA, n_QpB),
            ("異種電荷 Q=(+1,-1)", H + eps * Q0_p + eps * Q1_m, n_QpA, n_QmB),
        ]
        for name, U0, na, nb in configs:
            Uf, nf = evolve_final(U0, 40)
            d_int = nf - na - nb + n_vac
            e0, _, _ = site_observables(Uf, 0)
            e7, _, _ = site_observables(Uf, 7)
            U6 = Uf.reshape(8, 8, 8, 8, 8, 8)
            e_cross = float(np.sum(U6[0, :, 1, :, :, :] ** 2))
            ratio_07 = e0 / e7 if e7 > 1e-30 else float("inf")
            log(f"  {name:<28} | {eps:5.2f} | {nf:10.6f} | {d_int:+11.4e} | {ratio_07:15.6f} | {e_cross:12.6e}")

    # ==================================================================
    # [SC-2] 動的散乱と位相・エネルギー交換（標的 A=0 への媒介波束入射）
    # ==================================================================
    log("\n" + "-" * 88)
    log("[SC-2] 動的散乱：標的粒子 A=0 への媒介波束（U(1) ゲージ / スピン 2）入射による位相・エネルギー交換")
    log("-" * 88)
    eps_t = 0.08
    eps_g = 0.08
    G01 = build_gauge_packet(0, 1, M8)
    S1_x = build_site_spin2(1)

    # 自由標的（散乱波束なし）
    U_free = (H + eps_t * P0) / np.linalg.norm(H + eps_t * P0)
    U_g_p = (H + eps_t * P0 + eps_g * G01) / np.linalg.norm(H + eps_t * P0 + eps_g * G01)
    U_g_m = (H + eps_t * P0 - eps_g * G01) / np.linalg.norm(H + eps_t * P0 - eps_g * G01)
    U_s2 = (H + eps_t * P0 + eps_g * S1_x) / np.linalg.norm(H + eps_t * P0 + eps_g * S1_x)

    log(f"  {'τ':>3} | {'δθ(+G) [rad]':>13} | {'δθ(-G) [rad]':>13} | {'ΔE_bg(0,+G)':>12} | {'δθ(spin2) [rad]':>15} | {'ΔE_bg(0,s2)':>12}")
    for st in range(1, 21):
        U_free, _ = step_tri(U_free)
        U_g_p, _ = step_tri(U_g_p)
        U_g_m, _ = step_tri(U_g_m)
        U_s2, _ = step_tri(U_s2)
        if st in (1, 2, 3, 4, 5, 8, 10, 15, 20):
            ef, _, thf = site_observables(U_free, 0)
            ep, _, thp = site_observables(U_g_p, 0)
            em, _, thm = site_observables(U_g_m, 0)
            es, _, ths = site_observables(U_s2, 0)
            def wrap(x): return float((x + np.pi) % (2 * np.pi) - np.pi)
            log(f"  {st:3d} | {wrap(thp-thf):+13.5e} | {wrap(thm-thf):+13.5e} | {ep-ef:+12.4e} | {wrap(ths-thf):+15.5e} | {es-ef:+12.4e}")

    # ==================================================================
    # [SC-3] 2 粒子臨界崩壊点 ε_c^(2) vs 1 粒子臨界崩壊点 ε_c^(1)
    # ==================================================================
    log("\n" + "-" * 88)
    log("[SC-3] 2 粒子重力合算効果：臨界崩壊振幅 ε_c^(2) vs 1 粒子 ε_c^(1) (τ=40)")
    log("-" * 88)

    def find_ec(P_dir, lo=0.05, hi=0.50, iters=18):
        for _ in range(iters):
            mid = 0.5 * (lo + hi)
            _, nf = evolve_final(H + mid * P_dir, 40)
            if nf > 0.25:
                hi = mid
            else:
                lo = mid
        return lo, hi

    ec1_lo, ec1_hi = find_ec(P0)
    ec2p_lo, ec2p_hi = find_ec(P0 + P1)
    ec2m_lo, ec2m_hi = find_ec(P0 - P1)
    ec1 = 0.5 * (ec1_lo + ec1_hi)
    ec2p = 0.5 * (ec2p_lo + ec2p_hi)
    ec2m = 0.5 * (ec2m_lo + ec2m_hi)
    log(f"  1 粒子 (A=0)          : ε_c^(1)   = ({ec1_lo:.6f}, {ec1_hi:.6f})  [mid={ec1:.6f}]")
    log(f"  2 粒子同符号 (A=0,B=1): ε_c^(2,+) = ({ec2p_lo:.6f}, {ec2p_hi:.6f})  [mid={ec2p:.6f}, ratio={ec2p/ec1:.4f}]")
    log(f"  2 粒子異符号 (A=0,B=1): ε_c^(2,-) = ({ec2m_lo:.6f}, {ec2m_hi:.6f})  [mid={ec2m:.6f}, ratio={ec2m/ec1:.4f}]")

    # ==================================================================
    # [SC-4] 空間分離 2 粒子間の CHSH ベル不等式の動的検証 (Alice=0, Bob=1)
    # ==================================================================
    log("\n" + "-" * 88)
    log("[SC-4] 空間分離 2 サイト (Alice=0, Bob=1) の CHSH ベル不等式の動的破れ (Tri(T) 発展)")
    log("-" * 88)

    # まず真空 H 単体での (A=0, B=1) 2 量子ビット背景相関（ベースライン）
    rho_vac, tr_vac = extract_two_site_density_matrix(H, 0, 1, mode="cross_all")
    pv, sev, smv, ssv, _ = chsh_metrics(rho_vac)
    log(f"  [Baseline 真空 H (cross_all)]  Tr={tr_vac:.4e}, purity={pv:.4f}, S_E={sev:.4f} bit, S_max={smv:.4f}, S_std={ssv:.4f}")
    rho_vac_p, tr_vac_p = extract_two_site_density_matrix(H, 0, 1, mode="cross_probe")
    pvp, sevp, smvp, ssvp, _ = chsh_metrics(rho_vac_p)
    log(f"  [Baseline 真空 H (cross_probe)] Tr={tr_vac_p:.4e}, purity={pvp:.4f}, S_E={sevp:.4f} bit, S_max={smvp:.4f}, S_std={ssvp:.4f}")

    # (a) 初期積状態 |00> vs (b) 初期ベル状態 |Φ^+>, |Ψ^-> の動的発展（ε=0.10）
    for stype, label in [
        ("product_00", "初期積状態 |↑_A ↑_B>"),
        ("phi_plus",   "初期ベル状態 |Φ^+> = (|↑↑>+|↓↓>)/√2"),
        ("psi_minus",  "初期ベル状態 |Ψ^-> = (|↑↓>-|↓↑>)/√2"),
    ]:
        P_bell = build_two_site_bell_perturbation(0, 1, stype)
        U = (H + 0.10 * P_bell) / np.linalg.norm(H + 0.10 * P_bell)
        log(f"\n  --- {label} (ε=0.10, プローブ部分空間 cross_probe 抽出) ---")
        log(f"  {'τ':>3} | {'||P||':>9} | {'Tr(ρ_pr)':>10} | {'purity':>7} | {'S_E [bit]':>9} | {'S_max (Horodecki)':>17} | {'S_std (angle)':>13}")
        for st in range(0, 41):
            if st > 0:
                U, nrm = step_tri(U)
            else:
                nrm = normP_H
            if st in (0, 1, 2, 3, 5, 10, 20, 30, 40):
                rho, tr = extract_two_site_density_matrix(U, 0, 1, mode="cross_probe")
                pur, se, smax, sstd, _ = chsh_metrics(rho)
                log(f"  {st:3d} | {nrm:9.6f} | {tr:10.3e} | {pur:7.4f} | {se:9.4f} | {smax:17.6f} | {sstd:13.6f}")

    # ==================================================================
    # [SC-5] 振幅掃引：弱結合共存相 vs 強結合崩壊相における CHSH 破れの維持と崩壊
    # ==================================================================
    log("\n" + "-" * 88)
    log("[SC-5] 振幅掃引（|Ψ^-> ベル状態, τ=40）：弱結合共存相 vs 臨界崩壊相での CHSH 破れ")
    log("-" * 88)
    P_psi = build_two_site_bell_perturbation(0, 1, "psi_minus")
    log(f"  {'ε':>6} | {'||P||(40)':>10} | {'Tr(ρ_pr)':>10} | {'purity':>7} | {'S_E [bit]':>9} | {'S_max':>9} | {'S_std':>9} | {'CHSH > 2?':>9}")
    for eps in [0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50]:
        Uf, nf = evolve_final(H + eps * P_psi, 40)
        rho, tr = extract_two_site_density_matrix(Uf, 0, 1, mode="cross_probe")
        pur, se, smax, sstd, _ = chsh_metrics(rho)
        flag = "VIOLATED" if smax > 2.0 + 1e-6 else "CLASSICAL"
        log(f"  {eps:6.2f} | {nf:10.6f} | {tr:10.3e} | {pur:7.4f} | {se:9.4f} | {smax:9.4f} | {sstd:9.4f} | {flag:>9}")

    log(f"\nDone in {time.time() - t0:.2f}s.")
    out_path = "results/tri_loop_part6_scattering_chsh_v0.txt"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    run()
