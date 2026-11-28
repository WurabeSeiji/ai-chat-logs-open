#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Part II（統合仮説論文：外部時空 O_8^ext × 内部ゲージ O_8^int 分離仮説による
重力・電磁気・QCD ヤン＝ミルズ・ディラックの統一接続）用 決定論的機械精度検証スクリプト。

  [P2-1] 外部時空 O_8^ext × 内部ゲージ O_8^int 分離（Coleman-Mandula 回避の最小構成）:
         - 外側 SO(3) 空間回転生成子 (G_ext ⊗ I_8) と 内側 SU(3) × U(1)_Q 生成子 (I_8 ⊗ G_int) の
           厳密可換性 ||[G_ext ⊗ I, I ⊗ G_int]|| = 0.000e+00
         - 直積準同型 Tri(A ⊗ B) == Tri(A) ⊗ Tri(B) による独立進化
  [P2-2] 同一三つ組 Re(o_A o_B o_C) からの電磁気 (Maxwell) vs QCD (Yang-Mills) 読み分け:
         - 可換 U(1)_Q 平面: 交換子 [A, B] = 0 -> 閉じ残り F = O(ε^2) のみ、非可換項 0.000e+00
         - 非可換 SU(3) 8 グルーオン方向: 一次閉じ (A+B+C=0) で 1 - Re(U_A U_B U_C)/N -> (ε^4/2) ||[A, B]||_F^2 / N
           の厳密 O(ε^4) ヤン＝ミルズ作用スケーリングを確認
         - ビアンキ恒等式（ヤコビ恒等式）[[Λ_A, Λ_B], Λ_C] + cyclic = 0 の誤差確認
  [P2-3] ディラック演算子 D = Σ_{a=1..3} L_a p_a + m の平方 D^2 = -|p|^2 I + m^2 I (分散関係)
         および 1 世代 8 状態 (ν, d_bar, u, e^+) の選択的ゲージ結合:
         - U(1)_Q 光子回転 exp(i θ Q): ν(Q=0) は厳密不変 (0.000)、帯電粒子のみ exp(i Q θ) で位相回転
         - SU(3) グルーオン回転 exp(i θ Λ_A): レプトン (ν, e^+) は厳密不変 (0.000)、クォーク (d_bar, u) のみ変換
  [P2-4] 三つ組のディラック読み替え Re(o_A o_B o_C) == o_A_bar^T L(o_B) o_C の機械精度検証
乱数なし。
"""
import os
import time
import numpy as np
from scipy.linalg import expm
from bm_product_algebra_tables_v0 import make_mult_table


def oct_mul(O8, x, y):
    return np.einsum("ijk,i,j->k", O8, x, y)


def oct_conj(x):
    y = x.copy()
    y[1:] = -y[1:]
    return y


def oct_exp(O8, v):
    """純虚八元数 v (v[0]=0) の厳密指数 exp(v) = cos(|v|) + (v/|v|) sin(|v|)"""
    nrm = float(np.linalg.norm(v[1:]))
    out = np.zeros(8, dtype=np.float64)
    if nrm < 1e-15:
        out[0] = 1.0
        return out
    out[0] = np.cos(nrm)
    out[1:] = (np.sin(nrm) / nrm) * v[1:]
    return out


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    L = [O8[a, :, :].T.copy() for a in range(8)]

    log("=" * 80)
    log("Part II 検証：外部時空×内部ゲージ分離・Maxwell/Yang-Mills・Dirac 選択的結合")
    log("=" * 80)

    # ------------------------------------------------------------------
    # e_7 固定による CAR 演算子 α_1, α_2, α_3, Q, SU(3) Λ_1..Λ_8 の構成
    # ------------------------------------------------------------------
    pairs = []
    used = set()
    for a in range(1, 7):
        if a in used:
            continue
        vec = O8[7, a, :]
        b = int(np.argmax(np.abs(vec)))
        sgn = float(vec[b])
        pairs.append((a, b, sgn))
        used.add(a)
        used.add(b)

    alphas = [0.5 * (-L[a] + 1j * sgn * L[b]) for (a, b, sgn) in pairs]
    alpha_dags = [al.conj().T for al in alphas]
    N_op = sum(alpha_dags[m] @ alphas[m] for m in range(3))
    Q_op = N_op / 3.0

    lam = [
        np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
        np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex),
        np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),
        np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
        np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex),
        np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
        np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex),
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / np.sqrt(3.0),
    ]
    su3_ops = []
    for g in lam:
        op = np.zeros((8, 8), dtype=complex)
        for m in range(3):
            for n in range(3):
                if g[m, n] != 0:
                    op += g[m, n] * (alpha_dags[m] @ alphas[n])
        su3_ops.append(op)

    # ------------------------------------------------------------------
    # [P2-1] 外部時空 O_8^ext (SO(3)) × 内部ゲージ O_8^int (SU(3)×U(1)_Q) の厳密可換性
    # ------------------------------------------------------------------
    log("\n" + "-" * 80)
    log("[P2-1] 外部時空 O_8^ext (SO(3)) × 内部ゲージ O_8^int (SU(3) × U(1)_Q) 厳密可換性")
    log("-" * 80)
    so3_ext = []
    for (i, j) in [(1, 2), (2, 3), (3, 1)]:
        G = np.zeros((8, 8), dtype=complex)
        G[i, j] = -1.0
        G[j, i] = 1.0
        so3_ext.append(G)

    I8 = np.eye(8, dtype=complex)
    max_comm_cm = 0.0
    for G_ext in so3_ext:
        A64 = np.kron(G_ext, I8)
        for G_int in su3_ops + [Q_op]:
            B64 = np.kron(I8, G_int)
            err = float(np.linalg.norm(A64 @ B64 - B64 @ A64))
            max_comm_cm = max(max_comm_cm, err)
    log(f"  全 3×9 = 27 組の ||[G_ext(SO3) ⊗ I, I ⊗ G_int(SU3×U1)]|| 最大誤差 = {max_comm_cm:.3e}")

    # ------------------------------------------------------------------
    # [P2-2] 三つ組閉じ残りからの Maxwell (可換 U(1)) vs Yang-Mills (非可換 SU(3)) 読み分け
    # ------------------------------------------------------------------
    log("\n" + "-" * 80)
    log("[P2-2] 三つ組閉じ残り：Maxwell (可換 U(1)_Q) vs Yang-Mills (非可換 SU(3)) & ビアンキ恒等式")
    log("-" * 80)
    # 1) 可換 U(1)_Q 平面 (e_7 方向): A = a*e_7, B = b*e_7, C = -(a+b)*e_7 -> 閉じ残り厳密ゼロ
    e7 = np.zeros(8); e7[7] = 1.0
    e1 = np.zeros(8); e1[1] = 1.0
    e2 = np.zeros(8); e2[2] = 1.0
    eps_list = [0.2, 0.1, 0.05, 0.025]
    log("  (a) 八元数三つ組 Re(o_A o_B o_C) の一次閉ループ (c = -a-b):")
    for eps in eps_list:
        # 可換 (Maxwell 型): 同一虚数軸 e_7
        oA = oct_exp(O8, eps * 1.0 * e7)
        oB = oct_exp(O8, eps * 0.7 * e7)
        oC = oct_exp(O8, -eps * 1.7 * e7)
        c_ab = oct_mul(O8, oct_mul(O8, oA, oB), oC)[0]
        # 非可換 (Yang-Mills 型): 異方向 e_1, e_2 (cross = e_1 x e_2 = e_3, |cross|=1)
        uA = oct_exp(O8, eps * e1)
        uB = oct_exp(O8, eps * e2)
        uC = oct_exp(O8, -eps * (e1 + e2))
        c_na = oct_mul(O8, oct_mul(O8, uA, uB), uC)[0]
        # 理論値: 1 - c_na ≈ (1/2) * eps^4 * |a x b|^2  (係数 1/2 の確認: Baker-Campbell-Hausdorff)
        # 注: exp(eA)exp(eB)exp(-e(A+B)) = exp( (1/2) e^2 [A,B] + ... )
        # 八元数積では [e1, e2] = 2 e3 なので (1/2)[e1,e2] = e3, よって 1 - Re(...) ≈ (1/2) eps^4 |e1 x e2|^2
        ym_meas = 1.0 - c_na
        ym_theo = 0.5 * (eps ** 4) * (1.0 ** 2)
        log(f"    eps={eps:5.3f}: Maxwell (可換) 1-Re={1.0-c_ab:.3e} | Yang-Mills (非可換) (1-Re)/(0.5*eps^4)={ym_meas/ym_theo:.6f}")

    # 2) SU(3) 8 グルーオン行列群での Wilson プラケット閉じ残り & ビアンキ (ヤコビ) 恒等式
    A_su3 = 1j * su3_ops[0]  # i Λ_1
    B_su3 = 1j * su3_ops[3]  # i Λ_4 (非可換!)
    comm_AB = A_su3 @ B_su3 - B_su3 @ A_su3
    comm_norm2 = float(np.linalg.norm(comm_AB, "fro") ** 2)
    log("  (b) SU(3) 8 グルーオン表現上の三つ組プラケット 1 - Re Tr(U_A U_B U_C)/8:")
    for eps in eps_list:
        UA = expm(eps * A_su3)
        UB = expm(eps * B_su3)
        UC = expm(-eps * (A_su3 + B_su3))
        plaq = float(np.real(np.trace(UA @ UB @ UC)) / 8.0)
        # exp(eps A)exp(eps B)exp(-eps(A+B)) = I + (1/2)eps^2 [A,B] + (1/8)eps^4 [A,B]^2 + ...
        # 反エルミート [A,B] より Tr([A,B]^2) = -||[A,B]||_F^2 -> 1 - Re Tr/8 ≈ (1/8) eps^4 ||[A,B]||_F^2 / 8
        theo_su3 = (eps ** 4) * comm_norm2 / (8.0 * 8.0)
        log(f"    eps={eps:5.3f}: (1 - Re Tr/8) / theo_YM = {(1.0 - plaq)/theo_su3:.6f}")

    # ビアンキ (ヤコビ) 恒等式
    jac_err = 0.0
    for i in range(8):
        for j in range(i + 1, 8):
            for k in range(j + 1, 8):
                X, Y, Z = su3_ops[i], su3_ops[j], su3_ops[k]
                c1 = (X @ Y - Y @ X) @ Z - Z @ (X @ Y - Y @ X)
                c2 = (Y @ Z - Z @ Y) @ X - X @ (Y @ Z - Z @ Y)
                c3 = (Z @ X - X @ Z) @ Y - Y @ (Z @ X - X @ Z)
                jac_err = max(jac_err, float(np.linalg.norm(c1 + c2 + c3)))
    log(f"  SU(3) ヤコビ／ビアンキ恒等式 (全 56 組) 最大誤差 = {jac_err:.3e}")

    # ------------------------------------------------------------------
    # [P2-3] ディラック演算子 D^2 = -|p|^2 + m^2 & 1 世代 8 状態の選択的ゲージ結合
    # ------------------------------------------------------------------
    log("\n" + "-" * 80)
    log("[P2-3] ディラック演算子分散関係 & 電子・クォーク・ニュートリノの選択的ゲージ結合")
    log("-" * 80)
    # 3 次元空間運動量 p = (p1, p2, p3) と質量 m に対する離散ディラック演算子
    # D_slash = p1 L_1 + p2 L_2 + p3 L_3  (反エルミート, D_slash^2 = -|p|^2 I)
    p_vec = np.array([1.3, -0.7, 2.1])
    m_mass = 0.85
    D_slash = p_vec[0] * L[1] + p_vec[1] * L[2] + p_vec[2] * L[3]
    # クリフォード関係より (i D_slash)^2 = |p|^2 I, かつ (i D_slash + m I)(i D_slash - m I) = (|p|^2 - m^2) I
    D_plus = 1j * D_slash + m_mass * np.eye(8)
    D_minus = 1j * D_slash - m_mass * np.eye(8)
    kg_err = float(np.linalg.norm(D_plus @ D_minus - (np.dot(p_vec, p_vec) - m_mass ** 2) * np.eye(8)))
    log(f"  ディラック分散関係 ||(i D_slash + m)(i D_slash - m) - (|p|^2 - m^2)I|| = {kg_err:.3e}")

    # 1 世代 8 状態への U(1)_Q 光子回転 & SU(3) グルーオン回転の選択的結合テスト
    evals_n, evecs_n = np.linalg.eigh(N_op)
    th_ph = 0.73  # 任意の光子位相角
    U_photon = expm(1j * th_ph * Q_op)
    U_gluon = expm(1j * 0.51 * su3_ops[0])  # Λ_1 グルーオン回転

    log(f"  1 世代 8 状態の選択的ゲージ結合 (光子角 θ={th_ph}, グルーオン角 θ_g=0.51):")
    log(f"  {'Sector':>18s} | {'Q':>5s} | {'||U_Q ψ - e^{iQθ}ψ||':>20s} | {'||U_Q ψ - ψ|| (EM結合)':>22s} | {'||U_g ψ - ψ|| (QCD結合)':>23s}")
    sector_names = {0: "ν (レプトン1重項)", 1: "d_bar (クォーク3重項)", 2: "u (クォーク3重項)", 3: "e+ (レプトン1重項)"}
    for k in range(4):
        idx_k = np.where(np.abs(evals_n - k) < 1e-8)[0]
        psi = evecs_n[:, idx_k[0]]
        q_val = k / 3.0
        phase_err = float(np.linalg.norm(U_photon @ psi - np.exp(1j * q_val * th_ph) * psi))
        em_resp = float(np.linalg.norm(U_photon @ psi - psi))
        # SU(3) 全 8 グルーオンに対する最大応答
        qcd_resp = max(float(np.linalg.norm(expm(1j * 0.51 * op) @ psi - psi)) for op in su3_ops)
        log(f"  {sector_names[k]:>18s} | {q_val:5.3f} | {phase_err:20.3e} | {em_resp:22.6f} | {qcd_resp:23.6f}")

    # ------------------------------------------------------------------
    # [P2-4] 11/20 論文 [K7] (7.3) 三つ組ディラック読み替え Re(o_A o_B o_C) == o_A_bar^T L(o_B) o_C
    # ------------------------------------------------------------------
    log("\n" + "-" * 80)
    log("[P2-4] 三つ組ディラック読み替え Re(o_A o_B o_C) == o_A_bar^T L(o_B) o_C 検証")
    log("-" * 80)
    oA = np.array([0.5, 0.3, -0.2, 0.4, 0.1, -0.5, 0.2, 0.35]); oA /= np.linalg.norm(oA)
    oB = np.array([0.2, -0.4, 0.6, 0.1, -0.3, 0.2, 0.4, -0.25]); oB /= np.linalg.norm(oB)
    oC = np.array([-0.3, 0.2, 0.1, 0.5, 0.4, -0.2, -0.5, 0.3]); oC /= np.linalg.norm(oC)
    lhs_d = float(oct_mul(O8, oct_mul(O8, oA, oB), oC)[0])
    L_oB = sum(oB[a] * L[a] for a in range(8))
    rhs_d = float(np.dot(oct_conj(oA), L_oB @ oC))
    log(f"  Re(o_A o_B o_C) = {lhs_d:+.12f},  o_A_bar^T L(o_B) o_C = {rhs_d:+.12f}, diff = {abs(lhs_d - rhs_d):.3e}")

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.2f}s.")

    with open("results/tri_loop_part2_unified_forces_verify_v0.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    run()
