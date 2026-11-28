#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第 7 速報候補：八元数 Triality（スピン 1/2）、CAR フェルミオン代数（1:3:3:1 電荷スペクトル）、
               およびパリティ破れ（カイラル非対称性）の厳密検証スクリプト
  Part 1: Spin(8) Triality (28 次元) の厳密不変性 O8(Av x, y, z) + O8(x, As y, z) + O8(x, y, Ac z) = 0
          および空間 2π 回転 exp(2π A_v) = +I に対するスピノル符号反転 exp(2π A_s) = exp(2π A_c) = -I (スピン 1/2)
  Part 2: e_7 固定による O = C ⊕ C^3 分解、SU(3) (8 グルーオン) × U(1)_Q 分解、
          左乗算クリフォード代数 Cl_6 からの 3 つのフェルミオン生成消滅演算子 α_1, α_2, α_3 の構成、
          正準反交換関係 {α_a, α_b†} = δ_ab, {α_a, α_b} = 0 (パウリ排他則 α_a^2 = 0)、
          および電荷演算子 Q = N/3 の厳密スペクトル {0 (×1), 1/3 (×3), 2/3 (×3), 1 (×1)} = 1:3:3:1
  Part 3: H = O_8 ⊗ M_8 上でのカイラル非対称性（パリティ破れ）の直接力学検証：
          同回転 (J_7 = +J_in) vs 逆回転 (J_7 = -J_in)、および (I + σ_y) vs (I - σ_y) の τ=30 応答比較
乱数なし。
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


def run():
    t0 = time.time()
    os.makedirs("results", exist_ok=True)
    lines = []

    def log(s=""):
        print(s, flush=True)
        lines.append(s)

    O8 = np.array(make_mult_table(8), dtype=np.float64)
    T0 = tensor_prod_3d(O8, O8)

    log("=" * 80)
    log("八元数 Triality（スピン 1/2）・CAR フェルミオン電荷 (0, 1/3, 2/3, 1)・カイラル検証")
    log("=" * 80)

    # ------------------------------------------------------------------
    # Part 1: Spin(8) Triality (28 次元) と 2π 回転でのスピン 1/2 符号反転 (-I)
    # ------------------------------------------------------------------
    log("\n" + "-" * 80)
    log("[Part 1] Spin(8) Triality (28 次元) の厳密解と 2π 回転でのスピン 1/2 符号反転")
    log("-" * 80)

    # so(8) の 28 個の反対称基底 E_{ab} - E_{ba} (a < b)
    so8_basis = []
    so8_labels = []
    for a in range(8):
        for b in range(a + 1, 8):
            E = np.zeros((8, 8), dtype=np.float64)
            E[a, b] = 1.0
            E[b, a] = -1.0
            so8_basis.append(E)
            so8_labels.append((a, b))

    # 84 変数 (A_v, A_s, A_c) ∈ so(8)^3 に対する O8 不変条件の線形行列 (512 × 84)
    # O8[i,j,k] について sum_p (Av[i,p] O8[p,j,k] + As[j,p] O8[i,p,k] + Ac[k,p] O8[i,j,p]) = 0
    M_tri = np.zeros((512, 84), dtype=np.float64)
    for idx, E in enumerate(so8_basis):
        # 1st index (v)
        d1 = np.einsum("ip,pjk->ijk", E, O8).ravel()
        # 2nd index (s)
        d2 = np.einsum("jp,ipk->ijk", E, O8).ravel()
        # 3rd index (c)
        d3 = np.einsum("kp,ijp->ijk", E, O8).ravel()
        M_tri[:, idx] = d1
        M_tri[:, 28 + idx] = d2
        M_tri[:, 56 + idx] = d3

    _, s_tri, vt_tri = np.linalg.svd(M_tri, full_matrices=True)
    null_dim = int(np.sum(s_tri < 1e-10))
    # 84 列に対し rank = 56 なら null_dim = 28 (厳密に Spin(8)!)
    rank_tri = int(np.sum(s_tri > 1e-10))
    null_space = vt_tri[rank_tri:, :]  # (28, 84)
    log(f"  Triality 線形系 (512×84): rank = {rank_tri}, 厳密不変カーネル次元 = {84 - rank_tri} (理論値 Spin(8) = 28)")
    log(f"  最小非ゼロ特異値 = {s_tri[rank_tri-1]:.6f}, 最大ゼロ特異値 = {s_tri[rank_tri] if rank_tri < len(s_tri) else 0.0:.3e}")

    # 任意の A_v (28次元) に対して (A_s, A_c) を一意決定する写像を構成
    # null_space[:, :28] は 28×28 正則行列のはず
    V_part = null_space[:, :28]
    S_part = null_space[:, 28:56]
    C_part = null_space[:, 56:84]
    # 任意の係数 c_v (28,) に対し、w = V_part^{-1} c_v として c_s = S_part.T @ ... ではなく
    # 行ベクトル v_row = alpha @ V_part -> alpha = c_v @ inv(V_part)
    V_inv = np.linalg.inv(V_part)

    def get_triality_triple(c_v):
        alpha = c_v @ V_inv
        c_s = alpha @ S_part
        c_c = alpha @ C_part
        Av = sum(c_v[i] * so8_basis[i] for i in range(28))
        As = sum(c_s[i] * so8_basis[i] for i in range(28))
        Ac = sum(c_c[i] * so8_basis[i] for i in range(28))
        return Av, As, Ac

    # 空間回転 (e_1, e_2) 平面 (so8_labels で (1,2) のインデックス) で 2π 回転 & 4π 回転をテスト！
    idx_12 = so8_labels.index((1, 2))
    cv_12 = np.zeros(28)
    cv_12[idx_12] = 1.0
    Av, As, Ac = get_triality_triple(cv_12)
    inv_res = np.linalg.norm(
        np.einsum("ip,pjk->ijk", Av, O8)
        + np.einsum("jp,ipk->ijk", As, O8)
        + np.einsum("kp,ijp->ijk", Ac, O8)
    )
    Rv_2pi = expm(2 * np.pi * Av)
    Rs_2pi = expm(2 * np.pi * As)
    Rc_2pi = expm(2 * np.pi * Ac)
    Rs_4pi = expm(4 * np.pi * As)
    Rc_4pi = expm(4 * np.pi * Ac)

    log(f"  空間回転 (e1, e2) の Triality 微分不変残差 = {inv_res:.3e}")
    log(f"  2π 回転: ||R_v(2π) - (+I)|| = {np.linalg.norm(Rv_2pi - np.eye(8)):.3e}  (ベクトル: 2π で +1)")
    log(f"  2π 回転: ||R_s(2π) - (-I)|| = {np.linalg.norm(Rs_2pi + np.eye(8)):.3e}  (左手スピノル: 2π で -1！)")
    log(f"  2π 回転: ||R_c(2π) - (-I)|| = {np.linalg.norm(Rc_2pi + np.eye(8)):.3e}  (右手スピノル: 2π で -1！)")
    log(f"  4π 回転: ||R_s(4π) - (+I)|| = {np.linalg.norm(Rs_4pi - np.eye(8)):.3e}, ||R_c(4π) - (+I)|| = {np.linalg.norm(Rc_4pi - np.eye(8)):.3e}")

    # 有限回転 (θ = 0.7 rad) で O8(Rv x, Rs y, Rc z) == O8(x,y,z) を直接確認
    th_t = 0.7
    Rv_t, Rs_t, Rc_t = expm(th_t * Av), expm(th_t * As), expm(th_t * Ac)
    O8_rot = np.einsum("ip,jq,kr,pqr->ijk", Rv_t, Rs_t, Rc_t, O8, optimize=True)
    log(f"  有限 Triality 回転 (θ=0.7) 不変誤差 ||(Rv⊗Rs⊗Rc)O8 - O8|| / ||O8|| = {np.linalg.norm(O8_rot - O8)/np.linalg.norm(O8):.3e}")

    # ------------------------------------------------------------------
    # Part 2: e_7 固定 -> クリフォード代数 Cl_6 -> 3 つの CAR フェルミオン演算子 α_1, α_2, α_3
    #         -> 電荷 Q = N/3 の厳密スペクトル {0(×1), 1/3(×3), 2/3(×3), 1(×1)} = 1:3:3:1
    # ------------------------------------------------------------------
    log("\n" + "-" * 80)
    log("[Part 2] e_7 固定による Cl_6 フェルミオン代数 (CAR) と電荷 Q スペクトル (1:3:3:1)")
    log("-" * 80)

    # 左乗算演算子 L_a: (L_a v)_k = (e_a * v)_k = sum_j O8[a, j, k] v_j -> L_a[k, j] = O8[a, j, k]
    L = [O8[a, :, :].T.copy() for a in range(8)]
    # 確認: L_a^2 = -I (a=1..7), {L_a, L_b} = -2 δ_ab I
    cl_err = 0.0
    for a in range(1, 8):
        for b in range(1, 8):
            anti = L[a] @ L[b] + L[b] @ L[a]
            tgt = -2.0 * np.eye(8) if a == b else np.zeros((8, 8))
            cl_err = max(cl_err, float(np.linalg.norm(anti - tgt)))
    log(f"  左乗算クリフォード代数 {{L_a, L_b}} = -2 δ_ab I (a,b=1..7) 最大誤差 = {cl_err:.3e}")

    # e_7 固定による複素構造 J_7 = L_7 (J_7^2 = -I)
    J7 = L[7]
    # L_1..L_6 のうち、e_7 との積でペアになる 3 組を自動抽出 (e_a * e_7 = ±e_b)
    pairs = []
    used = set()
    for a in range(1, 7):
        if a in used:
            continue
        # e_7 * e_a の結果を探す
        vec = O8[7, a, :]
        b = int(np.argmax(np.abs(vec)))
        sgn = float(vec[b])
        pairs.append((a, b, sgn))
        used.add(a)
        used.add(b)
    log(f"  e_7 が定義する 3 つのカラー複素平面 (e_a, e_b = sgn * e_7*e_a): {pairs}")

    # 複素化空間 C^8 上で 3 つのフェルミオン消滅・生成演算子を構成 (Furey 2016 構成):
    # α_m = (1/2) * (-L_a + i * sgn * L_b)  (m = 1, 2, 3)
    # 確認: L_a^† = -L_a より α_m^† = (1/2) * (L_a + i * sgn * L_b)
    alphas = []
    alpha_dags = []
    for (a, b, sgn) in pairs:
        al = 0.5 * (-L[a] + 1j * sgn * L[b])
        al_dag = al.conj().T
        alphas.append(al)
        alpha_dags.append(al_dag)

    # CAR 検証: {α_m, α_n} = 0 (パウリ排他則 α_m^2 = 0 を含む), {α_m, α_n†} = δ_mn I
    car_aa_err = 0.0
    car_ad_err = 0.0
    for m in range(3):
        for n in range(3):
            e1 = float(np.linalg.norm(alphas[m] @ alphas[n] + alphas[n] @ alphas[m]))
            e2 = float(np.linalg.norm(alphas[m] @ alpha_dags[n] + alpha_dags[n] @ alphas[m] - (np.eye(8) if m == n else 0.0)))
            car_aa_err = max(car_aa_err, e1)
            car_ad_err = max(car_ad_err, e2)
    log(f"  フェルミオン反交換関係 {{α_m, α_n}} = 0 (パウリ排他則 α_m^2=0) 最大誤差 = {car_aa_err:.3e}")
    log(f"  フェルミオン反交換関係 {{α_m, α_n†}} = δ_mn I 最大誤差               = {car_ad_err:.3e}")

    # 数演算子 N = sum_m α_m† α_m と 電荷演算子 Q = N / 3 の固有値スペクトル
    N_op = sum(alpha_dags[m] @ alphas[m] for m in range(3))
    Q_op = N_op / 3.0
    evals_N = np.sort(np.linalg.eigvalsh(N_op))
    evals_Q = evals_N / 3.0
    log(f"  数演算子 N = Σ α_m†α_m の全 8 固有値 : {np.round(evals_N, 6).tolist()}")
    log(f"  電荷演算子 Q = N/3 の全 8 固有値     : {np.round(evals_Q, 6).tolist()}")

    # 縮退度カウント (1 : 3 : 3 : 1)
    deg = [int(np.sum(np.abs(evals_N - k) < 1e-8)) for k in range(4)]
    log(f"  電荷セクター縮退度 (Q = 0, 1/3, 2/3, 1): {deg[0]} : {deg[1]} : {deg[2]} : {deg[3]}  (厳密 1:3:3:1 = パスカル二項係数！)")

    # SU(3) 生成子 8 本の構成 (Λ_A = α† λ_A α) と Q との可換性 [SU(3), Q] = 0 の確認
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
    comm_q_err = 0.0
    for g in lam:
        op = np.zeros((8, 8), dtype=complex)
        for m in range(3):
            for n in range(3):
                if g[m, n] != 0:
                    op += g[m, n] * (alpha_dags[m] @ alphas[n])
        su3_ops.append(op)
        comm_q_err = max(comm_q_err, float(np.linalg.norm(op @ Q_op - Q_op @ op)))
    log(f"  SU(3) 8 グルーオン生成子と電荷 Q の可換性 ||[Λ_A, Q]|| 最大誤差 = {comm_q_err:.3e}")

    # 各電荷セクターにおける SU(3) カシミール固有値 (1重項=0, 3重項/反3重項=16/3)
    C2_su3 = sum(op @ op for op in su3_ops)
    # N_op の固有ベクトルで C2_su3 をブロック評価
    evals_n, evecs_n = np.linalg.eigh(N_op)
    for k in range(4):
        idx_k = np.where(np.abs(evals_n - k) < 1e-8)[0]
        V_k = evecs_n[:, idx_k]
        c2_block = V_k.conj().T @ C2_su3 @ V_k
        c2_eigs = np.round(np.linalg.eigvalsh(c2_block), 6).tolist()
        rep_name = "SU(3) 1重項 (レプトン型)" if abs(c2_eigs[0]) < 1e-6 else "SU(3) 3重項/反3重項 (クォーク型)"
        log(f"    N={k} (Q={k}/3, dim={len(idx_k)}): SU(3) Casimir C2 = {c2_eigs} -> {rep_name}")

    # ------------------------------------------------------------------
    # Part 3: H = O_8 ⊗ M_8 上でのカイラル非対称性（パリティ破れ）直接力学検証
    # ------------------------------------------------------------------
    log("\n" + "-" * 80)
    log("[Part 3] 第 1 層 H = O_8 ⊗ M_8 上でのカイラル非対称性（パリティ反転・左右巻き応答差）")
    log("-" * 80)
    U = T0.copy()
    U[0, 1, 0] += 1e-4
    U /= np.linalg.norm(U)
    for _ in range(200):
        U, _ = step_tri(U)
    H = U.copy()
    H_mat = H.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5).reshape(512, 512)
    _, _, vt_svd = np.linalg.svd(H_mat, full_matrices=False)
    M8 = vt_svd[0].reshape(8, 8, 8)

    # 1) H 自身の i <-> k 転置（空間パリティ反転 P）に対する対称／反対称成分比
    H_T = H.transpose(2, 1, 0)
    H_sym = 0.5 * (H + H_T)
    H_asym = 0.5 * (H - H_T)
    log(f"  背景 H のパリティ分解 (i<->k): E_sym = {np.sum(H_sym**2):.6f}, E_asym = {np.sum(H_asym**2):.6f}")

    # 2) 内側 b={2,3} チャネルで、実表現の左右カイラルプローブ (I + J) vs (I - J) [J = -i σ_y]
    #    および外側 O_8 のパリティ反転 O_8^T を掛けたプローブの τ=30 応答比較
    def make_chiral_probe(O_mat, inner_sign):
        # inner_sign = +1: M8 (右回り), -1: M8^T (逆回り = 複素共役 / 反対カイラリティ)
        P = np.zeros((64, 64, 64), dtype=np.float64)
        for A in range(8):
            for B in range(8):
                for C in range(8):
                    val = O_mat[A, B, C]
                    if val == 0.0:
                        continue
                    for a in [4, 5]:
                        for b in [2, 3]:
                            for c in [4, 5]:
                                mval = M8[a, b - 2, c] if inner_sign > 0 else M8[c, b - 2, a]
                                P[A * 8 + a, B * 8 + b, C * 8 + c] += val * mval
        return P / np.linalg.norm(P)

    O8_T = O8.transpose(2, 1, 0)  # パリティ反転した八元数
    chiral_tests = [
        ("O8   ⊗ M8   (同カイラル R,R)", O8, +1),
        ("O8   ⊗ M8^T (逆カイラル R,L)", O8, -1),
        ("O8^T ⊗ M8   (空間反転  L,R)", O8_T, +1),
        ("O8^T ⊗ M8^T (CP反転    L,L)", O8_T, -1),
    ]
    for label, Om, isgn in chiral_tests:
        P_ch = make_chiral_probe(Om, isgn)
        Ut = (H + 0.10 * P_ch)
        Ut /= np.linalg.norm(Ut)
        nP = 0.0
        for _ in range(30):
            Ut, nP = step_tri(Ut)
        T_pr = Ut.reshape(8, 8, 8, 8, 8, 8).transpose(0, 2, 4, 1, 3, 5)[:, :, :, :, 2:4, :]
        e_pr = float(np.sum(T_pr ** 2))
        nz = int(np.sum(np.abs(Ut) > 1e-12))
        log(f"  {label:30s}: tau=30 -> nz={nz:4d}, E_probe={e_pr:.6f}, ||P||={nP:.6f}")

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.2f}s.")

    with open("results/tri_loop_n64_triality_charge_v0.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    run()
