# -*- coding: utf-8 -*-
"""ChatGPT v1.1（2026-11-08）が v0.7 / exchange_core_v11 について述べた主張の確認。実験の走行ではない。
C1 命題3：ψ に共通の直交 U は総強度 ρ(ψ)=Σ_i|q_i(ψ)|² を各点で、各次数の全系パワー E_m を保存する（ChatGPT 主張 ~1e-15）。
C2 一回作用後の対の平面性は 0.849〜0.940 に落ちる（ChatGPT §11.1）。
C3 coefficients(odd=True, even=False) でも直流 c₀=φ/360 が残る（ChatGPT §7 移行表）。
C4 奇数倍音だけの対では read_pair（全次数 unwrap→mask）が Δ を誤読し、gcd=2 の半周の別解がある（ChatGPT 命題2・§7）。
C5 Gemini v1.0 の g^F=√Σ_m|Y_m|²/√(R_i²R_j²) は中心差に依らず一定、v0.7 の握手で使う |ΣY_m|/R_i² は相殺で減衰する（ChatGPT §6.2 表1）。
C6 Gemini v1.0 掲載の generator_natural（g^F×(π+|Φ|)）を v0.7 の条件に掛けると L の相対低下は 16.5〜27.7%（ChatGPT §11.1）。
"""
import json, math, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exchange_core_v11 as C
NG = 2048; psi = np.arange(NG) * 2 * math.pi / NG; dpsi = 2 * math.pi / NG; M = 100
theta = [97.0, 194.0, 58.0, 300.0]; phi0 = [0.0, 60.0, 120.0, 180.0]
axes = [C.unit(np.array(v, float)) for v in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1])]
Q = [C.make_pulse(psi, theta[i], phi0[i], 1.0, M, 1.0, axes[i]) for i in range(4)]
X = C.relation_matrix(Q, M); Xtot = X.sum(axis=2)
out = {}

def pipeline(Q, X):
    """check_core_v11.py と同じ手順で Φ, û を読み、G, U を作る"""
    Phi = [[0.0] * 4 for _ in range(4)]; Uax = [[None] * 4 for _ in range(4)]; Cax = [[0.0] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            if i == j: continue
            u, wv, c, plan = C.pair_frame(X[i, j]); Uax[i][j] = u; Cax[i][j] = c
    Uax, r_hat = C.orient_axes(Uax)
    for i in range(4):
        for j in range(4):
            if i == j: continue
            Phi[i][j] = C.read_pair(X[i, j], Uax[i][j], Cax[i][j])["Phi"]
    return Phi, Uax

Phi, Uax = pipeline(Q, X)
G = C.generator(4, Phi, Uax); U = C.simultaneous_unitary(G); Q1 = C.apply_unitary(U, Q)
# C1 命題3
rho0 = sum(np.sum(q * q, axis=1) for q in Q); rho1 = sum(np.sum(q * q, axis=1) for q in Q1)
out["C1_rho_pointwise_max_rel_change"] = float(np.max(np.abs(rho1 - rho0)) / np.max(rho0))
Em0 = C.relation_matrix(Q, M)[np.arange(4), np.arange(4), :, 0].sum(axis=0); Em1 = C.relation_matrix(Q1, M)[np.arange(4), np.arange(4), :, 0].sum(axis=0)
out["C1_Em_total_max_rel_change"] = float(np.max(np.abs(Em1 - Em0)) / np.max(Em0))
out["C1_per_wave_Em_changed_max_rel"] = float(max(np.max(np.abs(C.relation_matrix(Q1, M)[i, i, :, 0] - X[i, i, :, 0])) / np.max(X[i, i, :, 0]) for i in range(4)))
# C2 作用後の平面性
X1 = C.relation_matrix(Q1, M)
plan1 = sorted(float(C.pair_frame(X1[i, j])[3]) for i in range(4) for j in range(4) if i != j)
out["C2_planarity_after_min_max"] = [plan1[0], plan1[-1]]
# C3 奇数だけの直流
qodd = C.make_pulse(psi, 97.0, 0.0, 1.0, M, 1.0, axes[2], odd=True, even=False)
a, b = C.harmonic_coeffs(qodd); out["C3_odd_only_DC_norm"] = float(np.linalg.norm(a[0])); out["C3_odd_only_c0_expected"] = 1.0 / 360.0
# C4a 同軸の対（c=1、偶奇両方、Δ=30°、Φ=20°）：pair_frame の最小二乗 ŵ は Vec が n 一本の階数 1 のとき Vec·w=−Sc を m ごとに満たせず、偽の ŵ・c を返す（ChatGPT §5.1「階数 1 から n_i=n_j は結論できない」の裏側）
def read_one_pair(q1, q2, n_ref):
    Y = C.relation_matrix([q1, q2], M); u, wv, c, plan = C.pair_frame(Y[0, 1]); r = C.read_pair(Y[0, 1], u, c)
    sgn = 1.0 if np.dot(u, n_ref) >= 0 else -1.0
    return {"Phi_read_deg": sgn * r["Phi"] / C.DEG, "Delta_read_deg": sgn * r["Delta"] / C.DEG, "resid": r["resid"], "c_read": c, "planarity": plan, "k_w": float(np.linalg.norm(wv))}
out["C4a_coaxial_both_parity"] = dict(read_one_pair(C.make_pulse(psi, 130.0, 20.0, 1.0, M, 1.0, axes[2]), C.make_pulse(psi, 100.0, 0.0, 1.0, M, 1.0, axes[2]), axes[2]), Phi_true_deg=20.0, Delta_true_deg=30.0, c_true=1.0)
# C4b 軸の違う対（i, j、c=0）：偶奇両方は読める、奇数だけは read_pair（全次数 unwrap→mask）が壊れる（ChatGPT 命題2・§7）
out["C4b_skew_both_parity"] = dict(read_one_pair(C.make_pulse(psi, 130.0, 20.0, 1.0, M, 1.0, axes[0]), C.make_pulse(psi, 100.0, 0.0, 1.0, M, 1.0, axes[1]), C.unit(axes[0] + axes[1])), Phi_true_deg=20.0, Delta_true_deg=30.0, c_true=0.0)
out["C4b_skew_odd_only"] = dict(read_one_pair(C.make_pulse(psi, 130.0, 20.0, 1.0, M, 1.0, axes[0], odd=True, even=False), C.make_pulse(psi, 100.0, 0.0, 1.0, M, 1.0, axes[1], odd=True, even=False), C.unit(axes[0] + axes[1])), Phi_true_deg=20.0, Delta_true_deg=30.0, c_true=0.0)
# C4c 奇数だけで、偶数次数を外してから unwrap した場合（ChatGPT の「支持を先に選ぶ」）。gcd=2 なので (Φ,Δ) と (Φ+π·m0, Δ+π) の二候補
Yo = C.relation_matrix([C.make_pulse(psi, 130.0, 20.0, 1.0, M, 1.0, axes[0], odd=True, even=False), C.make_pulse(psi, 100.0, 0.0, 1.0, M, 1.0, axes[1], odd=True, even=False)], M)
u, wv, c, plan = C.pair_frame(Yo[0, 1]); sgn = 1.0 if np.dot(u, C.unit(axes[0] + axes[1])) >= 0 else -1.0
ms = np.arange(1, M + 1); Sc = Yo[0, 1, 1:, 0]; Vu = Yo[0, 1, 1:, 1:] @ u; s = math.sqrt(max((1 + c) / 2, 1e-12))
odd = (ms % 2 == 1); th = np.unwrap(np.arctan2(-s * Vu[odd], Sc[odd]), period=2 * math.pi, discont=math.pi)  # 奇数次数だけで隣接差 2Δ を unwrap
A_ = np.vstack([ms[odd], np.ones(int(odd.sum()))]).T; slope, icpt = np.linalg.lstsq(A_, th, rcond=None)[0]
out["C4c_odd_only_support_first"] = {"Phi_read_deg": sgn * ((icpt + math.pi) % (2 * math.pi) - math.pi) / C.DEG, "Delta_read_deg": sgn * (-slope) / C.DEG, "note": "奇数支持だけで unwrap。2Δ=60° の隣接差は unwrap 可だが、gcd=2 により Δ+180°, Φ+180° も同じ位相列を与える（別解）"}
# C5 g^F と g（同軸同形、幅 1°、M=100）
tab = []
for d in (0.0, 10.0, 90.0, 180.0):
    qa = C.make_pulse(psi, 100.0, 0.0, 1.0, M, 1.0, axes[2]); qc = C.make_pulse(psi, 100.0 + d, 0.0, 1.0, M, 1.0, axes[2])
    Y = C.relation_matrix([qa, qc], M); R2 = Y[0, 0, :, 0].sum() * Y[1, 1, :, 0].sum()
    gF = float(np.sqrt(np.sum(Y[0, 1] ** 2)) / math.sqrt(R2)); g = float(np.linalg.norm(Y[0, 1].sum(axis=0)) / math.sqrt(R2))
    tab.append({"center_diff_deg": d, "g_F_gemini": gF, "g_coherent_v07": g})
out["C5_gF_vs_g"] = tab
# C6 Gemini v1.0 掲載の generator_natural（忠実に写す：g^F × (π+|Φ|)）
def coupling_strengths(Y, R2):
    N = Y.shape[0]; g = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            g[i, j] = np.sqrt(np.sum(Y[i, j] ** 2)) / math.sqrt(max(R2[i] * R2[j], 1e-30))
    return g
def generator_natural(N, Phi, U_axis, g_coupling):
    G = np.zeros((4 * N, 4 * N))
    for i in range(N):
        for j in range(i + 1, N):
            phi_effective = g_coupling[i, j] * (math.pi + abs(Phi[i][j])); Ru = C.Rmat(np.concatenate([[0.0], U_axis[i][j]]))
            PD = np.zeros((N, N)); PD[i, i] = PD[j, j] = 0.5; PD[i, j] = PD[j, i] = -0.5
            G += phi_effective * np.kron(PD, Ru)
    return G
R2 = [X[i, i, :, 0].sum() for i in range(4)]; gF = coupling_strengths(X, R2)
Ug = C.simultaneous_unitary(generator_natural(4, Phi, Uax, gF)); Qg = C.apply_unitary(Ug, Q)
L0 = [C.read_self(X[i, i], Q[i], dpsi)["L"] for i in range(4)]; Lg = [C.read_self(C.relation_matrix(Qg, M)[i, i], Qg[i], dpsi)["L"] for i in range(4)]
out["C6_gemini_gF_offdiag"] = {f"{i+1}{j+1}": float(gF[i, j]) for i in range(4) for j in range(i + 1, 4)}
out["C6_gemini_L_rel_drop_percent"] = [100.0 * (1 - Lg[i] / L0[i]) for i in range(4)]
out["C6_gemini_U_orthogonality_err"] = float(np.max(np.abs(Ug.T @ Ug - np.eye(16))))
out["C6_gemini_phi_effective_rad"] = {f"{i+1}{j+1}": float(gF[i, j] * (math.pi + abs(Phi[i][j]))) for i in range(4) for j in range(i + 1, 4)}
os.makedirs("results", exist_ok=True)
json.dump(out, open(os.path.join("results", "check_chatgpt_claims_v1.json"), "w"), ensure_ascii=False, indent=2)
print(json.dumps(out, ensure_ascii=False, indent=2))
