# -*- coding: utf-8 -*-
"""中核関数の確認（N = 4、16 の関係）。実験の走行ではない。
K1 関係行列のエルミート性 X_ji = conj(X_ij)。K2 12 本の順序つき非対角関係の読出し（Φ・Δ・û・ŵ・c）が既知の初期値と一致し、Φ_ji = −Φ_ij、û_ji = û_ij。
K3 U = exp(A) の直交性、ブロックが右乗であること、全ノルム保存。K4 観測者ごとの位相の回転の実測と、16 の関係からの予測の一致。K5 自己関係の読み。"""
import json, math, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exchange_core_v11 as C
NG = 2048; psi = np.arange(NG) * 2 * math.pi / NG; dpsi = 2 * math.pi / NG; M = 100
theta = [97.0, 194.0, 58.0, 300.0]; phi0 = [0.0, 60.0, 120.0, 180.0]
axes = [np.array(v, float) for v in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1])]; axes = [C.unit(a) for a in axes]
Q = [C.make_pulse(psi, theta[i], phi0[i], 1.0, M, 1.0, axes[i]) for i in range(4)]
X = C.relation_matrix(Q, M); Xtot = X.sum(axis=2)   # 左共役の関係 Y（変数名は X のまま）
out = {}
# K1
herm = max(float(np.max(np.abs(X[i, j] - C.qconj(X[j, i])))) for i in range(4) for j in range(4)); out["K1_hermitian_max_err"] = herm
# K5 自己関係
out["K5_self"] = {f"{i+1}": C.read_self(X[i, i], Q[i], dpsi) for i in range(4)}
# K2 非対角（順序つき 12 本）
wrap = lambda d: (d + 180) % 360 - 180
K2 = []; Phi = [[0.0] * 4 for _ in range(4)]; Uax = [[None] * 4 for _ in range(4)]; Wax = [[None] * 4 for _ in range(4)]; Cax = [[0.0] * 4 for _ in range(4)]
for i in range(4):
    for j in range(4):
        if i == j: continue
        u, wv, c, plan = C.pair_frame(X[i, j]); Uax[i][j] = u; Wax[i][j] = wv; Cax[i][j] = (c, plan)
Uax, r_hat = C.orient_axes(Uax); out["orientation_reference_r_hat"] = r_hat.tolist()
for i in range(4):
    for j in range(4):
        if i == j: continue
        u = Uax[i][j]; wv = Wax[i][j]; c, plan = Cax[i][j]; r = C.read_pair(X[i, j], u, c)
        Phi[i][j] = r["Phi"]
        u_true = C.unit(axes[i] + axes[j]); u_true = u_true if np.dot(u_true, r_hat) >= 0 else -u_true
        w_true = C.unit(np.cross(axes[i], axes[j])) if np.linalg.norm(np.cross(axes[i], axes[j])) > 1e-12 else np.zeros(3)
        sgn = 1.0 if np.dot(C.unit(axes[i] + axes[j]), r_hat) >= 0 else -1.0   # û を反転した対は Φ, Δ の符号も反転して読まれる
        K2.append({"ij": f"{i+1}{j+1}", "Phi_read": r["Phi"] / C.DEG, "Phi_true": wrap(sgn * (phi0[i] - phi0[j])), "Delta_read": r["Delta"] / C.DEG, "Delta_true": wrap(sgn * (theta[i] - theta[j])),
                   "resid": r["resid"], "c_read": c, "c_true": float(axes[i] @ axes[j]), "u_err": float(np.linalg.norm(u - u_true)), "w_err": float(np.linalg.norm(wv - w_true)) if np.linalg.norm(w_true) > 0 else None, "planarity": plan})
out["K2_offdiag"] = K2
out["K2_max_err_Phi_deg"] = max(abs(wrap(e["Phi_read"] - e["Phi_true"])) for e in K2); out["K2_max_err_Delta_deg"] = max(abs(wrap(e["Delta_read"] - e["Delta_true"])) for e in K2)
out["K2_max_err_c"] = max(abs(e["c_read"] - e["c_true"]) for e in K2); out["K2_max_err_u"] = max(e["u_err"] for e in K2); out["K2_max_err_w"] = max(e["w_err"] for e in K2 if e["w_err"] is not None)
out["K2_order_consistency"] = {"max|Phi_ji + Phi_ij| deg": max(abs(wrap((Phi[i][j] + Phi[j][i]) / C.DEG)) for i in range(4) for j in range(4) if i != j),   # û を揃えたので Φ_ji = −Φ_ij
                               "max|u_ji - u_ij|": max(float(np.linalg.norm(Uax[i][j] - Uax[j][i])) for i in range(4) for j in range(4) if i != j)}
# K3 作用
G = C.generator(4, Phi, Uax); U = C.simultaneous_unitary(G)
out["K3_antisymmetric_G"] = float(np.max(np.abs(G + G.T))); out["K3_orthogonal_U"] = float(np.max(np.abs(U.T @ U - np.eye(16))))
H, berr = C.blocks_as_quaternions(U, 4); out["K3_blocks_are_right_mult_err"] = berr
Q1 = C.apply_unitary(U, Q); n0 = sum(float(np.sum(q * q)) * dpsi for q in Q); n1 = sum(float(np.sum(q * q)) * dpsi for q in Q1)
out["K3_norm_rel_drift"] = abs(n1 - n0) / n0
# K4 握手：実測 δ_i と 16 の関係からの予測
K4 = []
for i in range(4):
    am, axm = C.rotation_measured(Q[i], Q1[i], dpsi); ap, axp = C.handshake_prediction(Xtot, H, i)
    K4.append({"i": i + 1, "delta_measured_deg": am / C.DEG, "delta_predicted_deg": ap / C.DEG, "axis_err": float(np.linalg.norm(axm - axp)),
               "view_angle_terms |X_ij|/R_i^2": {f"{j+1}": float(np.linalg.norm(Xtot[i, j]) / Xtot[i, i, 0]) for j in range(4)}})
out["K4_handshake"] = K4; out["K4_max_err_deg"] = max(abs(e["delta_measured_deg"] - e["delta_predicted_deg"]) for e in K4)
out["gram_eigs_before"] = C.gram_eigs(Xtot).tolist(); out["gram_eigs_after"] = C.gram_eigs(C.relation_matrix(Q1, M).sum(axis=2)).tolist()
out["self_after"] = {f"{i+1}": C.read_self(C.relation_matrix(Q1, M)[i, i], Q1[i], dpsi) for i in range(4)}
json.dump(out, open(os.path.join("results", "check_core_v11.json"), "w"), ensure_ascii=False, indent=2)
print(json.dumps({k: v for k, v in out.items() if k not in ("K2_offdiag", "K4_handshake", "K5_self", "self_after")}, ensure_ascii=False, indent=2))
for e in K2: print(f"  {e['ij']}: Φ {e['Phi_read']:8.3f} (真 {e['Phi_true']:6.1f})  Δ {e['Delta_read']:8.3f} (真 {e['Delta_true']:6.1f})  c {e['c_read']:+.4f} ({e['c_true']:+.4f})  |û−真| {e['u_err']:.1e}  |ŵ−真| {'—' if e['w_err'] is None else format(e['w_err'], '.1e')}  resid {e['resid']:.1e}  平面性 {e['planarity']:.6f}")
for e in K4: print(f"  観測者 {e['i']}: δ 実測 {e['delta_measured_deg']:.4f}° 予測 {e['delta_predicted_deg']:.4f}°  軸差 {e['axis_err']:.1e}  視野角項 {e['view_angle_terms |X_ij|/R_i^2']}")
print("  自己関係（前）:", {k: {kk: round(vv, 5) for kk, vv in v.items()} for k, v in out["K5_self"].items()})
print("  自己関係（後）:", {k: {kk: round(vv, 5) for kk, vv in v.items()} for k, v in out["self_after"].items()})
