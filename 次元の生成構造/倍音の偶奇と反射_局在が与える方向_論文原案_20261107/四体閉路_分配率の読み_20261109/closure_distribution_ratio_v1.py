# -*- coding: utf-8 -*-
"""四体の閉路から束と分配率を読む（読みだけ。動力学は走らせない）。
要素 q_{i,m}(ψ)：v0.7 の 4 パルス（exchange_core_v11.make_pulse を import、再実装しない）を倍音ごとの関数に分けた 4×(M+1) 個。
 (1) 二要素の関係 Y_ab=⟨q̄_a q_b⟩ は巻数が等しい対にしか残らない（巻数が違う対は零）。
 (2) 束の閉路 Λ=⟨q̄_a q_b q̄_c q_d⟩、a∈層m、b=d∈層m+1、c∈層m+2。位相零で実正なら a,b,c は一つの束。要素の名前（どのパルスか）を使わず、閉路の位相だけで束を再構成する。
 (3) 距離の閉路 ⟨q̄_{a} q_{b} q̄_{b'} q_{a'}⟩、a,b∈層m、a',b'∈層m+1（a,a' 同束、b,b' 同束）。位相が −Δ（共通軸なら厳密）。
 (4) 束の分配率 Ω_ij=Σ_m Sc Y_{(i,m),(j,m)}/R_i²（符号つき）と、|Y| で取った場合の比較。
配置は v0.7 の K 確認と同じ（θ=97,194,58,300、φ₀=0,60,120,180、軸 x,y,z,(111)/√3、幅 1°、M=100、NG=2048）。対照として全軸 k（複素の場合）も走らせる。"""
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import exchange_core_v11 as C
NG = 2048; psi = np.arange(NG) * 2 * math.pi / NG; M = 100
AX_Q = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]; AX_K = [[0, 0, 1]] * 4
TH_INT = [97.0, 194.0, 58.0, 300.0]; PH_INT = [0.0, 60.0, 120.0, 180.0]                    # v0.7 と同じ整数度の配置
TH_GEN = [97.3, 194.71, 58.137, 300.5]; PH_GEN = [0.0, 60.7, 120.33, 180.9]               # 一般の角（整数度の偶然の閉路を避ける）
CONFIGS = {"quaternion_axes_xyz111_integer_deg": (AX_Q, TH_INT, PH_INT), "common_axis_k_integer_deg": (AX_K, TH_INT, PH_INT),
           "quaternion_axes_xyz111_generic_deg": (AX_Q, TH_GEN, PH_GEN), "common_axis_k_generic_deg": (AX_K, TH_GEN, PH_GEN)}
wrap = lambda d: (d + 180.0) % 360.0 - 180.0
def mean_q(x): return x.mean(axis=-2)
def ang_axis(q):
    v = q[1:]; nv = float(np.linalg.norm(v)); return math.degrees(math.atan2(nv, float(q[0]))), (v / nv if nv > 1e-300 else np.zeros(3))
def run(axes_in, theta, phi0):
    axes = [C.unit(np.array(a, float)) for a in axes_in]
    Q = [C.make_pulse(psi, theta[i], phi0[i], 1.0, M, 1.0, axes[i]) for i in range(4)]
    # 要素：倍音ごとの関数 q_{i,m}(ψ) = a_m cos mψ + b_m sin mψ
    E = np.zeros((4, M + 1, NG, 4))
    for i in range(4):
        a, b = C.harmonic_coeffs(Q[i])
        for m in range(M + 1): E[i, m] = a[m][None, :] * np.cos(m * psi)[:, None] + b[m][None, :] * np.sin(m * psi)[:, None]
    recon = max(float(np.max(np.abs(E[i].sum(axis=0) - Q[i]))) for i in range(4))
    out = {"reconstruction_max_err": recon}
    # (1) 二要素の関係：巻数が違う対は零
    same, diff = 0.0, 0.0; Yd = np.zeros((M + 1, 4, 4, 4))   # 同じ巻数のブロックだけ保存
    for m in range(M + 1):
        for mp in range(M + 1):
            blk = mean_q(C.qmul(C.qconj(E[:, m][:, None]), E[:, mp][None, :]))   # (4,4,4)
            if m == mp: Yd[m] = blk; same = max(same, float(np.max(np.linalg.norm(blk, axis=-1))))
            else: diff = max(diff, float(np.max(np.linalg.norm(blk, axis=-1))))
    out["Y_max_norm_same_winding"] = same; out["Y_max_norm_different_winding"] = diff
    # (2) 束の閉路：a∈層m, b=d∈層m+1, c∈層m+2 の全 4^3 組。名前を使わず位相零で結ぶ
    tol = 1e-8; links = []; phases_same = []; phases_mixed = []; n_coincident_mixed = 0
    for m in range(M - 1):
        A_ = E[:, m]; B_ = E[:, m + 1]; Cc = E[:, m + 2]
        for i in range(4):
            for j in range(4):
                ab = C.qmul(C.qconj(A_[i]), B_[j])                      # (NG,4)
                for k in range(4):
                    lam = mean_q(C.qmul(C.qmul(ab, C.qconj(Cc[k])), B_[j]))
                    a_deg, _ = ang_axis(lam); nrm = float(np.linalg.norm(lam))
                    if nrm < 1e-14: continue
                    rec = (m, i, j, k, a_deg, lam[0] > 0)
                    (phases_same if i == j == k else phases_mixed).append(a_deg)
                    if a_deg < tol and lam[0] > 0:
                        links.append(((i, m), (j, m + 1))); links.append(((j, m + 1), (k, m + 2)))
                        if not (i == j == k): n_coincident_mixed += 1
    # 連結成分 ＝ 束
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for u, v in links: parent[find(u)] = find(v)
    comps = {}
    for i in range(4):
        for m in range(M + 1): comps.setdefault(find((i, m)), []).append((i, m))
    bundles = sorted(comps.values(), key=lambda c: (c[0][0], c[0][1]))
    out["bundle_closure"] = {"n_quadruples": len(phases_same) + len(phases_mixed), "same_pulse_phase_max_deg": max(phases_same), "mixed_pulse_phase_min_deg": min(phases_mixed),
                             "n_mixed_quadruples_closing_by_coincidence": n_coincident_mixed, "n_bundles_found": len(bundles), "bundles_are_the_pulses": all(len({e[0] for e in c}) == 1 and len(c) == M + 1 for c in bundles)}
    # (3) 距離の閉路：⟨q̄_{i,m} q_{j,m} q̄_{j,m+1} q_{i,m+1}⟩ の位相と軸、真の Δ_ij と比較
    dist = []
    for i in range(4):
        for j in range(4):
            if i == j: continue
            angs = []; axs = []
            for m in range(M):
                lam = mean_q(C.qmul(C.qmul(C.qmul(C.qconj(E[i, m]), E[j, m]), C.qconj(E[j, m + 1])), E[i, m + 1]))
                a_deg, ax = ang_axis(lam); angs.append(a_deg); axs.append(ax)
            angs = np.array(angs); true = abs(wrap(theta[i] - theta[j]))
            dist.append({"ij": f"{i+1}{j+1}", "Delta_true_abs_deg": true, "loop_angle_mean_deg": float(angs.mean()), "loop_angle_std_deg": float(angs.std()),
                         "max_err_deg": float(np.max(np.abs(angs - true))), "axis_mean": np.mean(axs, axis=0).round(4).tolist(), "n_i_plus_n_j_unit": C.unit(axes[i] + axes[j]).round(4).tolist()})
    out["distance_closure"] = dist; out["distance_closure_max_err_deg"] = max(d["max_err_deg"] for d in dist)
    # (4) 束の分配率 Ω_ij（符号つき）と |Y| 版
    R2 = np.array([Yd[:, i, i, 0].sum() for i in range(4)])
    Om = np.array([[Yd[:, i, j, 0].sum() / R2[i] for j in range(4)] for i in range(4)])
    Om_abs = np.array([[np.linalg.norm(Yd[:, i, j].sum(axis=0)) / R2[i] for j in range(4)] for i in range(4)])
    layer_abs = np.array([[np.mean(np.abs(Yd[:, i, j, 0]) / np.sqrt(Yd[:, i, i, 0] * Yd[:, j, j, 0] + 1e-300)) for j in range(4)] for i in range(4)])
    out["Omega_signed_Sc_over_R2"] = Om.round(6).tolist(); out["Omega_abs_normY_over_R2"] = Om_abs.round(6).tolist()
    out["theta_deg"] = theta; out["phi0_deg"] = phi0; out["layerwise_mean_abs_Sc_ratio"] = layer_abs.round(4).tolist(); out["Delta_deg"] = [[wrap(theta[i] - theta[j]) for j in range(4)] for i in range(4)]
    return out
res = {name: run(*cfg) for name, cfg in CONFIGS.items()}
json.dump(res, open(os.path.join(HERE, "results", "closure_distribution_ratio_v1.json"), "w"), ensure_ascii=False, indent=2)
for name, o in res.items():
    print("==", name); print(" reconstruction", f"{o['reconstruction_max_err']:.1e}", "| Y same winding max", f"{o['Y_max_norm_same_winding']:.3e}", "| different winding max", f"{o['Y_max_norm_different_winding']:.1e}")
    b = o["bundle_closure"]; print(" bundle closure:", b)
    print(" distance closure max err deg", f"{o['distance_closure_max_err_deg']:.2e}")
    for d in o["distance_closure"]: print(f"   {d['ij']}: loop {d['loop_angle_mean_deg']:9.4f} ± {d['loop_angle_std_deg']:.1e}  true |Δ| {d['Delta_true_abs_deg']:6.1f}  axis {d['axis_mean']}  (n_i+n_j) {d['n_i_plus_n_j_unit']}")
    print(" Omega signed:"); print(np.array(o["Omega_signed_Sc_over_R2"])); print(" Omega |Y|:"); print(np.array(o["Omega_abs_normY_over_R2"])); print(" layerwise mean |Sc| ratio:"); print(np.array(o["layerwise_mean_abs_Sc_ratio"]))
