# -*- coding: utf-8 -*-
"""四元数 4 パルスの交換散乱と読出し（2026-11-08、木原氏の指示「読出しと相互作用による位相の回転を取り込んだ数値実験系」）。

骨格
  状態    ：閉じた弦 ψ ∈ [0, 2π) の上の四元数波形 q_i(ψ) = W + V_i i + V_j j + V_k k、N = 4。初期値は index_v9q4.html の 4 パルス
            q_i(ψ) = A_i [c₀ e^{n_i φ₀i} + Σ_m c_m e^{n_i (φ₀i + m(ψ − θ_i))}]、c_m = (2/π) sin(mπφ/360)/m、c₀ = φ/360（位相幅 φ、倍音 M、奇偶）
  ティック：交換 1 巡。対 (i, j) を順に、7/13・7/18 の二チャネル交換写像 U = P_S + e^{φ_a} P_D を四元数で書いた
            S = (q_i + q_j)/2、D = (q_i − q_j)/2、q_i′ = S + D g、q_j′ = S − D g、g = e^{n_ij φ_a}（単位四元数、右から掛ける）
            φ_a = π + Δ_F,ij。軸を共有する複素の場合は 7/13 の r a + t b、t a + r b に一致する（対照 C1）。正規化はしない（線形核、ノルム保存）
  Δ_F,ij  ：与えない。対の内部時計の差 Φ_ij を波形の交差項から読む。倍音 m ごとの交差四元数 X_m = ⟨q_{i,m} q̄_{j,m}⟩_ψ の縞の位相
            ϑ_m = Φ − mΔ を m にわたって直線に当て、切片 Φ_ij、傾き −Δ_ij。R_ij = sin²(Φ_ij/2)
  軸 n_ij ：各波の軸 n_i は波形のベクトル部の主方向（読出し）、対の軸は n_ij = unit(n_i + n_j)
  読出し  ：波ごとに ノルム・L・N_eff・偶奇の重み w_e, w_o・η・軸・交換 1 巡での位相の回転 δ_i（累積 Θ_i）。対ごとに Φ_ij・Δ_ij・R_ij。
            グラム行列の固有値と有効階数（閾値 0.05）。初期状態からの距離（再帰）。観測者 o の「一体」の読み |q_o + Σ_{j≠o} q_j|²
  対照    ：C1 交換写像が複素の r,t に一致（軸共有、Δ_F を直接与える）、C2 縞の読出しが既知の Φ, Δ を復元（軸が違う対も）、
            C3 対の順序を逆にした走行との最終状態の距離、C4 ノルム保存。C0 は control_original_v1.py（原本の表の再現）
"""
import argparse, csv, json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results"); FIG = os.path.join(HERE, "figs")
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
DEG = math.pi / 180.0

# ---------------- 四元数（配列の最後の軸が 4 成分 [w, x, y, z]）----------------
def qmul(p, q):
    pw, px, py, pz = p[..., 0], p[..., 1], p[..., 2], p[..., 3]
    qw, qx, qy, qz = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([pw * qw - px * qx - py * qy - pz * qz,
                     pw * qx + px * qw + py * qz - pz * qy,
                     pw * qy - px * qz + py * qw + pz * qx,
                     pw * qz + px * qy - py * qx + pz * qw], axis=-1)
def qconj(q): return q * np.array([1.0, -1.0, -1.0, -1.0])
def qexp(n, phi):                       # e^{n φ}、n は単位ベクトル (3,)
    return np.concatenate([[math.cos(phi)], math.sin(phi) * np.asarray(n, dtype=float)])
def unit(v):
    L = np.linalg.norm(v); return v / L if L > 1e-15 else np.array([0.0, 0.0, 1.0])

# ---------------- 波形 ----------------
def coefficients(phi_deg, M, odd=True, even=True):
    ms = [m for m in range(1, M + 1) if (m % 2 == 1 and odd) or (m % 2 == 0 and even)]
    cs = [(2.0 / math.pi) * math.sin(m * math.pi * phi_deg / 360.0) / m for m in ms]
    return np.array(ms), np.array(cs), phi_deg / 360.0
def make_pulse(psi, theta_deg, phi0_deg, phi_deg, M, A, n, odd=True, even=True):
    ms, cs, c0 = coefficients(phi_deg, M, odd, even)
    beta = phi0_deg * DEG + np.outer(ms, psi - theta_deg * DEG)                      # (len(ms), NG)
    re = c0 * math.cos(phi0_deg * DEG) + np.sum(cs[:, None] * np.cos(beta), axis=0)
    im = c0 * math.sin(phi0_deg * DEG) + np.sum(cs[:, None] * np.sin(beta), axis=0)
    q = np.zeros((len(psi), 4)); q[:, 0] = A * re; q[:, 1:] = A * np.outer(im, unit(n))
    return q
def norm2(q, dpsi): return float(np.sum(q * q) * dpsi)
def inner(p, q, dpsi): return np.sum(qmul(qconj(p), q), axis=0) * dpsi           # ⟨p, q⟩ = ∫ p̄ q、四元数

# ---------------- 読出し ----------------
def antipodal(q): return np.roll(q, q.shape[0] // 2, axis=0)
def parity(q, dpsi):
    eta = inner(q, antipodal(q), dpsi)[0] / norm2(q, dpsi)
    return (1 + eta) / 2, (1 - eta) / 2, eta
def axis_of(q, prev=None):
    V = q[:, 1:]; C = V.T @ V; w, U = np.linalg.eigh(C); n = U[:, -1]
    if prev is not None and np.dot(n, prev) < 0: n = -n
    elif prev is None and n[np.argmax(np.abs(n))] < 0: n = -n
    return unit(n), float(w[-1] / max(np.trace(C), 1e-300))
def localization(q, dpsi):
    p = np.sum(q * q, axis=1); p = p / np.sum(p); return float(np.sum(p * p))
def harmonic_coeffs(q):
    """各成分の実フーリエ係数。q_m(ψ) = a_m cos mψ + b_m sin mψ（a_m, b_m は四元数）。m = 0..NG/2"""
    NG = q.shape[0]; F = np.fft.rfft(q, axis=0) / NG            # (NG/2+1, 4)
    a = 2 * F.real; b = -2 * F.imag; a[0] /= 2
    return a, b
def harmonic_distribution(q, Mmax):
    a, b = harmonic_coeffs(q); pw = np.sum(a * a + b * b, axis=1) / 2; pw[0] *= 2
    pw = pw[:Mmax + 1]; tot = pw.sum(); return pw / tot if tot > 0 else pw
def n_eff(dist): return float(np.sum(np.arange(len(dist)) * dist))
def cross_harmonics(qi, qj, Mmax):
    """倍音 m ごとの交差四元数 X_m = ⟨q_{i,m} q̄_{j,m}⟩_ψ = (a_i ā_j + b_i b̄_j)/2"""
    ai, bi = harmonic_coeffs(qi); aj, bj = harmonic_coeffs(qj)
    X = (qmul(ai, qconj(aj)) + qmul(bi, qconj(bj))) / 2
    return X[:Mmax + 1]
def read_pair(qi, qj, ni, nj, Mmax):
    """縞の位相 ϑ_m = Φ − mΔ を倍音にわたって直線に当て、Φ_ij（切片）と Δ_ij（−傾き）を読む。軸が違う対は
    X_m = (c_i c_j/2)[(1 + n_i·n_j) cos ϑ + sin ϑ (n_i + n_j) − cos ϑ (n_i×n_j)] なので、û = unit(n_i + n_j)、s = √((1 + n_i·n_j)/2) で tan ϑ = s (Vec·û)/Sc"""
    X = cross_harmonics(qi, qj, Mmax); ms = np.arange(1, Mmax + 1); X = X[1:]
    dot = float(np.dot(ni, nj)); u = unit(ni + nj); s = math.sqrt(max((1 + dot) / 2, 1e-12))
    Sc = X[:, 0]; Vu = X[:, 1:] @ u
    th = np.unwrap(np.arctan2(s * Vu, Sc)); w = np.sqrt(Sc * Sc + Vu * Vu); w = w / (w.sum() + 1e-300)
    keep = w > 1e-3 * w.max()
    A_ = np.vstack([ms[keep], np.ones(keep.sum())]).T; W = np.sqrt(w[keep])
    slope, icpt = np.linalg.lstsq(A_ * W[:, None], th[keep] * W, rcond=None)[0]
    resid = float(np.sqrt(np.sum(w[keep] * (th[keep] - (slope * ms[keep] + icpt)) ** 2)))
    Phi = (icpt + math.pi) % (2 * math.pi) - math.pi; Delta = -slope
    return Phi, Delta, u, resid, X
def exchange(qi, qj, n_ij, Delta_F):
    """U = P_S + e^{n_ij φ_a} P_D、φ_a = π + Δ_F、右乗。軸共有の複素では r a + t b、t a + r b に一致"""
    g = qexp(n_ij, math.pi + Delta_F)
    S = (qi + qj) / 2; D = (qi - qj) / 2; Dg = qmul(D, g[None, :])
    return S + Dg, S - Dg
def rotation_of(q0, q1, dpsi, n_ref):
    """交換 1 巡で波がどれだけ回ったか：⟨q0, q1⟩ の偏角。符号は自分の軸 n_ref 方向"""
    ip = inner(q0, q1, dpsi); vec = ip[1:]; ang = math.atan2(np.linalg.norm(vec), ip[0])
    sgn = 1.0 if np.dot(vec, n_ref) >= 0 else -1.0
    return sgn * ang

# ---------------- 本体 ----------------
def run(args):
    NG = args.ng; psi = np.arange(NG) * 2 * math.pi / NG; dpsi = 2 * math.pi / NG
    theta = [97, 194, 58, 300]; phi = [1, 1, 1, 1]; M = [args.M] * 4; A = [1, 1, 1, 1]
    axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]]; phi0 = args.phi0
    Q = [make_pulse(psi, theta[i], phi0[i], phi[i], M[i], A[i], axes[i]) for i in range(4)]
    Q0 = [q.copy() for q in Q]; norm0 = sum(norm2(q, dpsi) for q in Q)
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    if args.order == "reverse": pairs = pairs[::-1]
    Mmax = args.M
    axis_prev = [None] * 4; Theta = [0.0] * 4
    rows = []; pair_rows = []
    def readout(t, Q, delta):
        nonlocal axis_prev
        rec = {"tick": t}
        G = np.zeros((4, 4)); nrm = [math.sqrt(norm2(q, dpsi)) for q in Q]
        for i in range(4):
            for j in range(4): G[i, j] = inner(Q[i], Q[j], dpsi)[0] / (nrm[i] * nrm[j])
        ev = np.sort(np.linalg.eigvalsh(G))
        for i, q in enumerate(Q):
            n, conc = axis_of(q, axis_prev[i]); axis_prev[i] = n
            we, wo, eta = parity(q, dpsi); dist = harmonic_distribution(q, Mmax)
            rec.update({f"norm2_{i+1}": norm2(q, dpsi), f"L_{i+1}": localization(q, dpsi), f"Neff_{i+1}": n_eff(dist),
                        f"we_{i+1}": we, f"wo_{i+1}": wo, f"eta_{i+1}": eta, f"nx_{i+1}": n[0], f"ny_{i+1}": n[1], f"nz_{i+1}": n[2], f"axisconc_{i+1}": conc,
                        f"delta_{i+1}": delta[i], f"Theta_{i+1}": Theta[i]})
        for k in range(4): rec[f"gram_ev{k+1}"] = ev[k]
        rec["eff_rank"] = int(np.sum(ev > args.thr)); rec["norm_total"] = sum(norm2(q, dpsi) for q in Q)
        d2 = sum(norm2(Q[i] - Q0[i], dpsi) for i in range(4)); rec["dist_initial"] = math.sqrt(d2 / norm0)
        o = args.observer - 1; others = sum(Q[j] for j in range(4) if j != o)
        rec["one_body_total"] = norm2(Q[o] + others, dpsi); rec["one_body_sum"] = norm2(Q[o], dpsi) + norm2(others, dpsi)
        rec["one_body_cross"] = rec["one_body_total"] - rec["one_body_sum"]
        rows.append(rec)
    delta = [0.0] * 4; readout(0, Q, delta)
    for t in range(1, args.ticks + 1):
        Qb = [q.copy() for q in Q]
        for (i, j) in pairs:
            ni, _ = axis_of(Q[i], axis_prev[i]); nj, _ = axis_of(Q[j], axis_prev[j])
            Phi, Delta, u, resid, X = read_pair(Q[i], Q[j], ni, nj, Mmax)
            Delta_F = Phi
            Q[i], Q[j] = exchange(Q[i], Q[j], u, Delta_F)
            pair_rows.append({"tick": t, "pair": f"{i+1}{j+1}", "Phi_deg": Phi / DEG, "Delta_deg": Delta / DEG, "R": math.sin(Phi / 2) ** 2,
                              "fit_resid": resid, "ux": u[0], "uy": u[1], "uz": u[2], "n_dot": float(np.dot(ni, nj))})
        delta = [rotation_of(Qb[i], Q[i], dpsi, axis_prev[i] if axis_prev[i] is not None else np.array(axes[i], float)) for i in range(4)]
        for i in range(4): Theta[i] += delta[i]
        readout(t, Q, delta)
    tag = args.tag
    with open(os.path.join(RES, f"timeseries_{tag}.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(os.path.join(RES, f"pairs_{tag}.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(pair_rows[0].keys())); w.writeheader(); w.writerows(pair_rows)
    summary = {"tag": tag, "ticks": args.ticks, "NG": NG, "M": args.M, "theta": theta, "phi0": phi0, "axes": axes, "order": args.order, "observer": args.observer,
               "norm_total_initial": rows[0]["norm_total"], "norm_total_final": rows[-1]["norm_total"],
               "norm_drift_max": max(abs(r["norm_total"] - rows[0]["norm_total"]) for r in rows),
               "dist_initial_min_after_0": min(r["dist_initial"] for r in rows[1:]), "tick_of_min": int(min(rows[1:], key=lambda r: r["dist_initial"])["tick"]),
               "eff_rank_values": sorted(set(r["eff_rank"] for r in rows)), "final": rows[-1]}
    json.dump(summary, open(os.path.join(RES, f"summary_{tag}.json"), "w"), ensure_ascii=False, indent=2)
    return rows, pair_rows, summary, Q, Q0, dpsi, psi

# ---------------- 対照 ----------------
def controls(args):
    out = {}
    NG = args.ng; psi = np.arange(NG) * 2 * math.pi / NG; dpsi = 2 * math.pi / NG
    # C1：軸共有・複素の極限で 7/13 の r a + t b、t a + r b に一致
    kax = [0, 0, 1]; qa = make_pulse(psi, 97, 0, 1, args.M, 1, kax); qb = make_pulse(psi, 194, 0, 1, args.M, 1, kax)
    worst = 0.0
    for Dl in [0.0, math.pi / 2, 2 * math.asin(math.sqrt(0.55)), math.pi]:
        tt = complex(math.cos(Dl / 2) * math.cos(Dl / 2), math.cos(Dl / 2) * math.sin(Dl / 2))           # t = e^{iΔ/2} cos(Δ/2)
        rr = -1j * complex(math.cos(Dl / 2), math.sin(Dl / 2)) * math.sin(Dl / 2)                      # r = −i e^{iΔ/2} sin(Δ/2)
        za = qa[:, 0] + 1j * qa[:, 3]; zb = qb[:, 0] + 1j * qb[:, 3]                                  # 軸 k の平面を複素数で
        a2 = rr * za + tt * zb; b2 = tt * za + rr * zb
        qa2, qb2 = exchange(qa, qb, np.array(kax, float), Dl)
        worst = max(worst, float(np.max(np.abs(qa2[:, 0] + 1j * qa2[:, 3] - a2))), float(np.max(np.abs(qb2[:, 0] + 1j * qb2[:, 3] - b2))),
                    float(np.max(np.abs(qa2[:, 1:3]))), float(np.max(np.abs(qb2[:, 1:3]))))
    out["C1_max_diff_vs_7_13_rt"] = worst
    # C2：縞の読出し（軸共有と軸が違う対）
    errs = []
    for (ni, nj) in [([0, 0, 1], [0, 0, 1]), ([1, 0, 0], [0, 1, 0]), ([1, 0, 0], [1, 1, 1])]:
        for (Phi, Dlt) in [(20, 30), (20, -30), (-70, 97), (150, -136)]:
            qi = make_pulse(psi, 100 + Dlt, 10 + Phi, 1, args.M, 1, ni); qj = make_pulse(psi, 100, 10, 1, args.M, 1, nj)
            P_, D_, u, resid, X = read_pair(qi, qj, unit(np.array(ni, float)), unit(np.array(nj, float)), args.M)
            ePhi = abs((P_ / DEG - Phi + 180) % 360 - 180); eD = abs((D_ / DEG - Dlt + 180) % 360 - 180)
            errs.append({"ni": ni, "nj": nj, "Phi": Phi, "Delta": Dlt, "Phi_read": P_ / DEG, "Delta_read": D_ / DEG, "err_Phi": ePhi, "err_Delta": eD, "resid": resid})
    out["C2_fringe_readout"] = errs; out["C2_max_err_deg"] = max(max(e["err_Phi"], e["err_Delta"]) for e in errs)
    # C4：ノルム保存（1 巡）
    Q = [make_pulse(psi, th, p0, 1, args.M, 1, ax) for th, p0, ax in zip([97, 194, 58, 300], args.phi0, [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]])]
    n0 = sum(norm2(q, dpsi) for q in Q)
    for (i, j) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]:
        ni, _ = axis_of(Q[i]); nj, _ = axis_of(Q[j]); Phi, Dl, u, resid, X = read_pair(Q[i], Q[j], ni, nj, args.M); Q[i], Q[j] = exchange(Q[i], Q[j], u, Phi)
    out["C4_norm_rel_drift_one_round"] = abs(sum(norm2(q, dpsi) for q in Q) - n0) / n0
    return out

def plots(rows, pair_rows, tag):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "Hiragino Sans"; plt.rcParams["axes.unicode_minus"] = False
    T = [r["tick"] for r in rows]; col = ["#1f77b4", "#9467bd", "#2ca02c", "#8c564b"]
    # 図1 Θ_i n_i の軌道
    fig = plt.figure(figsize=(6.5, 5.5)); ax = fig.add_subplot(111, projection="3d")
    for i in range(4):
        X = [r[f"Theta_{i+1}"] / DEG * r[f"nx_{i+1}"] for r in rows]; Y = [r[f"Theta_{i+1}"] / DEG * r[f"ny_{i+1}"] for r in rows]; Z = [r[f"Theta_{i+1}"] / DEG * r[f"nz_{i+1}"] for r in rows]
        ax.plot(X, Y, Z, color=col[i], lw=1.3, label=f"波 {i+1}"); ax.scatter([X[-1]], [Y[-1]], [Z[-1]], color=col[i], s=25)
    ax.set_xlabel("x（度）"); ax.set_ylabel("y（度）"); ax.set_zlabel("z（度）"); ax.legend(fontsize=8)
    ax.set_title("図1　交換で回った位相の累積 Θ_i を自分の軸 n_i に乗せた軌道 Θ_i n_i", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(FIG, f"fig1_theta_n_{tag}.svg")); fig.savefig(os.path.join(FIG, f"fig1_theta_n_{tag}.png"), dpi=110); plt.close(fig)
    # 図2 Φ_ij と R_ij
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
    for lab in sorted(set(p["pair"] for p in pair_rows)):
        pr = [p for p in pair_rows if p["pair"] == lab]
        axs[0].plot([p["tick"] for p in pr], [p["Phi_deg"] for p in pr], lw=1.0, label=lab); axs[1].plot([p["tick"] for p in pr], [p["R"] for p in pr], lw=1.0, label=lab)
    axs[0].set_xlabel("ティック"); axs[0].set_ylabel("読んだ内部時計の差 Φ_ij（度）"); axs[0].legend(fontsize=7, ncol=3); axs[1].set_xlabel("ティック"); axs[1].set_ylabel("R_ij = sin²(Φ_ij/2)"); axs[1].set_ylim(-0.02, 1.02)
    fig.suptitle("図2　対ごとに波形の交差項から読んだ Φ_ij と、それが決める R_ij（与えていない）", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(FIG, f"fig2_phi_R_{tag}.svg")); fig.savefig(os.path.join(FIG, f"fig2_phi_R_{tag}.png"), dpi=110); plt.close(fig)
    # 図3 L・N_eff・η
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.4))
    for i in range(4):
        axs[0].plot(T, [r[f"L_{i+1}"] for r in rows], color=col[i], lw=1.1, label=f"波 {i+1}"); axs[1].plot(T, [r[f"Neff_{i+1}"] for r in rows], color=col[i], lw=1.1); axs[2].plot(T, [r[f"eta_{i+1}"] for r in rows], color=col[i], lw=1.1)
    axs[0].set_ylabel("局在性 L = Σ p²"); axs[1].set_ylabel("有効倍音次数 N_eff"); axs[2].set_ylabel("偶奇 η = w_e − w_o"); axs[0].legend(fontsize=8)
    for ax in axs: ax.set_xlabel("ティック")
    fig.suptitle("図3　交換で移る局在性と倍音次数、偶奇の重み", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(FIG, f"fig3_L_Neff_eta_{tag}.svg")); fig.savefig(os.path.join(FIG, f"fig3_L_Neff_eta_{tag}.png"), dpi=110); plt.close(fig)
    # 図4 グラム行列・再帰・ノルム
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.4))
    for k in range(4): axs[0].plot(T, [r[f"gram_ev{k+1}"] for r in rows], lw=1.1, label=f"λ_{k+1}")
    axs[0].axhline(0.05, color="#999", lw=0.6, ls="--"); axs[0].set_ylabel("グラム行列の固有値"); axs[0].legend(fontsize=8)
    axs[1].plot(T, [r["dist_initial"] for r in rows], lw=1.1, color="#d62728"); axs[1].set_ylabel("初期状態からの距離 ‖Ψ_t − Ψ_0‖/‖Ψ_0‖")
    axs[2].plot(T, [r["norm_total"] for r in rows], lw=1.1, color="#1f77b4", label="全ノルム"); axs[2].plot(T, [r["one_body_cross"] for r in rows], lw=1.1, color="#2ca02c", label=f"観測者の一体の読みの交差項"); axs[2].legend(fontsize=8)
    for ax in axs: ax.set_xlabel("ティック")
    fig.suptitle("図4　有効階数（閾値 0.05）、再帰、ノルム保存と一体の読み", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(FIG, f"fig4_gram_recurrence_norm_{tag}.svg")); fig.savefig(os.path.join(FIG, f"fig4_gram_recurrence_norm_{tag}.png"), dpi=110); plt.close(fig)
    # 図5 位相の回転 δ_i
    fig, ax = plt.subplots(figsize=(7, 3.2))
    for i in range(4): ax.plot(T, [r[f"delta_{i+1}"] / DEG for r in rows], color=col[i], lw=1.1, label=f"波 {i+1}")
    ax.set_xlabel("ティック"); ax.set_ylabel("交換 1 巡の位相の回転 δ_i（度）"); ax.legend(fontsize=8); ax.set_title("図5　相互作用による位相の回転", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, f"fig5_delta_{tag}.svg")); fig.savefig(os.path.join(FIG, f"fig5_delta_{tag}.png"), dpi=110); plt.close(fig)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ticks", type=int, default=200); ap.add_argument("--ng", type=int, default=2048); ap.add_argument("--M", type=int, default=100)
    ap.add_argument("--phi0", type=float, nargs=4, default=[0.0, 60.0, 120.0, 180.0]); ap.add_argument("--order", choices=["forward", "reverse"], default="forward")
    ap.add_argument("--observer", type=int, default=3); ap.add_argument("--thr", type=float, default=0.05); ap.add_argument("--tag", default="main")
    ap.add_argument("--controls", action="store_true")
    args = ap.parse_args()
    if args.controls:
        out = controls(args); json.dump(out, open(os.path.join(RES, "controls.json"), "w"), ensure_ascii=False, indent=2)
        print(json.dumps({k: v for k, v in out.items() if k != "C2_fringe_readout"}, ensure_ascii=False, indent=2))
        for e in out["C2_fringe_readout"]: print(f"  C2 n_i={e['ni']} n_j={e['nj']} Φ={e['Phi']} Δ={e['Delta']} → 読み Φ={e['Phi_read']:.3f} Δ={e['Delta_read']:.3f} 誤差 {e['err_Phi']:.2e}/{e['err_Delta']:.2e} resid={e['resid']:.2e}")
    else:
        rows, pair_rows, summary, Q, Q0, dpsi, psi = run(args); plots(rows, pair_rows, args.tag)
        print(json.dumps({k: v for k, v in summary.items() if k != "final"}, ensure_ascii=False, indent=2))
