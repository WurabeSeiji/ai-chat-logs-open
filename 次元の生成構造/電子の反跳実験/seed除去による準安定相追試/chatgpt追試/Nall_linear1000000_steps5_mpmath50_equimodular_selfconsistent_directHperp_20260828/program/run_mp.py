# -*- coding: utf-8 -*-
"""多倍長（mpmath, 50 桁）で、等モジュラー自己無撞着親からの seedless 走行を刻み 2π/1000000・5 step だけ行う（N=5,8,10,16,20）。
 - 親 v：各 N の倍精度パッケージ（…linear100000_steps50…）の data/states_treatment.npz の Z[0] を mp に変換（倍精度で求めた親をそのまま入力とする。残差 r は倍精度水準のまま）
 - 相互作用：K_ij = Im(conj(v_i) v_j)（隣接辺）、更新：z ← exp(Δ·K) z を Taylor 級数（‖ΔK‖≈1e-6 なので 8 項で 1e-50 未満）で厳密に評価。正規化なし、seed なし
 - 読出し：H⊥ = ‖Z − p(p·Z) − q(q·Z)‖²（直交成分の直接計算、mp）、H_total = ‖Z‖²
 - 参考：親の自己無撞着残差 r を mp で再評価し、(r·τ)² と比較
出力：data/N{N}_mp_timeseries.csv, results/summary.json, figures/compare_N_L1000000_5_mp50.png"""
import os, json, numpy as np, mpmath as mp
mp.mp.dps = 50
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); P = os.path.dirname(H)
L = 1000000; STEPS = 5; DELTA = 2 * mp.pi / L; NS = (5, 8, 10, 16, 20)
def edges(n): return [(i, j) for i in range(n) for j in range(i + 1, n)]
def adjacency(n):
    E = edges(n); m = len(E); A = [[0] * m for _ in range(m)]
    for a in range(m):
        for b in range(m):
            if a != b and len(set(E[a]) & set(E[b])) == 1: A[a][b] = 1
    return A
def K_amp(A, v):
    m = len(v); K = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            if A[i][j]: K[i, j] = mp.im(mp.conj(v[i]) * v[j])
    return K
def expK_apply(K, z, delta, terms=10):
    """exp(delta K) z = Σ_k (delta K)^k z / k!"""
    out = z.copy(); term = z.copy()
    for k in range(1, terms + 1):
        term = (K * term) * (delta / k); out = out + term
    return out
def norm2(z): return mp.re(sum(mp.conj(z[i]) * z[i] for i in range(len(z))))
def residual(K, v):
    kv = mp.mpc(0, 1) * (K * v); mu = sum(mp.conj(v[i]) * kv[i] for i in range(len(v))) / norm2(v)
    d = kv - mu * v; return mp.sqrt(norm2(d) / norm2(v))
summary = {}
for N in NS:
    src = os.path.join(P, f"N{N}_linear100000_steps50_equimodular_selfconsistent_directHperp_20260828", "data", "states_treatment.npz")
    v0 = np.load(src)["Z"][0]; m = len(v0); A = adjacency(N)
    v = mp.matrix([mp.mpc(float(x.real), float(x.imag)) for x in v0])
    pr = mp.matrix([mp.mpf(float(x)) for x in v0.real]); pr = pr / mp.sqrt(norm2(pr))
    qi = mp.matrix([mp.mpf(float(x)) for x in v0.imag]); qi = qi - (sum(pr[i] * qi[i] for i in range(m))) * pr; qi = qi / mp.sqrt(norm2(qi))
    K0 = K_amp(A, v); r = residual(K0, v); Ht0 = norm2(v)
    rows = []; Z = v.copy()
    for t in range(STEPS + 1):
        pz = sum(pr[i] * Z[i] for i in range(m)); qz = sum(qi[i] * Z[i] for i in range(m))
        Zp = Z - pz * pr - qz * qi; hperp = norm2(Zp); htot = norm2(Z)
        rows.append((t, mp.nstr(hperp, 20), mp.nstr(htot, 25), mp.nstr(hperp / htot, 20), mp.nstr((r * t * DELTA) ** 2, 12)))
        if t < STEPS: Z = expK_apply(K_amp(A, Z), Z, DELTA)
    with open(os.path.join(H, "data", f"N{N}_mp_timeseries.csv"), "w") as f:
        f.write("step,H_perp,H_total,H_perp_over_H_total,pred_(r_tau)^2\n"); [f.write(",".join(map(str, rw)) + "\n") for rw in rows]
    summary[f"N{N}"] = {"M": m, "parent_residual_r_mp": mp.nstr(r, 12), "H_total_0": mp.nstr(Ht0, 25), "H_total_drift": mp.nstr(abs(norm2(Z) - Ht0), 6),
                        "f_step1": rows[1][3], "f_step5": rows[5][3], "pred_step1": rows[1][4], "pred_step5": rows[5][4]}
    print(f"N={N}: r={mp.nstr(r,6)}  f(step1)={rows[1][3]}  (rΔ)²={rows[1][4]}  f(step5)={rows[5][3]}  (rτ)²={rows[5][4]}  H_total drift={mp.nstr(abs(norm2(Z)-Ht0),3)}")
json.dump(summary, open(os.path.join(H, "results", "summary.json"), "w"), indent=1)
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Hiragino Sans", "font.size": 11}); COL = {5: "#d7263d", 8: "#e67e22", 10: "#2e8b57", 16: "#1f5fd8", 20: "#6c3483"}
fig, ax = plt.subplots(figsize=(11, 6.5))
for N in NS:
    rows = [l.strip().split(",") for l in open(os.path.join(H, "data", f"N{N}_mp_timeseries.csv"))][1:]
    st = [int(x[0]) for x in rows]; f = [float(x[3]) for x in rows]; pred = [float(x[4]) for x in rows]
    ax.semilogy(st[1:], f[1:], "o-", color=COL[N], label=f"N={N} 実測（mp 50 桁）"); ax.semilogy(st[1:], pred[1:], "--", color=COL[N], alpha=.6, label=f"N={N} 予言 (r·τ)²")
ax.set_xlabel("step"); ax.set_ylabel("H⊥ / H_total（mp 50 桁、直接読出し）"); ax.set_title("Δ=2π/1000000, 5 step, 多倍長 50 桁：倍精度の床 1e-31 の下"); ax.grid(True, which="both", alpha=.25); ax.legend(ncol=2, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(H, "figures", "compare_N_L1000000_5_mp50.png"), dpi=160)
