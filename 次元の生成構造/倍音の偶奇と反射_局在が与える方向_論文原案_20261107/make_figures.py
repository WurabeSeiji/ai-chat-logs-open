# -*- coding: utf-8 -*-
"""論文原案「倍音の偶奇が決める反射と、局在が与える方向」の図（SVG）と数値（results/values.json）を作る。
閉じた弦 ψ ∈ [0, 2π) の上の実数パルス。倍音は有限和（13 項、等振幅）。対蹠写像 P f(ψ) = f(ψ + π)。"""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Hiragino Sans"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["svg.fonttype"] = "none"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "図表", "figures")
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
PRE = os.path.join(RES, "preview"); os.makedirs(PRE, exist_ok=True)
def save(fig, name):
    fig.savefig(os.path.join(OUT, name)); fig.savefig(os.path.join(PRE, name.replace(".svg", ".png")), dpi=110)
NG = 4096
psi = np.arange(NG) * 2 * np.pi / NG
deg = 180 / np.pi
vals = {}

def harmonics(kind, N=13):
    if kind == "odd": return np.arange(1, 2 * N, 2)          # 1, 3, …, 25（13 項）
    if kind == "even": return np.arange(2, 2 * N + 1, 2)     # 2, 4, …, 26（13 項）
    return np.arange(1, N + 1)                               # 1, …, 13（13 項）
def pulse(ms, theta=0.0, c=None):
    c = np.ones(len(ms)) / len(ms) if c is None else c
    return sum(cm * np.cos(m * (psi - theta)) for m, cm in zip(ms, c))
def P(f): return np.roll(f, NG // 2)                          # 対蹠写像 ψ → ψ + π
def inner(f, g): return np.sum(np.conj(f) * g) * 2 * np.pi / NG
def norm2(f): return inner(f, f).real

# ---------- 図1：偶奇と対蹠点の像 ----------
u = (psi - np.pi) * deg                                       # −180 … 180
fig, axs = plt.subplots(1, 3, figsize=(12, 3.2))
for ax, kind, lab in zip(axs, ["odd", "even", "both"], ["奇数倍音だけ（m = 1, 3, …, 25）", "偶数倍音だけ（m = 2, 4, …, 26）", "両方（m = 1, …, 13）"]):
    f = np.roll(pulse(harmonics(kind), theta=0.0), NG // 2)   # 中心を u = 0 に
    ax.plot(u, f, lw=1.4, color="#1f77b4")
    ax.axvline(0, color="#999", lw=0.6); ax.axvline(180, color="#999", lw=0.6, ls="--"); ax.axvline(-180, color="#999", lw=0.6, ls="--")
    ax.set_xlim(-180, 180); ax.set_xlabel("u = ψ − θ（度）"); ax.set_title(lab, fontsize=10)
    ax.set_xticks([-180, -90, 0, 90, 180])
axs[0].set_ylabel("S(u)")
for ax, kind in zip(axs, ["odd", "even", "both"]):
    f = pulse(harmonics(kind)); r = inner(f, P(f)).real / norm2(f)
    ax.text(0.02, 0.92, f"〈S, PS〉/‖S‖² = {r:+.3f}", transform=ax.transAxes, fontsize=9)
    vals[f"fig1_eta_{kind}"] = float(r)
fig.suptitle("図1　対蹠写像 P: ψ → ψ + π のもとで、奇数倍音は符号が反転し（逆相の像）、偶数倍音は同相の像を持ち、両方なら一山で像がない", fontsize=10)
fig.tight_layout(); save(fig, "fig1_parity_antipodal_image.svg"); plt.close(fig)

# ---------- 図2：偶奇の釣り合いと交換の可視度 ----------
Ae = pulse(harmonics("even")); Ae = Ae / math.sqrt(norm2(Ae))
Ao = pulse(harmonics("odd"));  Ao = Ao / math.sqrt(norm2(Ao))
chis = np.linspace(0, np.pi / 2, 181)
eta_num = []
for chi in chis:
    A = math.cos(chi) * Ae + math.sin(chi) * Ao
    eta_num.append(inner(A, P(A)).real / norm2(A))
eta_num = np.array(eta_num); eta_th = np.cos(2 * chis)
vals["fig2_eta_max_err"] = float(np.max(np.abs(eta_num - eta_th)))
Delta = np.linspace(0, 2 * np.pi, 361)
fig, axs = plt.subplots(1, 2, figsize=(11, 3.6))
axs[0].plot(chis * deg, eta_num, lw=1.6, color="#1f77b4", label="η = 〈A, PA〉/‖A‖²（数値）")
axs[0].plot(chis * deg, 4 * np.cos(chis) ** 2 * np.sin(chis) ** 2, lw=1.6, color="#2ca02c", label="可視度 V = 4 w_e w_o = 1 − η²")
axs[0].set_xlabel("混合角 χ（度）　A = cos χ·A_e + sin χ·A_o"); axs[0].set_xlim(0, 90); axs[0].legend(fontsize=8); axs[0].set_title("偶奇の釣り合いが交換の可視度を決める", fontsize=10)
axs[0].axhline(0, color="#999", lw=0.5)
maxerr = 0.0
for chi, col, lab in [(0.0, "#d62728", "χ = 0（偶数だけ）"), (np.pi / 8, "#ff7f0e", "χ = 22.5°"), (np.pi / 4, "#2ca02c", "χ = 45°（釣り合い、一山）"), (3 * np.pi / 8, "#9467bd", "χ = 67.5°"), (np.pi / 2, "#1f77b4", "χ = 90°（奇数だけ）")]:
    A = math.cos(chi) * Ae + math.sin(chi) * Ao
    Se = (A + P(A)) / 2; Do = (A - P(A)) / 2
    R = []
    for d in Delta:
        Ap = Se + np.exp(1j * d) * Do
        R.append(abs(inner(P(A), Ap)) ** 2 / norm2(A) ** 2)
    R = np.array(R); we, wo = math.cos(chi) ** 2, math.sin(chi) ** 2
    Rth = (we - wo) ** 2 + 4 * we * wo * np.sin(Delta / 2) ** 2
    maxerr = max(maxerr, float(np.max(np.abs(R - Rth))))
    axs[1].plot(Delta * deg, R, lw=1.6, color=col, label=lab)
vals["fig2_R_max_err"] = maxerr
axs[1].set_xlabel("交換位相 Δ_F（度）"); axs[1].set_ylabel("R̃(Δ_F) = |〈PA, A′〉|²/‖A‖⁴"); axs[1].set_xlim(0, 360); axs[1].set_ylim(-0.02, 1.05); axs[1].legend(fontsize=8); axs[1].set_xticks([0, 90, 180, 270, 360])
axs[1].set_title("A′ = A_e + e^{iΔ_F} A_o。釣り合いで R̃ = R = sin²(Δ_F/2)、純偶奇では R̃ ≡ 1 で移る量 0（0° と 90°、22.5° と 67.5° は重なる）", fontsize=9)
fig.suptitle("図2　対蹠像への重なり R̃(Δ_F) = (w_e − w_o)² + 4 w_e w_o sin²(Δ_F/2)。数値（13 項の等振幅パルス）と式の最大差 %.1e" % maxerr, fontsize=10)
fig.tight_layout(); save(fig, "fig2_parity_balance_visibility.svg"); plt.close(fig)

# ---------- 図3：倍音ごとの交差項から方向（符号つき Δ と Φ）を読む ----------
ms = harmonics("both"); c = np.ones(len(ms)) / len(ms)
def fit_phase(Sc, Vc):
    ph = np.unwrap(np.arctan2(Vc, Sc)); A_ = np.vstack([ms, np.ones(len(ms))]).T
    slope, icpt = np.linalg.lstsq(A_, ph, rcond=None)[0]; return -slope, icpt   # ϑ_m = Φ − mΔ：傾き −Δ、切片 Φ
fig, axs = plt.subplots(1, 3, figsize=(13, 3.4))
cases = [(+30.0, 20.0, "#1f77b4", "Δ = +30°, Φ = 20°"), (-30.0, 20.0, "#d62728", "Δ = −30°, Φ = 20°")]
for D_, F_, col, lab in cases:
    ph = F_ / deg - ms * D_ / deg; Sc = c * c * np.cos(ph); Vc = c * c * np.sin(ph)   # q_i q̄_j の順：α_i − α_j = Φ − mΔ
    axs[0].plot(ms, Sc / (c * c), "o-", color=col, lw=1.2, ms=4, label=lab)
    axs[1].plot(ms, Vc / (c * c), "o-", color=col, lw=1.2, ms=4, label=lab)
    s, f0 = fit_phase(Sc, Vc); vals[f"fig3_fit_D{D_:+.0f}"] = {"Delta_deg": float(s * deg), "Phi_deg": float(f0 * deg)}
# Φ = 0 のときは Sc だけでは ±Δ が区別できない
for D_, col, lab in [(+30.0, "#1f77b4", "Δ = +30°, Φ = 0"), (-30.0, "#d62728", "Δ = −30°, Φ = 0")]:
    axs[2].plot(ms, np.cos(ms * D_ / deg), "o-" if D_ > 0 else "x--", color=col, lw=1.2, ms=5, label=lab + "（Sc）")
    axs[2].plot(ms, np.sin(-ms * D_ / deg), "s-" if D_ > 0 else "+--", color=col, lw=0.8, ms=4, alpha=0.6, label=lab + "（Vec）")
axs[0].set_title("Sc 部 cos(Φ − mΔ)", fontsize=10); axs[1].set_title("Vec 部 sin(Φ − mΔ)：符号が ±Δ で違う", fontsize=10); axs[2].set_title("Φ = 0：Sc は ±Δ で同じ、Vec で分かれる", fontsize=10)
for ax in axs: ax.set_xlabel("倍音 m"); ax.axhline(0, color="#999", lw=0.5); ax.legend(fontsize=7)
fig.suptitle("図3　対の交差項 q_i q̄_j を倍音ごとに並べると、縞の位相 Φ − mΔ の傾きが −Δ（符号つき）、切片が Φ。基本波 m = 1 だけでは決まらない", fontsize=10)
fig.tight_layout(); save(fig, "fig3_harmonic_fringes_direction.svg"); plt.close(fig)

# 四元数：軸が違う対の Vec 部の三成分（n_i, n_j, n_i×n_j）
fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
for ax, (ti, tj, f0i, f0j, lab) in zip(axs, [(40.0, 10.0, 20.0, 0.0, "θ_i − θ_j = +30°"), (-20.0, 10.0, 20.0, 0.0, "θ_i − θ_j = −30°")]):
    ai = (f0i - ms * ti) / deg; aj = (f0j - ms * tj) / deg
    a_ = np.sin(ai) * np.cos(aj); b_ = -np.cos(ai) * np.sin(aj); cc = -np.sin(ai) * np.sin(aj)
    ax.plot(ms, a_, "o-", ms=4, lw=1.2, label="n_i の係数  sin α_i cos α_j"); ax.plot(ms, b_, "s-", ms=4, lw=1.2, label="n_j の係数  −cos α_i sin α_j"); ax.plot(ms, cc, "^-", ms=4, lw=1.2, label="n_i×n_j の係数  −sin α_i sin α_j")
    ax.set_title(lab + "（φ₀i = 20°, φ₀j = 0, θ_j = 10°）", fontsize=10); ax.set_xlabel("倍音 m"); ax.axhline(0, color="#999", lw=0.5); ax.legend(fontsize=7)
fig.suptitle("図3b　軸が違う対の Vec 部の三つの係数列（n_i、n_j、n_i×n_j 方向）。α = φ₀ − mθ。三列から向きと裏表が決まる", fontsize=10)
fig.tight_layout(); save(fig, "fig3b_quaternion_vec_components.svg"); plt.close(fig)

# ---------- 図4：分解能は共有できる帯域で決まる ----------
Dgrid = np.linspace(-90, 90, 1801) / deg
def corr(Mi, Mj):
    Mc = min(Mi, Mj); m = np.arange(1, Mc + 1); C = sum(np.cos(k * Dgrid) for k in m) / (Mi * Mj); return C / C[len(C) // 2]
fig, ax = plt.subplots(figsize=(7, 3.4))
widths = {}
for Mi, Mj, col in [(13, 13, "#1f77b4"), (13, 39, "#2ca02c"), (39, 39, "#d62728")]:
    C = corr(Mi, Mj); ax.plot(Dgrid * deg, C, lw=1.5, color=col, label=f"M_i = {Mi}, M_j = {Mj}（共有帯域 {min(Mi, Mj)}）")
    half = Dgrid[np.where(C >= 0.5)[0]]; w = float((half.max() - half.min()) * deg); widths[f"{Mi}_{Mj}"] = w
ax.axhline(0.5, color="#999", lw=0.5, ls="--"); ax.set_xlabel("相対位置 Δ（度）"); ax.set_ylabel("相互相関 C(Δ)/C(0)"); ax.legend(fontsize=8); ax.set_xlim(-90, 90)
ax.set_title("半値全幅：" + "、".join(f"({k.replace('_', ',')}) {v:.1f}°" for k, v in widths.items()) + "。幅 × M_共有 ≈ 205°（∝ 1/M_共有）", fontsize=9)
fig.suptitle("図4　対の相互相関の幅は共有できる帯域 min(M_i, M_j) で決まる", fontsize=10)
fig.tight_layout(); save(fig, "fig4_resolution_common_band.svg"); plt.close(fig)
vals["fig4_halfwidth_deg"] = widths; vals["fig4_halfwidth_times_Mcommon_deg"] = {k: v * min(int(a) for a in k.split("_")) for k, v in widths.items()}

# ---------- 図5：凝縮体 = グラム行列の階数が落ちる ----------
def gram_eigs(thetas, M):
    m = np.arange(1, M + 1); n = len(thetas); G = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            G[i, j] = np.sum(np.cos(m * (thetas[i] - thetas[j]))) / M
    return np.sort(np.linalg.eigvalsh(G))
svals = np.linspace(0, 60, 241)
fig, axs = plt.subplots(1, 2, figsize=(11, 3.4))
thr = 0.05; res = {}
for ax, M in zip(axs, [13, 39]):
    E = np.array([gram_eigs(np.array([0.0, s / deg, 120 / deg, 240 / deg]), M) for s in svals])
    for k in range(4): ax.plot(svals, E[:, k], lw=1.4, label=f"λ_{k + 1}")
    ax.axhline(thr, color="#999", lw=0.6, ls="--"); ax.set_xlabel("波 1 と波 2 の間隔 s（度）"); ax.set_ylabel("グラム行列の固有値"); ax.set_title(f"M = {M}：4 波（θ = 0, s, 120°, 240°）", fontsize=10); ax.legend(fontsize=8)
    below = svals[E[:, 0] < thr]; s_res = float(below.max()) if len(below) else float("nan"); res[M] = s_res
    ax.axvline(s_res, color="#d62728", lw=0.8); ax.text(s_res + 1, 0.5, f"λ_min < {thr} は s ≦ {s_res:.1f}°", fontsize=8, color="#d62728")
vals["fig5_rank_drop_s_deg"] = res
fig.suptitle("図5　間隔が分解能より小さいと最小固有値が閾値の下に沈み、有効階数 4 → 3。二つの波が一つの凝縮体として読まれる", fontsize=10)
fig.tight_layout(); save(fig, "fig5_gram_rank_condensate.svg"); plt.close(fig)

with open(os.path.join(RES, "values.json"), "w", encoding="utf-8") as f: json.dump(vals, f, ensure_ascii=False, indent=2)
print(json.dumps(vals, ensure_ascii=False, indent=2))
