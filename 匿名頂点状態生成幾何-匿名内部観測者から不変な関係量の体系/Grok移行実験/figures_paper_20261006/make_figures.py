#!/usr/bin/env python3
"""論文用の図（4 系統：水素原子 e–p、電子–電子、陽子–陽子、水素分子 H–H）。
物理の計算は既存のプログラムをそのまま使う（import して関数を呼ぶか、保存済みの JSON を読む）。本スクリプトは描画と、描いた数値の data/*.json への保存だけ。
  水素原子      : ../exchange_engine_20261006/engine.py ＋ system_hydrogen.py（固定一式 ../exchange_rel_20261005 に回帰済み）、
                  層 1 の一原子の履歴は ../exchange_general_20261006/absorber_direction_recoil.py の関数（乱数列 20261006）
  ee・pp        : ../exchange_general_20261006/scattering_same_sign.py の mott / quadrupole_spectrum
  H₂            : ../exchange_engine_20261006/{audit_h2_*.json, validation_h2.json, radiative_association_h2.json}
出力：figures/fig01〜fig15 .png、data/*.json、figures_index.md
"""
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
ENGINE = HERE.parent / "exchange_engine_20261006"
GENERAL = HERE.parent / "exchange_general_20261006"
FIXED = HERE.parent / "exchange_rel_20261005"
FIG = HERE / "figures"; FIG.mkdir(exist_ok=True)
DATA = HERE / "data"; DATA.mkdir(exist_ok=True)
for p in (ENGINE, GENERAL, FIXED):
    sys.path.insert(0, str(p))
from matplotlib import font_manager  # noqa: E402
font_manager.fontManager.addfont = lambda path: None   # 固定一式の import 時のフォント登録を無効化（英語ラベルのみ使う）
import engine as EN  # noqa: E402
import system_hydrogen as SH  # noqa: E402
X = SH.X

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.dpi": 150, "legend.frameon": False})
KIND_COLOR = {"e1": "C0", "e2": "C2", "m1": "C1", "2ph": "C3", "gw": "C4", "e1_stim": "C0", "e2_stim": "C2"}
KIND_NAME = {"e1": "E1", "e2": "E2", "m1": "M1", "2ph": "2-photon", "gw": "GW"}
HC_EV_UM = 1.239841984   # eV·μm
INDEX = []


def jdump(name, obj):
    def conv(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, dict):
            return {str(k): conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [conv(v) for v in o]
        return o
    (DATA / name).write_text(json.dumps(conv(obj), indent=1), encoding="utf-8")


def save(fig, name, caption):
    fig.tight_layout()
    fig.savefig(FIG / (name + ".png"))
    plt.close(fig)
    INDEX.append((name, caption))
    print("saved", name)


# =====================================================================
# 系 1：水素原子（エンジン＋水素部品 = 固定一式）
# =====================================================================
states_H, channels_H, init_H = SH.build()
nH = len(states_H)
p0H = np.zeros(nH); p0H[init_H] = 1.0
R_H, _ = EN.generator(states_H, channels_H, 0.0)
tau_H, pinf_H = EN.exact(R_H, p0H)
out_H = -np.diag(R_H)
E1s = states_H[X.IDX[(1, 0, 1, -1)]]["E_eV"]
exc = np.array([s["E_eV"] - E1s for s in states_H])        # 励起エネルギー [eV]


def fig01_hydrogen_levels():
    """準位図（n, ℓ, j）と、光子数で重みづけた遷移（F で和）。縦軸は n（目盛に励起エネルギー）。"""
    agg = {}
    for a, b, kind, k, A in channels_H:
        nq = A * tau_H[a]
        if nq < 1e-4:
            continue
        sa, sb = states_H[a], states_H[b]
        key = ((sa["n"], sa["l"], sa["j"]), (sb["n"], sb["l"], sb["j"]), kind)
        agg[key] = agg.get(key, 0.0) + nq
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    def pos(n, l, j):
        return l + (0.18 if j > l else -0.18 if l > 0 else 0.0), n
    seen = set()
    for s in states_H:
        key = (s["n"], s["l"], s["j"])
        if key in seen:
            continue
        seen.add(key)
        x, y = pos(*key)
        ax.plot([x - 0.14, x + 0.14], [y, y], color="k", lw=1.4)
    for (ka, kb, kind), nq in sorted(agg.items(), key=lambda kv: kv[1]):
        xa, ya = pos(*ka); xb, yb = pos(*kb)
        ls = "--" if kind == "2ph" else ":" if kind in ("m1", "e2") else "-"
        ax.annotate("", xy=(xb, yb + 0.03), xytext=(xa, ya - 0.03),
                    arrowprops=dict(arrowstyle="-|>", color=KIND_COLOR[kind], lw=0.6 + 3.5 * np.sqrt(nq), ls=ls, alpha=0.85, shrinkA=0, shrinkB=0))
    for n in range(1, 6):
        e = 13.6 * (1 - 1 / n**2) if n > 1 else 0.0
        ax.text(4.75, n, "n=%d  %.2f eV" % (n, float(exc[[i for i, s in enumerate(states_H) if s["n"] == n and s["l"] == 0][0]])), va="center", fontsize=8)
    ax.set_xticks(range(5)); ax.set_xticklabels(["s", "p", "d", "f", "g"])
    ax.set_yticks(range(1, 6)); ax.set_ylabel("principal quantum number n (levels: j = l ± 1/2)")
    ax.set_xlim(-0.5, 5.6); ax.set_ylim(0.6, 5.4)
    for kind in ("e1", "2ph", "m1", "e2"):
        ax.plot([], [], color=KIND_COLOR[kind], ls="--" if kind == "2ph" else ":" if kind in ("m1", "e2") else "-", label=KIND_NAME[kind])
    ax.legend(loc="lower right", title="arrow width ∝ √(photons per cascade)")
    ax.set_title("H atom: 50 levels (n ≤ 5, hyperfine summed) and the cascade from 5p$_{3/2}$ F=1")
    save(fig, "fig01_H_levels_cascade", "水素原子：準位図と 5p₃/₂F1 からのカスケード（矢印の太さ ∝ √光子数、F は和）")
    jdump("fig01_H_cascade_flows.json", {"%s->%s:%s" % (ka, kb, kind): nq for (ka, kb, kind), nq in agg.items()})


def fig02_hydrogen_spectrum():
    """放出スペクトル：光子数（期待値）対 ħω、種別で色分け。下段：線幅 Γ = Γ_src + Γ_dst。"""
    rows = []
    for a, b, kind, k, A in channels_H:
        nq = A * tau_H[a]
        if nq <= 0:
            continue
        hw = states_H[a]["E_eV"] - states_H[b]["E_eV"]
        rows.append(dict(hw=hw, n=nq, kind=kind, G=out_H[a] + out_H[b], A=A, src=states_H[a]["name"], dst=states_H[b]["name"]))
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True)
    ax = axes[0]
    for kind in ("e1", "2ph", "e2", "m1", "gw"):
        rr = [r for r in rows if r["kind"] == kind]
        if not rr:
            continue
        ax.scatter([r["hw"] for r in rr], [r["n"] for r in rr], s=18, color=KIND_COLOR[kind], label=KIND_NAME[kind], zorder=3)
        ax.vlines([r["hw"] for r in rr], 1e-50, [r["n"] for r in rr], color=KIND_COLOR[kind], lw=0.6, alpha=0.6)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-50, 3)
    ax.set_ylabel("photons per cascade (expectation)")
    ax.legend(ncol=5, loc="upper center", bbox_to_anchor=(0.45, 1.0), fontsize=8)
    for name, key in (("Lyman-α", ((2, 1, 1, -1), (1, 0, 1, -1))), ("21 cm", ((1, 0, 1, 1), (1, 0, 1, -1)))):
        hw = states_H[X.IDX[key[0]]]["E_eV"] - states_H[X.IDX[key[1]]]["E_eV"]
        ax.annotate(name, xy=(hw, 1.0), xytext=(hw, 1.8), ha="center", fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.5))
    ax.set_title("H atom: emission record of the cascade (T = 0 absorber)")
    ax = axes[1]
    for kind in ("e1", "2ph", "e2", "m1", "gw"):
        rr = [r for r in rows if r["kind"] == kind]
        ax.scatter([r["hw"] for r in rr], [r["G"] for r in rr], s=14, color=KIND_COLOR[kind])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("photon energy ħω [eV]"); ax.set_ylabel("line width Γ = Γ$_{src}$ + Γ$_{dst}$ [s$^{-1}$]")
    ax2 = ax.twinx(); ax2.set_yscale("log"); ax2.set_ylim(EN.C / np.array(ax.get_ylim())[::-1]); ax2.set_ylabel("coherence length c/Γ [m]")
    ax2.spines["right"].set_visible(True)
    save(fig, "fig02_H_spectrum_widths", "水素原子：放出記録（光子数対 ħω、種別）と線幅 Γ・可干渉長 c/Γ")
    jdump("fig02_H_lines.json", rows)


def fig03_hydrogen_time():
    """時間発展（Radau）：殻 n ごとの占有と、種別ごとの累積放出。"""
    times, hist = EN.evolve(R_H, p0H, t_max=1e17, n_times=300)
    tplot = np.where(times > 0, times, 1e-13)
    cum = {}
    for a, b, kind, k, A in channels_H:
        hw = states_H[a]["E_eV"] - states_H[b]["E_eV"]
        flux = A * hist[:, a] * hw
        c = np.concatenate([[0.0], np.cumsum(0.5 * (flux[1:] + flux[:-1]) * np.diff(times))])
        cum[kind] = cum.get(kind, 0.0) + c
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    ax = axes[0]
    for n in range(1, 6):
        ax.plot(tplot, sum(hist[:, i] for i, s in enumerate(states_H) if s["n"] == n), lw=1.5, label="n=%d" % n)
    ax.set_xscale("log"); ax.set_xlabel("time [s]"); ax.set_ylabel("occupation"); ax.legend(ncol=2, fontsize=8)
    ax.set_title("occupation by shell (master equation, Radau)")
    ax = axes[1]
    for kind in ("e1", "2ph", "e2", "m1"):
        ax.plot(tplot, cum[kind], lw=1.5, color=KIND_COLOR[kind], label=KIND_NAME[kind])
    ax.plot(tplot, cum["gw"] * 1e47, lw=1.5, color=KIND_COLOR["gw"], label="GW ×10$^{47}$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-8, 30)
    ax.set_xlabel("time [s]"); ax.set_ylabel("cumulative emitted energy [eV]"); ax.legend(fontsize=8)
    ax.set_title("emission by kind (total 13.0545 eV)")
    save(fig, "fig03_H_time_evolution", "水素原子：殻ごとの占有の時間発展と種別ごとの累積放出")
    jdump("fig03_H_time.json", dict(times=times, cum={k: v for k, v in cum.items()}))


def fig04_hydrogen_layer2():
    """層 2：21 cm の二準位が背景 2.7255 K と平衡する過程、定常比の温度依存。"""
    i1, i0 = X.IDX[(1, 0, 1, 1)], X.IDX[(1, 0, 1, -1)]
    hw21 = states_H[i1]["E_eV"] - states_H[i0]["E_eV"]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    ax = axes[0]
    out = {}
    for T, ls in ((0.0, "--"), (2.7255, "-")):
        R, _ = EN.generator(states_H, channels_H, T)
        times, hist = EN.evolve(R, p0H, t_max=1e17, n_times=300)
        tplot = np.where(times > 0, times, 1e-13)
        ax.plot(tplot, hist[:, i1], ls=ls, color="C0", label="1s F=1, T=%g K" % T)
        ax.plot(tplot, hist[:, i0], ls=ls, color="C3", label="1s F=0, T=%g K" % T)
        out[str(T)] = dict(times=times, F1=hist[:, i1], F0=hist[:, i0])
    ax.set_xscale("log"); ax.set_xlabel("time [s]"); ax.set_ylabel("occupation"); ax.legend(fontsize=8)
    ax.axhline(0.745286, color="C0", lw=0.5, alpha=0.5); ax.axhline(0.254714, color="C3", lw=0.5, alpha=0.5)
    ax.set_title("layer 2: hyperfine ground doublet vs. absorber temperature")
    ax = axes[1]
    Ts = np.array([0.3, 0.5, 1.0, 2.0, 2.7255, 5.0, 10.0])
    ratios = []
    for T in Ts:
        R, _ = EN.generator(states_H, channels_H, T)
        _, hist = EN.evolve(R, p0H, t_max=1e17, n_times=3, rtol=1e-12)
        ratios.append(hist[-1][i1] / hist[-1][i0])
    Tc = np.linspace(0.2, 12, 300)
    ax.plot(Tc, 3 * np.exp(-hw21 / (EN.K_B_EV * Tc)), color="k", lw=1, label="3 e$^{-ħω/kT}$")
    ax.plot(Ts, ratios, "o", color="C1", label="engine steady state")
    ax.set_xlabel("absorber temperature T [K]"); ax.set_ylabel("p(F=1)/p(F=0) at t = 10$^{17}$ s"); ax.legend(fontsize=8)
    ax.set_title("21 cm spin temperature = absorber temperature")
    save(fig, "fig04_H_layer2_temperature", "水素原子 層 2：21 cm 二準位の平衡化（T = 0 と 2.7255 K）と定常比の温度依存（3e^{−ħω/kT} との一致）")
    out["ratio_vs_T"] = dict(T=Ts, ratio=np.array(ratios))
    jdump("fig04_H_layer2.json", out)


def fig05_hydrogen_layer1():
    """層 1：集団の Doppler 記録から v⃗ を復元（左・中）、一原子の 21 cm 事象列から蹴りの履歴を復元（右、absorber_direction_recoil の関数）。"""
    import absorber_direction_recoil as ADR
    rng = np.random.default_rng(7)
    v0 = np.array([0.0, 0.0, 370e3]); beta = 370e3 / EN.C
    ev = []
    for a in range(200):
        ev += [e for e in EN.simulate_atom(states_H, channels_H, 2.7255, v0, rng, init_H, n_events=8) if e["hw_lab"] > 1.0 and "2ph" not in e["kind"]]
    N = np.array([e["n_lab"] for e in ev]); d = np.array([e["delta"] for e in ev]); w = np.array([e["rel_width"] for e in ev])
    from scipy.optimize import least_squares
    b0 = np.linalg.lstsq(N, d, rcond=None)[0]
    sol = least_squares(lambda bv: (d - np.array([EN.doppler_delta(n, bv * EN.C) for n in N])) / w, x0=b0, loss="cauchy", f_scale=1.0)
    vfit = sol.x * EN.C
    npar = N[:, 2]
    model = np.array([EN.doppler_delta(n, vfit) for n in N])
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    ax = axes[0]
    ax.scatter(npar, d / beta, s=6, alpha=0.6, label="photons (Lyman/Balmer…, 200 atoms)")
    xx = np.linspace(-1, 1, 200)
    ax.plot(xx, [EN.doppler_delta(np.array([np.sqrt(1 - x * x), 0, x]), v0) / beta for x in xx], color="k", lw=1, label="exact Doppler, v⃗₀")
    ax.set_xlabel("n̂·v̂₀ (cos θ of photon direction)"); ax.set_ylabel("δ/β = (ω$_{lab}$/ω′ − 1)/β"); ax.legend(fontsize=8)
    ax.set_title("record (δ, n̂) of %d photons" % len(ev))
    ax = axes[1]
    res_ms = (d - model) * EN.C
    ax.hist(res_ms, bins=40, color="C0", alpha=0.8)
    ax.set_xlabel("residual (δ − model)·c [m/s]"); ax.set_ylabel("count")
    ax.set_title("fit: |v⃗_fit − v⃗₀| = %.2f m/s (first-order guess %.0f m/s)" % (np.linalg.norm(vfit - v0), np.linalg.norm(b0 * EN.C - v0)))
    ax = axes[2]
    T0 = 2.7255
    ch = ADR.build_channels(T0); tab = ADR.line_table(T0)
    rngA = np.random.default_rng(20261006)
    evA = ADR.simulate_atom(rngA, T0, np.array([0.0, 0.0, ADR.V0]), n_events=30, ch=ch)
    v0_fit, vs, nseq = ADR.reconstruct_history(evA, tab, T0, ch)
    seq_ev = [e for e in evA if "2ph" not in e["kind"] and e["hw_lab"] <= 1e-4 * ADR.E_SI]
    true_v = [np.linalg.norm(e["v_before"]) for e in seq_ev]
    v_true0 = seq_ev[0]["v_before"]                     # 21 cm 列の最初の事象の直前の速度（Lyman-α の蹴りの後）= 復元の対象
    rec_v = [np.linalg.norm(v) for v in vs] if vs is not None else []
    k = np.arange(len(true_v))
    ax.plot(k, (np.array(true_v) - true_v[0]) * 1e6, "o-", ms=3, lw=0.8, label="true |v⃗| before each 21 cm event")
    if rec_v:
        ax.plot(np.arange(len(rec_v)), (np.array(rec_v) - true_v[0]) * 1e6, "x--", ms=4, lw=0.8, label="reconstructed from record")
    ax.set_xlabel("21 cm event index (one atom, after the Lyman-α kick)"); ax.set_ylabel("|v⃗| − |v⃗| at first 21 cm event [μm/s]"); ax.legend(fontsize=8)
    err0 = np.linalg.norm(v0_fit - v_true0) if v0_fit is not None else float("nan")
    ax.set_title("one atom, 21 cm kicks (1.9 μm/s each): |v⃗₀ fit − true| = %.1e m/s, %d events" % (err0, nseq), fontsize=9)
    save(fig, "fig05_H_layer1_doppler_recoil", "水素原子 層 1：集団 200 原子の Doppler 記録と厳密模型での v⃗ 復元、一原子の 21 cm 事象列からの蹴り履歴の復元")
    jdump("fig05_H_layer1.json", dict(npar=npar, delta=d, residual_ms=res_ms, vfit=vfit, v_lin=b0 * EN.C, true_v=true_v, rec_v=rec_v, v0_fit=v0_fit, v_true0=v_true0, err0=err0))


def fig06_hydrogen_layer3():
    """層 3：m_F = 0 の 5p→1s F0 の角分布（π 光）、2.7255 K での 1s F=1 の整列 vs 背景の四重極 a₂。"""
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    ax = axes[0]
    th = np.linspace(0, np.pi, 400)
    ax.plot(np.degrees(th), [EN.pattern(1, 0, np.cos(t)) * 4 * np.pi for t in th], label="Δm = 0 (π): 3/2·sin²θ")
    ax.plot(np.degrees(th), [EN.pattern(1, 1, np.cos(t)) * 4 * np.pi for t in th], label="Δm = ±1 (σ): 3/4·(1+cos²θ)")
    ax.axhline(1.0, color="k", lw=0.8, ls="--", label="isotropic")
    ax.set_xlabel("θ from the mark axis ẑ [deg]"); ax.set_ylabel("4π·P(θ)"); ax.legend(fontsize=8)
    ax.set_title("E1 angular patterns; 5p$_{3/2}$F1 m=0 → 1s F0 is pure π")
    ax = axes[1]
    beta = 370e3 / EN.C
    i1 = X.IDX[(1, 0, 1, 1)]
    Rm0, mb, edm, _ = EN.generator_m(states_H, channels_H, 0.0)
    midx = {s: k for k, s in enumerate(mb)}
    p0m = np.zeros(len(mb)); p0m[midx[(init_H, 0)]] = 1.0
    a2s = [0.0, 1e-6, 4e-6, 1e-5]
    al = []
    for a2 in a2s:
        Rm, mb, edm, _ = EN.generator_m(states_H, channels_H, 2.7255, beta, a2)
        _, histm = EN.evolve(Rm, p0m, n_times=2, rtol=1e-12)
        pm = {mm: histm[-1][midx[(i1, mm)]] for mm in (-1, 0, 1)}
        al.append((pm[1] + pm[-1] - 2 * pm[0]) / sum(pm.values()))
    ax.plot([a / 1e-6 for a in a2s], al, "o-", color="C1")
    ax.set_xlabel("background quadrupole a₂ [10$^{-6}$] (dipole β = 1.23×10$^{-3}$ fixed)"); ax.set_ylabel("alignment (p₊+p₋−2p₀) of 1s F=1")
    ax.set_title("steady alignment of 1s F=1 at 2.7255 K vs. background quadrupole a₂", fontsize=9)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    save(fig, "fig06_H_layer3_alignment", "水素原子 層 3：Δm の型の角分布（m_F=0 の 5p→1s F0 は純 π）と 2.7255 K の定常整列 vs 背景の四重極 a₂")
    jdump("fig06_H_layer3.json", dict(a2=a2s, alignment=al))


# =====================================================================
# 系 2・3：電子–電子、陽子–陽子（scattering_same_sign の関数）
# =====================================================================
def pair_params(SS, pair):
    T = SS.T
    p = T.PAIRS[pair]; c = T.pair_couplings(p)
    mu = c["mu"]; q2 = p["q_a"] * p["q_b"] / 9.0
    v_c = np.sqrt(2.0 * SS.E_REL / (mu * SS.ME_EV))
    eta = q2 * SS.ALPHA / v_c
    k = mu * v_c / SS.ALPHA
    a_c = eta / k; b = SS.L_HBAR / k; e = np.sqrt(1.0 + (b / a_c) ** 2)
    return dict(mu=mu, q2=q2, v_c=v_c, eta=eta, k=k, a_c=a_c, b=b, e=e, name=p["name"])


def fig07_mott():
    import scattering_same_sign as SS      # import で scattering_same_sign.md が同じ内容で再生成される
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    th = np.radians(np.linspace(5, 175, 1200))
    out = {}
    for ax, pair in zip(axes, ("ee", "pp")):
        pp = pair_params(SS, pair)
        m, r = SS.mott(th, pp["eta"], pp["k"])
        f = (SS.A0_M * 1e9) ** 2
        ax.plot(np.degrees(th), r * f, color="0.5", lw=1, label="distinguishable (both counted)")
        ax.plot(np.degrees(th), m * f, color="C0", lw=1, label="Mott (identical spin-½, unpolarized)")
        ax.set_yscale("log"); ax.set_xlabel("c.m. scattering angle θ [deg]"); ax.set_ylabel("dσ/dΩ [nm²/sr]")
        ax.set_title("%s: η = %.2f, E = %.3f eV, k = %.3f/a₀" % (pair, pp["eta"], SS.E_REL, pp["k"])); ax.legend(fontsize=8)
        out[pair] = dict(theta_deg=np.degrees(th), mott=m * f, distinguishable=r * f, params={k: v for k, v in pp.items() if k != "name"})
    save(fig, "fig07_ee_pp_mott", "電子–電子・陽子–陽子：Mott 断面積（同種スピン ½、90° で 1/2）と区別できる粒子の古典値。E = |E(5p₃/₂F1)|")
    jdump("fig07_mott.json", out)


def fig08_bremsstrahlung():
    import scattering_same_sign as SS
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    out = {}
    for ax, pair in zip(axes, ("ee", "pp")):
        pp = pair_params(SS, pair)
        om, Sq, P_int, pars = SS.quadrupole_spectrum(pp["e"])
        pref = (pp["q2"] * 27.211386245988 / pp["a_c"]) * (pp["a_c"] * SS.A0_M / (pp["v_c"] * X.C_SI)) * pp["v_c"] ** 5 / (180 * np.pi)
        unit_w = pp["v_c"] * X.C_SI / (pp["a_c"] * SS.A0_M)
        hw = om * unit_w * EN.HBAR_EVS
        dEdw = pref * Sq
        E_rad = np.trapezoid(dEdw, om * unit_w)                 # ∫dE/dω dω [eV]
        hw_char = unit_w * EN.HBAR_EVS                          # ħ(v/a)：特徴的な光子エネルギー
        Nph = E_rad / hw_char                                   # scattering_same_sign.md と同じ定義（E_rad / ħ(v/a)）
        mask = dEdw > 1e-9 * dEdw.max()                         # 数値積分の床（1e-9 以下）は描かない
        ax.plot(hw[mask], dEdw[mask], color="C2", lw=1.2)
        ax.axvline(hw_char, color="0.6", lw=0.6, ls="--"); ax.text(hw_char, dEdw.max() * 0.5, " ħv/a", fontsize=8, color="0.4")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("ħω [eV]"); ax.set_ylabel("dE/dω [eV·s]")
        ax.set_title("%s: quadrupole bremsstrahlung, b = L/μv, e = %.3f\nE_rad = %.2e eV, photons E_rad/ħ(v/a) = %.1e, Parseval %.4f" % (pair, pp["e"], E_rad, Nph, pars / P_int), fontsize=9)
        out[pair] = dict(hw_eV=hw, dEdw=dEdw, E_rad_eV=E_rad, hw_char_eV=hw_char, photons=Nph, e=pp["e"], parseval_ratio=pars / P_int)
    save(fig, "fig08_ee_pp_bremsstrahlung", "電子–電子・陽子–陽子：四重極制動放射のスペクトル（古典、双極子は恒等的に 0）と一回の衝突の光子期待値")
    jdump("fig08_bremsstrahlung.json", out)


# =====================================================================
# 系 4：H₂（保存済み JSON）
# =====================================================================
def load(name):
    p = ENGINE / name
    if not p.exists():
        raise SystemExit("missing %s — run ../exchange_engine_20261006/run_all.sh first" % p)
    return json.loads(p.read_text())


def fig09_h2_levels():
    A13 = load("audit_h2_v1J3.json"); A141 = load("audit_h2_v14J1.json")
    st = A13["states"]
    E00 = [s["E_eV"] for s in st if s["label"] == [0, 0]][0]
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for s in st:
        ax.plot([s["J"] - 0.35, s["J"] + 0.35], [s["E_eV"] - E00] * 2, color="0.6" if s["I"] == 0 else "0.3", lw=0.8)
    for A, color, lab in ((A13, "C0", "from (1,3)"), (A141, "C3", "from (14,1)")):
        for ln in A["lines"]:
            if ln["n"] < 0.02:
                continue
            a, b = st[ln["src"]], st[ln["dst"]]
            ax.annotate("", xy=(b["J"], b["E_eV"] - E00), xytext=(a["J"], a["E_eV"] - E00),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.4 + 3.0 * np.sqrt(ln["n"]), alpha=0.7, shrinkA=0, shrinkB=0))
        ax.plot([], [], color=color, label=lab)
    for v in (0, 1, 2, 5, 10, 14):
        s0 = [s for s in st if s["label"] == [v, 0]][0]
        ax.text(-1.2, s0["E_eV"] - E00, "v=%d" % v, fontsize=8, va="center")
    ax.set_xlabel("rotational quantum number J (even: para I=0, odd: ortho I=1)"); ax.set_ylabel("energy above (0,0) [eV]  (dissociation at %.4f eV)" % (-E00))
    ax.set_xlim(-1.8, 31.8); ax.legend(title="E2 cascades (arrow ∝ √photons, n > 0.02)", fontsize=8)
    ax.set_title("H₂ X¹Σg⁺: 301 bound rovibrational levels (BO + adiabatic, sinc-DVR)")
    save(fig, "fig09_H2_levels_cascades", "H₂：301 の振動回転準位（横軸 J、縦軸 (0,0) からの励起）と (1,3)・(14,1) からの E2 カスケード")


def fig10_h2_spectra():
    fig, axes = plt.subplots(3, 1, figsize=(7.5, 7.5), sharex=True)
    for ax, name in zip(axes, ("audit_h2_v1J3.json", "audit_h2_v3J5.json", "audit_h2_v14J1.json")):
        A = load(name)
        for kind in ("e2", "m1", "gw"):
            rr = [ln for ln in A["lines"] if ln["kind"] == kind and ln["n"] > 1e-9]
            if not rr:
                continue
            lam = [HC_EV_UM / ln["hw_eV"] for ln in rr]
            ax.vlines(lam, 1e-9, [ln["n"] for ln in rr], color=KIND_COLOR[kind], lw=0.8, alpha=0.7)
            ax.scatter(lam, [ln["n"] for ln in rr], s=10, color=KIND_COLOR[kind], label=KIND_NAME[kind])
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(1e-9, 3)
        ax.set_ylabel("photons per cascade")
        ax.set_title("start (v,J) = (%d,%d): %.3f eV emitted, %.2f E2 + %.3f M1 photons" % (*A["init"], sum(A["fluxes"].values()), A["photons"].get("e2", 0), A["photons"].get("m1", 0)), fontsize=9)
        ax.legend(fontsize=8, loc="upper right")
    axes[-1].set_xlabel("wavelength λ [μm]")
    save(fig, "fig10_H2_cascade_spectra", "H₂：カスケードの放出記録（光子数対波長）、開始 (1,3)・(3,5)・(14,1)")


def fig11_h2_validation():
    V = load("validation_h2.json")
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    ax = axes[0]
    pr = [p for p in V["pairs"] if p["Aq_cds"] > 0]
    dv = np.array([p["u"][0] - p["l"][0] for p in pr])
    xc = np.array([p["Aq_cds"] for p in pr]); ym = np.array([p["Aq_mine"] * (p["sigma_cds"] / p["sigma_mine"]) ** 5 for p in pr])
    sc = ax.scatter(xc, ym, c=np.clip(dv, 0, 6), s=5, cmap="viridis", alpha=0.7)
    ax.plot([1e-16, 1e-5], [1e-16, 1e-5], color="k", lw=0.6)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("A$_{E2}$ Roueff+2019 (CDS) [s$^{-1}$]"); ax.set_ylabel("A$_{E2}$ this work (σ⁵ rescaled to CDS σ)")
    cb = fig.colorbar(sc, ax=ax); cb.set_label("Δv")
    ax.set_title("%d E2 lines, median |rel. diff.| = %.1e" % (len(pr), np.median(np.abs(ym / xc - 1))), fontsize=9)
    ax = axes[1]
    rel = ym / xc - 1
    ax.hist(np.log10(np.abs(rel) + 1e-12), bins=50, color="C2", alpha=0.8)
    ax.set_xlabel("log₁₀ |A$_{this}$/A$_{CDS}$ − 1|"); ax.set_ylabel("transitions")
    ax.set_title("Θ(R): WSD 1998 (this work) vs PK (CDS)", fontsize=9)
    ax = axes[2]
    lv = [l for l in V["levels"] if l["E_cds_cm"] is not None]
    vs = np.array([l["v"] for l in lv]); dd = np.array([l["E_mine_cm"] - l["E_cds_cm"] for l in lv]); Js = np.array([l["J"] for l in lv])
    sc = ax.scatter(vs + 0.02 * Js, dd, c=Js, s=8, cmap="plasma")
    cb = fig.colorbar(sc, ax=ax); cb.set_label("J")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("vibrational quantum number v"); ax.set_ylabel("E(this, BO+adiabatic) − E(CDS, full) [cm$^{-1}$]")
    ax.set_title("levels: nonadiabatic + rel. + QED not included", fontsize=9)
    save(fig, "fig11_H2_validation", "H₂ の検証：A_E2 の CDS（Roueff+2019）との比較（4669 本）、相対差の分布、準位エネルギーの差（非断熱・相対論・QED 抜きの分）")


def fig12_h2_radiative_association():
    Rj = load("radiative_association_h2.json")
    CM = 219474.6313632; KB = 0.695034800
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    ax = axes[0]
    for J in range(0, 5):
        tab = Rj["tabs"][str(J)]
        E = np.array([r["E_h"] for r in tab]) * CM
        P = np.array([2 * np.pi * r["rho"] * 2.4188843265857e-17 * (r["a_e2"] + r["a_m1"]) for r in tab])
        ax.plot(E, P, "o-", ms=3, lw=0.8, label="J = %d" % J)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("collision energy E [cm$^{-1}$]"); ax.set_ylabel("capture probability per pass P$_J$(E)")
    ax.legend(fontsize=8); ax.set_title("P$_J$(E): box continuum → bound (E2+M1); (14,4) resonance at J = 4", fontsize=9)
    ax = axes[1]
    rows = Rj["sigma_rows"]
    E = np.array([r["E_cm"] for r in rows]); sig = np.array([r["sigma_cm2"] for r in rows]); kE = np.array([r["k_cm3s"] for r in rows])
    ax.plot(E, sig, "o-", color="C0", label="σ(E) [cm²]")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("E [cm$^{-1}$]"); ax.set_ylabel("σ [cm²]")
    ax2 = ax.twinx(); ax2.plot(E, kE, "s--", color="C3", label="k(E) = σv [cm³ s$^{-1}$]"); ax2.set_yscale("log"); ax2.set_ylabel("k(E) [cm³ s$^{-1}$]")
    ax2.spines["right"].set_visible(True)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels(); ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="center left")
    ax.set_title("H + H → H₂ + γ: σ(E), k(E) (singlet ¼ × spin weights)", fontsize=9)
    ax = axes[2]
    Ef = np.array(Rj["kE_fine"]["E_h"]) * CM; kf = np.array(Rj["kE_fine"]["k_cm3s"])
    ax.plot(Ef / KB, kf, color="C3", lw=1, label="k(E) vs E/k_B")
    Ts = [r["T"] for r in Rj["kT_rows"]]; kT = [r["k_cm3s"] for r in Rj["kT_rows"]]
    ax.plot(Ts, kT, "o", color="k", label="Maxwell–Boltzmann k(T)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("T or E/k$_B$ [K]"); ax.set_ylabel("rate coefficient [cm³ s$^{-1}$]"); ax.legend(fontsize=8)
    ax.set_title("k(T) = 2–4 × 10$^{-29}$ cm³ s$^{-1}$; spikes = shape resonances\n(quasi-bound (v,J) behind the centrifugal barrier, width not resolved by the box)", fontsize=8)
    save(fig, "fig12_H2_radiative_association", "H + H → H₂ + 光子：通過 1 回あたりの捕獲確率 P_J(E)、σ(E) と k(E)、熱平均 k(T)。k(E) の鋭いピークは遠心力障壁の内側の準束縛 (v,J)（形状共鳴）で、幅は箱の離散化では解像されない")


def fig13_h2_layer1():
    A = load("audit_h2_v1J3.json")
    ev = A["layer1"]["events"]; v0 = np.array(A["layer1"]["v0"]); vfit = np.array(A["layer1"]["vfit"]); beta = np.linalg.norm(v0) / EN.C
    N = np.array([e["n"] for e in ev]); d = np.array([e["delta"] for e in ev])
    model = np.array([EN.doppler_delta(n, vfit) for n in N])
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    ax = axes[0]
    ax.scatter(N[:, 2], d / beta, s=6, alpha=0.6, label="H₂ photons (%d, 200 molecules)" % len(ev))
    xx = np.linspace(-1, 1, 200)
    ax.plot(xx, [EN.doppler_delta(np.array([np.sqrt(1 - x * x), 0, x]), v0) / beta for x in xx], color="k", lw=1, label="exact Doppler, v⃗₀")
    ax.set_xlabel("n̂·v̂₀"); ax.set_ylabel("δ/β"); ax.legend(fontsize=8); ax.set_title("layer 1 record for H₂ cascade from (1,3)")
    ax = axes[1]
    ax.hist((d - model) * EN.C, bins=40, color="C2", alpha=0.8)
    ax.set_xlabel("residual (δ − model)·c [m/s]"); ax.set_ylabel("count")
    ax.set_title("|v⃗_fit − v⃗₀| = %.3f m/s (first-order %.0f m/s); recoil %.2f–%.2f m/s" % (np.linalg.norm(vfit - v0), np.linalg.norm(np.array(A["layer1"]["v_lin"]) - v0), min(e["dv"] for e in ev), max(e["dv"] for e in ev)))
    save(fig, "fig13_H2_layer1_doppler", "H₂ 層 1：Doppler 記録と v⃗ の復元（線幅 1e-22 で水素原子より高精度）")


def fig14_h2_layer3():
    A = load("audit_h2_v1J3.json"); B = load("audit_h2_v14J1.json")
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    for ax, Aj, ttl in ((axes[0], A, "(1,3) start, m = 0"), (axes[1], B, "(14,1) start, m = 0")):
        items = sorted(Aj["layer3"]["lines"].items(), key=lambda kv: -sum(kv[1].values()))[:7]
        names = [k.replace("->", "→") for k, _ in items]
        tot = [sum(v.values()) for _, v in items]
        n0 = [v.get("0", 0) / t for (_, v), t in zip(items, tot)]; np_ = [v.get("1", 0) / t for (_, v), t in zip(items, tot)]; nm = [v.get("-1", 0) / t for (_, v), t in zip(items, tot)]
        y = np.arange(len(items))
        ax.barh(y, n0, color="C0", label="Δm = 0"); ax.barh(y, np_, left=n0, color="C1", label="Δm = +1"); ax.barh(y, nm, left=np.array(n0) + np.array(np_), color="C2", label="Δm = −1")
        ax.set_yticks(y); ax.set_yticklabels(names, fontsize=7); ax.invert_yaxis(); ax.set_xlabel("fraction of photons by Δm")
        fm = Aj["layer3"]["final_m"]
        ax.set_title("%s\nfinal (0,1): p(m=−1,0,+1) = %s, alignment %+.3f" % (ttl, "/".join("%.3f" % fm[k] for k in ("-1", "0", "1")), Aj["layer3"]["align"]), fontsize=8)
        ax.legend(fontsize=7, loc="lower right")
    save(fig, "fig14_H2_layer3_polarization", "H₂ 層 3：初期 m = 0 からの各線の Δm 分配（偏光）と終状態 (0,1) の m 分布（背景は ⟨n̄⟩ ≈ 0 で緩和しない）")


def fig15_h2_gravity():
    A = load("audit_h2_v1J3.json"); st = A["states"]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.8))
    ax = axes[0]
    for J in (0, 1, 2, 3, 10, 20):
        pts = sorted([(s["label"][0], s["E_grav_cm"]) for s in st if s["J"] == J])
        ax.plot([p[0] for p in pts], [-p[1] for p in pts], "o-", ms=3, lw=0.8, label="J = %d" % J)
    ax.set_xlabel("v"); ax.set_ylabel("−⟨−G m_H²/R⟩ = G m_H²⟨1/R⟩ [cm$^{-1}$]"); ax.legend(fontsize=8)
    ax.set_title("inter-atomic gravity ⟨−Gm_H²/R⟩ per level (separate account)", fontsize=9)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax = axes[1]
    import system_h2 as S2                                  # GW の辺は光子数が 1e-35 で JSON の閾値以下なので、部品から直接取る
    st2, ch2, _, _ = S2.build()
    e2 = {(a, b): A for a, b, kind, k, A in ch2 if kind == "e2"}
    gw = {(a, b): A for a, b, kind, k, A in ch2 if kind == "gw"}
    keys = [k for k in e2 if k in gw and e2[k] > 0]
    hw = np.array([st2[k[0]]["E_eV"] - st2[k[1]]["E_eV"] for k in keys]); ratio = np.array([gw[k] / e2[k] for k in keys])
    dv = np.array([st2[k[0]]["label"][0] - st2[k[1]]["label"][0] for k in keys])
    sc = ax.scatter(hw, ratio, s=6, c=np.clip(dv, 0, 6), cmap="viridis", alpha=0.7)
    cb = fig.colorbar(sc, ax=ax); cb.set_label("Δv")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("ħω [eV]"); ax.set_ylabel("A$_{GW}$ / A$_{E2}$ per line")
    ax.set_title("GW/E2 per line = 4(α$_{Gp}$/α)|⟨q_m⟩|²/|⟨Θ⟩|² (%d lines)" % len(keys), fontsize=9)
    save(fig, "fig15_H2_gravity_accounting", "H₂：準位ごとの原子間重力 ⟨−Gm_H²/R⟩（別勘定）と、線ごとの重力波／E2 の率の比（⟨Θ⟩ が相殺で小さい Δv の大きい線で比が上がる）")
    jdump("fig15_H2_gravity.json", dict(hw_eV=hw, ratio=ratio, dv=dv))


def main():
    fig01_hydrogen_levels(); fig02_hydrogen_spectrum(); fig03_hydrogen_time()
    fig04_hydrogen_layer2(); fig05_hydrogen_layer1(); fig06_hydrogen_layer3()
    fig07_mott(); fig08_bremsstrahlung()
    fig09_h2_levels(); fig10_h2_spectra(); fig11_h2_validation(); fig12_h2_radiative_association()
    fig13_h2_layer1(); fig14_h2_layer3(); fig15_h2_gravity()
    lines = ["# 図の一覧（figures/）", "", "| 図 | 内容 |", "|---|---|"] + ["| %s.png | %s |" % (n, c) for n, c in INDEX]
    (HERE / "figures_index.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
