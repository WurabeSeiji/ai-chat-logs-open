#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重力論文図生成（日英）: 保存済みJSONのみから描画（走行なし・決定論）
figG1 閉塞監査 / figG2 ω場 / figG3 質量因子化 / figG4 γ(t)=t/L /
figG5 二体普遍引力 / figG6 双線形分光 / figG7 巡回保存"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

JA_FONT = None
for f in ("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
          "/System/Library/Fonts/Hiragino Sans GB.ttc"):
    if Path(f).exists():
        font_manager.fontManager.addfont(f)
        JA_FONT = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams["axes.unicode_minus"] = False
HERE = Path(__file__).resolve().parent
BAT = HERE.parent / "周期表追試_局所θ_v1"

T = {
"ja": dict(
 g1_title="図G1  閉塞監査: 重力は状態側に入れられない・両エンジンは同一不変量を保存",
 g1_labels=["素の二体力学", "V(x)注入\n(標準理論式)", "ゲージ側操作", "N体 N=5\n(Cayley)", "N体 N=8\n(Cayley)"],
 g1_y="閉塞ドリフト（相対）",
 g2_title="図G2  時計場 ω(x): 実測 vs 読出し核からの解析予言（力学なし）",
 g2_x="質量中心からの距離 dx", g2_y="ω(x)", g2_l=["実測（T=200）", "解析予言"],
 g3_title="図G3  質量生成の因子化 m(f,v)=y(f)·g(v)（ヒッグス機構の構造）",
 g3a_x="海振幅 v（凝縮値）", g3a_y="g(v)=m(f,v)/m(f,0.2)",
 g3a_g="v^-0.92", g3b_x="解析予言 m_pred", g3b_y="動力学実測 m",
 g3b_note="15点: 比0.983・CV1.2%",
 g4_title="図G4  空間ゲージの線形成長: γ(t)=t/L（k̇=∇ω 運動量則）",
 g4_x="時刻 t（ステップ）", g4_y="γ = 〈|δτ_x|〉/〈|δτ_t|〉",
 g4_l=["実測", "予言 t/L（L=48.7・核幾何から）"],
 g5_title="図G5  広帯域海の二体時計結合: 全10ペア×全分離で引力",
 g5_x="ペア", g5_y="位相平均結合 Ē",
 g5_note="直交対比 1.32/1.13 >1 ⇒ σも重力荷（不変量方向）",
 g6_title="図G6  双線形署名: 時計場は場の差周波数のみを拾う（double copy）",
 g6a_t="(a) 場のスペクトル（線形量）", g6b_t="(b) 時計場のスペクトル",
 g6_x="角周波数 ω/step", g6_y="スペクトル強度",
 g6_l1="場の2線 ν₁,ν₂", g6_l2="差 |ν₂−ν₁|（時計主線・誤差0.0%）",
 g7_title="図G7  巡回巻き保存: 局所化の「漏れ」＝Nyquist折返し（大域系と同一署名）",
 g7a_y="Q_wind（整数持ち上げ）", g7b_y="η端パワー比", g7_x="衝突数 j",
 sfx="_ja"),
"en": dict(
 g1_title="Fig. G1  Closure audit: gravity cannot enter the state side; both engines preserve one invariant",
 g1_labels=["bare two-body", "V(x) injection\n(standard-theory)", "gauge-side op", "N-body N=5\n(Cayley)", "N-body N=8\n(Cayley)"],
 g1_y="closure drift (relative)",
 g2_title="Fig. G2  The clock field ω(x): measured vs analytic prediction from the readout kernel (no dynamics)",
 g2_x="distance from mass center dx", g2_y="ω(x)", g2_l=["measured (T=200)", "analytic prediction"],
 g3_title="Fig. G3  Factorization of mass generation m(f,v)=y(f)·g(v) (the Higgs-mechanism structure)",
 g3a_x="sea amplitude v (condensate value)", g3a_y="g(v)=m(f,v)/m(f,0.2)",
 g3a_g="v^-0.92", g3b_x="analytic prediction m_pred", g3b_y="measured m (dynamics)",
 g3b_note="15 points: ratio 0.983, CV 1.2%",
 g4_title="Fig. G4  Linear growth of the spatial gauge: γ(t)=t/L (the momentum law k̇=∇ω)",
 g4_x="time t (steps)", g4_y="γ = ⟨|δτ_x|⟩/⟨|δτ_t|⟩",
 g4_l=["measured", "prediction t/L (L=48.7 from kernel geometry)"],
 g5_title="Fig. G5  Two-body clock coupling in the broadband sea: attraction for all 10 pairs, all separations",
 g5_x="pair", g5_y="phase-averaged coupling Ē",
 g5_note="orthogonal-pair ratios 1.32/1.13 > 1 ⇒ σ also gravitates (invariant direction)",
 g6_title="Fig. G6  The bilinear signature: the clock field picks up only difference frequencies (double copy)",
 g6a_t="(a) spectrum of the field (linear quantity)", g6b_t="(b) spectrum of the clock field",
 g6_x="angular frequency ω/step", g6_y="spectral power",
 g6_l1="two field lines ν₁,ν₂", g6_l2="difference |ν₂−ν₁| (clock main line, 0.0% error)",
 g7_title="Fig. G7  Cyclic winding conservation: the localization \"leak\" = Nyquist fold-back (same signature as global)",
 g7a_y="Q_wind (integer lift)", g7b_y="η edge power fraction", g7_x="collision j",
 sfx="_en"),
}

d1 = json.loads((HERE / "result_g1_closure_and_engines_v1.json").read_text())
d3 = json.loads((HERE / "result_g3_higgs_factorization_v1.json").read_text())
d4 = json.loads((HERE / "result_g4_momentum_law_v1.json").read_text())
d5 = json.loads((HERE / "result_g5_twobody_broadband_v1.json").read_text())
d6 = json.loads((HERE / "result_g6_bilinear_spectroscopy_v1.json").read_text())
d7 = json.loads((BAT / "pre_aliasing_result_v8__smooth.json").read_text())

for lang, t in T.items():
    plt.rcParams["font.family"] = JA_FONT if lang == "ja" else "DejaVu Sans"

    # G1
    vals = [d1["bare"]["rel_drift"], d1["V_inject"]["rel_drift"], 1e-16,
            d1["Nbody_N5"]["zero_square_drift"], d1["Nbody_N8"]["zero_square_drift"]]
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    cols = ["tab:green", "tab:red", "tab:green", "tab:blue", "tab:blue"]
    ax.bar(range(5), vals, color=cols, alpha=0.85)
    ax.set_yscale("log"); ax.set_ylim(1e-17, 10)
    ax.set_xticks(range(5)); ax.set_xticklabels(t["g1_labels"], fontsize=8.5)
    ax.set_ylabel(t["g1_y"]); ax.set_title(t["g1_title"], fontsize=10.5)
    ax.axhline(1e-12, color="k", lw=0.7, ls=":")
    for i, v in enumerate(vals):
        ax.text(i, v * 2, f"{v:.0e}" if v > 1e-16 else "0", ha="center", fontsize=8.5)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout(); fig.savefig(HERE / f"fig_g1_closure{t['sfx']}_v1.png", dpi=150); plt.close(fig)

    # G2
    om_m = np.array(d4["omega_profile_meas"]); om_p = np.array(d4["omega_profile_pred"])
    n = len(om_m); c = d4["center"]
    xs = np.arange(n)
    dx = ((xs - c + n // 2) % n) - n // 2
    order = np.argsort(dx)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(dx[order], om_m[order], lw=1.6, color="tab:red", label=t["g2_l"][0])
    ax.plot(dx[order], om_p[order], lw=1.2, ls="--", color="k", label=t["g2_l"][1])
    ax.set_xlim(-100, 100)
    ax.set_xlabel(t["g2_x"]); ax.set_ylabel(t["g2_y"])
    ax.set_title(t["g2_title"], fontsize=10.5)
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(HERE / f"fig_g2_omega_field{t['sfx']}_v1.png", dpi=150); plt.close(fig)

    # G3
    M = np.array(d3["M_meas"]); P = np.array(d3["M_pred"])
    vs = d3["vs"]; fs = d3["fs"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    mks = ["o", "s", "^"]
    for i, f in enumerate(fs):
        axes[0].loglog(vs, M[i] / M[i][2], mks[i] + "-", ms=6, label=f"f={f}")
    gx = np.array([min(vs), max(vs)])
    axes[0].loglog(gx, (gx / 0.2) ** d3["g_exponent"], "k--", lw=1, label=t["g3a_g"])
    axes[0].set_xlabel(t["g3a_x"]); axes[0].set_ylabel(t["g3a_y"])
    axes[0].legend(fontsize=8.5); axes[0].grid(alpha=0.3, which="both")
    axes[1].plot(P.flatten(), M.flatten(), "o", ms=7, color="tab:red")
    lim = [0, max(M.max(), P.max()) * 1.1]
    axes[1].plot(lim, lim, "k--", lw=1)
    axes[1].set_xlabel(t["g3b_x"]); axes[1].set_ylabel(t["g3b_y"])
    axes[1].text(0.05, 0.9, t["g3b_note"], transform=axes[1].transAxes, fontsize=9,
                 bbox=dict(boxstyle="round", fc="#f2f2f2"))
    axes[1].grid(alpha=0.3)
    fig.suptitle(t["g3_title"], fontsize=11)
    fig.tight_layout(); fig.savefig(HERE / f"fig_g3_higgs{t['sfx']}_v1.png", dpi=150); plt.close(fig)

    # G4
    ts = sorted(int(s) for s in d4["gamma_t"])
    gs = [d4["gamma_t"][str(s)] for s in ts]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(ts, gs, "o", ms=9, color="tab:red", label=t["g4_l"][0])
    tt_ = np.array([0, max(ts) * 1.1])
    ax.plot(tt_, tt_ * d4["invL"], "k--", lw=1.2, label=t["g4_l"][1])
    ax.set_xlabel(t["g4_x"]); ax.set_ylabel(t["g4_y"])
    ax.set_title(t["g4_title"], fontsize=10.5)
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(HERE / f"fig_g4_momentum_law{t['sfx']}_v1.png", dpi=150); plt.close(fig)

    # G5
    keys = list(d5["pairs"].keys())
    Eb = [d5["pairs"][k]["E_bar"] for k in keys]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.bar(range(len(keys)), Eb, color="tab:blue", alpha=0.85)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(keys))); ax.set_xticklabels(keys, fontsize=8, rotation=30)
    ax.set_xlabel(t["g5_x"]); ax.set_ylabel(t["g5_y"])
    ax.set_title(t["g5_title"], fontsize=10.5)
    ax.text(0.02, 0.08, t["g5_note"], transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle="round", fc="#f2f2f2"))
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout(); fig.savefig(HERE / f"fig_g5_twobody{t['sfx']}_v1.png", dpi=150); plt.close(fig)

    # G6
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    frf = np.array(d6["spec_field"]["fr"]); Sf = np.array(d6["spec_field"]["S"])
    frc = np.array(d6["spec_clock"]["fr"]); Sc = np.array(d6["spec_clock"]["S"])
    axes[0].semilogy(frf, Sf + 1e-12, lw=1.0, color="tab:blue")
    for nu in d6["nu_field"]:
        axes[0].axvline(nu, color="tab:red", lw=1.0, ls="--")
    axes[0].set_title(t["g6a_t"], fontsize=10)
    axes[0].set_xlim(0, 0.5)
    axes[1].semilogy(frc, Sc + 1e-12, lw=1.0, color="tab:purple")
    axes[1].axvline(d6["delta_nu"], color="tab:red", lw=1.2, ls="--")
    for nu in d6["nu_field"]:
        axes[1].axvline(nu, color="gray", lw=0.8, ls=":")
    axes[1].set_title(t["g6b_t"], fontsize=10)
    axes[1].set_xlim(0, 0.5)
    for a_ in axes:
        a_.set_xlabel(t["g6_x"]); a_.set_ylabel(t["g6_y"]); a_.grid(alpha=0.3)
    handles = [plt.Line2D([], [], color="tab:red", ls="--", label=t["g6_l2"]),
               plt.Line2D([], [], color="gray", ls=":", label=t["g6_l1"])]
    axes[1].legend(handles=handles, fontsize=8)
    fig.suptitle(t["g6_title"], fontsize=11)
    fig.tight_layout(); fig.savefig(HERE / f"fig_g6_bilinear{t['sfx']}_v1.png", dpi=150); plt.close(fig)

    # G7
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    for name, col in (("D", "tab:green"), ("S1", "tab:red")):
        cc = d7["cases"][name]
        axes[0].plot(cc["j"], cc["Q_wind"], lw=1.5, color=col, label=name)
        axes[1].semilogy(cc["j"], np.maximum(cc["edge_frac"], 1e-20), lw=1.5, color=col,
                         label=f"{name} (corr={cc['corr_dQ_edge']:+.2f})")
    axes[0].set_ylabel(t["g7a_y"]); axes[1].set_ylabel(t["g7b_y"])
    for a_ in axes:
        a_.set_xlabel(t["g7_x"]); a_.legend(fontsize=9); a_.grid(alpha=0.3)
    fig.suptitle(t["g7_title"], fontsize=11)
    fig.tight_layout(); fig.savefig(HERE / f"fig_g7_cyclic{t['sfx']}_v1.png", dpi=150); plt.close(fig)

print("figs G1-G7 (ja/en) done")
