#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper figures, English edition: figs P1-P10 drawn ONLY from saved JSONs (no runs, deterministic).
Same data as make_paper_figures_periodic_v1/v2.py; labels translated for the EN edition of the paper."""
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
HERE = Path(__file__).resolve().parent

# Fig P1: clock universality
d = json.loads((HERE / "番地走査_v1" / "periodic_address_scan_result_v1.json").read_text())
Ns = [r["N"] for r in d["scan"]]
cl = [r["clock_over_pi72"] / 5 for r in d["scan"]]
d2 = json.loads((HERE / "census_longwindow_result_v2.json").read_text())
Ns2 = [r["N"] for r in d2["scan"]]
cl2 = [r["clock_over_pi72_step"] for r in d2["scan"]]
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.axhline(1.0, color="k", lw=0.8, ls="--", label="ω=π/72 (one clock cycle = 144 steps)")
ax.plot(Ns, cl, "o", ms=5, color="tab:blue", label="short window T=4000 (N=4..144)")
ax.plot(Ns2, cl2, "s", ms=6, color="tab:red", label="long window T=42000 (N=5..16, ±0.1%)")
ax.set_xscale("log"); ax.set_xlabel("N (relational wavenumber = effective energy)")
ax.set_ylabel("ω_clock / (π/72)")
ax.set_title("Fig. P1  Universality of the collective clock: ω=π/72/step across N=4–144")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p1_clock_universality_en_v1.png", dpi=150); plt.close(fig)

# Fig P2: stable species vs resonances
d1 = json.loads((HERE / "pre_mode_census_result_v1.json").read_text())
fig, ax = plt.subplots(figsize=(7, 4.2))
for dd, mk, col, lab in ((d1, "o", "tab:orange", "short window T=4000 (resonances visible)"),
                          (d2, "s", "tab:blue", "long window T=42000 (only stable species remain)")):
    xs, ys = [], []
    for r in dd["scan"]:
        for pl in r["planes"]:
            xs.append(r["N"]); ys.append(pl["rho"])
    ax.plot(xs, ys, mk, ms=4, alpha=0.6, color=col, label=lab)
ax.axhline(1.0, color="k", lw=0.8, ls="--")
ax.set_xlabel("N"); ax.set_ylabel("rotation number ρ = f_mode / f_clock")
ax.set_title("Fig. P2  Short-window sidebands (resonances) contract to 1/1 in the long window:\nthe only stable species is the massless ground species")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p2_stable_vs_resonance_en_v1.png", dpi=150); plt.close(fig)

# Fig P3: charged lifetime and walk
d3 = json.loads((HERE / "pre_charged_stability_result_v3.json").read_text())
w = d3["windows"]; qs = d3["q_series"]
fig, ax = plt.subplots(figsize=(7, 4.2))
for tag, col in (("+1", "tab:red"), ("+3", "tab:purple"), ("0", "tab:gray"), ("-1", "tab:blue")):
    ax.plot(w, qs[tag], lw=1.5, color=col, label=f"winding q={tag}")
ax.set_xlabel("collision j"); ax.set_ylabel("fermionic power weight")
ax.set_title("Fig. P3  Metastability of charged species: q=+1 holds for τ≈13,400,\nthen sum-rule walk to the partner +3")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p3_charged_lifetime_walk_en_v1.png", dpi=150); plt.close(fig)

# Fig P4: readout rectification
d5 = json.loads((HERE / "pre_readout_rectification_result_v5.json").read_text())
fig, ax = plt.subplots(figsize=(7, 4.2))
Js = [3, 4, 5, 6]; width = 0.35
for i, a in enumerate(d5["analyses"][:2]):
    vals = [a["folds"][str(J)]["q1_concentration_last"] or 0.0 for J in Js]
    ax.bar(np.arange(len(Js)) + (i - 0.5) * width, vals, width,
           label=a["label"].replace("_", " ").replace("海", "sea").replace("種", " species"),
           color=["tab:red", "tab:blue"][i], alpha=0.85)
ax.set_xticks(range(len(Js))); ax.set_xticklabels([f"J={J}" for J in Js])
ax.set_ylabel("|q|=1 concentration (final window, charged class)")
ax.set_title("Fig. P4  Readout rectification: only the divide-by-3 observation clock\nreads all charged content at |q|=1")
ax.axhline(1.0, color="k", lw=0.8, ls="--")
ax.legend(fontsize=9); ax.grid(alpha=0.3, axis="y")
fig.tight_layout(); fig.savefig(HERE / "fig_p4_rectification_en_v1.png", dpi=150); plt.close(fig)

# Fig P5: ledger identity
d6 = json.loads((HERE / "pre_signed_charge_result_v6.json").read_text())
c = d6["cases"]["D_v3orig"]
Q3 = np.array(c["Q3_series"]); Qw = np.array(c["Q_wind_series"])
W = (Qw - Q3) / 3
jw = np.arange(1, len(Q3) + 1) * d6["J_WIN"]
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(jw, Q3, lw=2, color="tab:red", label="readable charge Q3 (mod-3 folded)")
ax.plot(jw, 3 * (W - W[0]) * (-1) + Q3[0], lw=1.2, ls="--", color="k",
        label="Q3(0) − 3ΔW (ledger prediction)")
ax2 = ax.twinx()
ax2.plot(jw, Qw, lw=1, color="tab:green", alpha=0.7, label="Q_wind (exactly conserved)")
ax2.set_ylabel("Q_wind", color="tab:green")
ax.set_xlabel("collision j"); ax.set_ylabel("Q3")
ax.set_title("Fig. P5  Ledger identity ΔQ3=−3ΔW (accuracy 7e-10):\nreadable charge is carried into neutral composites")
ax.legend(fontsize=9, loc="upper right"); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p5_ledger_en_v1.png", dpi=150); plt.close(fig)

# Fig P6: cyclic conservation and Nyquist wrap
d8 = json.loads((HERE / "pre_aliasing_result_v8.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
for name, col in (("D", "tab:green"), ("S1", "tab:red")):
    cc = d8["cases"][name]
    axes[0].plot(cc["j"], cc["Q_wind"], lw=1.5, color=col, label=name)
    axes[1].semilogy(cc["j"], np.maximum(cc["edge_frac"], 1e-20), lw=1.5, color=col,
                     label=f"{name} (corr={cc['corr_dQ_edge']:+.2f})")
axes[0].set_xlabel("collision j"); axes[0].set_ylabel("Q_wind (integer lift)")
axes[0].set_title("(a) integer winding charge: D = exactly conserved in band")
axes[1].set_xlabel("collision j"); axes[1].set_ylabel("η edge power fraction (|m| ≥ ne/2−4)")
axes[1].set_title("(b) accumulation at the Nyquist edge correlates with breaking")
for a_ in axes: a_.legend(fontsize=9); a_.grid(alpha=0.3)
fig.suptitle("Fig. P6  Winding conservation is cyclic mod ne: apparent breaking = register wrap-around", y=1.00)
fig.tight_layout(); fig.savefig(HERE / "fig_p6_cyclic_conservation_en_v1.png", dpi=150); plt.close(fig)

# Fig P7: divisor-class theorem
d = json.loads((HERE / "pre_v10b_longcoupling_result_v1.json").read_text())
rows = [r for r in d["rows"] if r["settle"] == 4000 and r["m"] >= 1]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
colors = {1: "tab:red", 2: "tab:blue", 4: "tab:green"}
for r in rows:
    g = math.gcd(r["m"], 16)
    for ax, key, lab in ((axes[0], "mass2", "compensated mass²"), (axes[1], "S", "polarization S"),
                          (axes[2], "retention", "band retention")):
        ax.plot(r["m"], r[key], "o", ms=9, color=colors[g])
for ax, key, lab in ((axes[0], "mass2", "compensated mass²"), (axes[1], "S", "polarization S"),
                      (axes[2], "retention", "band retention")):
    ax.set_xlabel("winding m"); ax.set_ylabel(lab); ax.grid(alpha=0.3)
axes[0].set_title("odd {1,3,5,7} match to 5 digits")
handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=f"gcd(m,16)={g}")
           for g, c in colors.items()]
axes[2].legend(handles=handles, fontsize=8)
fig.suptitle("Fig. P7  Divisor-class theorem: in-sea properties depend only on gcd(m,16) (settle=4000)")
fig.tight_layout(); fig.savefig(HERE / "fig_p7_divisor_class_en_v1.png", dpi=150); plt.close(fig)

# Fig P8: Z2 phase quantization
d = json.loads((HERE / "pre_covering_degree_result_v13b.json").read_text())
fig, ax = plt.subplots(figsize=(7, 4.2))
for name, lab, col, mk in (("帯電census(D)", "charged census (D)", "tab:red", "o"),
                            ("中性m=0束", "neutral m=0 bundle", "tab:blue", "s")):
    pk = d["cases"][name]["peaks"]
    ax.plot([p["k"] for p in pk], [p["phase"] / np.pi for p in pk], mk + "-",
            color=col, label=lab, ms=7)
ax.axhline(0, color="k", lw=0.6, ls="--"); ax.axhline(1, color="k", lw=0.6, ls="--")
ax.axhline(-1, color="k", lw=0.6, ls="--")
ax.set_xlabel("observable-recurrence peak k"); ax.set_ylabel("amplitude phase Φ/π")
ax.set_title("Fig. P8  Z₂ phase quantization: charged species locked to Φ∈{0,π}\n(double cover); the neutral bundle drifts")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p8_z2_quantization_en_v1.png", dpi=150); plt.close(fig)

# Fig P9: confinement
d = json.loads((HERE / "pre_confinement_result_v15.json").read_text())
fig, ax = plt.subplots(figsize=(7, 4.2))
for name, lab, col in (("quark型m=+2+海", "quark type m=+2 + sea", "tab:red"),
                        ("electron型m=+3+海", "electron type m=+3 + sea", "tab:blue"),
                        ("quark型m=+2孤立", "quark type m=+2 isolated", "tab:gray")):
    c = d["cases"][name]
    ax.plot(c["windows"], c["f_read"], lw=1.8, color=col, label=lab)
ax.set_xlabel("collision j"); ax.set_ylabel("readable power fraction f_read (m≡0 mod 3)")
ax.set_title("Fig. P9  Confinement = mod-3 readability: quark type strictly unreadable in isolation;\nin the sea, readable only through hadronization")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p9_confinement_en_v1.png", dpi=150); plt.close(fig)

# Fig P10: neutrino row + cross table
d17 = json.loads((HERE / "pre_neutrino_row_result_v17b.json").read_text())
d14 = json.loads((HERE / "pre_spin_statistics_cross_result_v14.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
axes[0].plot(range(1, d17["n_peaks"] + 1), d17["phis_over_pi"], "o-",
             color="tab:purple", ms=7)
axes[0].axhline(0, color="k", lw=0.6, ls="--"); axes[0].axhline(1, color="k", lw=0.6, ls="--")
axes[0].axhline(-1, color="k", lw=0.6, ls="--")
axes[0].set_xlabel("recurrence peak"); axes[0].set_ylabel("Φ/π")
axes[0].set_title(f"(a) ν candidate (pure m=0, F band): Qz2={d17['Qz2']:.2f}, cover 2")
labels = []; vals = []
for key, r in d14["cells"].items():
    labels.append(key.replace("偶(F分類)", "F band").replace("奇(B分類)", "B band").replace("|", "\n"))
    vals.append(2 if "2" in r["verdict"] else (1 if r["verdict"] == "被覆度1" else 0))
axes[1].bar(range(len(vals)), vals, color=["tab:red" if v == 2 else "tab:blue" for v in vals])
axes[1].set_xticks(range(len(labels))); axes[1].set_xticklabels(labels, fontsize=8)
axes[1].set_yticks([1, 2]); axes[1].set_yticklabels(["cover 1", "cover 2"])
axes[1].set_title("(b) cross table: covering degree depends only on χ parity")
fig.suptitle("Fig. P10  The neutrino row and the spin-statistics correspondence")
fig.tight_layout(); fig.savefig(HERE / "fig_p10_nu_spinstat_en_v1.png", dpi=150); plt.close(fig)
print("EN figures P1-P10 done")
