#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文図生成 v1: 保存済みJSONのみから描画（走行なし・決定論）"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for f in ("/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
          "/System/Library/Fonts/Hiragino Sans GB.ttc"):
    if Path(f).exists():
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name()
        break
plt.rcParams["axes.unicode_minus"] = False
HERE = Path(__file__).resolve().parent

# 図P1: 時計普遍性
d = json.loads((HERE / "番地走査_v1" / "periodic_address_scan_result_v1.json").read_text())
Ns = [r["N"] for r in d["scan"]]
cl = [r["clock_over_pi72"] / 5 for r in d["scan"]]
d2 = json.loads((HERE / "census_longwindow_result_v2.json").read_text())
Ns2 = [r["N"] for r in d2["scan"]]
cl2 = [r["clock_over_pi72_step"] for r in d2["scan"]]
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.axhline(1.0, color="k", lw=0.8, ls="--", label="ω=π/72（時計一周=144ステップ）")
ax.plot(Ns, cl, "o", ms=5, color="tab:blue", label="短窓 T=4000（N=4..144）")
ax.plot(Ns2, cl2, "s", ms=6, color="tab:red", label="長窓 T=42000（N=5..16、±0.1%）")
ax.set_xscale("log"); ax.set_xlabel("N（関係波数＝実効エネルギー）")
ax.set_ylabel("ω_clock / (π/72)")
ax.set_title("図P1  集団時計の普遍性：ω=π/72/step は N=4〜144 の全域で成立")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p1_clock_universality_v1.png", dpi=150); plt.close(fig)

# 図P2: 安定種と共鳴（v1短窓ρ vs v2長窓ρ）
d1 = json.loads((HERE / "pre_mode_census_result_v1.json").read_text())
fig, ax = plt.subplots(figsize=(7, 4.2))
for dd, mk, col, lab in ((d1, "o", "tab:orange", "短窓 T=4000（共鳴が見える）"),
                          (d2, "s", "tab:blue", "長窓 T=42000（安定種のみ残る）")):
    xs, ys = [], []
    for r in dd["scan"]:
        for pl in r["planes"]:
            xs.append(r["N"]); ys.append(pl["rho"])
    ax.plot(xs, ys, mk, ms=4, alpha=0.6, color=col, label=lab)
ax.axhline(1.0, color="k", lw=0.8, ls="--")
ax.set_xlabel("N"); ax.set_ylabel("回転数 ρ = f_モード / f_時計")
ax.set_title("図P2  短窓の側帯（共鳴）は長窓で 1/1 へ収縮：安定種は無質量基底種のみ")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p2_stable_vs_resonance_v1.png", dpi=150); plt.close(fig)

# 図P3: 帯電種の寿命とウォーク（v3）
d3 = json.loads((HERE / "pre_charged_stability_result_v3.json").read_text())
w = d3["windows"]; qs = d3["q_series"]
fig, ax = plt.subplots(figsize=(7, 4.2))
for tag, col in (("+1", "tab:red"), ("+3", "tab:purple"), ("0", "tab:gray"), ("-1", "tab:blue")):
    ax.plot(w, qs[tag], lw=1.5, color=col, label=f"巻き数 q={tag}")
ax.set_xlabel("衝突数 j"); ax.set_ylabel("フェルミオン的パワー重み")
ax.set_title(f"図P3  帯電種の準安定性：q=+1 は τ≈13,400 で保持、相棒 +3 へ和則ウォーク")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p3_charged_lifetime_walk_v1.png", dpi=150); plt.close(fig)

# 図P4: 読出し整流（v5）
d5 = json.loads((HERE / "pre_readout_rectification_result_v5.json").read_text())
fig, ax = plt.subplots(figsize=(7, 4.2))
Js = [3, 4, 5, 6]; width = 0.35
for i, a in enumerate(d5["analyses"][:2]):
    vals = [a["folds"][str(J)]["q1_concentration_last"] or 0.0 for J in Js]
    ax.bar(np.arange(len(Js)) + (i - 0.5) * width, vals, width,
           label=a["label"].replace("_", " "),
           color=["tab:red", "tab:blue"][i], alpha=0.85)
ax.set_xticks(range(len(Js))); ax.set_xticklabels([f"J={J}" for J in Js])
ax.set_ylabel("|q|=1 集中度（最終窓・荷電類）")
ax.set_title("図P4  読出し整流：分母3の観測時計だけが全荷電内容を |q|=1 に読む")
ax.axhline(1.0, color="k", lw=0.8, ls="--")
ax.legend(fontsize=9); ax.grid(alpha=0.3, axis="y")
fig.tight_layout(); fig.savefig(HERE / "fig_p4_rectification_v1.png", dpi=150); plt.close(fig)

# 図P5: 簿記恒等式（v6 D系列）
d6 = json.loads((HERE / "pre_signed_charge_result_v6.json").read_text())
c = d6["cases"]["D_v3orig"]
Q3 = np.array(c["Q3_series"]); Qw = np.array(c["Q_wind_series"])
W = (Qw - Q3) / 3
jw = np.arange(1, len(Q3) + 1) * d6["J_WIN"]
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(jw, Q3, lw=2, color="tab:red", label="読める電荷 Q3（mod3折返し）")
ax.plot(jw, 3 * (W - W[0]) * (-1) + Q3[0], lw=1.2, ls="--", color="k",
        label="Q3(0) − 3ΔW（簿記の予言）")
ax2 = ax.twinx()
ax2.plot(jw, Qw, lw=1, color="tab:green", alpha=0.7, label="Q_wind（厳密保存）")
ax2.set_ylabel("Q_wind", color="tab:green")
ax.set_xlabel("衝突数 j"); ax.set_ylabel("Q3")
ax.set_title("図P5  簿記恒等式 ΔQ3=−3ΔW（精度 7e-10）：読める電荷は中性複合へ持ち込まれる")
ax.legend(fontsize=9, loc="upper right"); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p5_ledger_v1.png", dpi=150); plt.close(fig)

# 図P6: 巡回保存とNyquist折返し（v8）
d8 = json.loads((HERE / "pre_aliasing_result_v8.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
for name, col in (("D", "tab:green"), ("S1", "tab:red")):
    cc = d8["cases"][name]
    axes[0].plot(cc["j"], cc["Q_wind"], lw=1.5, color=col, label=name)
    axes[1].semilogy(cc["j"], np.maximum(cc["edge_frac"], 1e-20), lw=1.5, color=col,
                     label=f"{name}（corr={cc['corr_dQ_edge']:+.2f}）")
axes[0].set_xlabel("衝突数 j"); axes[0].set_ylabel("Q_wind（整数持ち上げ）")
axes[0].set_title("(a) 整数巻き数電荷：D=帯域内で厳密保存")
axes[1].set_xlabel("衝突数 j"); axes[1].set_ylabel("η端パワー比（|m| >= ne/2-4）")
axes[1].set_title("(b) Nyquist端への蓄積が破れと相関")
for a_ in axes: a_.legend(fontsize=9); a_.grid(alpha=0.3)
fig.suptitle("図P6  巻き数保存は mod ne の巡回保存：見かけの破れ＝レジスタ折返し", y=1.00)
fig.tight_layout(); fig.savefig(HERE / "fig_p6_cyclic_conservation_v1.png", dpi=150); plt.close(fig)

print("図6枚生成完了: fig_p1〜p6")
