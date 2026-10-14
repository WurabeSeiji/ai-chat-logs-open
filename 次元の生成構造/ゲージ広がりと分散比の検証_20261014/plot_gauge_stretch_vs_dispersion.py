# -*- coding: utf-8 -*-
"""図化：ゲージの広がり s と ω/k = sqrt(1+1/κ_N^2) の比較"""
import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"; plt.rcParams["svg.fonttype"] = "none"
R = json.load(open("results_gauge_stretch_vs_dispersion.json", encoding="utf-8"))["results"]
N = [r["N"] for r in R]
fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
for ax, key, title in zip(axs, ["S1", "S3a", "S3b"],
        ["円周 cos(nφ/2)（κ=N/2）", "S³ cos(nψ)（κ=N）", "S³ U_n/(n+1)（κ=N+1）"]):
    ax.plot(N, [r[key]["s_int"] for r in R], "o-", color="#c0392b", label="s（ピーク区間）")
    ax.plot(N, [r[key]["s_pt"] for r in R], "s--", color="#e67e22", label="s（ピーク点）")
    ax.plot(N, [r[key]["w_over_k"] for r in R], "^-", color="#1f77b4", label="ω/k = √(1+1/κ²)")
    ax.axhline(np.sqrt(2), color="0.6", lw=0.8, ls=":"); ax.text(N[-1], np.sqrt(2)+0.04, "√2", ha="right", color="0.4")
    ax.set_title(title, fontsize=11); ax.set_xlabel("最高次 N（全次数）"); ax.grid(alpha=0.3)
axs[0].set_ylabel("比"); axs[0].legend(fontsize=9)
fig.suptitle("ゲージの広がり s と ω/k の比較（R² = Σx_n²、本当の二乗）", fontsize=12)
fig.tight_layout(); fig.savefig("fig_gauge_stretch_vs_dispersion.png", dpi=160); fig.savefig("fig_gauge_stretch_vs_dispersion.svg")
print("ok")
