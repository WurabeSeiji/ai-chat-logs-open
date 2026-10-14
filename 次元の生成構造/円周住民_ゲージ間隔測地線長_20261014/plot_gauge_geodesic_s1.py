# -*- coding: utf-8 -*-
"""図化：ゲージ間隔に対する測地線長の不等間隔（円周・二次元）"""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["svg.fonttype"] = "none"

res = json.load(open("results_gauge_geodesic_s1.json", encoding="utf-8"))["results"]
dat = np.load("data_gauge_geodesic_s1.npz")
th = dat["theta"]
LAB = {"all_N2":"全次数 N=2","all_N4":"全次数 N=4","all_N18":"全次数 N=18",
       "odd_N17":"奇数次 N=17","even_N18":"偶数次 N=18"}
COL = {"all_N2":"#8c6bb1","all_N4":"#6a51a3","all_N18":"#c0392b",
       "odd_N17":"#1f77b4","even_N18":"#2ca02c"}
deg = lambda x: np.degrees(np.asarray(x))

# fig1: |R(θ)| とゲージ目盛
fig, axs = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
for ax, k in zip(axs, ["all_N18", "odd_N17", "even_N18"]):
    ax.plot(deg(th), dat[k], color=COL[k], lw=1.6)
    for m in res[k]["marks"]:
        ax.axvline(deg(m), color="0.6", lw=0.6, ls=":")
    ax.set_ylabel("|R(θ)| = R/i")
    ax.set_title(LAB[k] + "　点線＝ゲージ目盛（最高次成分の節）", fontsize=10, loc="left")
axs[-1].set_xlabel("位相 θ − θ0 [度]")
axs[-1].set_xlim(-180, 180)
fig.suptitle("図1　R² = Σx_n²（本当の二乗）から決まる振幅 R の位相依存", fontsize=12)
fig.tight_layout(); fig.savefig("fig1_R_vs_phase.png", dpi=160); fig.savefig("fig1_R_vs_phase.svg")

# fig2: 各ゲージ区間の測地線長比
order = ["all_N2", "all_N4", "all_N18", "odd_N17", "even_N18"]
fig, axs = plt.subplots(len(order), 1, figsize=(9, 11), sharex=True)
for ax, k in zip(axs, order):
    r = res[k]; w = 360.0/r["N"]
    ax.bar(deg(r["centers"]), r["ratio"], width=w*0.85, color=COL[k], alpha=0.8)
    ax.axhline(1.0, color="k", lw=0.7)
    ax.set_ylim(0.5, 1.35); ax.set_ylabel("Δs/平均Δs")
    ax.set_title("%s　最大/最小 = %.3f" % (LAB[k], r["max_over_min"]), fontsize=10, loc="left")
axs[-1].set_xlabel("ゲージ区間の中心位相 θ − θ0 [度]（区間幅は位相で等間隔 360°/N）")
axs[-1].set_xlim(-200, 260)
fig.suptitle("図2　同じ一目盛りに対する測地線長の比（比なので共通因子 i は消える）", fontsize=12)
fig.tight_layout(); fig.savefig("fig2_gauge_interval_lengths.png", dpi=160); fig.savefig("fig2_gauge_interval_lengths.svg")

# fig3: 住民の物差し（位相の等間隔目盛 vs 測地線長での位置）
fig, ax = plt.subplots(figsize=(10, 5))
for j, k in enumerate(order):
    r = res[k]; N = r["N"]
    s = np.concatenate([[0], np.cumsum(r["ds_over_i"])]); s = s/s[-1]
    p = np.linspace(0, 1, N+1)
    y0 = j*1.0
    ax.vlines(p, y0+0.05, y0+0.35, color="0.55", lw=1)
    ax.vlines(s, y0-0.35, y0-0.05, color=COL[k], lw=1.4)
    for a, b in zip(p, s):
        ax.plot([a, b], [y0+0.05, y0-0.05], color="0.8", lw=0.6)
    ax.text(-0.02, y0, LAB[k], ha="right", va="center", fontsize=10)
ax.set_xlim(-0.25, 1.02); ax.set_ylim(-0.6, len(order)-0.4)
ax.set_yticks([]); ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xticklabels(["一周の始点", "", "半周", "", "一周の終点"])
ax.set_xlabel("上：位相で等間隔のゲージ目盛　／　下：同じ目盛の測地線長での位置（全長で規格化）")
for sp in ["top", "right", "left"]: ax.spines[sp].set_visible(False)
ax.set_title("図3　住民の物差し：同じ一目盛りが場所によって違う長さになる", fontsize=12)
fig.tight_layout(); fig.savefig("fig3_resident_ruler.png", dpi=160); fig.savefig("fig3_resident_ruler.svg")
print("done")
