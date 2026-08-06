#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""追加図5枚（fig_s8〜fig_s12）: 4D並進・質量/電荷/スピン/分類の分布

データ: translation_4d_readout_result_v1.json / distribution_readouts_result_v1.json
（いずれも対照検定つき・保存済み）
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = ["Hiragino Sans", "Arial Unicode MS", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False
HERE = Path(__file__).resolve().parent
tr = json.loads((HERE / "translation_4d_readout_result_v1.json").read_text())
di = json.loads((HERE / "distribution_readouts_result_v1.json").read_text())

# ============ fig_s8 4D並進（xyz＋τ） ============
ts = np.array([r["t"] for r in tr["rows"]])
U = np.array(tr["unwrapped"])
VEL = tr["VEL"]; x0 = tr["x0"]
fig = plt.figure(figsize=(11.5, 9))
ax3 = fig.add_subplot(2, 2, 1, projection="3d")
tt = np.linspace(0, ts[-1], 200)
pred = np.array([[x0[a] + VEL[a] * t for a in range(3)] for t in tt])
ax3.plot(pred[:, 0], pred[:, 1], pred[:, 2], "k--", lw=1, label="予言（等速直線）")
sc = ax3.scatter(U[:, 0], U[:, 1], U[:, 2], c=ts, cmap="viridis", s=55, zorder=5)
for i in range(len(ts)):
    if i % 3 == 0:
        ax3.text(U[i, 0], U[i, 1], U[i, 2] + 0.6, f"τ={ts[i]}", fontsize=7)
ax3.set_xlabel("x"); ax3.set_ylabel("y"); ax3.set_zlabel("z")
ax3.set_title("(a) xyz読出しの3次元軌道（等間隔の点=等速並進）\n"
              "τ=0〜900・10点全て偏差0.000・PR一定3.39", fontsize=10)
ax3.legend(fontsize=8)
labels = ["x", "y", "z"]
for a in range(3):
    ax = fig.add_subplot(2, 2, 2 + a)
    ax.plot(tt, x0[a] + VEL[a] * tt, "k--", lw=1, label=f"予言 傾き={VEL[a]}")
    ax.plot(ts, U[:, a], "o", color="tab:red", ms=7, label="実測")
    ax.set_xlabel("τ（固有時間 = R′時計の刻み, step）")
    ax.set_ylabel(f"{labels[a]}（展開表示, セル）")
    ax.set_title(f"({chr(98+a)}) {labels[a]}–τ 面への投影：直線＝等速並進", fontsize=10)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
fig.suptitle("xyz＋τ（R′固有時間）の4次元読出し：物質パケットの等速並進", fontsize=13)
fig.tight_layout()
fig.savefig(HERE / "fig_s8_translation4d_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s8 saved")

# ============ fig_s9 質量分布 ============
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9), sharey=True)
for ax, key, name in ((axes[0], "early", "成長期（j=400–500）"),
                       (axes[1], "late", "熱平衡期（j=3000–3200）")):
    mh = di[key]["mass_hist"]
    e = np.array(mh["edges"]); c = (e[:-1] + e[1:]) / 2; w = e[1] - e[0]
    ax.bar(c - w * 0.25, np.array(mh["ferm"]) / max(sum(mh["ferm"]), 1e-300), w * 0.3,
           color="tab:red", label="フェルミオン的（偶bin≥4）")
    ax.bar(c + w * 0.05, np.array(mh["bos"]) / max(sum(mh["bos"]), 1e-300), w * 0.3,
           color="tab:blue", label="ボゾン的（奇bin）")
    mc_ = di["charged"]["mass_hist"]
    ax.bar(c + w * 0.35, np.array(mc_["ferm"]) / max(sum(mc_["ferm"]), 1e-300), w * 0.3,
           color="tab:green", alpha=0.8, label="帯電構造（参考・両パネル同一）")
    ax.set_xlabel("質量²/T²（モード別 非コヒーレンス）")
    ax.set_title(name, fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
axes[0].set_ylabel("パワー重み割合")
fig.suptitle("質量として読み出せる量の分布：海のモードは両時期とも最大非コヒーレンス（質量²/T²≈1）に集中\n"
             "——海は「重い」。軽い（光的）状態は海には現れず、コヒーレントな構造に限られる（図4の固有モード 1.2e-7 参照）", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "fig_s9_mass_distribution_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s9 saved")

# ============ fig_s10 電荷分布 ============
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
ms_ = sorted(int(m) for m in di["early"]["charge_distribution"])
for key, name, col in (("early", "成長期", "tab:orange"), ("late", "熱平衡期", "tab:gray")):
    vals = [di[key]["charge_distribution"][str(m)] for m in ms_]
    tot = max(sum(vals), 1e-300)
    axes[0].plot(ms_, np.array(vals) / tot, "o-", color=col, label=name)
axes[0].set_xlabel("電荷（巻き数 m）"); axes[0].set_ylabel("パワー割合")
axes[0].set_title("(a) 対称な海：電荷は多数の値に分布\n（±1に安定しない——オープン課題）", fontsize=10)
axes[0].legend(fontsize=9); axes[0].grid(alpha=0.3)
vals_c = [max(di["charged"]["charge_distribution"][str(m)], 1e-12) for m in ms_]
axes[1].bar([str(m) for m in ms_], vals_c,
            color=["tab:red" if m == 1 else ("tab:green" if m == 3 else "tab:gray") for m in ms_])
axes[1].set_yscale("log"); axes[1].set_ylim(1e-12, 3)
axes[1].set_xlabel("電荷（巻き数 m）"); axes[1].set_ylabel("パワー（対数）")
axes[1].set_title("(b) 帯電した構造（census型）：種 m=1（赤）と\n和則相棒 m=3（緑）——複数の電荷値が共存", fontsize=10)
fig.suptitle("電荷として読み出せる量の分布——±1への安定化は未解決（オープン課題として明記）", fontsize=11)
fig.tight_layout()
fig.savefig(HERE / "fig_s10_charge_distribution_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s10 saved")

# ============ fig_s11 ボゾン/フェルミオン分布とスピン分布 ============
fig, axes = plt.subplots(2, 2, figsize=(11, 7.2))
for ax, key, name in ((axes[0, 0], "early", "成長期"), (axes[0, 1], "late", "熱平衡期")):
    kk = np.array(di[key]["k_signed"]); Pk = np.array(di[key]["P_k"])
    fm = np.array(di[key]["ferm_mask"]); bm = np.array(di[key]["bos_mask"])
    o = np.argsort(kk)
    ax.semilogy(kk[o][bm[o]], np.maximum(Pk[o][bm[o]], 1e-12), ".", color="tab:blue",
                ms=4, label=f"ボゾン的 総和 {di[key]['class_power']['bosonic']:.1f}")
    ax.semilogy(kk[o][fm[o]], np.maximum(Pk[o][fm[o]], 1e-12), ".", color="tab:red",
                ms=4, label=f"フェルミオン的 総和 {di[key]['class_power']['fermionic']:.1f}")
    ax.set_xlim(-80, 80)
    ax.set_xlabel("生bin k（符号つき）"); ax.set_ylabel("パワー P(k)")
    ax.set_title(f"(a{'１' if key=='early' else '２'}) {name}：クラス別パワースペクトル", fontsize=10)
    ax.legend(fontsize=8)
for ax, hkey, xlabel, ttl in ((axes[1, 0], "smag_hist", "|S|/T（コヒーレンス度）",
                                "(b) スピン的な量の大きさの分布"),
                               (axes[1, 1], "sz_hist", "s_z = Z/T",
                                "(c) スピンz成分の分布")):
    names = {"early": "成長", "late": "熱平衡", "charged": "帯電構造"}
    for key, ls, alpha in (("early", "-", 0.9), ("late", "--", 0.6), ("charged", ":", 0.9)):
        hh = di[key][hkey]
        e = np.array(hh["edges"]); c = (e[:-1] + e[1:]) / 2
        ax.plot(c, np.array(hh["ferm"]) / max(sum(hh["ferm"]), 1e-300), ls, color="tab:red",
                alpha=alpha, lw=2 if key == "charged" else 1.4,
                label=f"フェルミオン的（{names[key]}）")
        if key != "charged":
            ax.plot(c, np.array(hh["bos"]) / max(sum(hh["bos"]), 1e-300), ls, color="tab:blue",
                    alpha=alpha, label=f"ボゾン的（{names[key]}）")
    ax.set_xlabel(xlabel); ax.set_ylabel("パワー重み割合")
    ax.set_title(ttl + "（探索的読出し＝GramのBlochベクトル）\n海はほぼ無偏極（|S|≈0.03）・帯電種は約5倍偏極（|S|≈0.13）だが完全偏極には遠い", fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
fig.suptitle("ボゾン的な波とフェルミオン的な波の分布と、それぞれのスピン的な量の分布", fontsize=12)
fig.tight_layout()
fig.savefig(HERE / "fig_s11_class_spin_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s11 saved")

# ============ fig_s12 粒子的/反粒子的/灰色 ============
fig, axes = plt.subplots(1, 2, figsize=(10, 3.9))
for key, name, col in (("early", "対称な海・成長期", "tab:orange"),
                        ("late", "対称な海・熱平衡期", "tab:gray"),
                        ("charged", "帯電した構造", "tab:red")):
    bh = di[key]["balance_hist"]
    e = np.array(bh["edges"]); c = (e[:-1] + e[1:]) / 2
    h = np.array(bh["h"], float)
    axes[0].plot(c, h / max(h.sum(), 1e-300), "o-", color=col, label=name, ms=4)
axes[0].axvspan(-0.5, 0.5, color="#eeeeee", zorder=0)
axes[0].annotate("灰色領域", (0, axes[0].get_ylim()[1] * 0.02), ha="center", fontsize=9)
axes[0].set_xlabel("符号均衡 balance = (P₊−P₋)/(P₊+P₋)")
axes[0].set_ylabel("パワー重み割合")
axes[0].set_title("(a) フェルミオン的な波の粒子/反粒子バランス分布", fontsize=10)
axes[0].legend(fontsize=8); axes[0].grid(alpha=0.3)
cats = ["粒子的\n(b>0.5)", "灰色\n(|b|≤0.5)", "反粒子的\n(b<−0.5)"]
x = np.arange(3); w = 0.26
for i, (key, name, col) in enumerate((("early", "対称な海・成長期", "tab:orange"),
                                        ("late", "対称な海・熱平衡期", "tab:gray"),
                                        ("charged", "帯電した構造", "tab:red"))):
    fr = di[key]["fractions"]
    axes[1].bar(x + (i - 1) * w, [fr["particle"], fr["gray"], fr["antiparticle"]], w,
                color=col, label=name)
axes[1].set_xticks(x); axes[1].set_xticklabels(cats)
axes[1].set_ylabel("パワー割合")
axes[1].set_title("(b) 三分類の割合：対称な海は灰色、帯電種は粒子的", fontsize=10)
axes[1].legend(fontsize=8)
fig.suptitle("フェルミオン的な波の粒子的・反粒子的・灰色状態の分布", fontsize=12)
fig.tight_layout()
fig.savefig(HERE / "fig_s12_particle_gray_v1.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("fig_s12 saved")
print("ALL ADDENDA FIGURES DONE")
