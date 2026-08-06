#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文図生成 v2: 仮説提案版の追加図（fig_p7-p10、保存済みJSONのみから描画）"""
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

# 図P7: 約数類定理（v10b）
d = json.loads((HERE / "pre_v10b_longcoupling_result_v1.json").read_text())
rows = [r for r in d["rows"] if r["settle"] == 4000 and r["m"] >= 1]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
colors = {1: "tab:red", 2: "tab:blue", 4: "tab:green"}
import math
for r in rows:
    g = math.gcd(r["m"], 16)
    for ax, key, lab in ((axes[0], "mass2", "質量²(補償)"), (axes[1], "S", "偏極 S"),
                          (axes[2], "retention", "帯保持率")):
        ax.plot(r["m"], r[key], "o", ms=9, color=colors[g])
for ax, key, lab in ((axes[0], "mass2", "質量²(補償)"), (axes[1], "S", "偏極 S"),
                      (axes[2], "retention", "帯保持率")):
    ax.set_xlabel("巻き数 m"); ax.set_ylabel(lab); ax.grid(alpha=0.3)
axes[0].set_title("奇数{1,3,5,7}が5桁一致")
handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=f"gcd(m,16)={g}")
           for g, c in colors.items()]
axes[2].legend(handles=handles, fontsize=8)
fig.suptitle("図P7 約数類定理: 海中の種の性質は gcd(m,16) のみに依存（settle=4000）")
fig.tight_layout(); fig.savefig(HERE / "fig_p7_divisor_class_v1.png", dpi=150); plt.close(fig)

# 図P8: Z₂位相量子化（v13b）
d = json.loads((HERE / "pre_covering_degree_result_v13b.json").read_text())
fig, ax = plt.subplots(figsize=(7, 4.2))
for name, col, mk in (("帯電census(D)", "tab:red", "o"), ("中性m=0束", "tab:blue", "s")):
    pk = d["cases"][name]["peaks"]
    ax.plot([p["k"] for p in pk], [p["phase"] / np.pi for p in pk], mk + "-",
            color=col, label=name, ms=7)
ax.axhline(0, color="k", lw=0.6, ls="--"); ax.axhline(1, color="k", lw=0.6, ls="--")
ax.axhline(-1, color="k", lw=0.6, ls="--")
ax.set_xlabel("観測回帰ピーク番号 k"); ax.set_ylabel("振幅位相 Φ/π")
ax.set_title("図P8 Z₂位相量子化: 帯電種は Φ∈{0,π} に厳密量子化（二重被覆）・中性束は連続")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p8_z2_quantization_v1.png", dpi=150); plt.close(fig)

# 図P9: 閉じ込め検定（v15）
d = json.loads((HERE / "pre_confinement_result_v15.json").read_text())
fig, ax = plt.subplots(figsize=(7, 4.2))
for name, col in (("quark型m=+2+海", "tab:red"), ("electron型m=+3+海", "tab:blue"),
                    ("quark型m=+2孤立", "tab:gray")):
    c = d["cases"][name]
    ax.plot(c["windows"], c["f_read"], lw=1.8, color=col, label=name)
ax.set_xlabel("衝突数 j"); ax.set_ylabel("可読パワー分率 f_read（m≡0 mod 3）")
ax.set_title("図P9 閉じ込め=mod3可読性: クォーク型は孤立で厳密非可読・海中でハドロン化分のみ可読")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(HERE / "fig_p9_confinement_v1.png", dpi=150); plt.close(fig)

# 図P10: ν行（v17b）＋スピン統計クロス表（v14）
d17 = json.loads((HERE / "pre_neutrino_row_result_v17b.json").read_text())
d14 = json.loads((HERE / "pre_spin_statistics_cross_result_v14.json").read_text())
fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
axes[0].plot(range(1, d17["n_peaks"] + 1), d17["phis_over_pi"], "o-",
             color="tab:purple", ms=7)
axes[0].axhline(0, color="k", lw=0.6, ls="--"); axes[0].axhline(1, color="k", lw=0.6, ls="--")
axes[0].axhline(-1, color="k", lw=0.6, ls="--")
axes[0].set_xlabel("回帰ピーク番号"); axes[0].set_ylabel("Φ/π")
axes[0].set_title(f"(a) ν候補（純m=0・F帯）: Qz2={d17['Qz2']:.2f}・被覆度2")
labels = []; vals = []
for key, r in d14["cells"].items():
    labels.append(key.replace("偶(F分類)", "F帯").replace("奇(B分類)", "B帯"))
    vals.append(2 if "2" in r["verdict"] else (1 if r["verdict"] == "被覆度1" else 0))
axes[1].bar(range(len(vals)), vals, color=["tab:red" if v == 2 else "tab:blue" for v in vals])
axes[1].set_xticks(range(len(labels))); axes[1].set_xticklabels(labels, fontsize=8)
axes[1].set_yticks([1, 2]); axes[1].set_yticklabels(["被覆度1", "被覆度2"])
axes[1].set_title("(b) クロス表: 被覆度はχパリティのみに依存")
fig.suptitle("図P10 ν行の成立とスピン統計対応")
fig.tight_layout(); fig.savefig(HERE / "fig_p10_nu_spinstat_v1.png", dpi=150); plt.close(fig)
print("fig_p7〜p10 生成完了")
