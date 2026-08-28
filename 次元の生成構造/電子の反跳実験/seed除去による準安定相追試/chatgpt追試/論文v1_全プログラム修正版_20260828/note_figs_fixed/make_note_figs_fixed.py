"""note 図の修正版再生成：原本 make_note_figs.py から図1・図3・図6 の節をそのまま抜き出し、入力パスだけ引数化。図2・図4（数値直書き）、図5（生成元不明の入力）、図7（幾何図）は対象外。"""
# -*- coding: utf-8 -*-
"""note 記事用の一般向け図（日本語ラベル）。データは N5_dynamics_followup zip の CSV。"""
import os, sys, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = sys.argv[2]  # 修正版 followup の CSV フォルダ（パス変更のみ）
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Hiragino Sans", "font.size": 12, "axes.titlesize": 15, "axes.titleweight": "bold",
                     "axes.labelsize": 13, "mathtext.fontset": "dejavusans", "axes.spines.top": False, "axes.spines.right": False})
RED, BLUE, GREEN, GRAY, NAVY = "#d7263d", "#1f5fd8", "#2e8b57", "#777777", "#12245e"

def read(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for k in rows[0]:
        vals = [r[k] for r in rows]
        try:
            out[k] = np.array([float(v) for v in vals])
        except ValueError:
            out[k] = vals
    return out

# ---------------------------------------------------------------- 図1 移送
d = read(os.path.join(HERE, "pump_depletion_timeseries.csv"))
fig, ax = plt.subplots(figsize=(10, 5.6), dpi=170)
ax.plot(d["step"], d["H_parallel"], color=BLUE, lw=2.4, label="もとの平面にあった成分 H∥")
ax.plot(d["step"], d["H_perp"], color=RED, lw=2.4, label="新しい方向へ移った成分 H⊥")
ax.plot(d["step"], d["H_total"], color=GREEN, lw=2.0, ls="--", label="合計 H∥ + H⊥（厳密に一定）")
ax.axvspan(0, 449, color="#fff1f1", zorder=0); ax.text(225, 0.5, "急拡大\n（〜step 449）", ha="center", color=RED, fontsize=12)
ax.text(2600, 0.72, "H⊥ = 0.678", color=RED, fontsize=13); ax.text(2600, 0.27, "H∥ = 0.322", color=BLUE, fontsize=13)
ax.set_xlim(0, 5000); ax.set_ylim(-0.02, 1.08); ax.set_xlabel("ステップ"); ax.set_ylabel("二乗振幅（全体を 1 とする）")
ax.set_title("急拡大の正体は「移送」──合計は増えていない")
ax.legend(loc="center right", frameon=False); fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig1_transfer.png")); plt.close(fig)

# ---------------------------------------------------------------- 図3 Floquet 単位円
fl = read(os.path.join(HERE, "floquet_spectrum.csv"))
m = fl["fd_eps"] == 3e-6
re_, im_, mod = fl["eig_re"][m], fl["eig_im"][m], fl["modulus"][m]
fig, ax = plt.subplots(figsize=(8.6, 7.6), dpi=170)
th = np.linspace(0, 2 * np.pi, 400); ax.plot(np.cos(th), np.sin(th), color=GRAY, lw=1.2, ls="--", label="単位円（ここに乗れば「そのまま」）")
ax.axhline(0, color="#cccccc", lw=0.8); ax.axvline(0, color="#cccccc", lw=0.8)
stable = mod < 1 - 1e-6; neutral = np.abs(mod - 1) <= 1e-6; unst = mod > 1 + 1e-6
ax.plot(re_[stable], im_[stable], "o", color=BLUE, ms=9, label="縮む方向（|μ| < 1）")
ax.plot(re_[neutral], im_[neutral], "o", color=GRAY, ms=9, mfc="white", label="変わらない方向（|μ| = 1）")
ax.plot(re_[unst], im_[unst], "o", color=RED, ms=13, label="伸びる方向（|μ| > 1）＝倒れる方向")
ax.annotate("μ₁ = 1.0901（2 重）\n最速で倒れる 2 次元", (1.0901, 0), xytext=(0.55, 0.55), fontsize=12, color=RED,
            arrowprops=dict(arrowstyle="->", color=RED))
ax.annotate("μ₂ = 1.0526（2 重）", (1.0526, 0), xytext=(0.35, -0.55), fontsize=11, color=RED, arrowprops=dict(arrowstyle="->", color=RED))
ax.set_aspect("equal"); ax.set_xlim(-1.25, 1.35); ax.set_ylim(-1.2, 1.2)
ax.set_xlabel("実部"); ax.set_ylabel("虚部"); ax.set_title("倒れかけの波の「倒れ方」を、波そのものから計算する")
ax.legend(loc="lower left", fontsize=10, frameon=False); fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig3_floquet_circle.png")); plt.close(fig)

# ---------------------------------------------------------------- 図6 seed 掃引
sd = read(os.path.join(HERE, "N5_moduli_seed_sweep.csv"))
fig, ax = plt.subplots(figsize=(10, 5.2), dpi=170)
ax.axhline(0, color=GRAY, lw=1, ls="--")
ax.bar(sd["seed"], sd["relative_phase_mod_pi_rad"], color=BLUE, width=0.55)
for s_, v in zip(sd["seed"], sd["relative_phase_mod_pi_rad"]):
    ax.text(s_, v + (0.012 if v >= 0 else -0.03), f"{v:+.3f}", ha="center", fontsize=10, color="#333")
ax.text(0.98, 0.95, "8 回とも：4 群（3+3+2+2）は同じ、振幅も全て 0.1\nしかし二つの距離族の「相対位相」だけは走行ごとに違う", transform=ax.transAxes, va="top", ha="right", fontsize=12,
        bbox=dict(boxstyle="round,pad=0.5", fc="#f5f6fc", ec=NAVY))
ax.set_xlabel("乱数の種（親状態を作る乱数）"); ax.set_ylabel("二距離族の相対位相 [rad]"); ax.set_ylim(-0.16, 0.42)
ax.set_title("決まっているものと、決まっていないもの"); fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig6_seed_phase.png")); plt.close(fig)

