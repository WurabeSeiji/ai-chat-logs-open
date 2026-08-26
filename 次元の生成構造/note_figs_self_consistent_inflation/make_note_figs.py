# -*- coding: utf-8 -*-
"""note 記事用の一般向け図（日本語ラベル）。データは N5_dynamics_followup zip の CSV。"""
import os, sys, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
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

# ---------------------------------------------------------------- 図2 onset vs residual
eps = np.array([3.8728081613e-7, 1.8154355031e-9, 5.0849521984e-11, 2.3845913288e-13]); on = np.array([72, 134, 176, 238])
x = -np.log(eps); a, b = 11.616225, -99.563139
fig, ax = plt.subplots(figsize=(10, 5.6), dpi=170)
xx = np.linspace(13, 30, 10); ax.plot(xx, a * xx + b, color=GRAY, lw=1.5, ls="--", label=f"直線あてはめ：傾き {a:.3f}（R² = 0.99999）")
ax.plot(x, on, "o", color=RED, ms=11, label="実測（4 本の走行）", zorder=3)
for xi, yi, e in zip(x, on, eps):
    ax.annotate(f"最初のズレ {e:.1e}\n→ {int(yi)} step で始動", (xi, yi), textcoords="offset points", xytext=(12, -28), fontsize=10, color="#333")
ax.text(13.5, 225, "予言：傾きは 1/ln μ₁ = 11.593\n（次の図の固有値から）", fontsize=12, color=NAVY,
        bbox=dict(boxstyle="round,pad=0.5", fc="#eef2ff", ec=NAVY))
ax.set_xlabel("最初のズレの小ささ  −ln ε  （右ほどズレが小さい）"); ax.set_ylabel("急拡大が始まるステップ")
ax.set_title("ズレが小さいほど、始まりは遅れる──しかも対数で正確に"); ax.set_xlim(13, 30); ax.set_ylim(50, 260)
ax.legend(loc="lower right", frameon=False); fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig2_onset_residual.png")); plt.close(fig)

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

# ---------------------------------------------------------------- 図4 三重整合
fig, ax = plt.subplots(figsize=(11, 6.2), dpi=170); ax.axis("off"); ax.set_xlim(0, 11); ax.set_ylim(0, 6.2)
def box(x, y, w, h, title, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15,rounding_size=0.25", fc="white", ec=color, lw=2.2))
    ax.text(x + w / 2, y + h - 0.45, title, ha="center", va="center", fontsize=13, weight="bold", color=color)
    ax.text(x + w / 2, y + h / 2 - 0.35, body, ha="center", va="center", fontsize=12, color="#222", linespacing=1.6)
box(0.3, 3.7, 3.2, 2.2, "① 時間発展をそのまま測る", "急拡大の成長率\n0.172513 / step", BLUE)
box(0.3, 0.4, 3.2, 2.2, "② 最初のズレを変えて測る", "始動時刻の傾き\n11.616", GREEN)
box(7.3, 2.0, 3.4, 2.4, "③ 倒れる速さを計算する", "線形化の固有値\nμ₁ = 1.090086569", RED)
ax.add_patch(FancyBboxPatch((4.2, 1.9), 2.6, 2.6, boxstyle="round,pad=0.15,rounding_size=0.3", fc="#fff8e1", ec=NAVY, lw=2.5))
ax.text(5.5, 3.85, "一つの数", ha="center", fontsize=14, weight="bold", color=NAVY)
ax.text(5.5, 3.05, "2 ln μ₁ = 0.172514", ha="center", fontsize=13, color=BLUE)
ax.text(5.5, 2.45, "1 / ln μ₁ = 11.593", ha="center", fontsize=13, color=GREEN)
for (x0, y0, x1, y1, c) in [(3.5, 4.8, 4.25, 3.6, BLUE), (3.5, 1.5, 4.25, 2.7, GREEN), (7.3, 3.2, 6.8, 3.2, RED)]:
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=22, color=c, lw=2.2))
ax.text(5.5, 0.55, "独立な三つの測定が、同じ固有値で説明される（一致 0.0006 % と 0.2 %）", ha="center", fontsize=12.5, color="#333")
ax.set_title("三重整合──始まりの機構は「倒れかけの澄んだ波の線形不安定性」", fontsize=15, weight="bold", color=NAVY)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig4_triple.png")); plt.close(fig)

# ---------------------------------------------------------------- 図5 エントロピー
e = read(os.path.join(HERE, "N5_spectral_entropy_timeseries.csv"))
fig, ax = plt.subplots(figsize=(10, 5.4), dpi=170)
ax.plot(e["step"], e["entropy_over_lnM"], color=NAVY, lw=2.4)
ax.axhline(1.0, color=GRAY, lw=1, ls="--"); ax.text(4950, 1.0008, "完全な均等分配 = 1", ha="right", fontsize=11, color=GRAY)
ax.annotate("急拡大の最中に\nいったん下がる（step 375）", (375, 0.97322), xytext=(1200, 0.978), fontsize=11, arrowprops=dict(arrowstyle="->", color="#333"))
ax.annotate("step 5000 で 1.000000", (5000, 1.0), xytext=(3300, 0.9905), fontsize=11, arrowprops=dict(arrowstyle="->", color="#333"))
ax.set_xlim(0, 5000); ax.set_ylim(0.972, 1.003); ax.set_xlabel("ステップ"); ax.set_ylabel("振幅のばらつきの指標  S / ln M")
ax.set_title("止まった後に起きること──10 本の波の振幅は完全に均等になる")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig5_entropy.png")); plt.close(fig)

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

# ---------------------------------------------------------------- 図7 K5 と四角錐（日本語）
CLASS = {"A+": ["12", "13", "45"], "A-": ["14", "15", "23"], "B+": ["24", "35"], "B-": ["25", "34"]}
EDGE = {e_: c for c, es in CLASS.items() for e_ in es}
STYLE = {"A+": (RED, "-"), "A-": (RED, "--"), "B+": (BLUE, "-"), "B-": (BLUE, "--")}
K5 = {"1": (0.0, 1.0), "2": (-0.95, 0.31), "3": (-0.59, -0.81), "4": (0.59, -0.81), "5": (0.95, 0.31)}
PYR = {"1": (0, 0, 1.25), "2": (-1, -1, 0), "4": (1, -1, 0), "3": (1, 1, 0), "5": (-1, 1, 0)}
def proj(p): x, y, z = p; return (x * 0.9 + 0.55 * y, z * 0.85 + 0.35 * y - 0.075 * x)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.8), dpi=170)
for e_, c in EDGE.items():
    col, ls = STYLE[c]; p, q = K5[e_[0]], K5[e_[1]]
    a1.plot([p[0], q[0]], [p[1], q[1]], color=col, ls=ls, lw=2.6)
    a1.text((p[0] + q[0]) / 2 * 1.1, (p[1] + q[1]) / 2 * 1.1, e_, color=col, fontsize=10, ha="center", va="center", bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))
    P, Q = proj(PYR[e_[0]]), proj(PYR[e_[1]]); diag = e_ in ("23", "45")
    a2.plot([P[0], Q[0]], [P[1], Q[1]], color=col, ls=ls, lw=2.6 if not diag else 1.8, alpha=1 if not diag else 0.7)
for v, (x, y) in K5.items():
    a1.add_patch(plt.Circle((x, y), 0.12, fc="white", ec="black", lw=1.4, zorder=4)); a1.text(x, y, v, ha="center", va="center", fontsize=12, zorder=5)
for v, pt in PYR.items():
    x, y = proj(pt); a2.add_patch(plt.Circle((x, y), 0.12, fc=RED if v == "1" else BLUE, ec="black", lw=0.8, zorder=4)); a2.text(x, y, v, ha="center", va="center", fontsize=11, color="white", weight="bold", zorder=5)
from matplotlib.lines import Line2D
h = [Line2D([0], [0], color=STYLE[c][0], ls=STYLE[c][1], lw=2.6, label=lbl) for c, lbl in [("A+", "A+（3 本）"), ("A-", "A−（3 本）＝ A+ の符号反転"), ("B+", "B+（2 本）"), ("B-", "B−（2 本）＝ B+ の符号反転")]]
a1.legend(handles=h, loc="lower center", fontsize=10, frameon=False, ncol=2, bbox_to_anchor=(0.5, -0.02))
for a in (a1, a2): a.set_aspect("equal"); a.axis("off")
a1.set_xlim(-1.3, 1.3); a1.set_ylim(-1.5, 1.25); a2.set_xlim(-1.9, 1.9); a2.set_ylim(-1.35, 1.6)
a1.set_title("5 体の 10 本の関係は 4 つの群に分かれる", fontsize=13); a2.set_title("同じ 10 本を四角錐として読む\n（8 本の外郭 ＋ 2 本の底面対角線）", fontsize=13)
fig.suptitle("N = 5 の落ち着いた形：3 + 3 + 2 + 2", fontsize=15, weight="bold", color=NAVY)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig7_k5_pyramid.png")); plt.close(fig)
print("done:", sorted(os.listdir(OUT)))
