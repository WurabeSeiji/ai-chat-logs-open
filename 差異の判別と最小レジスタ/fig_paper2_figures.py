#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第二篇「全射影系における読出しの不在と零二乗和形式への移行」論文用の図（1枚2パネル）
  (a) 実数の全射影系: 同心球面族（葉層）＋葉上の位相座標——同じ光線は同じθ、θ=0は任意
  (b) 零二乗和形式: 錐 Σy²=0（断面図）——実数解は原点のみ、非自明解は複素化で開く
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_LEAF = "#1f77b4"
C_RAY = "#999999"
C_PT = "#d62728"
C_CONE = "#9467bd"
C_PHASE = "#ff7f0e"

fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

# ---------- (a) 全射影系 ----------
ax = axes[0]
t = np.linspace(0, 2 * np.pi, 400)
radii = [0.5, 1.0, 1.6, 2.3]
labels = [r"$R_0$", r"$R_1$", r"$R_2$", r"$R_3$"]
for R, lb in zip(radii, labels):
    ax.plot(R * np.cos(t), R * np.sin(t), color=C_LEAF, lw=1.6)
    ax.annotate(lb, xy=(R * np.cos(0.35), R * np.sin(0.35)),
                fontsize=9, color=C_LEAF)
for ang in np.linspace(0, 2 * np.pi, 12, endpoint=False):
    ax.plot([0, 2.55 * np.cos(ang)], [0, 2.55 * np.sin(ang)],
            color=C_RAY, lw=0.7, alpha=0.8, zorder=0)

theta_ref = 0.0
ax.plot([0, 2.55 * np.cos(theta_ref)], [0, 2.55 * np.sin(theta_ref)],
        color="black", lw=1.0, ls="--", alpha=0.8, zorder=1)
ax.text(1.05, -0.22, r"$\theta=0$ (arbitrary)", fontsize=8, color="black")

ang0 = 1.05
ax.plot([0, 2.55 * np.cos(ang0)], [0, 2.55 * np.sin(ang0)],
        color=C_PHASE, lw=1.8, alpha=0.95, zorder=2)
arc_r = 0.42
arc_t = np.linspace(theta_ref, ang0, 80)
ax.plot(arc_r * np.cos(arc_t), arc_r * np.sin(arc_t),
        color=C_PHASE, lw=1.4)
ax.text(0.37 * np.cos(ang0 / 2), 0.37 * np.sin(ang0 / 2) + 0.08,
        r"$\theta$", fontsize=11, color=C_PHASE)

for R in radii:
    ax.plot([R * np.cos(ang0)], [R * np.sin(ang0)], "o", color=C_PT, ms=6)
ax.annotate(r"$x(R,\theta)=R\,\mathbf{u}(\theta)$" + "\nsame ray = same " + r"$\theta$",
            xy=(radii[2] * np.cos(ang0), radii[2] * np.sin(ang0)),
            xytext=(-2.45, 2.05), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.9))
ax.annotate(r"$(R_0,\theta)\mapsto(R_1,\theta)$" + "\nonly " + r"$R$" + " scales",
            xy=(radii[1] * np.cos(ang0), radii[1] * np.sin(ang0)),
            xytext=(0.9, -2.25), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.9))
ax.plot([0], [0], "k.", ms=6)
ax.set_xlim(-2.7, 2.7)
ax.set_ylim(-2.7, 2.7)
ax.set_aspect("equal")
ax.set_title("(a) Bounded angular phase on each leaf\nsame ray = same θ; origin arbitrary")
ax.axis("off")

# ---------- (b) 零二乗和形式（錐の断面） ----------
ax = axes[1]
# 錐 x² + y² − r² = 0 の (x, r) 断面: r = ±x の二直線
xx = np.linspace(-2.4, 2.4, 200)
ax.plot(xx, np.abs(xx), color=C_CONE, lw=2.2, label=r"$\sum_k y_k^2 = 0$ (cone, section)")
ax.plot(xx, -np.abs(xx), color=C_CONE, lw=2.2)
ax.fill_between(xx, np.abs(xx), 2.6, color=C_CONE, alpha=0.06)
ax.fill_between(xx, -np.abs(xx), -2.6, color=C_CONE, alpha=0.06)
# 実数解 = 原点のみ
ax.plot([0], [0], "o", color=C_PT, ms=9, zorder=5)
ax.annotate("only real solution:\nthe origin ([2] App. C)",
            xy=(0, 0), xytext=(-2.3, -1.9), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.9))
# 旧形式（右辺 R²）＝葉 1 枚を対比で薄く
th = np.linspace(0, 2 * np.pi, 200)
ax.plot(1.5 * np.cos(th), 1.5 * np.sin(th), color=C_LEAF, lw=1.2, ls="--",
        label=r"$\sum_n x_n^2 = R^2$ (one leaf, real solutions)")
ax.annotate("all-axis symmetry forces\nthe right-hand side to zero",
            xy=(1.15, 1.15), xytext=(0.15, 2.15), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.9))
ax.set_xlim(-2.6, 2.6)
ax.set_ylim(-2.6, 2.6)
ax.set_aspect("equal")
ax.set_title("(b) Null quadratic form: real solutions vanish,\nnon-trivial solutions require complexification")
ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
ax.axis("off")

fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(f"fig_paper2_transition.{ext}", dpi=200, bbox_inches="tight")
print("done: fig_paper2_transition (png+svg)")
