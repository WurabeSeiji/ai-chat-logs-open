#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第一篇「共役複素対の二つの積と射影の階梯」論文用の図（3枚）
  fig_paper1_ladder   : 図1 射影の階梯（螺旋→余弦波→直線）＋端面図
  fig_paper1_oddharm  : 図2 奇数倍音の脈動螺旋と局在波
  fig_paper1_fiber    : 図3 読出しのファイバー階層（Z2 / Z / S1）
PNG と SVG を両方出力。ラベルは英語（英語版と共用）。
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "."

C_Z = "#1f77b4"      # Z helix
C_ZBAR = "#ff7f0e"   # conjugate mirror
C_SUM = "#2ca02c"    # partial projection (cosine)
C_NORM = "#d62728"   # zero projection (line)
C_W = "#9467bd"      # pulsating spiral
C_LOC = "#e6a817"    # localized wave

# ---------------- 図1 ----------------

def fig1():
    fig = plt.figure(figsize=(11, 4.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    th = np.linspace(0, 4 * np.pi, 800)
    # 軸配置: x=θ, y=Im(奥行き), z=Re(縦) —— 読める量(実数値)を縦に立てる
    ax.plot(th, np.sin(th), np.cos(th), color=C_Z, lw=1.8,
            label=r"$Z=e^{i\theta}$  (identity: helix)")
    ax.plot(th, -np.sin(th), np.cos(th), color=C_ZBAR, lw=1.2, alpha=0.85,
            label=r"$\bar{Z}=e^{-i\theta}$  (mirror helix)")
    ax.plot(th, 0 * th, 2 * np.cos(th), color=C_SUM, lw=2.4,
            label=r"$Z+\bar{Z}=2\cos\theta$  (partial projection)")
    ax.plot(th, 0 * th, 1 + 0 * th, color=C_NORM, lw=2.4,
            label=r"$Z\bar{Z}=1$  (zero projection)")
    # Im=0 の鉛直面をうっすら
    Xp, Zp = np.meshgrid(np.linspace(0, 4 * np.pi, 2), np.linspace(-2.2, 2.2, 2))
    ax.plot_surface(Xp, 0 * Xp, Zp, alpha=0.06, color="gray")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel("Im")
    ax.set_zlabel("Re")
    ax.set_yticks([-1, 0, 1])
    ax.tick_params(axis="y", labelsize=7)
    ax.set_title("(a) Ladder of projections: helix / cosine / line")
    ax.view_init(elev=16, azim=-72)
    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.9)
    ax.set_box_aspect((2.6, 1, 1))

    ax2 = fig.add_subplot(1, 2, 2)
    t = np.linspace(0, 2 * np.pi, 400)
    ax2.plot(np.cos(t), np.sin(t), color=C_Z, lw=2.0,
             label=r"end view of helix: circle $a^2+b^2=\rho^2$")
    ax2.plot([-2, 2], [0, 0], color=C_SUM, lw=2.6, alpha=0.9,
             label=r"shadow of $Z+\bar{Z}$: segment $[-2\rho,\,2\rho]$")
    ax2.plot([1], [0], "o", color=C_NORM, ms=9,
             label=r"$Z\bar{Z}=\rho^2$: a single point")
    for ang in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        ax2.plot([np.cos(ang)], [np.sin(ang)], ".", color=C_Z, ms=4)
    ax2.annotate("all $\\theta$ are solutions",
                 xy=(np.cos(2.2), np.sin(2.2)), xytext=(-2.05, 1.55),
                 fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=0.8))
    ax2.set_xlim(-2.3, 2.3)
    ax2.set_ylim(-1.8, 1.8)
    ax2.set_aspect("equal")
    ax2.set_xlabel("Re")
    ax2.set_ylabel("Im")
    ax2.set_title("(b) End view along $\\theta$: the Pythagorean circle")
    ax2.legend(loc="lower left", fontsize=7.5, framealpha=0.9)
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{OUT}/fig_paper1_ladder.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ---------------- 図2 ----------------

def fig2():
    fig = plt.figure(figsize=(11, 4.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    th = np.linspace(0, 2 * np.pi, 1200)
    W = (1 + 2 * np.cos(2 * th)) * np.exp(1j * 3 * th)
    ax.plot(th, W.imag, W.real, color=C_W, lw=1.8,
            label=r"$W=(1+2\cos 2\theta)\,e^{i3\theta}$")
    r_env = np.abs(1 + 2 * np.cos(2 * th))
    ax.plot(th, 0 * th, r_env, color="gray", lw=0.8, ls="--", alpha=0.7,
            label=r"radius $|1+2\cos 2\theta|$")
    nodes = [np.arccos(-0.5) / 2, (2 * np.pi - np.arccos(-0.5)) / 2]
    for nd in nodes:
        ax.scatter([nd], [0], [0], color="k", s=14)
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel("Im")
    ax.set_zlabel("Re")
    ax.set_yticks([-3, 0, 3])
    ax.tick_params(axis="y", labelsize=7)
    ax.set_title("(a) $Z_1+Z_3+Z_5$: pulsating helix (winding of the central mode)")
    ax.view_init(elev=16, azim=-72)
    ax.legend(loc="upper left", fontsize=7.5, framealpha=0.9)
    ax.set_box_aspect((2.2, 1, 1))

    ax2 = fig.add_subplot(1, 2, 2)
    t = np.linspace(-np.pi / 2 + 1e-4, np.pi / 2 - 1e-4, 1200)
    S = np.sin(6 * t) / np.sin(t)
    ax2.plot(t, S, color=C_LOC, lw=2.0,
             label=r"$W+\bar{W}=\sin 6\theta/\sin\theta$  (localized wave)")
    ax2.plot(t, (1 + 2 * np.cos(2 * t)) ** 2, color=C_W, lw=1.4, ls="-.",
             label=r"$|W|^2=(1+2\cos 2\theta)^2$  (diagonal + cross)")
    ax2.axhline(3, color=C_NORM, lw=1.6, ls="--",
                label=r"$\sum_n Z_n\bar{Z}_n = 3$  (diagonal only: Parseval)")
    ax2.set_xlabel(r"$\theta$")
    ax2.set_title("(b) What each readout retains")
    ax2.legend(loc="upper right", fontsize=7.5, framealpha=0.9)
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{OUT}/fig_paper1_oddharm.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

# ---------------- 図3 ----------------

def fig3():
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    # (a) y = x^2 : Z2 fiber
    ax = axes[0]
    x = np.linspace(-2, 2, 300)
    ax.plot(x, x ** 2, color=C_Z, lw=2)
    y0 = 2.2
    x0 = np.sqrt(y0)
    ax.axhline(y0, color="gray", lw=0.9, ls="--")
    ax.plot([x0, -x0], [y0, y0], "o", color=C_NORM, ms=8)
    ax.annotate(r"$+x_0$", xy=(x0, y0), xytext=(x0 + 0.1, y0 + 0.25), fontsize=10)
    ax.annotate(r"$-x_0$", xy=(-x0, y0), xytext=(-x0 - 0.55, y0 + 0.25), fontsize=10)
    ax.set_title(r"(a) $y=x^2$ : fiber $\{\pm x\}\cong\mathbb{Z}_2$" + "\n(1 bit lost)")
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-0.3, 4.4)
    ax.grid(alpha=0.25)
    # (b) R -> S1 : Z fiber
    ax = axes[1]
    t = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(t), np.sin(t), color=C_Z, lw=2)
    ph = 0.9
    for k in range(4):
        yy = 1.5 + 0.5 * k
        ax.plot([np.cos(ph) * 0.2 + 1.6], [yy], "o", color=C_NORM, ms=6)
        ax.annotate(rf"$\theta_0+{2*k}\pi$" if k else r"$\theta_0$",
                    xy=(1.75, yy - 0.06), fontsize=8.5)
    ax.annotate("", xy=(np.cos(ph), np.sin(ph)), xytext=(1.6, 1.45),
                arrowprops=dict(arrowstyle="->", lw=1.0, color="gray"))
    ax.plot([np.cos(ph)], [np.sin(ph)], "o", color=C_SUM, ms=8)
    ax.set_title(r"(b) $\theta\ \mathrm{mod}\ 2\pi$ : fiber $\cong\mathbb{Z}$" + "\n(winding number lost)")
    ax.set_xlim(-1.6, 3.2)
    ax.set_ylim(-1.6, 3.4)
    ax.set_aspect("equal")
    ax.axis("off")
    # (c) ZZbar : S1 fiber
    ax = axes[2]
    ax.plot(np.cos(t), np.sin(t), color=C_Z, lw=2)
    for ang in np.linspace(0, 2 * np.pi, 14, endpoint=False):
        ax.annotate("", xy=(0.12 * np.cos(ang), 0.12 * np.sin(ang)),
                    xytext=(0.92 * np.cos(ang), 0.92 * np.sin(ang)),
                    arrowprops=dict(arrowstyle="->", lw=0.6, color="gray", alpha=0.7))
    ax.plot([0], [0], "o", color=C_NORM, ms=10)
    ax.annotate(r"$\rho^2$", xy=(0.08, -0.28), fontsize=11)
    ax.set_title(r"(c) $Z\mapsto Z\bar{Z}$ : fiber $\cong S^1$" + "\n(entire phase lost)")
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{OUT}/fig_paper1_fiber.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    print("done: fig_paper1_ladder / fig_paper1_oddharm / fig_paper1_fiber (png+svg)")
