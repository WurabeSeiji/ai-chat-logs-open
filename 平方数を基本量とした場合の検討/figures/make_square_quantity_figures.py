#!/usr/bin/env python3
"""Generate figures for the square-quantity readout paper.

All labels are in English so the same figures can be used in the
English manuscript. The script writes PNG and SVG versions.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


OUT = Path(__file__).resolve().parent


def save(fig, stem):
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def fig01_2d_square_map():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    ax_s, ax_p = axes

    theta = np.linspace(0, np.pi / 2, 240)
    x = np.cos(theta)
    y = np.sin(theta)
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(theta)))

    ax_s.plot(x, y, color="black", lw=1.6)
    sample_t = np.linspace(0.08, np.pi / 2 - 0.08, 7)
    sx = np.cos(sample_t)
    sy = np.sin(sample_t)
    sample_colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(sample_t)))
    ax_s.scatter(sx, sy, c=sample_colors, s=58, zorder=3, edgecolor="white", linewidth=0.8)
    ax_s.plot([0, 1], [0, 0], color="0.55", lw=1)
    ax_s.plot([0, 0], [0, 1], color="0.55", lw=1)
    ax_s.set_aspect("equal")
    ax_s.set_xlim(-0.06, 1.08)
    ax_s.set_ylim(-0.06, 1.08)
    ax_s.set_xlabel(r"$x$")
    ax_s.set_ylabel(r"$y$")
    ax_s.set_title(r"Root side: positive quadrant of $x^2+y^2=1$")
    ax_s.text(0.34, 0.78, "curved arc", fontsize=10)

    X = x * x
    Y = y * y
    ax_p.plot(X, Y, color="black", lw=1.6)
    ax_p.scatter(sx * sx, sy * sy, c=sample_colors, s=58, zorder=3, edgecolor="white", linewidth=0.8)
    ax_p.plot([0, 1], [1, 0], color="black", lw=1.6)
    ax_p.fill_between([0, 1], [1, 0], [0, 0], color="tab:blue", alpha=0.08)
    ax_p.set_aspect("equal")
    ax_p.set_xlim(-0.06, 1.08)
    ax_p.set_ylim(-0.06, 1.08)
    ax_p.set_xlabel(r"$X=x^2$")
    ax_p.set_ylabel(r"$Y=y^2$")
    ax_p.set_title(r"Square-quantity side: $X+Y=1$")
    ax_p.text(0.33, 0.58, "straight line", fontsize=10)

    fig.suptitle("Fig. 1. Squaring sends a curved positive arc to a linear simplex edge", fontsize=13)
    fig.text(
        0.5,
        0.02,
        r"Positive root is one-to-one: $(x,y)\leftrightarrow(X,Y)=(x^2,y^2)$, with $x=\sqrt{X}$, $y=\sqrt{Y}$.",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    save(fig, "fig01_2d_square_map")


def fig02_3d_square_map():
    fig = plt.figure(figsize=(12, 5.8))
    ax_s = fig.add_subplot(1, 2, 1, projection="3d")
    ax_p = fig.add_subplot(1, 2, 2, projection="3d")

    phi = np.linspace(0, np.pi / 2, 36)
    theta = np.linspace(0, np.pi / 2, 36)
    Phi, Theta = np.meshgrid(phi, theta)
    x = np.cos(Phi) * np.sin(Theta)
    y = np.sin(Phi) * np.sin(Theta)
    z = np.cos(Theta)
    ax_s.plot_surface(x, y, z, color="#4c78a8", alpha=0.34, linewidth=0, antialiased=True, shade=False)

    grid = np.linspace(0, np.pi / 2, 8)
    for t in grid:
        ph = np.linspace(0, np.pi / 2, 100)
        ax_s.plot(np.cos(ph) * np.sin(t), np.sin(ph) * np.sin(t), np.cos(t), color="0.35", lw=0.7)
        th = np.linspace(0, np.pi / 2, 100)
        ax_s.plot(np.cos(t) * np.sin(th), np.sin(t) * np.sin(th), np.cos(th), color="0.35", lw=0.7)

    samples = np.array(
        [
            [0.64, 0.25, 0.11],
            [0.50, 0.30, 0.20],
            [0.34, 0.18, 0.48],
            [0.18, 0.58, 0.24],
            [0.10, 0.25, 0.65],
        ]
    )
    sample_colors = plt.cm.viridis(np.linspace(0.10, 0.90, len(samples)))
    roots = np.sqrt(samples)
    ax_s.scatter(
        roots[:, 0],
        roots[:, 1],
        roots[:, 2],
        c=sample_colors,
        s=48,
        edgecolor="white",
        linewidth=0.8,
        zorder=5,
    )
    ax_s.scatter([1, 0, 0], [0, 1, 0], [0, 0, 1], color="black", s=32)
    ax_s.text(1.05, 0, 0, "x", fontsize=10)
    ax_s.text(0, 1.05, 0, "y", fontsize=10)
    ax_s.text(0, 0, 1.05, "z", fontsize=10)
    ax_s.set_title(r"Root side: octant of $x^2+y^2+z^2=1$")
    ax_s.set_xlim(-0.04, 1.08)
    ax_s.set_ylim(-0.04, 1.08)
    ax_s.set_zlim(-0.04, 1.08)
    ax_s.set_box_aspect((1, 1, 1))
    ax_s.view_init(elev=24, azim=38)
    ax_s.set_axis_off()

    tri = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    ax_p.add_collection3d(Poly3DCollection([tri], facecolor="#f58518", alpha=0.24, edgecolor="black", linewidth=1.3))
    for i, j in ((0, 1), (1, 2), (2, 0)):
        ax_p.plot([tri[i, 0], tri[j, 0]], [tri[i, 1], tri[j, 1]], [tri[i, 2], tri[j, 2]], color="black", lw=1.4)
    for t in np.linspace(0.1, 0.9, 5):
        ax_p.plot([t, t], [0, 1 - t], [1 - t, 0], color="0.45", lw=0.7)
        ax_p.plot([0, 1 - t], [t, t], [1 - t, 0], color="0.45", lw=0.7)
        ax_p.plot([0, 1 - t], [1 - t, 0], [t, t], color="0.45", lw=0.7)
    ax_p.scatter(
        samples[:, 0],
        samples[:, 1],
        samples[:, 2],
        c=sample_colors,
        s=48,
        edgecolor="white",
        linewidth=0.8,
        zorder=5,
    )
    ax_p.scatter([1, 0, 0], [0, 1, 0], [0, 0, 1], color="black", s=32)
    ax_p.text(1.06, 0, 0, "X", fontsize=10)
    ax_p.text(0, 1.06, 0, "Y", fontsize=10)
    ax_p.text(0, 0, 1.06, "Z", fontsize=10)
    ax_p.set_title(r"Square-quantity side: $X+Y+Z=1$")
    ax_p.set_xlim(-0.04, 1.08)
    ax_p.set_ylim(-0.04, 1.08)
    ax_p.set_zlim(-0.04, 1.08)
    ax_p.set_box_aspect((1, 1, 0.75))
    ax_p.view_init(elev=24, azim=38)
    ax_p.set_axis_off()

    fig.suptitle("Fig. 2. In 3D, the positive spherical octant reads as the standard simplex", fontsize=13)
    fig.text(
        0.5,
        0.02,
        r"The same map $(X,Y,Z)=(x^2,y^2,z^2)$ changes a spherical constraint into an affine plane; matched colors mark exact corresponding points.",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    save(fig, "fig02_3d_square_map")


def fig03_motion_readouts():
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4))
    T = np.linspace(0, 1, 300)
    t = np.sqrt(T)
    A = 0.72
    B = 1.2

    axes[0, 0].plot(T, A * T, color="#4c78a8", lw=2.4)
    axes[0, 0].set_title(r"P side: linear readout $X=A\,T$")
    axes[0, 0].set_xlabel(r"$T=t^2$")
    axes[0, 0].set_ylabel(r"$X=x^2$")
    axes[0, 0].grid(alpha=0.28)
    axes[0, 0].text(0.08, 0.58, f"A = {A:.2f}", fontsize=10)

    axes[0, 1].plot(t, np.sqrt(A) * t, color="#4c78a8", lw=2.4)
    axes[0, 1].set_title(r"Root side: $x=\sqrt{A}\,t$")
    axes[0, 1].set_xlabel(r"$t$")
    axes[0, 1].set_ylabel(r"$x$")
    axes[0, 1].grid(alpha=0.28)
    axes[0, 1].text(0.08, 0.78, "same uniform-motion form", fontsize=10)

    axes[1, 0].plot(T, 0.5 * B * T * T, color="#f58518", lw=2.4)
    axes[1, 0].set_title(r"P side: quadratic readout $X=\frac{1}{2}B\,T^2$")
    axes[1, 0].set_xlabel(r"$T=t^2$")
    axes[1, 0].set_ylabel(r"$X=x^2$")
    axes[1, 0].grid(alpha=0.28)
    axes[1, 0].text(0.08, 0.43, f"B = {B:.2f}", fontsize=10)

    axes[1, 1].plot(t, np.sqrt(B / 2.0) * t * t, color="#f58518", lw=2.4)
    axes[1, 1].set_title(r"Root side: $x=\sqrt{B/2}\,t^2$")
    axes[1, 1].set_xlabel(r"$t$")
    axes[1, 1].set_ylabel(r"$x$")
    axes[1, 1].grid(alpha=0.28)
    axes[1, 1].text(0.08, 0.58, "same constant-acceleration form", fontsize=10)

    for ax in axes.flat:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    fig.suptitle("Fig. 3. Elementary motion forms are recovered by square-root readout", fontsize=13)
    fig.text(
        0.5,
        0.02,
        "No physical law is assumed here; these are algebraic readouts of equations placed on the square-quantity side.",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    save(fig, "fig03_motion_readouts")


def fig04_quadratic_readings():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4))

    t = np.linspace(0, 2 * np.pi, 500)
    axes[0].plot(np.cos(t), np.sin(t), color="#4c78a8", lw=2.2)
    axes[0].set_aspect("equal")
    axes[0].set_title(r"Elliptic reading" + "\n" + r"$x^2+y^2=1$")
    axes[0].set_xlabel(r"$x$")
    axes[0].set_ylabel(r"$y$")
    axes[0].grid(alpha=0.22)

    y = np.linspace(-1.05, 1.05, 220)
    axes[1].plot(y, y, color="#54a24b", lw=2.2)
    axes[1].plot(y, -y, color="#54a24b", lw=2.2)
    axes[1].set_aspect("equal")
    axes[1].set_title(r"Double-cone reading" + "\n" + r"$x^2-y^2=0$")
    axes[1].set_xlabel(r"$x$")
    axes[1].set_ylabel(r"$y$")
    axes[1].grid(alpha=0.22)

    u = np.linspace(-1.4, 1.4, 500)
    axes[2].plot(np.cosh(u), np.sinh(u), color="#e45756", lw=2.2)
    axes[2].plot(-np.cosh(u), np.sinh(u), color="#e45756", lw=2.2)
    axes[2].set_aspect("equal")
    axes[2].set_title(r"Hyperbolic reading" + "\n" + r"$x^2-y^2=1$")
    axes[2].set_xlabel(r"$x$")
    axes[2].set_ylabel(r"$y$")
    axes[2].grid(alpha=0.22)

    for ax in axes:
        ax.axhline(0, color="0.25", lw=0.8)
        ax.axvline(0, color="0.25", lw=0.8)

    fig.suptitle("Fig. 4. Different root-side geometries can share a simple quadratic-form origin", fontsize=13)
    fig.text(
        0.5,
        0.01,
        "This figure is only geometric: the signs and right-hand side choose an elliptic, conic, or hyperbolic reading.",
        ha="center",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.92))
    save(fig, "fig04_quadratic_readings")


def theta_of_R(R):
    return np.arccos(-(np.tan(1.0 / (2.0 * R)) ** 2))


def ks_of_R(R):
    theta = theta_of_R(R)
    return R * R * (4.0 * theta - 2.0 * np.pi)


def fig05_area_coefficient():
    Rmin = 2.0 / np.pi + 1e-4
    R = np.linspace(Rmin, 8.0, 900)
    ks = ks_of_R(R)

    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.plot(R, ks, color="#b279a2", lw=2.4)
    ax.axhline(1.0, color="0.25", lw=1.0, linestyle="--", label="flat area = 1")
    for rv in (1.0, 2.0, 3.0, 5.0):
        ax.scatter([rv], [ks_of_R(rv)], color="#b279a2", s=42, zorder=3)
        ax.text(rv + 0.06, ks_of_R(rv) + 0.012, f"R={rv:g}", fontsize=9)
    ax.set_xlim(0.55, 8.0)
    ax.set_ylim(0.96, min(1.55, np.nanmax(ks) * 1.02))
    ax.set_xlabel(r"curvature radius $R$")
    ax.set_ylabel(r"area correction coefficient $k_s(R)$")
    ax.set_title("Fig. 5. Exact area coefficient for a unit geodesic square")
    ax.grid(alpha=0.28)
    ax.legend(loc="upper right")
    ax.text(
        0.62,
        1.43,
        r"$k_s(R)=R^2\,[4\arccos(-\tan^2(1/(2R)))-2\pi]$" + "\n"
        r"domain shown: $R>2/\pi$ for a positive-curvature unit square",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="0.75"),
    )
    fig.tight_layout()
    save(fig, "fig05_area_coefficient_ks")


def main():
    fig01_2d_square_map()
    fig02_3d_square_map()
    fig03_motion_readouts()
    fig04_quadratic_readings()
    fig05_area_coefficient()
    print("Generated square-quantity readout figures in", OUT)


if __name__ == "__main__":
    main()
