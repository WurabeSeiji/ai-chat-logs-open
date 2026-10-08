#!/usr/bin/env python3
"""第十二思考実験の図 1〜3（模式図、データなし）。英語ラベル、SVG と PNG を出力する。

図 1：光円錐と中心投影は同じ直角三角形。
    (a) (r, t) 平面の光円錐。傾き v の世界線は円錐の内側。
    (b) 固有時間の軸 τ を足すと r² + τ² + (it)² = 0。同じ世界線が円錐の表面に乗る。
    (c) 中心投影 ℓ² = |x|² + R_c²。r ↔ |x|、τ = √(R² + ΣQᵢ²) ↔ R_c、t ↔ ℓ。sin θ = v。
図 2：観測者ごとに読める軸と、隠れた軸の回転対称性の次元 (8−k)(7−k)/2。
図 3：隠れた 5 方向の 3+2 分割。S(U(3)×U(2)) ⊂ SU(5)。
出力：fig01_lightcone_central_projection.{svg,png}、fig02_observers_symmetry.{svg,png}、fig03_split_3_2.{svg,png}
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle, Arc  # noqa: E402

HERE = Path(__file__).resolve().parent
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "figure.dpi": 150, "savefig.dpi": 150,
                     "svg.hashsalt": "te12", "svg.fonttype": "path"})
C_R, C_TAU, C_T = "C0", "C2", "C3"   # r ↔ |x|, τ ↔ R_c, t ↔ ℓ
V = 0.45                             # 一般の速度（α とは置かない）
TAU = np.sqrt(1 - V ** 2)


def save(fig, stem):
    for ext in ("svg", "png"):
        fig.savefig(HERE / f"{stem}.{ext}", metadata={"Date": None} if ext == "svg" else None)
    print("saved", stem + ".svg/.png")


def angle_mark(ax, center, r, theta1, theta2, label, color="k", lab_r=1.35):
    ax.add_patch(Arc(center, 2 * r, 2 * r, theta1=theta1, theta2=theta2, color=color, lw=0.9))
    m = np.deg2rad((theta1 + theta2) / 2)
    ax.text(center[0] + lab_r * r * np.cos(m), center[1] + lab_r * r * np.sin(m), label, ha="center", va="center", color=color)


# ---------------- 図 1 ----------------
fig = plt.figure(figsize=(13.5, 4.9))

# (a) light cone in (r, t)
ax = fig.add_subplot(1, 3, 1)
s = np.linspace(0, 1.25, 2)
ax.fill_between([-1.25, 0, 1.25], [1.25, 0, 1.25], 1.25, color="0.93")
ax.plot(s, s, color="k", lw=1.4); ax.plot(-s, s, color="k", lw=1.4)
ax.plot(V * s, s, color=C_T, lw=2.2)
ax.plot([V], [1], "o", color=C_T, ms=4)
ax.plot([0, V], [1, 1], color=C_R, lw=2.2)
ax.text(V / 2, 1.04, "r = v t", color=C_R, ha="center", va="bottom")
ax.text(V + 0.05, 0.96, "event (r, t)", fontsize=8, va="top")
ax.text(0.9, 0.72, "light cone\nr = t", fontsize=8, ha="left")
ax.text(V * 0.62 + 0.04, 0.62, "worldline\nv < 1", fontsize=8, color=C_T, ha="left")
ax.text(0, -0.08, "O", ha="center", va="top")
ax.set_xlim(-1.3, 1.3); ax.set_ylim(-0.15, 1.3); ax.set_aspect("equal")
ax.set_xlabel("r"); ax.set_ylabel("t")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("(a) light cone in (r, t):\nthe worldline lies inside", fontsize=10)

# (b) add the τ axis
ax = fig.add_subplot(1, 3, 2, projection="3d")
th = np.linspace(0, 2 * np.pi, 160)
tt = np.linspace(0, 1.25, 9)
TH, TT = np.meshgrid(th, tt)
ax.plot_wireframe(TT * np.cos(TH), TT * np.sin(TH), TT, color="0.80", lw=0.4, rstride=1, cstride=16)
ax.plot(np.cos(th), np.sin(th), np.ones_like(th), color="0.45", lw=1.0)
P = np.array([V, TAU, 1.0])
sl = np.linspace(0, 1.2, 2)
ax.plot(P[0] * sl, P[1] * sl, P[2] * sl, color=C_T, lw=2.2)
ax.plot([0, 0], [0, 0], [0, 1], color="0.3", lw=0.9, ls="--")
ax.plot([0, 0], [0, TAU], [1, 1], color=C_TAU, lw=2.4)
ax.plot([0, V], [TAU, TAU], [1, 1], color=C_R, lw=2.4)
ax.plot([0, V], [0, TAU], [1, 1], color=C_T, lw=2.4, ls=(0, (4, 2)))
ax.scatter([V], [TAU], [1], color=C_T, s=18)
ax.text(0.02, TAU / 2, 1.08, "$\\tau$", color=C_TAU, fontsize=11)
ax.text(V / 2, TAU + 0.06, 1.06, "r", color=C_R, fontsize=11)
ax.text(V / 2 + 0.05, TAU / 2 - 0.25, 0.93, "t", color=C_T, fontsize=11)
ax.text(0.05, -0.05, -0.08, "O", fontsize=9)
ax.set_xlabel("r"); ax.set_ylabel("$\\tau$"); ax.set_zlabel("t", labelpad=-2)
ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.25); ax.set_zlim(0, 1.3)
ax.view_init(elev=22, azim=-62)
ax.set_box_aspect((1, 1, 0.75))
ax.set_title("(b) add the proper-time axis $\\tau$:\n$r^2+\\tau^2+(it)^2=0$, the worldline lies on the cone", fontsize=10, pad=0)

# (c) central projection
ax = fig.add_subplot(1, 3, 3)
ax.add_patch(plt.Circle((0, 0), TAU, fill=False, color="0.45", lw=1.0))
ax.plot([-1.0, 1.0], [TAU, TAU], color="0.25", lw=1.0)
ax.text(-1.0, TAU + 0.03, "tangent space (x, y, z)", fontsize=8, va="bottom")
ax.plot([0, 0], [0, TAU], color=C_TAU, lw=2.4)
ax.plot([0, V], [TAU, TAU], color=C_R, lw=2.4)
ax.plot([0, V * 1.12], [0, TAU * 1.12], color=C_T, lw=2.4)
proj = TAU * np.array([V, TAU])
ax.plot(*proj, "o", color="k", ms=4)
ax.plot([V], [TAU], "o", color=C_T, ms=4)
ax.text(-0.05, TAU / 2, "$R_c$", color=C_TAU, ha="right", va="center", fontsize=10)
ax.text(V / 2, TAU + 0.035, "$|x|$", color=C_R, ha="center", va="bottom", fontsize=10)
ax.text(V / 2 + 0.07, TAU / 2 - 0.04, "$\\ell$", color=C_T, ha="left", va="center", fontsize=11)
ax.text(V + 0.04, TAU - 0.02, "P", va="top")
ax.text(proj[0] + 0.04, proj[1] - 0.05, "$\\pi(P)=R_c\\,P/\\ell$", fontsize=8, va="top")
ax.text(0.02, -0.06, "O", va="top")
angle_mark(ax, (0, 0), 0.22, np.rad2deg(np.arctan2(TAU, V)), 90, "$\\theta$", lab_r=1.45)
ax.text(-0.62, -0.5, "sphere of radius $R_c$", fontsize=8, color="0.35")
ax.text(-1.02, -0.98, "$r \\leftrightarrow |x|$      $\\tau=\\sqrt{R^2+\\sum Q_i^2} \\leftrightarrow R_c$      $t \\leftrightarrow \\ell$\n$\\sin\\theta = |x|/\\ell = r/t = v$", fontsize=8.4, va="top")
ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.3, 1.15); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("(c) central projection: $\\ell^2=|x|^2+R_c^2$,\nthe same right triangle", fontsize=10)

fig.subplots_adjust(left=0.04, right=0.99, top=0.86, bottom=0.08, wspace=0.12)
save(fig, "fig01_lightcone_central_projection")
plt.close(fig)

# ---------------- 図 2 ----------------
axes_names = ["x", "y", "z", "t", "R", "Q₁", "Q₂", "Q₃"]
observers = [("reads x, y, z", {0, 1, 2}), ("reads Q₁, Q₂, Q₃", {5, 6, 7}), ("reads t only", {3})]
fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.0, 3.9), gridspec_kw={"width_ratios": [1.25, 1]})
for row, (lab, readable) in enumerate(observers):
    y = len(observers) - 1 - row
    for j in range(8):
        on = j in readable
        ax.add_patch(Rectangle((j, y), 0.92, 0.82, color="C0" if on else "0.85"))
        ax.text(j + 0.46, y + 0.41, axes_names[j], ha="center", va="center", color="w" if on else "0.35", fontsize=10)
    k = len(readable)
    ax.text(-0.15, y + 0.41, lab, ha="right", va="center", fontsize=9)
    ax.text(8.15, y + 0.41, f"k = {k}: hidden {8 - k} seen only as ρ_hidden", ha="left", va="center", fontsize=8.4)
ax.add_patch(Rectangle((0, -0.75), 0.5, 0.4, color="C0")); ax.text(0.6, -0.55, "directly readable", va="center", fontsize=8)
ax.add_patch(Rectangle((3.6, -0.75), 0.5, 0.4, color="0.85")); ax.text(4.2, -0.55, "hidden (indistinguishable)", va="center", fontsize=8)
ax.set_xlim(-2.9, 12.6); ax.set_ylim(-0.9, 3.0); ax.axis("off")
ax.set_title("the same eight unnamed directions, three observation maps", fontsize=10)

ks = np.arange(1, 8)
dims = (8 - ks) * (7 - ks) // 2
cols = ["C3" if k == 1 else ("C0" if k == 3 else "0.6") for k in ks]
bx.bar(ks, dims, color=cols, width=0.65)
for k, d in zip(ks, dims):
    bx.text(k, d + 0.4, str(d), ha="center", va="bottom", fontsize=8)
bx.text(1.45, 19.0, "t only: O(7)", color="C3", fontsize=8.4, va="center")
bx.text(3.35, 10.2, "x, y, z: O(5) type", color="C0", fontsize=8.4, va="center")
bx.set_xticks(ks); bx.set_xlabel("number of directly readable directions k")
bx.set_ylabel("dimension of hidden rotations")
bx.set_ylim(0, 24.5)
bx.spines[["top", "right"]].set_visible(False)
bx.set_title("dim O(8−k) = (8−k)(7−k)/2\n(O(7−k, 1) when t is hidden: same dimension)", fontsize=10)
fig.tight_layout()
save(fig, "fig02_observers_symmetry")
plt.close(fig)

# ---------------- 図 3 ----------------
names = ["Q₁", "Q₂", "Q₃", "t", "R"]
fig, ax = plt.subplots(figsize=(6.4, 4.6))
for i in range(5):
    for j in range(5):
        blk3 = i < 3 and j < 3
        blk2 = i >= 3 and j >= 3
        col = "C0" if blk3 else ("C1" if blk2 else "0.92")
        ax.add_patch(Rectangle((j, 4 - i), 0.94, 0.94, color=col, alpha=0.85 if (blk3 or blk2) else 1.0))
for j, n in enumerate(names):
    ax.text(j + 0.47, 5.12, n, ha="center", va="bottom", fontsize=10)
    ax.text(-0.12, 4 - j + 0.47, n, ha="right", va="center", fontsize=10)
ax.text(1.47, 2.5 + 0.47, "U(3)", ha="center", va="center", color="w", fontsize=13, weight="bold")
ax.text(4.0, 0.5 + 0.47, "U(2)", ha="center", va="center", color="w", fontsize=12, weight="bold")
ax.text(4.0, 3.5, "mixing\n3 ↔ 2\ndistinguishable", ha="center", va="center", fontsize=7.6, color="0.4")
ax.text(1.5, 0.95, "mixing 3 ↔ 2: distinguishable", ha="center", va="center", fontsize=7.6, color="0.4")
ax.text(-0.9, -0.35, "common phase / det = 1:  S(U(3)×U(2)) ⊂ SU(5)\nS(U(3)×U(2)) ≅ (SU(3)×SU(2)×U(1)) / ℤ₆", fontsize=8.6, va="top")
ax.set_xlim(-0.95, 5.3); ax.set_ylim(-1.3, 5.6); ax.set_aspect("equal"); ax.axis("off")
ax.set_title("3+2 split of the five hidden directions", fontsize=10)
fig.tight_layout()
save(fig, "fig03_split_3_2")
plt.close(fig)
