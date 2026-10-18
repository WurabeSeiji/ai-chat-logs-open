# -*- coding: utf-8 -*-
"""四元数による状態空間から時空への写像：論文用の静的な図（PNG と PDF）を作る。

図 1  写像の連鎖 {zₙ} → {qₙ} → (1＋3) と、各段の全投射・部分投射
図 2  四元円 S³ を実部 q₀＝cos θ で切った球面（θ＝0°, 30°, 45°, 60°, 90°）と子午線
図 3  北極の接平面への中心投影（q₃＝0 の切断の 3 次元図と、q₂＝q₃＝0 の断面）
図 4  光円錐と本当の二乗（q₀ と |v| の平面）
図 5  正十六胞体 {±1, ±i, ±j, ±k}：頂点と Σq₀² : Σ|v|² ＝ 1 : 3
図 6  最も対称な形 q ＝ ±e^(π n/4)（n は正二十面体）：像 u ＝ n と Σq² ＝ 24 cos 2θ

点はすべて四元数の式から直接計算する（近似なし）。
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

import sys
EN = len(sys.argv) > 1 and sys.argv[1] == "en"
L = (lambda ja, en: en if EN else ja)  # 英語版は引数 en。計算と描画は同じで、文字だけを英語にする
plt.rcParams["font.family"] = ["Hiragino Sans", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

C_S3, C_PROJ, C_GEO, C_GRID, C_INK, C_MUTED = "#2a56c6", "#c4470f", "#0e7a6e", "#b9c2ce", "#18212c", "#596577"


def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return np.array([a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
                     a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0])


def icosa_dirs():
    p = (1 + 5 ** 0.5) / 2
    raw = [(0, 1, p), (0, -1, p), (0, 1, -p), (0, -1, -p), (1, p, 0), (-1, p, 0), (1, -p, 0), (-1, -p, 0),
           (p, 0, 1), (-p, 0, 1), (p, 0, -1), (-p, 0, -1)]
    n = np.array(raw, float)
    return n / np.linalg.norm(n, axis=1)[:, None]


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}{'_en' if EN else ''}.{ext}", dpi=200, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("saved", name + ("_en" if EN else ""))


def wire_sphere(ax, r, color, alpha, n_lat=7, n_lon=10, center=(0, 0, 0)):
    t = np.linspace(0, 2*np.pi, 120)
    for i in range(1, n_lat):
        ph = np.pi*i/n_lat
        ax.plot(center[0] + r*np.sin(ph)*np.cos(t), center[1] + r*np.sin(ph)*np.sin(t), center[2] + r*np.cos(ph)*np.ones_like(t),
                color=color, lw=0.6, alpha=alpha)
    s = np.linspace(0, np.pi, 80)
    for j in range(n_lon):
        a = 2*np.pi*j/n_lon
        ax.plot(center[0] + r*np.sin(s)*np.cos(a), center[1] + r*np.sin(s)*np.sin(a), center[2] + r*np.cos(s),
                color=color, lw=0.6, alpha=alpha)


def clean3d(ax, lim, labels):
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel(labels[0], labelpad=-8); ax.set_ylabel(labels[1], labelpad=-8); ax.set_zlabel(labels[2], labelpad=-8)
    ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
    ax.view_init(elev=18, azim=-58)


# ---------------------------------------------------------------- 図 1
def fig1():
    fig, ax = plt.subplots(figsize=(11, 3.4))
    ax.set_xlim(0, 11); ax.set_ylim(0, 3.4); ax.axis("off")
    boxes = [(0.3, L("{zₙ}  n＝1, …, 2N", "{zₙ}  n = 1, …, 2N"), L("複素数（実数 4N 個）", "complex numbers (4N reals)")),
             (4.15, L("{qₙ}  n＝1, …, N", "{qₙ}  n = 1, …, N"), L("四元数 qₙ ∈ ℍ（実数 4N 個）", "quaternions qₙ ∈ ℍ (4N reals)")),
             (8.0, L("(1＋3) の読み出し", "(1+3) readout"), L("N 個の点（t, xyz）", "N points (t, xyz)"))]
    for x, top, sub in boxes:
        ax.add_patch(FancyBboxPatch((x, 1.45), 2.7, 1.25, boxstyle="round,pad=0.05", fc="white", ec=C_INK, lw=1.2))
        ax.text(x + 1.35, 2.25, top, ha="center", va="center", fontsize=12.5)
        ax.text(x + 1.35, 1.78, sub, ha="center", va="center", fontsize=10, color=C_MUTED)
    for x0, x1, a, b, col in [(3.05, 4.1, L("二つずつ組にする  q ＝ z₁ ＋ z₂ j", "pair them  q = z₁ + z₂ j"), L("全投射（逆に戻せる）", "full projection (invertible)"), C_GEO),
                              (6.9, 7.95, L("中心投影  u ＝ v/q₀ と実部の側", "central projection  u = v/q₀ and the real part"), L("部分投射（|q|・±・赤道が逃げる）", "partial projection (|q|, ±, equator escape)"), C_PROJ)]:
        ax.annotate("", xy=(x1, 2.08), xytext=(x0, 2.08), arrowprops=dict(arrowstyle="-|>", lw=1.8, color=col))
        xm = (x0 + x1)/2
        ax.text(xm, 1.05, a, ha="center", va="center", fontsize=9.5)
        ax.text(xm, 0.62, b, ha="center", va="center", fontsize=9.5, color=col)
    ax.text(5.5, 3.15, L("条件：Σqₙ ＝ 0、Σqₙ² ＝ 0（無名性・自己無撞着・区別可能性の公理）", "Conditions: Σqₙ = 0, Σqₙ² = 0 (anonymity, self-consistency, distinguishability axioms)"), ha="center", fontsize=10, color=C_MUTED)
    save(fig, "fig1_mapping_chain")


# ---------------------------------------------------------------- 図 2
def fig2():
    thetas = [0, 30, 45, 60, 90]
    fig = plt.figure(figsize=(14, 5.6))
    for k, td in enumerate(thetas):
        th = np.deg2rad(td)
        ax = fig.add_subplot(2, 5, k + 1, projection="3d")
        wire_sphere(ax, 1, C_GRID, 0.6)
        r = np.sin(th)
        if r > 1e-6:
            wire_sphere(ax, r, C_S3, 0.95, 8, 12)
        else:
            ax.scatter([0], [0], [0], color=C_S3, s=25)
        clean3d(ax, 1.05, ["q₁", "q₂", "q₃"])
        ax.set_title(L(f"θ ＝ {td}°   q₀ ＝ {np.cos(th):.3f}\n切断球の半径 sin θ ＝ {r:.3f}", f"θ = {td}°   q₀ = {np.cos(th):.3f}\nslice radius sin θ = {r:.3f}"), fontsize=10)
        ax2 = fig.add_subplot(2, 5, k + 6)
        t = np.linspace(0, 2*np.pi, 300)
        ax2.plot(np.cos(t), np.sin(t), color=C_MUTED, lw=1)
        ax2.axhline(0, color=C_GRID, lw=0.8); ax2.axvline(0, color=C_GRID, lw=0.8)
        ax2.plot([-r, r], [np.cos(th)]*2, color=C_S3, lw=2.5)
        ax2.scatter([-r, r], [np.cos(th)]*2, color=C_S3, s=18, zorder=3)
        ax2.set_aspect("equal"); ax2.set_xlim(-1.25, 1.25); ax2.set_ylim(-1.25, 1.25)
        ax2.set_xticks([]); ax2.set_yticks([])
        ax2.set_xlabel(L("|v| の方向（子午線）", "direction of |v| (meridian)"), fontsize=9); ax2.set_ylabel("q₀", fontsize=9)
        for s in ax2.spines.values():
            s.set_visible(False)
    fig.suptitle(L("図 2　四元円 S³：q₀²＋q₁²＋q₂²＋q₃²＝1 を実部 q₀＝cos θ で切ると、(q₁, q₂, q₃) の中の半径 sin θ の球面 S²（上）。下は子午線の上での切断の位置", "Figure 2. Cutting the quaternion circle S³: q₀²+q₁²+q₂²+q₃²=1 at real part q₀=cos θ gives a sphere S² of radius sin θ in (q₁, q₂, q₃) (top). Bottom: position of the cut on the meridian"),
                 fontsize=11, y=1.0)
    save(fig, "fig2_s3_latitude_slices")


# ---------------------------------------------------------------- 図 3
def fig3():
    fig = plt.figure(figsize=(12.5, 5.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    th = np.deg2rad(40)
    wire_sphere(ax, 1, C_GRID, 0.55)
    g = np.linspace(-1.6, 1.6, 9)
    for x in g:
        ax.plot([x, x], [-1.6, 1.6], [1, 1], color=C_GRID, lw=0.6)
        ax.plot([-1.6, 1.6], [x, x], [1, 1], color=C_GRID, lw=0.6)
    a = np.linspace(0, 2*np.pi, 200)
    ax.plot(np.sin(th)*np.cos(a), np.sin(th)*np.sin(a), np.cos(th)*np.ones_like(a), color=C_S3, lw=2)
    T = np.tan(th)
    ax.plot(T*np.cos(a), T*np.sin(a), np.ones_like(a), color=C_PROJ, lw=2)
    for k in range(8):
        b = 2*np.pi*k/8
        ax.plot([0, T*np.cos(b)], [0, T*np.sin(b)], [0, 1], color=C_PROJ, lw=0.8, ls="--")
        ax.scatter([np.sin(th)*np.cos(b)], [np.sin(th)*np.sin(b)], [np.cos(th)], color=C_S3, s=12)
    ax.scatter([0], [0], [1], color=C_INK, s=20); ax.text(0, 0, 1.12, L("北極 1", "north pole 1"), fontsize=9)
    ax.scatter([0], [0], [0], color=C_MUTED, s=10)
    clean3d(ax, 1.6, ["q₁", "q₂", "q₀"])
    ax.set_title(L("q₃＝0 の 3 次元（q₀, q₁, q₂）：θ＝40° の緯線（青）と\n接平面 q₀＝1 の上の像（橙、半径 tan θ）", "3D slice q₃=0 (q₀, q₁, q₂): latitude circle at θ=40° (blue) and\nits image on the tangent plane q₀=1 (orange, radius tan θ)"), fontsize=10)
    ax2 = fig.add_subplot(1, 2, 2)
    t = np.linspace(0, 2*np.pi, 400)
    ax2.plot(np.cos(t), np.sin(t), color=C_MUTED, lw=1)
    ax2.axhline(1, color=C_GRID, lw=1.2); ax2.text(-1.0, 1.06, L("接線 q₀＝1", "tangent q₀=1"), fontsize=9, color=C_MUTED, ha="center")
    ax2.axhline(0, color=C_GRID, lw=0.6); ax2.axvline(0, color=C_GRID, lw=0.6)
    for td, col in [(20, C_S3), (40, C_S3), (60, C_S3), (120, C_PROJ)]:
        th = np.deg2rad(td)
        p = np.array([np.sin(th), np.cos(th)])
        u = np.tan(th)
        if np.cos(th) > 0:
            ax2.plot([0, u], [0, 1], color=col, lw=1, ls="--")
        else:
            ax2.plot([p[0], 0, u], [p[1], 0, 1], color=col, lw=1, ls="--")
        ax2.scatter([p[0]], [p[1]], color=col, s=22, zorder=3)
        ax2.scatter([u], [1], color=C_PROJ, s=22, zorder=3)
        ax2.text(p[0] + 0.07, p[1] - 0.04, L(f"θ＝{td}°", f"θ={td}°"), fontsize=8.5, color=col, ha="left", va="top")
        ax2.text(u, 1.1, f"{u:+.2f}", fontsize=8.5, color=C_PROJ, ha="center")
    ax2.set_aspect("equal"); ax2.set_xlim(-2.3, 2.3); ax2.set_ylim(-1.3, 1.45)
    ax2.set_xlabel(L("q₁（q₂＝q₃＝0 の断面）", "q₁ (cross-section q₂=q₃=0)")); ax2.set_ylabel("q₀")
    ax2.set_title(L("断面：角 θ の点は接線の上の tan θ に写る。θ＝120° の点は\n原点の反対側を通り、−q と同じ点（tan 120° ＝ −1.73）に写る", "Cross-section: a point at angle θ maps to tan θ on the tangent. The point at θ=120°\npasses through the origin and maps to the same point as −q (tan 120° = −1.73)"), fontsize=10)
    fig.suptitle(L("図 3　北極の接平面への中心投影（グノモン射影）：u ＝ v/q₀ ＝ n tan θ", "Figure 3. Central projection onto the tangent plane at the north pole (gnomonic projection): u = v/q₀ = n tan θ"), fontsize=11.5)
    save(fig, "fig3_gnomonic_projection")


# ---------------------------------------------------------------- 図 4
def fig4():
    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    t = np.linspace(-np.pi/2, np.pi/2, 300)
    ax.plot(np.cos(t), np.sin(t), color=C_MUTED, lw=1)
    ax.plot([0, 1.3], [0, 1.3], color=C_INK, lw=1, ls="--"); ax.plot([0, 1.3], [0, -1.3], color=C_INK, lw=1, ls="--")
    ax.fill_between([0, 1.3], [0, 1.3], [1.3, 1.3], color=C_S3, alpha=0.06)
    ax.fill_between([0, 1.3], [0, -1.3], [-1.3, -1.3], color=C_S3, alpha=0.06)
    ax.text(0.25, 1.15, L("円錐の内側\n（q₀² ＞ |v|²）", "inside the cone\n(q₀² > |v|²)"), fontsize=9, color=C_S3)
    ax.text(1.0, 0.12, L("円錐の外側\n（q₀² ＜ |v|²）", "outside the cone\n(q₀² < |v|²)"), fontsize=9, color=C_PROJ)
    ax.text(0.95, 1.02, L("q₀ ＝ |v|", "q₀ = |v|"), fontsize=9); ax.text(0.95, -1.08, L("q₀ ＝ −|v|", "q₀ = −|v|"), fontsize=9)
    for td in (30, 45, 60):
        th = np.deg2rad(td)
        q = (np.sin(th), np.cos(th))
        q2 = (abs(np.sin(2*th)), np.cos(2*th))
        ax.scatter(*q, color=C_S3, s=40, zorder=4); ax.text(q[0] - 0.04, q[1], L(f"q（θ＝{td}°）", f"q (θ={td}°)"), fontsize=9, color=C_S3, ha="right", va="center")
        ax.scatter(*q2, color=C_PROJ, s=30, zorder=4, marker="s"); ax.text(q2[0] + 0.04, q2[1] - 0.07, L(f"q²（2θ＝{2*td}°）", f"q² (2θ={2*td}°)"), fontsize=9, color=C_PROJ)
        ax.annotate("", xy=q2, xytext=q, arrowprops=dict(arrowstyle="->", color=C_MUTED, lw=0.8))
    r0 = 2 ** -0.5
    ax.scatter([r0, r0], [r0, -r0], color=C_GEO, s=60, facecolors="none", lw=1.5, zorder=5)
    ax.text(r0 + 0.05, -r0 - 0.1, L("緑の円：四元円と円錐の交わり（θ＝45°, 135°）", "green circles: quaternion circle meets the cone (θ=45°, 135°)"), fontsize=8.5, color=C_GEO, ha="center", va="top")
    ax.axhline(0, color=C_GRID, lw=0.6); ax.axvline(0, color=C_GRID, lw=0.6)
    ax.set_aspect("equal"); ax.set_xlim(-0.3, 1.5); ax.set_ylim(-1.3, 1.3)
    ax.set_xlabel(L("|v|（虚部の大きさ）", "|v| (size of the imaginary part)")); ax.set_ylabel(L("q₀（実部）", "q₀ (real part)"))
    ax.set_title(L("図 4　光円錐と本当の二乗：q² ＝ cos 2θ ＋ n sin 2θ、実部 q₀² − |v|² ＝ cos 2θ。\nθ＝45° で q は円錐に乗り、q² は実部 0（赤道）に来る", "Figure 4. Light cone and the true square: q² = cos 2θ + n sin 2θ, real part q₀² − |v|² = cos 2θ.\nAt θ=45°, q lies on the cone and q² comes to real part 0 (the equator)"), fontsize=10.5)
    save(fig, "fig4_light_cone_true_square")


# ---------------------------------------------------------------- 図 5
def fig5():
    E = np.eye(4)
    V = np.array([E[0], -E[0], E[1], -E[1], E[2], -E[2], E[3], -E[3]])
    names = ["1", "−1", "i", "−i", "j", "−j", "k", "−k"]
    fig = plt.figure(figsize=(14, 4.8))
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    for p in range(8):
        for q in range(p + 1, 8):
            if abs(V[p] @ V[q]) < 1e-12:
                s = np.linspace(0, np.pi/2, 40)
                arc = np.outer(np.cos(s), V[p]) + np.outer(np.sin(s), V[q])
                ax.plot(arc[:, 1], arc[:, 2], arc[:, 3], color=C_MUTED, lw=0.8)
    for v, nm in zip(V, names):
        ax.scatter([v[1]], [v[2]], [v[3]], color=C_S3, s=30, alpha=0.4 + 0.6*(v[0] + 1)/2)
        if nm in ("1", "−1"):
            continue
        ax.text(v[1] + 0.06, v[2] + 0.06, v[3] + 0.06, nm, fontsize=9)
    ax.text(0.08, -0.25, -0.12, "±1", fontsize=9)
    clean3d(ax, 1.1, ["q₁", "q₂", "q₃"])
    ax.set_title(L("頂点と 24 本の辺（大円の弧）の虚部\n±1 は原点、±i, ±j, ±k は正八面体", "vertices (imaginary parts) and 24 edges\n±1 at the origin; ±i, ±j, ±k: octahedron"), fontsize=10)
    ax2 = fig.add_subplot(1, 3, 2)
    t = np.linspace(-np.pi/2, np.pi/2, 300)
    ax2.plot(np.cos(t), np.sin(t), color=C_MUTED, lw=1)
    ax2.plot([0, 1.3], [0, 1.3], color=C_INK, lw=1, ls="--"); ax2.plot([0, 1.3], [0, -1.3], color=C_INK, lw=1, ls="--")
    for v in V:
        q2 = qmul(v, v)
        ax2.scatter(np.linalg.norm(v[1:]), v[0], color=C_S3, s=50, zorder=3)
        ax2.scatter(np.linalg.norm(q2[1:]), q2[0], color=C_PROJ, s=30, marker="s", zorder=4)
    ax2.text(0.05, 1.05, L("±1：円錐の内側", "±1: inside the cone"), fontsize=9, color=C_S3); ax2.text(1.03, 0.08, L("±i, ±j, ±k：外側", "±i, ±j, ±k: outside"), fontsize=9, color=C_S3)
    ax2.text(0.05, 0.85, L("q²：(±1)² ＝ 1", "q²: (±1)² = 1"), fontsize=9, color=C_PROJ); ax2.text(0.8, -0.85, L("q²：(±i)² ＝ −1 など", "q²: (±i)² = −1, etc."), fontsize=9, color=C_PROJ)
    ax2.set_aspect("equal"); ax2.set_xlim(-0.05, 1.45); ax2.set_ylim(-1.3, 1.3)
    ax2.set_xlabel("|v|"); ax2.set_ylabel("q₀"); ax2.set_title(L("光円錐の上での位置（青：q、橙：q²）", "light cone (blue: q, orange: q²)"), fontsize=10)
    ax3 = fig.add_subplot(1, 3, 3)
    st2 = np.sum(V[:, 0]**2); sv2 = np.sum(V[:, 1:]**2)
    S2 = np.sum([qmul(v, v) for v in V], axis=0)
    ax3.barh([1, 0], [st2, sv2], color=[C_S3, C_PROJ])
    ax3.set_yticks([1, 0]); ax3.set_yticklabels(["Σq₀²", "Σ|v|²"])
    for y, val in [(1, st2), (0, sv2)]:
        ax3.text(val + 0.1, y, f"{val:.0f}", va="center")
    ax3.set_xlim(0, 7.2)
    ax3.set_title(L(f"二次モーメントは 2I（等方的な 13 本を満たす）。\nΣq₀² : Σ|v|² ＝ 1 : 3 で Σq² ＝ {S2[0]:+.0f}（≠ 0）", f"second moment is 2I (the 13 isotropic conditions hold).\nΣq₀² : Σ|v|² = 1 : 3, so Σq² = {S2[0]:+.0f} (≠ 0)"), fontsize=10)
    for s in ("top", "right"):
        ax3.spines[s].set_visible(False)
    fig.suptitle(L("図 5　正十六胞体 {±1, ±i, ±j, ±k}（四元数群 Q₈）：等方的な配置は Σq² ＝ 0 を満たさない", "Figure 5. The 16-cell {±1, ±i, ±j, ±k} (quaternion group Q₈): the isotropic configuration does not satisfy Σq² = 0"), fontsize=11.5, y=1.04)
    save(fig, "fig5_16cell_isotropic_vs_true_square")


# ---------------------------------------------------------------- 図 6
def fig6():
    n = icosa_dirs()
    fig = plt.figure(figsize=(13, 5.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    wire_sphere(ax, 1, C_INK, 0.25, 7, 10)
    d = n @ n.T
    best = np.max(d[d < 1 - 1e-9])
    for a in range(12):
        for b in range(a + 1, 12):
            if abs(d[a, b] - best) < 1e-9:
                ax.plot(*np.c_[n[a], n[b]], color=C_PROJ, lw=1.1)
    ax.scatter(n[:, 0], n[:, 1], n[:, 2], color=C_PROJ, s=30)
    ax.scatter([0], [0], [0], color=C_INK, s=15)
    clean3d(ax, 1.15, ["u₁", "u₂", "u₃"])
    ax.set_title(L("接超平面 q₀＝1 の中の像 u ＝ v/q₀ ＝ n（24 個の q が 12 点に写る）。\nすべて光円錐の像 |u|＝1 の上にあり、辺は直線（正二十面体）", "image u = v/q₀ = n in the tangent hyperplane q₀=1 (24 q map to 12 points).\nAll lie on |u|=1, the image of the light cone; edges are straight (icosahedron)"), fontsize=10)
    ax2 = fig.add_subplot(1, 2, 2)
    th = np.linspace(0, np.pi, 721)
    qs = []
    for t_ in th:
        S = np.zeros(4)
        for m in n:
            for sg in (1, -1):
                q = sg*np.r_[np.cos(t_), np.sin(t_)*m]
                S += qmul(q, q)
        qs.append(S)
    qs = np.array(qs)
    ax2.plot(np.rad2deg(th), qs[:, 0], color=C_S3, lw=2, label=L("Σq² の実部（計算）", "real part of Σq² (computed)"))
    ax2.plot(np.rad2deg(th), 24*np.cos(2*th), color=C_PROJ, lw=1, ls="--", label="24 cos 2θ")
    ax2.plot(np.rad2deg(th), np.linalg.norm(qs[:, 1:], axis=1), color=C_GEO, lw=1.5, label=L("|Σq² のベクトル部|", "|vector part of Σq²|"))
    ax2.axhline(0, color=C_GRID, lw=0.8)
    for z in (45, 135):
        ax2.axvline(z, color=C_MUTED, lw=0.8, ls=":")
        ax2.text(z, 21, f"{z}°", ha="center", fontsize=9)
    ax2.set_xlim(0, 180); ax2.set_xticks(range(0, 181, 45))
    ax2.set_xlabel(L("全体の緯度 θ（度）", "overall latitude θ (degrees)")); ax2.set_ylabel(L("Σq²（24 個の和）", "Σq² (sum over 24)"))
    ax2.legend(fontsize=9, loc="lower left", frameon=False)
    ax2.set_title(L("q ＝ ±(cos θ ＋ n sin θ) の 24 個で、Σq² ＝ 24 cos 2θ。\n0 になるのは θ＝45°（q ＝ ±e^(π n/4)）と 135° だけ", "For the 24 quaternions q = ±(cos θ + n sin θ), Σq² = 24 cos 2θ.\nIt vanishes only at θ=45° (q = ±e^(π n/4)) and 135°"), fontsize=10)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    fig.suptitle(L("図 6　最も対称な形 q ＝ ±e^(π n/4) ＝ ±(1 ＋ n)/√2（n は正二十面体の 12 方向、N ＝ 24）", "Figure 6. The most symmetric form q = ±e^(π n/4) = ±(1 + n)/√2 (n the 12 icosahedral directions, N = 24)"), fontsize=11.5)
    save(fig, "fig6_most_symmetric_icosahedral")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
