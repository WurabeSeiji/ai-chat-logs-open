# -*- coding: utf-8 -*-
"""観測は放射と吸収の同時成立である ― 論文用の静的な図（PNG・SVG・PDF）を作る。

図 1  視野角と大きさ（P1）：二体の関係点 X で交わる二本の見る線、共有される視野角 θ_AB、
      大きさ D_A・D_B と距離 d＝D/θ。同じ視野角なら相手が大きいほど遠い、視野角が狭いほど遠い、
      自分自身は二本の線が平行で θ_AA＝0（関係点なし）
図 2  三体：三つの視野角と大きさの比は読め、全体のスケール λ は読めない（λ 倍の写しで角は同じ）
図 3  倍音（P3）：有限和 N＝13 で厳密に計算したパルス形と自己相関 E(θ)＝Σ w_m cos mθ。単音、矩形波（奇数倍音 m≤13、c_m∝sin(mπ/2)/m）、
      等振幅の局在波（m＝1…13、奇数・偶数）。三角波は矩形波の M→∞ の極限として参考に点線で重ねる
図 4  ツィレルソン限界（3.4）：単位ベクトル a, a′, b, b′ と b±b′、S(δ)＝|3cos nδ − cos 3nδ| の 2 と 2√2
図 5  二者相関（§4）：E(θ)＝−cos nθ（n＝1, 2, 4）と最大の破れを与える角の刻み π/(4n)

点はすべて式から直接計算する。乱数は使わない。
幾何の約束（図 1・2）：二つの円板の内共通接線を「見る線」とし、その交点 X を関係点と呼ぶ。
  交点の角が視野角 θ_AB（A から見ても B から見ても同じ）、X から各円板の中心までが d＝(D/2)/sin(θ/2)、
  小角で d＝D/θ。自分自身では二本の線が平行になり θ＝0、関係点は無限遠に行く。
実行:  python make_figures.py [en]     （en で英語のラベル。計算と描画は同じで文字だけが違う）
出力:  figures/figN_*.png / .svg / .pdf（英語版は figN_*_en）
依存:  numpy, matplotlib
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc, FancyArrowPatch, Wedge

EN = len(sys.argv) > 1 and sys.argv[1] == "en"
L = (lambda ja, en: en if EN else ja)
SUF = "_en" if EN else ""
plt.rcParams["font.family"] = ["Hiragino Sans", "DejaVu Sans"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 10
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

C_A, C_B, C_C, C_INK, C_MUTED, C_GRID, C_LINE = "#c4470f", "#2a56c6", "#0e7a6e", "#18212c", "#596577", "#c9d0d8", "#7a8699"


def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(OUT / f"{name}{SUF}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved", name + SUF)


# ---------------------------------------------------------------- 幾何の小道具
def inner_center(ca, ra, cb, rb):
    """内相似中心（内共通接線の交点）"""
    ca, cb = np.asarray(ca, float), np.asarray(cb, float)
    return (rb * ca + ra * cb) / (ra + rb)


def half_angle(ca, ra, cb, rb):
    """内共通接線の交角の半分：sin(θ/2) ＝ (r_a ＋ r_b)/|c_b − c_a|"""
    Lab = np.hypot(*(np.asarray(cb, float) - np.asarray(ca, float)))
    return np.arcsin((ra + rb) / Lab)


def tangent_pts(c, r, x):
    """外点 x から円 (c, r) への接点（＋側, −側）"""
    c, x = np.asarray(c, float), np.asarray(x, float)
    v = x - c
    d = np.hypot(*v)
    beta = np.arctan2(v[1], v[0])
    gam = np.arccos(r / d)
    return [c + r * np.array([np.cos(beta + s * gam), np.sin(beta + s * gam)]) for s in (+1, -1)]


def extend(p, q, before=0.0, after=0.0):
    """線分 p→q を両側へ伸ばした端点"""
    p, q = np.asarray(p, float), np.asarray(q, float)
    u = (q - p) / np.hypot(*(q - p))
    return p - before * u, q + after * u


def draw_pair(ax, ca, ra, cb, rb, col_a, col_b, lab_a, lab_b, theta_label, ext=0.6, arc_r=None, show_X=True, lw=1.4, arc_both=True, label_gap=0.3, show_X_label=True, label_perp=False, wedge=False):
    """二体と内共通接線、関係点 X と視野角の弧。戻り値 (X, θ)"""
    X = inner_center(ca, ra, cb, rb)
    h = half_angle(ca, ra, cb, rb)
    ta = tangent_pts(ca, ra, X)
    tb = tangent_pts(cb, rb, X)
    # 線 1：A の＋側 → B の−側、線 2：A の−側 → B の＋側
    for pa, pb in ((ta[0], tb[1]), (ta[1], tb[0])):
        p, q = extend(pa, pb, ext, ext)
        ax.plot([p[0], q[0]], [p[1], q[1]], color=C_LINE, lw=lw, zorder=1)
    for c, r, col, lab in ((ca, ra, col_a, lab_a), (cb, rb, col_b, lab_b)):
        ax.add_patch(Circle(c, r, facecolor=col, alpha=0.18, edgecolor=col, lw=1.6, zorder=2))
        ax.plot([c[0]], [c[1]], marker="+", color=col, ms=7, mew=1.2, zorder=3)
        ax.text(c[0], c[1] + r + 0.12, lab, ha="center", va="bottom", color=col, fontsize=11, fontweight="bold")
    if show_X:
        ax.plot([X[0]], [X[1]], "o", color=C_INK, ms=4.5, zorder=4)
        if show_X_label:
            ax.text(X[0], X[1] - 0.16, "X", ha="center", va="top", color=C_INK, fontsize=10)
    # 視野角の弧（B の側と A の側、同じ角）
    base = np.degrees(np.arctan2(cb[1] - ca[1], cb[0] - ca[0]))
    if arc_r is None:
        arc_r = 0.55 * min(np.hypot(*(X - np.asarray(ca))), np.hypot(*(X - np.asarray(cb))))
    sides = ((base, +1), (base + 180, -1)) if arc_both else ((base, +1),)
    for side, s in sides:
        if wedge:
            ax.add_patch(Wedge(X, arc_r, side - np.degrees(h), side + np.degrees(h), facecolor=C_INK, alpha=0.16, edgecolor=C_INK, lw=0.8, zorder=3))
        else:
            ax.add_patch(Arc(X, 2 * arc_r, 2 * arc_r, angle=0, theta1=side - np.degrees(h), theta2=side + np.degrees(h), color=C_INK, lw=1.2, zorder=3))
    ang = np.radians(base) + (np.pi / 2 if label_perp else 0.0)
    ax.text(X[0] + (arc_r + label_gap) * np.cos(ang), X[1] + (arc_r + label_gap) * np.sin(ang), theta_label, ha="center", va="center", color=C_INK, fontsize=11,
            bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.92), zorder=6)
    return X, 2 * h


def dim_line(ax, p, q, y, text, col=C_MUTED, dy=0.12, fontsize=9.5):
    """寸法線（水平）"""
    ax.annotate("", xy=(q, y), xytext=(p, y), arrowprops=dict(arrowstyle="<->", color=col, lw=1.0, shrinkA=0, shrinkB=0))
    ax.text((p + q) / 2, y - dy, text, ha="center", va="top", color=col, fontsize=fontsize)


def diam_bar(ax, c, r, x_off, text, col, fontsize=9.5):
    """直径の縦の寸法線"""
    x = c[0] + x_off
    ax.annotate("", xy=(x, c[1] + r), xytext=(x, c[1] - r), arrowprops=dict(arrowstyle="<->", color=col, lw=1.0, shrinkA=0, shrinkB=0))
    ax.text(x + (0.12 if x_off > 0 else -0.12), c[1], text, ha="left" if x_off > 0 else "right", va="center", color=col, fontsize=fontsize)


def tidy(ax, xlim, ylim, title=None):
    ax.set_aspect("equal")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=11, loc="left", color=C_INK)


# ================================================================ 図 1  視野角と大きさ（P1）
def fig1():
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 8.2))
    # (a) 二体の関係
    ax = axs[0, 0]
    ca, ra, cb, rb = np.array([0.0, 0.0]), 0.55, np.array([5.6, 0.0]), 1.0
    X, th = draw_pair(ax, ca, ra, cb, rb, C_A, C_B, "A", "B", r"$\theta_{AB}$")
    diam_bar(ax, ca, ra, -0.95, r"$D_A$", C_A)
    diam_bar(ax, cb, rb, 1.4, r"$D_B$", C_B)
    dim_line(ax, X[0], ca[0], -1.55, r"$d_{B\to A}=D_A/\theta_{AB}$", C_A)
    dim_line(ax, X[0], cb[0], -1.55, r"$d_{A\to B}=D_B/\theta_{AB}$", C_B)
    ax.text(X[0], 2.05, L("関係点 X：二本の見る線の交点\n視野角 $\\theta_{AB}$ は A から見ても B から見ても同じ",
                          "relation point X: crossing of the two lines of sight\n$\\theta_{AB}$ is the same seen from A and from B"),
            ha="center", va="top", fontsize=9.5, color=C_INK)
    ax.text(2.8, -2.35, L("$\\theta_{AB}$ ＝ %.1f°、$d_{A\\to B}/d_{B\\to A}=D_B/D_A$ ＝ %.2f（式どおり。厳密には $d=D/(2\\sin(\\theta/2))$、小角で $D/\\theta$）" % (np.degrees(th), rb / ra),
                          "$\\theta_{AB}$ = %.1f°, $d_{A\\to B}/d_{B\\to A}=D_B/D_A$ = %.2f (as the formula says; exactly $d=D/(2\\sin(\\theta/2))$, $D/\\theta$ for small angles)" % (np.degrees(th), rb / ra)),
            ha="center", va="top", fontsize=9, color=C_MUTED)
    tidy(ax, (-1.9, 7.5), (-2.9, 2.75), L("(a) 二体の関係：共有される視野角と、相手の大きさ", "(a) two bodies: shared viewing angle and the partner's size"))

    # (b) 同じ視野角なら相手が大きいほど遠い
    ax = axs[0, 1]
    h = np.radians(9.0)
    xmax = 9.6
    for s in (+1, -1):
        ax.plot([0, xmax], [0, s * xmax * np.tan(h)], color=C_LINE, lw=1.4)
    ax.plot([0], [0], "o", color=C_INK, ms=4.5)
    ax.text(0, -0.2, "X", ha="center", va="top", fontsize=10)
    ax.add_patch(Arc((0, 0), 2.2, 2.2, theta1=-np.degrees(h), theta2=np.degrees(h), color=C_INK, lw=1.2))
    ax.text(1.4, 0, r"$\theta$", ha="center", va="center", fontsize=11)
    r0 = 0.33
    for k, col in ((1, C_B), (2, C_C), (3, C_A)):
        r = k * r0
        x = r / np.sin(h)
        ax.add_patch(Circle((x, 0), r, facecolor=col, alpha=0.18, edgecolor=col, lw=1.6))
        ax.plot([x], [0], marker="+", color=col, ms=7, mew=1.2)
        ax.text(x, r + 0.12, L("D", "D") if k == 1 else f"{k}D", ha="center", va="bottom", color=col, fontsize=10.5, fontweight="bold")
        dim_line(ax, 0, x, -1.45 - 0.42 * (k - 1), (r"$d=D/\theta$" if k == 1 else rf"${k}d$"), col, dy=0.1, fontsize=9)
    ax.text(4.8, 1.9, L("同じ視野角 θ、相手が大きいほど遠い", "same viewing angle θ: a larger partner is farther"), ha="center", va="center", fontsize=9.5, color=C_INK)
    tidy(ax, (-0.8, 10.3), (-2.9, 2.3), L("(b) 同じ視野角なら相手が大きいほど遠い", "(b) at the same viewing angle, a larger partner is farther"))

    # (c) 視野角が狭いほど遠い
    ax = axs[1, 0]
    r0 = 0.42
    rows = ((18.0, C_A, 1.45), (10.0, C_C, 0.0), (5.5, C_B, -1.45))
    for deg, col, y0 in rows:
        h = np.radians(deg / 2)
        x = r0 / np.sin(h)
        xe = x + r0 / np.tan(h) * 0.35 + 0.6
        for s in (+1, -1):
            ax.plot([0, xe], [y0, y0 + s * xe * np.tan(h)], color=C_LINE, lw=1.3)
        ax.plot([0], [y0], "o", color=C_INK, ms=4)
        ax.add_patch(Arc((0, y0), 1.6, 1.6, theta1=-deg / 2, theta2=deg / 2, color=C_INK, lw=1.1))
        ax.add_patch(Circle((x, y0), r0, facecolor=col, alpha=0.18, edgecolor=col, lw=1.6))
        ax.plot([x], [y0], marker="+", color=col, ms=7, mew=1.2)
        ax.text(x, y0 + r0 + 0.1, "D", ha="center", va="bottom", color=col, fontsize=10.5, fontweight="bold")
        ax.text(-0.25, y0, "θ = %g°" % deg, ha="right", va="center", fontsize=9.5, color=C_INK)
        dim_line(ax, 0, x, y0 - 0.62, r"$d=D/\theta$" + ("" if deg == 18 else ""), col, dy=0.08, fontsize=8.5)
    ax.text(4.9, 2.55, L("同じ大きさ D、視野角が狭いほど遠い", "same size D: a narrower viewing angle is farther"), ha="center", va="center", fontsize=9.5, color=C_INK)
    tidy(ax, (-1.9, 10.3), (-2.5, 2.9), L("(c) 視野角が狭いほど遠い", "(c) a narrower viewing angle is farther"))

    # (d) 自分自身：平行、θ_AA ＝ 0
    ax = axs[1, 1]
    r = 0.9
    ax.add_patch(Circle((2.2, 0), r, facecolor=C_A, alpha=0.18, edgecolor=C_A, lw=1.6))
    ax.plot([2.2], [0], marker="+", color=C_A, ms=7, mew=1.2)
    ax.text(2.2, r + 0.12, "A", ha="center", va="bottom", color=C_A, fontsize=11, fontweight="bold")
    for s in (+1, -1):
        ax.plot([-0.8, 7.2], [s * r, s * r], color=C_LINE, lw=1.4)
    ax.annotate("", xy=(7.4, 0), xytext=(5.2, 0), arrowprops=dict(arrowstyle="->", color=C_MUTED, lw=1.0))
    ax.text(6.3, 0.15, L("関係点は無限遠", "relation point at infinity"), ha="center", va="bottom", fontsize=9, color=C_MUTED)
    ax.text(4.4, -1.55, L("見る線と見られる線が平行：$\\theta_{AA}$ ＝ 0\n交点がないので距離は定義できない（自分自身は観測できない）",
                          "the seeing and the seen lines are parallel: $\\theta_{AA}$ = 0\nno crossing, so no distance (one cannot observe oneself)"),
            ha="center", va="top", fontsize=9.5, color=C_INK)
    tidy(ax, (-1.2, 7.9), (-2.9, 2.3), L("(d) 自分自身：二本の線が平行で視野角 0", "(d) oneself: the two lines are parallel, viewing angle 0"))
    fig.suptitle(L("図 1　視野角と大きさ（P1）：読めるのは視野角と大きさの比だけ", "Fig. 1  Viewing angle and size (P1): only the angle and size ratios are readable"), fontsize=12, color=C_INK, y=0.995)
    fig.tight_layout()
    save(fig, "fig1_viewing_angle_and_size")


# ================================================================ 図 2  三体とスケールのゲージ
def fig2():
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.6))
    bodies = [(np.array([0.0, 0.0]), 0.5, C_A, "A"), (np.array([6.2, 0.4]), 0.95, C_B, "B"), (np.array([2.4, 4.4]), 0.35, C_C, "C")]
    pairs = [(0, 1, r"$\theta_{AB}$"), (1, 2, r"$\theta_{BC}$"), (2, 0, r"$\theta_{CA}$")]
    angles = {}
    for ax, lam, ttl in ((axs[0], 1.0, L("(a) 三体 A, B, C", "(a) three bodies A, B, C")), (axs[1], 1.6, L("(b) 全体を λ ＝ 1.6 倍した写し", "(b) the same scaled by λ = 1.6"))):
        for (i, j, lab) in pairs:
            ci, ri, coli, li = bodies[i]
            cj, rj, colj, lj = bodies[j]
            X, th = draw_pair(ax, lam * ci, lam * ri, lam * cj, lam * rj, coli, colj, "", "", lab, ext=0.3, arc_r=0.9 * lam, show_X=True, arc_both=False, label_gap=-0.55 * lam, show_X_label=False, label_perp=True, wedge=True)
            angles[(lab, lam)] = np.degrees(th)
        for c, r, col, lab in bodies:
            ax.text(lam * c[0], lam * (c[1] + r) + 0.14, lab, ha="center", va="bottom", color=col, fontsize=11, fontweight="bold")
        txt = "\n".join("%s ＝ %.1f°" % (lab, angles[(lab, lam)]) for (_, _, lab) in pairs)
        ax.text(9.1, 6.2, txt, ha="left", va="top", fontsize=10.5, color=C_INK, linespacing=1.5)
        ax.text(9.1, 3.9, L("大きさの比\n$D_A : D_B : D_C$ ＝ 1 : %.1f : %.1f" % (bodies[1][1] / bodies[0][1], bodies[2][1] / bodies[0][1]),
                           "size ratios\n$D_A : D_B : D_C$ = 1 : %.1f : %.1f" % (bodies[1][1] / bodies[0][1], bodies[2][1] / bodies[0][1])), ha="left", va="top", fontsize=10, color=C_INK)
        tidy(ax, (-1.6, 13.2), (-1.8, 8.4), ttl)
    fig.suptitle(L("図 2　三体：三つの視野角と大きさの比は読め、全体のスケール λ は読めない（角は同じ）",
                   "Fig. 2  Three bodies: the three viewing angles and the size ratios are readable; the overall scale λ is not (angles coincide)"), fontsize=12, color=C_INK, y=0.995)
    fig.tight_layout()
    save(fig, "fig2_three_bodies_gauge")


# ================================================================ 図 3  倍音とパルス形の自己相関（P3）
def coeffs(phi, M):
    """矩形パルス（幅 φ 度）のフーリエ係数：c_0 ＝ φ/360、c_m ＝ (2/(mπ)) sin(mπφ/360)"""
    m = np.arange(1, M + 1)
    return phi / 360.0, (2.0 / (m * np.pi)) * np.sin(m * np.pi * phi / 360.0), m


def fig3():
    """有限和 N ＝ 13 で厳密に計算する。矩形波：奇数倍音 m ≤ N、c_m ∝ sin(mπ/2)/m（符号は交互、弦の半分で +1・残りで −1 の矩形）。局在波：等振幅、m ＝ 1…N（奇数・偶数）"""
    N = 13
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.8))
    x = np.linspace(-180, 180, 2881)
    m_all = np.arange(1, N + 1)
    m_odd = m_all[m_all % 2 == 1]

    def partial(ms, cs, xx):
        """Σ_m c_m cos(m xx)（xx は度）"""
        return (cs[:, None] * np.cos(np.radians(ms[:, None] * xx[None, :]))).sum(0)

    cases = [
        (L("単音（m ＝ 1）", "single tone (m = 1)"), C_B, np.array([1]), np.array([1.0])),
        (L("矩形波（奇数倍音 m ≤ %d、c_m ∝ sin(mπ/2)/m）" % N, "square wave (odd harmonics m ≤ %d, c_m ∝ sin(mπ/2)/m)" % N), C_A, m_odd, np.sin(m_odd * np.pi / 2) / m_odd),
        (L("局在波（等振幅、m ＝ 1…%d、奇数・偶数）" % N, "localized wave (equal amplitudes, m = 1…%d, odd and even)" % N), C_C, m_all, np.ones(N)),
    ]
    ax = axs[0]
    for lab, col, ms, cs in cases:
        f = partial(ms, cs, x)
        ax.plot(x, f / np.abs(f).max(), color=col, lw=1.6, label=lab)
    ax.axhline(0, color=C_GRID, lw=0.8)
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xlabel(L("弦の上の角 ψ − θ₀ [°]", "angle on the string ψ − θ₀ [deg]"))
    ax.set_ylabel(L("パルス f ＝ Σ c_m cos m(ψ−θ₀)（最大 1 に正規化）", "pulse f = Σ c_m cos m(ψ−θ₀) (scaled to max 1)"))
    ax.set_title(L("(a) パルス形（有限和、N ＝ %d）" % N, "(a) pulse shapes (finite sums, N = %d)" % N), loc="left", fontsize=11)
    ax.legend(fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1, frameon=False)
    ax.grid(color=C_GRID, lw=0.6)
    ax = axs[1]
    report = {}
    for lab, col, ms, cs in cases:
        w = cs ** 2 / (cs ** 2).sum()
        E = partial(ms, w, x)
        ax.plot(x, E, color=col, lw=1.6, label=lab)
        report[lab] = (E.min(), partial(ms, w, np.array([180.0]))[0])
    tri = 1 - 2 * np.abs(np.radians(x)) / np.pi
    ax.plot(x, tri, color=C_INK, lw=0.9, ls=(0, (2, 2)), label=L("三角波 1 − 2|θ|/π（矩形波の M → ∞ の極限、参考）", "triangle wave 1 − 2|θ|/π (square wave as M → ∞, reference)"))
    ax.axhline(0, color=C_GRID, lw=0.8)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xlabel(L("位相差 θ [°]", "phase difference θ [deg]"))
    ax.set_ylabel(L("E(θ) ＝ Σ w_m cos mθ、w_m ＝ c_m²/Σc_m²", "E(θ) = Σ w_m cos mθ,  w_m = c_m²/Σc_m²"))
    ax.set_title(L("(b) 自己相関（同じ軸、有限和、N ＝ %d）" % N, "(b) autocorrelation (same axis, finite sums, N = %d)" % N), loc="left", fontsize=11)
    ax.legend(fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1, frameon=False)
    ax.grid(color=C_GRID, lw=0.6)
    fig.suptitle(L("図 3　倍音（P3）：相関はパルス形の自己相関。単音は cos θ、矩形波（奇数倍音）は三角波に近づき、等振幅の局在波は θ ＝ 0 の峰に集まる",
                   "Fig. 3  Harmonics (P3): the correlation is the pulse's autocorrelation (single tone, square wave, equal-amplitude localized wave; N = 13)"), fontsize=11, color=C_INK, y=1.0)
    fig.tight_layout()
    save(fig, "fig3_harmonics_autocorrelation")
    for lab, (emin, epi) in report.items():
        print("  E min = %.4f, E(180°) = %.4f  <- %s" % (emin, epi, lab))


# ================================================================ 図 4  ツィレルソン限界（3.4）
def fig4():
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.0), gridspec_kw=dict(width_ratios=[1, 1.25]))
    ax = axs[0]
    ax.add_patch(Circle((0, 0), 1, fill=False, color=C_GRID, lw=1.0))
    vecs = {"a": (0.0, C_A), "a'": (90.0, C_A), "b": (45.0, C_B), "b'": (-45.0, C_B)}
    for name, (deg, col) in vecs.items():
        v = np.array([np.cos(np.radians(deg)), np.sin(np.radians(deg))])
        ax.add_patch(FancyArrowPatch((0, 0), v, arrowstyle="-|>", mutation_scale=14, color=col, lw=1.8))
        ax.text(1.13 * v[0], 1.13 * v[1], r"$\mathbf{%s}$" % name.replace("'", r"'"), ha="center", va="center", color=col, fontsize=12)
    b, bp = np.array([np.cos(np.radians(45)), np.sin(np.radians(45))]), np.array([np.cos(np.radians(-45)), np.sin(np.radians(-45))])
    for v, lab, dx in ((b + bp, r"$\mathbf{b}+\mathbf{b}'$", 0.0), (b - bp, r"$\mathbf{b}-\mathbf{b}'$", 0.0)):
        ax.add_patch(FancyArrowPatch((0, 0), v, arrowstyle="-|>", mutation_scale=14, color=C_C, lw=1.6, ls=(0, (3, 2))))
        ax.text(v[0] + 0.1, v[1] + 0.1, lab + r"  $(|\cdot|=\sqrt{2})$", ha="left" if v[0] > 0.5 else "center", va="bottom", color=C_C, fontsize=10)
    ax.text(-1.45, -1.25, L("E ＝ −a·b。S ＝ a·(b＋b′) ＋ a′·(b−b′) ＝ √2 ＋ √2 ＝ 2√2\n等号：b ⊥ b′、a ∥ b＋b′、a′ ∥ b−b′（nα ＝ 0, nα′ ＝ 90°, nβ ＝ ±45°）",
                            "E = −a·b.  S = a·(b+b′) + a′·(b−b′) = √2 + √2 = 2√2\nequality: b ⊥ b′, a ∥ b+b′, a′ ∥ b−b′ (nα = 0, nα′ = 90°, nβ = ±45°)"),
            ha="left", va="top", fontsize=9, color=C_INK)
    tidy(ax, (-1.5, 1.9), (-1.85, 1.95), L("(a) 読出し平面の単位ベクトル（n 倍した角）", "(a) unit vectors in the readout plane (angles times n)"))
    ax = axs[1]
    d = np.linspace(0, 90, 1801)
    for n, col in ((1, C_B), (2, C_C), (4, C_A)):
        S = np.abs(3 * np.cos(np.radians(n * d)) - np.cos(np.radians(3 * n * d)))
        ax.plot(d, S, color=col, lw=1.6, label=L("n ＝ %d（最大は δ ＝ %g°）" % (n, 45 / n), "n = %d (max at δ = %g°)" % (n, 45 / n)))
    ax.axhline(2, color=C_INK, lw=1.0, ls=(0, (4, 2)))
    ax.text(89, 2.03, L("ベル（局所隠れ変数）2", "Bell (local hidden variables) 2"), ha="right", va="bottom", fontsize=9, color=C_INK)
    ax.axhline(2 * np.sqrt(2), color=C_INK, lw=1.0, ls=(0, (4, 2)))
    ax.text(89, 2 * np.sqrt(2) + 0.03, L("ツィレルソン 2√2 ≈ 2.828", "Tsirelson 2√2 ≈ 2.828"), ha="right", va="bottom", fontsize=9, color=C_INK)
    ax.axhline(4, color=C_MUTED, lw=0.9, ls=(0, (1, 2)))
    ax.text(89, 4.03, L("PR 箱 4（内積では表せない）", "PR box 4 (not an inner product)"), ha="right", va="bottom", fontsize=9, color=C_MUTED)
    ax.set_xlim(0, 90)
    ax.set_ylim(0, 4.4)
    ax.set_xlabel(L("設定角の刻み δ [°]（α ＝ 0, α′ ＝ 2δ, β ＝ δ, β′ ＝ −δ）", "setting step δ [deg] (α = 0, α′ = 2δ, β = δ, β′ = −δ)"))
    ax.set_ylabel(L("CHSH 量 S(δ) ＝ |3cos nδ − cos 3nδ|", "CHSH S(δ) = |3cos nδ − cos 3nδ|"))
    ax.set_title(L("(b) S は 2√2 を超えず、最大の刻みは 45°/n", "(b) S never exceeds 2√2; the maximizing step is 45°/n"), loc="left", fontsize=11)
    ax.legend(fontsize=9, loc="upper left", bbox_to_anchor=(0.01, 0.90))
    ax.grid(color=C_GRID, lw=0.6)
    fig.suptitle(L("図 4　ツィレルソン限界（3.4）：相関が長さ 1 の向きの内積なら、CHSH 量は 2√2 を超えない",
                   "Fig. 4  Tsirelson bound (3.4): if the correlation is an inner product of unit vectors, CHSH never exceeds 2√2"), fontsize=12, color=C_INK, y=1.0)
    fig.tight_layout()
    save(fig, "fig4_tsirelson_bound")


# ================================================================ 図 5  二者相関 E ＝ −cos nθ と最適角（§4）
def fig5():
    fig, ax = plt.subplots(figsize=(11.5, 4.4))
    th = np.linspace(0, 360, 2881)
    for n, col, lab in ((1, C_B, L("スピン½ n ＝ 1（周期 360°）", "spin-½  n = 1 (period 360°)")), (2, C_C, L("光子 n ＝ 2（周期 180°）", "photon  n = 2 (period 180°)")), (4, C_A, L("重力波 n ＝ 4（周期 90°）", "gravitational wave  n = 4 (period 90°)"))):
        E = -np.cos(np.radians(n * th))
        ax.plot(th, E, color=col, lw=1.6, label=lab)
        for k, mk in ((45.0 / n, "o"), (135.0 / n, "s")):
            ax.plot([k], [-np.cos(np.radians(n * k))], marker=mk, color=col, ms=6, mec=C_INK, mew=0.6, zorder=5)
    ax.axhline(0, color=C_GRID, lw=0.8)
    ax.set_xlim(0, 360)
    ax.set_ylim(-1.15, 1.15)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.set_xlabel(L("相対角 θ ＝ α − β [°]", "relative angle θ = α − β [deg]"))
    ax.set_ylabel(L("E(θ) ＝ −cos nθ", "E(θ) = −cos nθ"))
    ax.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=False)
    ax.grid(color=C_GRID, lw=0.6)
    ax.text(182, 1.02, L("○：最大の破れを与える差 |α−β| ＝ π/(4n) ＝ 45°（n＝1）、22.5°（n＝2）、11.25°（n＝4）（三つの対）　□：3π/(4n)（残りの一つ）",
                         "○: |α−β| = π/(4n) = 45° (n=1), 22.5° (n=2), 11.25° (n=4) at maximal violation (three pairs)   □: 3π/(4n) (the remaining pair)"), ha="center", va="bottom", fontsize=9, color=C_INK)
    ax.set_title(L("図 5　二者相関 E(θ) ＝ −cos nθ（§4）：上限は n によらず 2√2、最適角の刻みは 1/n で縮む",
                   "Fig. 5  Two-party correlation E(θ) = −cos nθ (§4): the bound 2√2 is n-independent; the optimal step shrinks as 1/n"), fontsize=12, color=C_INK, loc="left")
    fig.tight_layout()
    save(fig, "fig5_correlation_cos_n_theta")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5()
