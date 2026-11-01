# -*- coding: utf-8 -*-
"""EPR パラドックスはパラドックスではない（補遺）― 論文用の静的な図（PNG・SVG・PDF）。
図 1  閉包の対：q_B = −q_A は対蹠、位置と速度が同時に逆向き。時間を進めても対蹠。A を軸に落としても B は変わらない
図 2  順序への非依存：S から A・B への null の関係、A と B の読出し事象は空間的に離れ、枠によって順序が入れ替わる（ローレンツ変換）
点はすべて式から直接計算する。乱数は使わない。
実行:  python make_figures.py [en]     出力: figures/figN_*.png/.svg/.pdf（英語版は _en）"""
from pathlib import Path
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
EN = len(sys.argv) > 1 and sys.argv[1] == "en"; L = (lambda ja, en: en if EN else ja); SUF = "_en" if EN else ""
plt.rcParams["font.family"] = ["Hiragino Sans", "DejaVu Sans"]; plt.rcParams["pdf.fonttype"] = 42; plt.rcParams["svg.fonttype"] = "none"; plt.rcParams["axes.unicode_minus"] = False; plt.rcParams["font.size"] = 10
OUT = Path(__file__).resolve().parent / "figures"; OUT.mkdir(exist_ok=True)
C_A, C_B, C_C, C_INK, C_MUTED, C_GRID, C_LINE = "#c4470f", "#2a56c6", "#0e7a6e", "#18212c", "#596577", "#c9d0d8", "#7a8699"
def save(fig, name):
    for ext in ("png", "svg", "pdf"): fig.savefig(OUT / f"{name}{SUF}.{ext}", dpi=200 if ext == "png" else None, bbox_inches="tight", facecolor="white")
    plt.close(fig); print("saved", name + SUF)
def arrow(ax, p, q, col, lw=1.8, ls="-"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=13, color=col, lw=lw, linestyle=ls, zorder=4))
def tidy(ax, xl, yl, title):
    ax.set_aspect("equal"); ax.set_xlim(*xl); ax.set_ylim(*yl); ax.axis("off"); ax.set_title(title, fontsize=11, loc="left", color=C_INK)

def fig1():
    A, w = 1.0, 1.0
    fig, axs = plt.subplots(1, 3, figsize=(12.0, 4.3))
    # (a)
    ax = axs[0]; th = 0.6
    ax.add_patch(Circle((0, 0), A, fill=False, color=C_GRID, lw=1.0))
    xa = np.array([np.cos(th), np.sin(th)]) * A; xb = -xa; va = w * np.array([-np.sin(th), np.cos(th)]) * A; vb = -va
    ax.plot([0], [0], "o", color=C_INK, ms=5); ax.text(0.0, -0.17, "S", ha="center", va="top", fontsize=11)
    arrow(ax, (0, 0), xa, C_A); arrow(ax, (0, 0), xb, C_B)
    ax.text(*(xa * 1.22), r"$\boldsymbol{x}_A$", color=C_A, ha="center", va="center", fontsize=11); ax.text(*(xb * 1.22), r"$\boldsymbol{x}_B=-\boldsymbol{x}_A$", color=C_B, ha="center", va="center", fontsize=11)
    arrow(ax, xa, xa + 0.55 * va, C_A, ls=(0, (3, 2))); arrow(ax, xb, xb + 0.55 * vb, C_B, ls=(0, (3, 2)))
    ax.text(*(xa + 0.75 * va), r"$\boldsymbol{v}_A$", color=C_A, fontsize=10, ha="center"); ax.text(*(xb + 0.75 * vb), r"$\boldsymbol{v}_B=-\boldsymbol{v}_A$", color=C_B, fontsize=10, ha="center")
    ax.plot(*xa, "o", color=C_A, ms=7); ax.plot(*xb, "o", color=C_B, ms=7)
    ax.text(0, -1.55, L("閉包 $q_A+q_B=0$ ⇒ 位置の和と速度の和が同時にゼロ", "closure $q_A+q_B=0$ ⇒ positions and velocities sum to zero at once"), ha="center", fontsize=9.5)
    tidy(ax, (-1.75, 1.75), (-1.8, 1.5), L("(a) 対蹠の対：$\\boldsymbol{x}_B=-\\boldsymbol{x}_A$、$\\boldsymbol{p}_B=-\\boldsymbol{p}_A$", "(a) antipodal pair: $\\boldsymbol{x}_B=-\\boldsymbol{x}_A$, $\\boldsymbol{p}_B=-\\boldsymbol{p}_A$"))
    # (b)
    ax = axs[1]
    ax.add_patch(Circle((0, 0), A, fill=False, color=C_GRID, lw=1.0)); ax.plot([0], [0], "o", color=C_INK, ms=5)
    for k, tt in enumerate((0.0, 0.9, 1.8)):
        ang = th + w * tt; xa = A * np.array([np.cos(ang), np.sin(ang)]); xb = -xa; al = 0.35 + 0.3 * k
        ax.plot(*xa, "o", color=C_A, ms=7, alpha=al); ax.plot(*xb, "o", color=C_B, ms=7, alpha=al)
        ax.plot([xa[0], xb[0]], [xa[1], xb[1]], color=C_LINE, lw=1.0, alpha=al)
        ax.text(*(xa * 1.2), "t%d" % k, color=C_A, fontsize=9, ha="center", alpha=al); ax.text(*(xb * 1.2), "t%d" % k, color=C_B, fontsize=9, ha="center", alpha=al)
    ax.text(0, -1.55, L("同じ角速度で回り、常に対蹠（$q_B=-q_A$ は全時刻で成立）", "both turn at the same rate and stay antipodal ($q_B=-q_A$ at all times)"), ha="center", fontsize=9.5)
    tidy(ax, (-1.75, 1.75), (-1.8, 1.5), L("(b) 時間を進めても閉包は保たれる", "(b) the closure is kept as time runs"))
    # (c)
    ax = axs[2]; ang = th; xa = A * np.array([np.cos(ang), np.sin(ang)]); xb = -xa; nA = np.array([1.0, 0.0])
    ax.add_patch(Circle((0, 0), A, fill=False, color=C_GRID, lw=1.0)); ax.plot([0], [0], "o", color=C_INK, ms=5)
    ax.plot([-1.5, 1.5], [0, 0], color=C_INK, lw=1.0, ls=(0, (4, 2))); ax.text(1.55, 0.0, r"$\boldsymbol{n}_A$", va="center", fontsize=11)
    ax.plot(*xa, "o", color=C_A, ms=7, alpha=0.4); arrow(ax, xa, (np.dot(xa, nA), 0.0), C_A, ls=(0, (2, 2)), lw=1.2)
    ax.plot([np.dot(xa, nA)], [0], "s", color=C_A, ms=8); ax.text(np.dot(xa, nA), 0.17, L("A：＋に落ちる", "A: falls to +"), color=C_A, ha="center", fontsize=9.5)
    ax.plot(*xb, "o", color=C_B, ms=7); ax.text(*(xb * 1.25), L("B：何も変わらない", "B: nothing changes"), color=C_B, ha="center", fontsize=9.5)
    ax.text(0, -1.55, L("A を軸に対して破壊的に読む（P4、既知）。B の状態には作用しない", "A is read destructively against its axis (P4, known); B is untouched"), ha="center", fontsize=9.5)
    tidy(ax, (-1.75, 1.95), (-1.8, 1.5), L("(c) A の読出しは B に出来事を起こさない", "(c) reading A is no event at B"))
    fig.suptitle(L("図 1　閉包の対：一つの閉包が位置と運動量を同時に決め、読出しは一端を読むだけ", "Fig. 1  The closed pair: one closure fixes position and momentum at once; a readout reads one end"), fontsize=12, color=C_INK, y=1.0)
    fig.tight_layout(); save(fig, "fig1_closure_pair")

def fig2():
    # 事象：S=(0,0) から null 線で A=(−d, d), B=(d, d)。読出し事象 EA=(−d, tA), EB=(d, tB)。枠 2 は速度 β で boost。
    d, tA, tB = 1.0, 1.4, 1.7; beta = 0.6; g = 1 / np.sqrt(1 - beta ** 2)
    def boost(x, t): return g * (x - beta * t), g * (t - beta * x)
    fig, axs = plt.subplots(1, 2, figsize=(11.0, 4.6))
    for ax, bb, ttl in ((axs[0], 0.0, L("(a) 枠 1：A の読出しが先（t_A < t_B）", "(a) frame 1: A is read first ($t_A<t_B$)")), (axs[1], beta, L("(b) 枠 2（速度 0.6c）：B の読出しが先", "(b) frame 2 (speed 0.6c): B is read first"))):
        def T(x, t):
            if bb == 0: return x, t
            return boost(x, t)
        pts = {"S": T(0, 0), "A": T(-d, d), "B": T(d, d), "EA": T(-d, tA), "EB": T(d, tB)}
        for a, b_ in (("S", "A"), ("S", "B")):
            ax.plot([pts[a][0], pts[b_][0]], [pts[a][1], pts[b_][1]], color=C_C, lw=1.6)
        for a, b_, col in (("A", "EA", C_A), ("B", "EB", C_B)):
            ax.plot([pts[a][0], pts[b_][0]], [pts[a][1], pts[b_][1]], color=col, lw=1.6)
        ax.plot(*pts["S"], "o", color=C_INK, ms=6); ax.text(pts["S"][0], pts["S"][1] - 0.15, "S", ha="center", va="top", fontsize=11)
        ax.plot(*pts["EA"], "s", color=C_A, ms=8); ax.plot(*pts["EB"], "s", color=C_B, ms=8)
        ax.text(pts["EA"][0] - 0.1, pts["EA"][1], L("A の読出し", "readout at A"), color=C_A, ha="right", va="center", fontsize=10)
        ax.text(pts["EB"][0] + 0.1, pts["EB"][1], L("B の読出し", "readout at B"), color=C_B, ha="left", va="center", fontsize=10)
        for name, col in (("EA", C_A), ("EB", C_B)):
            ax.axhline(pts[name][1], color=col, lw=0.8, ls=(0, (2, 3)))
        ax.text(-2.4, pts["EA"][1], "t = %.2f" % pts["EA"][1], color=C_A, fontsize=9, va="bottom"); ax.text(-2.4, pts["EB"][1], "t = %.2f" % pts["EB"][1], color=C_B, fontsize=9, va="bottom")
        ax.text(0.0, -0.75, L("S–A、S–B は時空間隔ゼロの一回の関係（P2）\n二つの読出しは空間的に離れ、順序は枠で変わる。記録の突き合わせは同じ", "S–A and S–B are one relation of zero interval (P2)\nthe two readouts are spacelike; their order depends on the frame, the compared records do not"), ha="center", va="top", fontsize=9.5, color=C_INK)
        ax.set_xlim(-2.6, 2.6); ax.set_ylim(-1.5, 3.4); ax.set_xlabel("x"); ax.set_ylabel("t"); ax.set_title(ttl, loc="left", fontsize=11); ax.grid(color=C_GRID, lw=0.6); ax.set_aspect("equal")
    fig.suptitle(L("図 2　順序への非依存：A が先か B が先かは枠に依存し、相関は依存しない", "Fig. 2  Order independence: which readout comes first depends on the frame; the correlation does not"), fontsize=12, color=C_INK, y=1.0)
    fig.tight_layout(); save(fig, "fig2_order_independence")

if __name__ == "__main__":
    fig1(); fig2()
