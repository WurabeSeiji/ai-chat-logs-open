# -*- coding: utf-8 -*-
"""note 記事 第 1 節の図：二つの閉塞を満たす点の集まりと、満たさない点の集まり（SVG と PNG）。"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sys
EN = len(sys.argv) > 1 and sys.argv[1] == "en"
L = (lambda ja, en: en if EN else ja)

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["svg.fonttype"] = "none"
OUT = Path(__file__).resolve().parent

z = np.exp(2j * np.pi * np.arange(5) / 5)          # 1 の 5 乗根：Σz＝0、Σz²＝0
zp = z ** 2                                         # 二乗の層 z′＝z²
w = np.exp(1j * np.deg2rad([10, 30, 50, 70, 90]))  # 片寄った集まり：閉塞を満たさない

panels = [(z, L("(a) z の層", "(a) layer z"), "#1f5fa8"), (zp, L("(b) 二乗の層 z′＝z²", "(b) squared layer z′ = z²"), "#b5532a"), (w, L("(c) 片寄った集まり", "(c) a lopsided set"), "#666666")]

fig, axes = plt.subplots(1, 3, figsize=(12, 4.6))
for ax, (pts, title, col) in zip(axes, panels):
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color="#cccccc", lw=0.8, zorder=0)
    ax.axhline(0, color="#999999", lw=0.6, zorder=0)
    ax.axvline(0, color="#999999", lw=0.6, zorder=0)
    for k, p in enumerate(pts):
        ax.plot([0, p.real], [0, p.imag], color=col, lw=0.6, alpha=0.4, zorder=1)
        ax.scatter(p.real, p.imag, s=70, color=col, zorder=3)
        ax.annotate(f"{k + 1}", (p.real, p.imag), xytext=(1.16 * p.real, 1.16 * p.imag),
                    ha="center", va="center", fontsize=10, color=col)
    c = pts.mean()
    ax.scatter(c.real, c.imag, marker="x", s=110, color="black", lw=2, zorder=4)
    ax.annotate(L("重心", "centroid"), (c.real, c.imag), xytext=(c.real - 0.42, c.imag - 0.2), fontsize=10)
    s1, s2 = pts.sum(), (pts ** 2).sum()
    a, b = pts.real, pts.imag
    ax.set_title(title, fontsize=13)
    ax.text(0, -1.62,
            L(f"Σz ＝ {abs(s1):.2f}　 Σz² ＝ {abs(s2):.2f}\n縦と横の広がり Σa²：Σb² ＝ {np.sum(a**2):.2f}：{np.sum(b**2):.2f}\n縦と横の相関 Σab ＝ {np.sum(a*b):.2f}",
              f"|Σz| = {abs(s1):.2f}    |Σz²| = {abs(s2):.2f}\nspread  Σa² : Σb² = {np.sum(a**2):.2f} : {np.sum(b**2):.2f}\ncorrelation  Σab = {np.sum(a*b):.2f}"),
            ha="center", va="top", fontsize=10, linespacing=1.6)
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(1.42, 0.06, L("実", "Re"), fontsize=9, color="#777777", ha="right")
    ax.text(0.06, 1.38, L("虚", "Im"), fontsize=9, color="#777777")

fig.text(0.5, 0.005, L("(a) と (b) は、点の番号の並びが変わるだけで、同じ正五角形になる。(c) は重心が原点から外れ、Σz も Σz² も 0 にならない。",
         "(a) and (b) are the same regular pentagon; only the numbering of the points changes. In (c) the centroid is off the origin and neither Σz nor Σz² is 0."),
         ha="center", fontsize=10.5)
fig.subplots_adjust(left=0.02, right=0.98, top=0.9, bottom=0.25, wspace=0.08)
for ext in ("svg", "png"):
    fig.savefig(OUT / f"fig1_two_closures_one_definition{'_en' if EN else ''}.{ext}", dpi=200, facecolor="white")
print("ok")
