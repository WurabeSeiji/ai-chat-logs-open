# -*- coding: utf-8 -*-
"""概念図 v0.1：次数を逆（共役）に読む。
すべて閉じた円周の上の整数次の余弦の有限和。描いた式は README を参照。
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Hiragino Sans"
plt.rcParams["axes.unicode_minus"] = False

EN = len(sys.argv) > 1 and sys.argv[1] == "en"   # 英語版の図（文字だけを替える）
def L(ja, en):
    return en if EN else ja

M = 10            # 最初の読みの最大次数（円周を M 分割）
K = 10            # 逆の読みで大きい円周に入る区画の数
CELL = 360.0 / K  # 一区画の角
A_BASE = 10.0     # (e) で足す最低位のベース波の振幅（結晶の歯の高さ M に合わせた）
psi = np.linspace(0.0, 360.0, 14401)
rad = np.deg2rad

C_READ = "#1f6fb2"   # 読める側
C_INT = "#555555"    # 内部（結晶）
C_RES = "#c0392b"    # 分解能

fig = plt.figure(figsize=(10.5, 14.5))
gs = fig.add_gridspec(6, 1, height_ratios=[1.5, 1.0, 1.0, 1.0, 0.8, 0.8], hspace=0.62)

# ---------- (a) 次数の梯子 ----------
ax = fig.add_subplot(gs[0])
ax.set_xscale("log")
ax.set_xlim(30, 0.07)
ax.set_ylim(-0.1, 2.5)
ax.set_yticks([])
ax.set_xticks([10, 5, 2, 1, 0.5, 0.2, 0.1])
ax.set_xticklabels(["10", "5", "2", "1", "1/2", "1/5", "1/10"])
ax.set_xlabel(L("波長（ベース波に見える波の波長を 1 とする）　　← 粗い　　　細かい →", "wavelength (the wavelength of the apparent base wave is 1)        ← coarse          fine →"))
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
# 上段：最初の読み
y0 = 1.45
ax.plot([13, 0.07], [y0, y0], color="#bbbbbb", lw=0.8)
for m in range(1, M + 1):
    ax.plot([1.0 / m] * 2, [y0, y0 + 0.5], color=C_READ, lw=2)
ax.plot([1.0 / M] * 2, [y0, y0 + 0.62], color=C_RES, lw=3)
ax.text(1.0, y0 + 0.56, L("ベース波（1 次）", "base wave (order 1)"), ha="center", va="bottom", fontsize=9, color=C_READ)
ax.text(1.0 / M, y0 + 0.68, L("分解能（10 次）", "resolution (order 10)"), ha="center", va="bottom", fontsize=9, color=C_RES)
ax.text(0.32, y0 + 0.62, L("読める：次数 1〜10", "readable: orders 1–10"), ha="center", va="bottom", fontsize=9, color=C_READ)
ax.text(28, y0 + 0.22, L("最初の読み", "first\nreading"), ha="left", va="center", fontsize=11, weight="bold")
ax.text(4.2, y0 + 0.22, L("（何もない）", "(nothing)"), ha="center", va="center", fontsize=9, color="#999999")
# 下段：逆に読む
y1 = 0.15
ax.plot([13, 0.07], [y1, y1], color="#bbbbbb", lw=0.8)
for n in range(1, K):
    ax.plot([float(K) / n] * 2, [y1, y1 + 0.5], color=C_READ, lw=2)
for m in range(2, M + 1):
    ax.plot([1.0 / m] * 2, [y1, y1 + 0.5], color=C_INT, lw=2)
ax.plot([1.0, 1.0], [y1, y1 + 0.62], color=C_RES, lw=3)
ax.annotate("", xy=(0.072, y1 + 0.25), xytext=(0.095, y1 + 0.25),
            arrowprops=dict(arrowstyle="->", color=C_INT, lw=1.5))
ax.text(1.0, y1 + 0.68, L("観測者の分解能\n（ベース波に見える波）", "observer's resolution\n(apparent base wave)"), ha="center", va="bottom", fontsize=9, color=C_RES)
ax.text(4.6, y1 + 0.62, L("読める：分解能より粗い波の中心位相", "readable: centre phase of\nwaves coarser than the resolution"), ha="center", va="bottom", fontsize=9, color=C_READ)
ax.text(0.27, y1 + 0.62, L("内部：硬い結晶（次数は無限）", "internal: rigid crystal\n(orders without end)"), ha="center", va="bottom", fontsize=9, color=C_INT)
ax.text(28, y1 + 0.22, L("逆に読む", "reverse\nreading"), ha="left", va="center", fontsize=11, weight="bold")
ax.set_title(L("(a) 次数の梯子：同じ構造の二つの読み", "(a) Ladder of orders: two readings of the same structure"), loc="left", fontsize=11)

def deg_axis(ax, cells=True, label=False):
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, CELL))
    if cells:
        for k in range(K + 1):
            ax.axvline(k * CELL, color="#cccccc", lw=0.7, zorder=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if label:
        ax.set_xlabel(L("円周の上の角 ψ（度）", "angle ψ on the circle (degrees)"))

# ---------- (b) 最初の読み ----------
ax = fig.add_subplot(gs[1])
for m in range(1, M + 1):
    ax.plot(psi, np.cos(m * rad(psi - 180.0)), color="#9db8d0", lw=0.6)
f1 = sum(np.cos(m * rad(psi - 180.0)) for m in range(1, M + 1))
ax.plot(psi, f1, color=C_READ, lw=1.8)
deg_axis(ax)
ax.set_ylabel(L("振幅", "amplitude"))
ax.set_title(L("(b) 最初の読み：一つの円周の倍音 1〜10 次（細線）とその和（太線）。10 次が円周を 10 分割する", "(b) First reading: harmonics 1–10 of one circle (thin) and their sum (thick). Order 10 divides the circle into 10"),
             loc="left", fontsize=10)

# ---------- (c) 逆に読む：分解能より上 ----------
ax = fig.add_subplot(gs[2])
g = sum(np.cos(K * j * rad(psi - CELL / 2)) for j in range(1, M + 1))
ax.plot(psi, g, color=C_INT, lw=0.9, label=L("硬い結晶：大きい円周の 10, 20, …, 100 次の和", "rigid crystal: sum of orders 10, 20, …, 100 of the large circle"))
ax.plot(psi, np.cos(K * rad(psi - CELL / 2)), color=C_RES, lw=2.0,
        label=L("分解能でローパス（10 次まで）：ベース波に見える", "low-passed at the resolution (up to order 10): looks like the base wave"))
deg_axis(ax)
ax.set_ylabel(L("振幅", "amplitude"))
ax.legend(loc="upper right", fontsize=8, frameon=False, ncol=2, bbox_to_anchor=(1.0, 1.02))
ax.set_ylim(-3.2, 13.5)
ax.set_title(L("(c) 逆に読む・分解能より上：(b) の構造は一区画の中に入り、10 区画に同じ形が並ぶ", "(c) Reverse reading, above the resolution: the structure of (b) fits in one cell; the same shape fills all 10 cells"), loc="left", fontsize=10)

# ---------- (d) 逆に読む：分解能より下 ----------
ax = fig.add_subplot(gs[3])
thetas = [90.0, 102.0, 114.0, 126.0]
cols = ["#1f6fb2", "#4a9bd6", "#e08a2e", "#b85c00"]
for th, c in zip(thetas, cols):
    ax.plot(psi, np.cos(rad(psi - th)), color=c, lw=1.5, label=L("θ₀＝%d°", "θ₀ = %d°") % th)
    ax.plot([th], [1.12], marker="v", color=c, ms=7)
    k = int(th // CELL)
    ax.axvspan(k * CELL, (k + 1) * CELL, color=c, alpha=0.07, lw=0)
deg_axis(ax)
ax.legend(loc="upper right", fontsize=8, frameon=False, ncol=4)
ax.set_ylim(-1.15, 1.75)
ax.set_ylabel(L("振幅", "amplitude"))
ax.set_title(L("(d) 逆に読む・分解能より下：最も低位のベース波。中心位相 θ₀（▼）は連続に動く。読めるのは θ₀ のある区画だけ", "(d) Reverse reading, below the resolution: the coarsest base wave. Its centre phase θ₀ (▼) moves\ncontinuously; only the cell that contains θ₀ is readable"),
             loc="left", fontsize=10)

# ---------- (e) 重ね合わせのパワー ----------
for i, (th, c) in enumerate([(102.0, cols[1]), (114.0, cols[2])]):
    ax = fig.add_subplot(gs[4 + i])
    p = (g + A_BASE * np.cos(rad(psi - th))) ** 2
    ax.plot(psi, p, color=c, lw=0.9)
    imax = int(np.argmax(p))
    ax.plot([psi[imax]], [p[imax]], marker="v", color=C_RES, ms=8)
    k = int(th // CELL)
    ax.axvspan(k * CELL, (k + 1) * CELL, color=c, alpha=0.10, lw=0)
    ax.axvline(th, color=c, lw=1.0, ls="--")
    deg_axis(ax, label=(i == 1))
    ax.set_ylabel(L("パワー", "power"))
    ax.set_ylim(0, 470)
    if i == 0:
        ax.set_title(L("(e) パワー（硬い結晶＋最も低位のベース波）²：θ₀ を 102° から 114° へ動かすと、最大の歯（▼）は隣の区画へ一段で移る", "(e) Power (rigid crystal + coarsest base wave)²: moving θ₀ from 102° to 114° shifts the largest tooth (▼)\nto the next cell in one step"),
                     loc="left", fontsize=10)
    ax.text(358, 420, L("θ₀＝%d°（破線）", "θ₀ = %d° (dashed)") % th, ha="right", va="top", fontsize=9, color=c)

fig.savefig("figs/fig_reverse_reading_v0.1%s.png" % ("_en" if EN else ""), dpi=150, bbox_inches="tight")
print("saved")
# 確認の数値
for th in thetas:
    p = (g + A_BASE * np.cos(rad(psi - th))) ** 2
    print("theta0=%.0f  max tooth at psi=%.2f  cell=%d  cell_of_theta0=%d"
          % (th, psi[int(np.argmax(p))], int(psi[int(np.argmax(p))] // CELL), int(th // CELL)))
