#!/usr/bin/env python3
"""B中心読み直しの三円図 v1（加速度論文 v5 用の作図案）

親閉鎖の円（半径R）の同一円周上にある A・B 二体について、
未来位相位置中心の回転の代わりに、B を中心とする
  R'  : B の固有半径（子閉鎖）
  R'' : B から位相差θにある A までの関係距離（弦）
の二つの円を重ね、三つの円（R, R', R''）で幾何を読む。
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Arc
import matplotlib

matplotlib.rcParams["font.family"] = "Hiragino Sans"

BG = "#faf7f0"
INK = "#2b3a42"
GRAY = "#8ba0ad"
TEAL = "#106e6e"
ORANGE = "#c05b3c"
BLUE = "#3a6ea5"

R = 1.0
theta = np.deg2rad(70.0)          # 位相差
angB = np.deg2rad(-55.0)          # Bの位置角
angA = angB - theta               # Aの位置角
O = np.array([0.0, 0.0])
A = R * np.array([np.cos(angA), np.sin(angA)])
B = R * np.array([np.cos(angB), np.sin(angB)])
Rpp = float(np.linalg.norm(A - B))   # R'' = 2R sin(θ/2)
Rp = 0.24                            # R'（子閉鎖・図示用）

fig, ax = plt.subplots(figsize=(12.5, 11), facecolor=BG)
ax.set_facecolor(BG)

# 親円
ax.add_patch(Circle(O, R, fill=False, ec=GRAY, lw=2.2))
# B中心 R'' 円（Aを厳密に通る）
ax.add_patch(Circle(B, Rpp, fill=False, ec=BLUE, lw=1.8, ls=(0, (6, 4))))
# B中心 R' 円（子閉鎖）
ax.add_patch(Circle(B, Rp, fill=False, ec=ORANGE, lw=2.0))

# 中心と半径線
ax.plot(*O, "o", color=INK, ms=5)
ax.annotate("O（親閉鎖の中心）", O, xytext=(0.03, 0.06), textcoords="data", fontsize=12, color=INK)
angR = np.deg2rad(148)
Rend = R * np.array([np.cos(angR), np.sin(angR)])
ax.plot([O[0], Rend[0]], [O[1], Rend[1]], color=GRAY, lw=1.4)
ax.annotate("R", 0.55 * Rend, fontsize=15, color=INK, ha="center",
            xytext=0.55 * Rend + np.array([-0.05, 0.05]), textcoords="data")

# O→A, O→B の関係線と位相差θ
for P in (A, B):
    ax.plot([O[0], P[0]], [O[1], P[1]], color=GRAY, lw=1.0, ls=":")
arc_r = 0.30
ax.add_patch(Arc(O, 2*arc_r, 2*arc_r, angle=0,
                 theta1=np.rad2deg(angA), theta2=np.rad2deg(angB), color=INK, lw=1.4))
mid = (angA + angB) / 2
ax.annotate("θ（位相差）", (arc_r + 0.13) * np.array([np.cos(mid), np.sin(mid)]) + np.array([0.16, 0.0]),
            fontsize=13, color=INK, ha="center")

# A, B 点
ax.plot(*A, "o", color=TEAL, ms=13)
ax.plot(*B, "o", color=ORANGE, ms=13)
ax.annotate("A", A + np.array([-0.11, -0.08]), fontsize=17, color=TEAL, weight="bold")
ax.annotate("B", B + np.array([0.06, -0.10]), fontsize=17, color=ORANGE, weight="bold")

# B→A 弦 = R''
ax.plot([B[0], A[0]], [B[1], A[1]], color=BLUE, lw=2.2)
mAB = 0.5 * (A + B)
ax.annotate("R″ = 2R sin(θ/2)", mAB + np.array([0.05, 0.07]), fontsize=14, color=BLUE, ha="center")

# R' ラベル（Bから子閉鎖円へ）
dirp = np.array([np.cos(np.deg2rad(20)), np.sin(np.deg2rad(20))])
ax.plot([B[0], B[0] + Rp * dirp[0]], [B[1], B[1] + Rp * dirp[1]], color=ORANGE, lw=1.6)
ax.annotate("R′", B + (Rp + 0.07) * dirp, fontsize=15, color=ORANGE, ha="center")

# 加速度矢印：Aの向心加速度（→O）と、Bが読む方向（A→B）
def arrow(p, q, color, lw=2.6):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=22,
                                 color=color, lw=lw, shrinkA=0, shrinkB=0))

arrow(A, A + 0.34 * (O - A) / np.linalg.norm(O - A), INK)
ax.annotate("向心加速度 a = Rω²\n（親閉鎖上の厳密な恒等式）",
            A + 0.46 * (O - A) / np.linalg.norm(O - A) + np.array([-0.72, 0.05]),
            fontsize=11.5, color=INK)

uAB = (B - A) / np.linalg.norm(B - A)
arrow(A, A + 0.30 * uAB, BLUE)
ax.annotate("B が読む加速度成分（等価原理で重力加速度）",
            A + np.array([-0.28, -0.22]), fontsize=11.5, color=BLUE)

# 説明ボックス
box1 = ("A² + B² + (i R)² = 0（親閉鎖）\n"
        "R″ = 2R sin(θ/2)：B から A への関係距離\n"
        "R′：B の固有半径（B 自身が子閉鎖 Σx_m² = R_m² を私有）")
ax.text(-1.58, 1.30, box1, fontsize=12.5, color=INK, va="top",
        bbox=dict(boxstyle="round,pad=0.55", fc="white", ec=GRAY, lw=1.2))

box2 = ("作図規則\n"
        "1. 三つの円はすべて閉鎖から読まれる（背景座標なし）\n"
        "2. 回転中心は未来位相位置ではなく、親の中心 O と\n"
        "    同一円周上の実在 B のみ\n"
        "3. B の読みは R′ と R″ の比だけで書く（比だけが物理）")
ax.text(0.42, 1.30, box2, fontsize=11.5, color=INK, va="top",
        bbox=dict(boxstyle="round,pad=0.55", fc="#eef3f7", ec=GRAY, lw=1.2))

box3 = ("図の主張\n"
        "同一円周上の A は、B から見ると距離 R″ に静止して\n"
        "共回転している。その向心加速度を B の固有スケール R′\n"
        "で読み直すと、遠心力と釣り合う重力加速度が R′/R″ の\n"
        "比から出る——等価原理により、B の内部観測者には\n"
        "これが重力として読まれる。")
ax.text(0.42, -1.02, box3, fontsize=11.5, color=INK, va="top",
        bbox=dict(boxstyle="round,pad=0.55", fc="white", ec=GRAY, lw=1.2))

ax.set_title("AB二体閉鎖における B 中心読み直しの幾何図（三つの円 R・R′・R″）",
             fontsize=17, color=INK, pad=18)
ax.text(0, 1.52, "未来位相位置中心の回転を、同一円周上の B を中心とする R′・R″ の関係読みに置き換える",
        fontsize=12.5, color=GRAY, ha="center")

ax.set_xlim(-1.65, 1.95)
ax.set_ylim(-1.75, 1.62)
ax.set_aspect("equal")
ax.axis("off")
fig.savefig("fig_b_centered_three_circles_v1.png", dpi=160, bbox_inches="tight", facecolor=BG)
print("saved; R''=", round(Rpp, 6), " (=2R sin(θ/2) =", round(2*R*np.sin(theta/2), 6), ")")
