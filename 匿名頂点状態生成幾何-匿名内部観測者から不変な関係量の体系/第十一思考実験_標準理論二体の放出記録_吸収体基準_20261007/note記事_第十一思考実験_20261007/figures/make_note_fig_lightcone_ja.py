#!/usr/bin/env python3
"""note 記事用の図（一般向け、日本語ラベル、模式図でデータなし）。

左：空間で見た絵。中心に二体（陽子と電子）、遠方に放射を完全に吸収する球殻。
    吸収体の静止系で見ると球、二体が動いて見える系の同時刻面で切ると楕円 r(θ) = t₀/(1 − β cos θ)。
右：時空で見た絵。二体の世界線から出た各放出事象が、自分の未来光円錐の表面（null 面）に沿って
    未来の無限遠（吸収体）に着地し、記録（時刻・方向・振動数・偏光）になる。
    過去の無限遠からは背景放射 2.7255 K が入ってくる。
出力：note_fig_lightcone_absorber_ja.png
論文の図 16（figures/make_fig16_lightcone.py）を一般向けに描き直したもの。
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
plt.rcParams.update({"font.family": "Hiragino Sans", "font.size": 11, "figure.dpi": 150,
                     "savefig.dpi": 200, "axes.unicode_minus": False})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 6.4))

# ------------------------------------------------------------------ 左：空間の絵
ax1.set_aspect("equal")
ax1.axis("off")
ax1.set_xlim(-1.5, 2.05)
ax1.set_ylim(-1.66, 1.42)
ax1.set_anchor("N")
ax1.set_title("空間で見る：二体と、無限遠の完全吸収体", fontsize=13, pad=10)

# 吸収体（静止系で球）
ax1.add_patch(Circle((0, 0), 1.0, fill=False, color="C0", lw=3.0))
# 二体が速度 β で動いて見える系の同時刻面で切った断面（楕円、二体は焦点）
beta, t0 = 0.30, 1.0
a = t0 / (1 - beta**2)
b = a * np.sqrt(1 - beta**2)
ax1.add_patch(Ellipse((a * beta, 0), 2 * a, 2 * b, fill=False, color="C3", lw=2.6, ls="--"))

# 二体：陽子と電子
ax1.add_patch(Circle((0, 0), 0.14, fill=False, color="0.4", lw=0.9, ls=":"))
ax1.scatter([0], [0], s=160, color="0.15", zorder=6)
ax1.scatter([0.14 * np.cos(0.9)], [0.14 * np.sin(0.9)], s=40, color="C2", zorder=6)
ax1.text(0.02, -0.21, "陽子", fontsize=10, ha="center", va="top")
ax1.text(0.17, 0.17, "電子", fontsize=10, ha="left", va="bottom", color="C2")

# 出ていく放射
for ang in (0.35, 1.25, 2.3, 3.4, 4.3, 5.3):
    ax1.add_patch(FancyArrowPatch((0.25 * np.cos(ang), 0.25 * np.sin(ang)),
                                  (0.90 * np.cos(ang), 0.90 * np.sin(ang)),
                                  arrowstyle="-|>", mutation_scale=14, lw=1.3, color="0.45"))
ax1.text(-0.97, 0.62, "放射（電磁波・重力波）\nが出ていき、吸収体で\n吸収される", fontsize=10,
         color="0.3", ha="right", va="center")

# 入ってくる背景放射（波線）
for ang in (0.0, 1.9, 4.6):
    s = np.linspace(0.98, 0.42, 60)
    wig = 0.025 * np.sin(2 * np.pi * (s - 0.42) / 0.14)
    x = s * np.cos(ang) - wig * np.sin(ang)
    y = s * np.sin(ang) + wig * np.cos(ang)
    ax1.plot(x, y, color="C1", lw=1.3)
    ax1.add_patch(FancyArrowPatch((x[-2], y[-2]), (x[-1], y[-1]), arrowstyle="-|>",
                                  mutation_scale=12, lw=1.0, color="C1"))
ax1.text(1.27, -0.62, "背景放射（2.7255 K）\nが入ってくる", fontsize=10, color="C1", ha="left", va="center")

ax1.text(0.0, 1.07, "無限遠の完全吸収体\n（吸収体の静止系で見ると球）", fontsize=10.5, color="C0",
         ha="center", va="bottom")
ax1.text(a * beta + 0.0, -b - 0.06, "二体が速度 v で動いて見える系では楕円\n（Doppler 偏移と光行差の幾何）",
         fontsize=10.5, color="C3", ha="center", va="top")
ax1.text(0.2, -1.63, "二体の運動は、放射と吸収の記録だけから読み出す。座標格子や計量は状態に持たない。",
         fontsize=9.5, color="0.3", ha="center", va="bottom")

# ------------------------------------------------------------------ 右：時空の絵
ax2.set_aspect("equal")
ax2.axis("off")
ax2.set_xlim(-0.80, 1.70)
ax2.set_ylim(-1.34, 1.30)
ax2.set_anchor("N")
ax2.set_title("時空で見る：各放出事象は自分の光円錐で吸収体に着地する", fontsize=13, pad=10)

# 世界線
ax2.plot([0, 0], [-1.0, 1.0], color="k", lw=2.6)
ax2.text(-0.03, 1.03, "二体の世界線（水素原子など）", fontsize=10, ha="right", va="bottom")

# 未来の無限遠（吸収体）と過去の無限遠（背景）
ax2.plot([1.0, 0.0], [0.0, 1.0], color="C0", lw=3.2)
ax2.annotate("未来の無限遠 ＝ 吸収体\n（未来光円錐の表面、null 面）", xy=(0.24, 0.76), xytext=(0.50, 0.98),
             fontsize=10.5, color="C0", ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color="C0", lw=0.8))
ax2.plot([1.0, 0.0], [0.0, -1.0], color="C1", lw=3.2)
ax2.annotate("過去の無限遠 ＝ 入ってくる背景放射\n（最終散乱面、2.7255 K）", xy=(0.24, -0.76), xytext=(0.50, -0.98),
             fontsize=10.5, color="C1", ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color="C1", lw=0.8))

# 放出事象と記録
events = ((-0.78, "放出 1", "u₁"), (-0.36, "放出 2", "u₂"), (0.08, "放出 3", "u₃"))
for te, lab, u in events:
    s = (1 - te) / 2
    ax2.plot([0, s], [te, te + s], color="C0", lw=1.1)
    ax2.scatter([0], [te], color="k", s=22, zorder=5)
    ax2.scatter([s], [te + s], color="C0", s=30, zorder=5)
    ax2.text(-0.04, te, lab, fontsize=10, ha="right", va="center")
    ax2.text(s + 0.04, te + s - 0.01, "記録（時刻 %s、方向、振動数、偏光）" % u, fontsize=9.5, va="center", ha="left")

# 入射する背景の光線
for ta in (-0.9, -0.56, -0.14):
    s = (1 + ta) / 2
    ax2.plot([s, 0], [ta - s, ta], color="C1", lw=0.9, ls="--")

# 事象の間
ax2.annotate("", xy=(-0.12, 0.64), xytext=(-0.12, 0.26), arrowprops=dict(arrowstyle="-|>", lw=1.2, color="0.3"))
ax2.text(-0.17, 0.45, "放出と放出の間：\n読み出せる変化はない", fontsize=10, color="0.3", ha="right", va="center")

# 軸の向き
ax2.annotate("", xy=(-0.62, -0.78), xytext=(-0.62, -1.08), arrowprops=dict(arrowstyle="-|>", lw=1.0, color="0.4"))
ax2.annotate("", xy=(-0.32, -1.08), xytext=(-0.62, -1.08), arrowprops=dict(arrowstyle="-|>", lw=1.0, color="0.4"))
ax2.text(-0.60, -0.74, "時間", fontsize=9.5, color="0.4", ha="center", va="bottom")
ax2.text(-0.28, -1.08, "空間", fontsize=9.5, color="0.4", ha="left", va="center")

ax2.text(-0.78, -1.32,
         "光は 45 度の線（光円錐の表面）を進む。記録はその線に沿って変わらずに運ばれ、\n"
         "Doppler 偏移は線の両端の時計の違いとして現れる。",
         fontsize=9.5, color="0.3", ha="left", va="bottom")

fig.tight_layout()
out = HERE / "note_fig_lightcone_absorber_ja.png"
fig.savefig(out)
print("saved", out)
