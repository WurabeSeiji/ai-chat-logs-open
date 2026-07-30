#!/usr/bin/env python3
"""
同振幅の奇数倍音の重ね合わせによるシャープな波形の生成デモ。

f(x) = Σ_{n=1,3,5,...,17} sin(nx)

各奇数倍音は x = π/2 で位相が交互に ±1 となるが、
同振幅で重ねると波形はパルス状に集中し、鋭いピークが立つ。
（Fourier 級数の Dirichlet 核と同型の集中現象）

出力: odd_harmonics_sharp_wave.svg / .png（本スクリプトと同じフォルダ）
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# macOS の日本語フォント
plt.rcParams["font.family"] = "Hiragino Sans"
plt.rcParams["axes.unicode_minus"] = False

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
N_MAX = 17  # 最高次の奇数倍音

x = np.linspace(0, 2 * np.pi, 4000)
odd_orders = np.arange(1, N_MAX + 1, 2)  # 1, 3, 5, ..., 17
harmonics = [np.sin(n * x) for n in odd_orders]
total = np.sum(harmonics, axis=0)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# --- 上段: 個々の奇数倍音 ---
cmap = plt.get_cmap("viridis")
for i, (n, h) in enumerate(zip(odd_orders, harmonics)):
    ax1.plot(x, h, lw=1.0, color=cmap(i / (len(odd_orders) - 1)),
             label=f"sin({n}x)")
ax1.set_title(f"同振幅の奇数倍音（1〜{N_MAX}次、{len(odd_orders)}成分）")
ax1.set_ylabel("振幅")
ax1.legend(loc="upper right", fontsize=8, ncol=3)
ax1.grid(alpha=0.3)

# --- 下段: 重ね合わせ ---
ax2.plot(x, total, lw=1.8, color="crimson")
ax2.axhline(0, color="gray", lw=0.5)
ax2.set_title(f"重ね合わせ  f(x) = Σ sin(nx)  (n = 1, 3, ..., {N_MAX})")
ax2.set_xlabel("x [rad]")
ax2.set_ylabel("振幅")
ax2.grid(alpha=0.3)

peak_idx = np.argmax(total)
ax2.annotate(f"鋭いピーク（最大値 {total[peak_idx]:.2f}）",
             xy=(x[peak_idx], total[peak_idx]),
             xytext=(x[peak_idx] + 0.7, total[peak_idx] - 0.5),
             arrowprops=dict(arrowstyle="->", color="black"))

ax2.set_xticks(np.arange(0, 2.1 * np.pi, np.pi / 2))
ax2.set_xticklabels(["0", "π/2", "π", "3π/2", "2π"])

fig.suptitle("古典波動論：同振幅奇数倍音の重ね合わせによるシャープな波形", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])

for ext in ("svg", "png"):
    path = os.path.join(OUT_DIR, f"odd_harmonics_sharp_wave.{ext}")
    fig.savefig(path, dpi=150)
    print(f"saved: {path}")
