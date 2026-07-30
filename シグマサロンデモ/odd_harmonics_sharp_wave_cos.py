#!/usr/bin/env python3
"""
同振幅の奇数倍音（cos 版）の重ね合わせによるシャープな波形の生成デモ。

f(x) = Σ_{n=1,3,5,...,17} cos(nx)

cos では x = 0 で全成分が +1 に揃うため、ピークは x = 0（および 2π）に立ち、
ピーク値はちょうど成分数 N = 9 になる。x = π では交互に ∓1 で -9 の負ピーク。
閉形式: Σ_{k=0}^{N-1} cos((2k+1)x) = sin(2Nx) / (2 sin x)

出力: odd_harmonics_sharp_wave_cos.svg / .png（本スクリプトと同じフォルダ）
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

x = np.linspace(-np.pi, np.pi, 4000)
odd_orders = np.arange(1, N_MAX + 1, 2)  # 1, 3, 5, ..., 17
harmonics = [np.cos(n * x) for n in odd_orders]
total = np.sum(harmonics, axis=0)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# --- 上段: 個々の奇数倍音 ---
cmap = plt.get_cmap("viridis")
for i, (n, h) in enumerate(zip(odd_orders, harmonics)):
    ax1.plot(x, h, lw=1.0, color=cmap(i / (len(odd_orders) - 1)),
             label=f"cos({n}x)")
ax1.set_title(f"同振幅の奇数倍音・cos 版（1〜{N_MAX}次、{len(odd_orders)}成分）")
ax1.set_ylabel("振幅")
ax1.legend(loc="upper right", fontsize=8, ncol=3)
ax1.grid(alpha=0.3)

# --- 下段: 重ね合わせ ---
ax2.plot(x, total, lw=1.8, color="crimson")
ax2.axhline(0, color="gray", lw=0.5)
ax2.set_title(f"重ね合わせ  f(x) = Σ cos(nx)  (n = 1, 3, ..., {N_MAX})")
ax2.set_xlabel("x [rad]")
ax2.set_ylabel("振幅")
ax2.grid(alpha=0.3)

peak_idx = np.argmax(total)
ax2.annotate(f"x=0 で全成分同位相（最大値 {total[peak_idx]:.2f} = 成分数）",
             xy=(x[peak_idx], total[peak_idx]),
             xytext=(x[peak_idx] + 0.5, total[peak_idx] - 0.8),
             arrowprops=dict(arrowstyle="->", color="black"))

ax2.set_xticks(np.arange(-np.pi, 1.1 * np.pi, np.pi / 2))
ax2.set_xticklabels(["-π", "-π/2", "0", "π/2", "π"])

fig.suptitle("古典波動論：同振幅奇数倍音（cos）の重ね合わせによるシャープな波形", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])

for ext in ("svg", "png"):
    path = os.path.join(OUT_DIR, f"odd_harmonics_sharp_wave_cos.{ext}")
    fig.savefig(path, dpi=150)
    print(f"saved: {path}")
