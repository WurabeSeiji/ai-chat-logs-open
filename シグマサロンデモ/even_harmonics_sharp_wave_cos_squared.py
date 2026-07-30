#!/usr/bin/env python3
"""
同振幅の偶数倍音（cos 版）の重ね合わせ波形と、その振幅の2乗（強度分布）のデモ。

f(x)  = Σ_{n=2,4,6,...,16} cos(nx)
      = (sin(17x) / sin x − 1) / 2

奇数倍音版との対比:
- 偶数倍音のみの和は周期 π（基本波の半周期）を持つ。
- x = ±π でも全成分が +1 に揃うため、ピークは全て正（奇数版は ±π で −N）。
- ピーク値は成分数 N = 8、2乗の強度ピークは N² = 64。

出力: even_harmonics_sharp_wave_cos_squared.svg / .png（本スクリプトと同じフォルダ）
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
N_MAX = 17  # この次数以下の偶数倍音を使う

x = np.linspace(-np.pi, np.pi, 4000)
even_orders = np.arange(2, N_MAX + 1, 2)  # 2, 4, 6, ..., 16
total = np.sum([np.cos(n * x) for n in even_orders], axis=0)
intensity = total ** 2

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# --- 上段: 重ね合わせ波形（振幅） ---
ax1.plot(x, total, lw=1.8, color="crimson")
ax1.axhline(0, color="gray", lw=0.5)
ax1.set_title(f"重ね合わせ  f(x) = Σ cos(nx)  (n = 2, 4, ..., {even_orders[-1]}、{len(even_orders)}成分)")
ax1.set_ylabel("振幅")
ax1.grid(alpha=0.3)

# --- 下段: 振幅の2乗（強度） ---
ax2.plot(x, intensity, lw=1.8, color="darkblue")
ax2.fill_between(x, intensity, alpha=0.2, color="darkblue")
ax2.set_title("振幅の2乗  f(x)²  ── 強度分布")
ax2.set_xlabel("x [rad]")
ax2.set_ylabel("強度")
ax2.grid(alpha=0.3)

peak_idx = np.argmax(intensity)
ax2.annotate(f"最大値 {intensity[peak_idx]:.2f} = 成分数²（8² = 64）",
             xy=(x[peak_idx], intensity[peak_idx]),
             xytext=(x[peak_idx] + 0.5, intensity[peak_idx] * 0.85),
             arrowprops=dict(arrowstyle="->", color="black"))

ax2.set_xticks(np.arange(-np.pi, 1.1 * np.pi, np.pi / 2))
ax2.set_xticklabels(["-π", "-π/2", "0", "π/2", "π"])

fig.suptitle("古典波動論：同振幅偶数倍音（cos）の重ね合わせと振幅の2乗", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])

for ext in ("svg", "png"):
    path = os.path.join(OUT_DIR, f"even_harmonics_sharp_wave_cos_squared.{ext}")
    fig.savefig(path, dpi=150)
    print(f"saved: {path}")
