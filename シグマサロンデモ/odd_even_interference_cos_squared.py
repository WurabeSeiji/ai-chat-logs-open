#!/usr/bin/env python3
"""
奇数倍音と偶数倍音の干渉デモ（cos²・強度分布の重ね描き）。

f_odd(x)  = Σ cos(nx)  (n = 1, 3, ..., 17;  9成分)   → 強度ピーク 81  （x = 0, ±π）
f_even(x) = Σ cos(nx)  (n = 2, 4, ..., 16;  8成分)   → 強度ピーク 64  （x = 0, ±π）
f_sum(x)  = f_odd + f_even = Σ_{n=1}^{17} cos(nx)     → 強度ピーク 289 = 17²（x = 0 のみ）

振幅レベルで足すと、x = 0 では 9 + 8 = 17 の完全コヒーレント加算、
x = ±π では −9 + 8 = −1 のほぼ完全な相殺が起きる。
偶数倍音が作った π 周期の「複製ピーク」を奇数倍音が消し、局在は x = 0 に一意化される。

出力: odd_even_interference_cos_squared.svg / .png（本スクリプトと同じフォルダ）
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
N_MAX = 17

x = np.linspace(-np.pi, np.pi, 8000)
odd_orders = np.arange(1, N_MAX + 1, 2)   # 1, 3, ..., 17（9成分）
even_orders = np.arange(2, N_MAX + 1, 2)  # 2, 4, ..., 16（8成分）

f_odd = np.sum([np.cos(n * x) for n in odd_orders], axis=0)
f_even = np.sum([np.cos(n * x) for n in even_orders], axis=0)
f_sum = f_odd + f_even  # 振幅で足してから2乗する（干渉）

I_odd, I_even, I_sum = f_odd**2, f_even**2, f_sum**2

fig, ax = plt.subplots(figsize=(11, 6.5))

ax.plot(x, I_odd, lw=1.5, color="green",
        label=f"奇数倍音² (n=1,3,...,17)  ピーク {I_odd.max():.0f}")
ax.plot(x, I_even, lw=1.5, color="red",
        label=f"偶数倍音² (n=2,4,...,16)  ピーク {I_even.max():.0f}")
ax.plot(x, I_sum, lw=2.0, color="blue",
        label=f"合成 (奇+偶)²  ピーク {I_sum.max():.0f} = 17²")
ax.fill_between(x, I_sum, alpha=0.12, color="blue")

ax.annotate("x = 0：完全コヒーレント加算\n(9 + 8)² = 289",
            xy=(0, I_sum.max()), xytext=(0.55, 255),
            arrowprops=dict(arrowstyle="->", color="black"))
ax.annotate("x = ±π：ほぼ完全な相殺\n(−9 + 8)² = 1",
            xy=(np.pi, 1), xytext=(1.55, 90),
            arrowprops=dict(arrowstyle="->", color="black"))

ax.set_xlabel("x [rad]")
ax.set_ylabel("強度  f(x)²")
ax.set_title("奇数倍音と偶数倍音の干渉：強度分布 cos² の比較（振幅で加算 → 2乗）")
ax.set_xticks(np.arange(-np.pi, 1.1 * np.pi, np.pi / 2))
ax.set_xticklabels(["-π", "-π/2", "0", "π/2", "π"])
ax.grid(alpha=0.3)
ax.legend(loc="upper left", fontsize=10)

fig.tight_layout()

for ext in ("svg", "png"):
    path = os.path.join(OUT_DIR, f"odd_even_interference_cos_squared.{ext}")
    fig.savefig(path, dpi=150)
    print(f"saved: {path}")

# 数値確認
i0 = np.argmin(np.abs(x))
ipi = len(x) - 1
print(f"x=0 : f_odd={f_odd[i0]:.3f}, f_even={f_even[i0]:.3f}, "
      f"f_sum={f_sum[i0]:.3f}, I_sum={I_sum[i0]:.3f}")
print(f"x=π : f_odd={f_odd[ipi]:.3f}, f_even={f_even[ipi]:.3f}, "
      f"f_sum={f_sum[ipi]:.3f}, I_sum={I_sum[ipi]:.3f}")
