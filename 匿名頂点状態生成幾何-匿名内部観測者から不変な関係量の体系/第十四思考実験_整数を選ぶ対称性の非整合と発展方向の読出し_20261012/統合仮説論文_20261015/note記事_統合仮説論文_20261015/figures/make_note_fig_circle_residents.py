"""note 用の説明図：円周の住民の波（論文 §2 の設定と §2.7 の式）。
上段：円周（xy 平面）の上に虚数振幅を垂直に描く。n=0 は ±iA の二状態、n=1 は一周で符号が反転し二周で戻る、n=2 は一周で戻る。
下段：等振幅の倍音の和。奇数次 N=17 は二周に一度、偶数次 N=18 は一周に一度、両方 N=18 は一周目の終わりで打ち消し合い二周目に戻る。
f_odd = sum_{k} cos((2k+1) phi/2), f_even = sum_{k} cos(k phi), f_all = f_odd + f_even（iA を単位とする）。"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

phi = np.linspace(0, 4 * np.pi, 4001)
lap1 = phi <= 2 * np.pi
fig = plt.figure(figsize=(12.5, 7.4), constrained_layout=True)
# upper row: circle with imaginary amplitude drawn vertically
cases = ((0, "n = 0 : two states +iA and -iA"), (1, "n = 1 : flips sign after one lap, returns after two"), (2, "n = 2 : returns after one lap"))
for col, (n, title) in enumerate(cases):
    ax = fig.add_subplot(2, 3, col + 1, projection="3d")
    ax.plot(np.cos(phi[lap1]), np.sin(phi[lap1]), 0 * phi[lap1], color="0.6", lw=1)
    if n == 0:
        ax.plot(np.cos(phi[lap1]), np.sin(phi[lap1]), np.ones(lap1.sum()), color="tab:blue", lw=2, label="+iA")
        ax.plot(np.cos(phi[lap1]), np.sin(phi[lap1]), -np.ones(lap1.sum()), color="tab:orange", lw=2, ls="--", label="-iA")
        ax.legend(loc="upper left", fontsize=8)
    else:
        amp = np.cos(n * phi / 2)
        ax.plot(np.cos(phi[lap1]), np.sin(phi[lap1]), amp[lap1], color="tab:blue", lw=2, label="lap 1")
        if n % 2 == 1:
            ax.plot(np.cos(phi[~lap1]), np.sin(phi[~lap1]), amp[~lap1], color="tab:orange", lw=2, ls="--", label="lap 2")
        ax.legend(loc="upper left", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_zlim(-1.2, 1.2); ax.set_zlabel("amplitude / A", fontsize=8, labelpad=2)
    ax.set_xticks([]); ax.set_yticks([]); ax.view_init(elev=28, azim=-60)
# lower row: sums of equal-amplitude harmonics over two laps
f_odd = np.sum([np.cos((2 * k + 1) * phi / 2) for k in range(0, 9)], axis=0)      # n = 1,3,...,17
f_even = np.sum([np.cos(k * phi) for k in range(1, 10)], axis=0)                  # n = 2,4,...,18
sums = ((f_odd, "odd n = 1..17 : peak every two laps"), (f_even, "even n = 2..18 : peak every lap"), (f_odd + f_even, "odd + even n = 1..18 : zero at end of lap 1"))
for col, (f, title) in enumerate(sums):
    ax = fig.add_subplot(2, 3, 4 + col)
    ax.plot(phi / np.pi, f, color="tab:blue", lw=1.6)
    ax.axvline(2, color="0.5", lw=1, ls=":"); ax.text(2.03, ax.get_ylim()[1] * 0.9, "one lap", fontsize=8, color="0.4")
    ax.axhline(0, color="0.8", lw=0.8)
    ax.set_xlabel("phase / pi  (0 to 4 pi = two laps)"); ax.set_ylabel("imaginary amplitude / A")
    ax.set_title(title, fontsize=10)
fig.savefig(Path(__file__).with_name("fig00_circle_residents_waves.png"), dpi=150)
print("saved")
