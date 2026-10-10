"""note 用の説明図：基本波の整数倍音 N 本を位相をそろえて重ねると、山が鋭くなる。
|Psi_N|^2 = [sin(N phi/2) / (N sin(phi/2))]^2 を N = 1, 4, 16 で描く（論文 §2 の式）。"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

phi = np.linspace(-np.pi, np.pi, 4001)
fig, ax = plt.subplots(figsize=(8, 4.2), constrained_layout=True)
for n, c in ((1, "tab:gray"), (4, "tab:blue"), (16, "tab:orange")):
    s = np.sum([np.exp(1j * h * phi) for h in range(1, n + 1)], axis=0) / n
    ax.plot(phi / np.pi, np.abs(s) ** 2, color=c, lw=2, label=("N = 1 (single tone)" if n == 1 else f"N = {n} harmonics"))
ax.set_xlabel("phase / pi  (one period of the fundamental)")
ax.set_ylabel("|Psi_N|^2")
ax.set_title("Aligning N harmonics: the peak width shrinks as 1/N")
ax.legend()
fig.savefig(Path(__file__).with_name("fig00_harmonic_peak_sharpening.png"), dpi=150)
