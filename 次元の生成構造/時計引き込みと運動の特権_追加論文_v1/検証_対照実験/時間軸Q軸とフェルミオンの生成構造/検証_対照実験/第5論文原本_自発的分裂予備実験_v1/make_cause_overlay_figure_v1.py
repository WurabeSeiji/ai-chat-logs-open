#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""拡大と停止の原因を1枚に重ねる図：f(τ) と σ₂/σ₁ の N双対スケール偏差。

run_cause_instrumented_v1.py が出力した cause_*.csv（tau, f, sigma2_over_sigma1）
を読み、N=5,40,300 を縦3段に並べる。各段で
  左軸(log): 休眠フラクション f(τ)  … 拡大の結果
  右軸(lin): ε/scale = 2(N-1)·(1/2 - σ₂/σ₁) … スペクトル偏差（N双対スケール）
を重ね、閾値交差 τ を縦線で示す。増幅中は ε/scale が親値にロックし、
準安定転移で動くことを可視化する。
"""

import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(HERE, "cause_instrumented_result_v1")
N_VALUES = (5, 40, 300)
DELTA = 1e-15
SEED = 0


def load(n):
    t, f, s = [], [], []
    tag = f"N{n:05d}_delta{DELTA:.0e}_seed{SEED}"
    with open(os.path.join(RESULT_DIR, f"cause_{tag}.csv")) as fh:
        for r in csv.DictReader(fh):
            t.append(int(r["tau"]))
            f.append(float(r["f"]))
            s.append(float(r["sigma2_over_sigma1"]))
    return np.array(t), np.array(f), np.array(s)


def main():
    fig, axes = plt.subplots(len(N_VALUES), 1, figsize=(9.0, 10.0))
    for ax, n in zip(axes, N_VALUES):
        t, f, s = load(n)
        scale = 1.0 / (2 * (n - 1))
        eps_scaled = (0.5 - s) / scale
        cross = int(t[np.argmax(f > 0.05)])

        # 左軸: f(τ) 対数
        ax.semilogy(t, f, color="#1f77b4", lw=1.4, label=r"$f(\tau)$ (dormant, left)")
        ax.set_ylabel(r"$f(\tau)$", color="#1f77b4")
        ax.tick_params(axis="y", labelcolor="#1f77b4")
        ax.set_ylim(1e-31, 2.0)
        ax.axvline(cross, color="0.4", ls=":", lw=1.0)
        ax.text(cross, 3e-31, f"  crossing τ={cross}", fontsize=8, color="0.3")

        # 右軸: ε/scale = 2(N-1)(1/2 - σ2/σ1) 線形
        ax2 = ax.twinx()
        ax2.plot(t, eps_scaled, color="#d62728", lw=1.3,
                 label=r"$\varepsilon/$scale$=2(N-1)(1/2-\sigma_2/\sigma_1)$ (right)")
        ax2.axhline(1.0, color="#d62728", ls="--", lw=0.7, alpha=0.6)
        ax2.set_ylabel(r"$\varepsilon \cdot 2(N-1)$", color="#d62728")
        ax2.tick_params(axis="y", labelcolor="#d62728")
        ax2.set_ylim(0.0, max(2.3, np.nanmax(eps_scaled) * 1.1))

        ax.set_title(f"N={n} (M={n*(n-1)//2:,}): expansion result $f$ vs spectral-deviation "
                     f"order parameter", fontsize=10)
        ax.grid(alpha=0.25, which="both")
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc="center right", fontsize=8)
    axes[-1].set_xlabel(r"$\tau$")
    fig.suptitle(r"Cause overlay: dormant fraction $f$ and $\sigma_2/\sigma_1$ deviation "
                 r"on the $N$-dual scale ($\delta=10^{-15}$)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = os.path.join(RESULT_DIR, "cause_overlay_f_sigma_ratio_v1.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
