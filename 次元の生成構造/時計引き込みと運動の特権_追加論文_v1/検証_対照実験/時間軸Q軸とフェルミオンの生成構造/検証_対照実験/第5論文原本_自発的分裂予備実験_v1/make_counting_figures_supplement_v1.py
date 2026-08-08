#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
第6論文用 補足図の生成:
  (1) 機構図: 振幅の熱化と振動数の帯域端凝集の同時進行（主張2の「原因」）
  (2) U2改良版: 凝集前（中間時刻）と凝集後（終時刻）の規約格子の対比
      ——中間時刻では規約依存の豊かさ、終時刻では凝集による平坦化
"""

import math
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
import run_n_scaling_lowrank_v1 as lr
from run_counting_ceiling_v1 import weighted_spectrum, n_read, OMEGA, RESULT_DIR


def expansion_with_pr(n, steps=4000, sub=25, seed=0):
    sys_lr = lr.LowRankSystem(n)
    rng = np.random.default_rng(70260722 + seed)
    Z = lr.zero_closure_generic(rng, sys_lr.m)
    taus, spectra, prs = [], [], []
    for t in range(steps + 1):
        if t % sub == 0:
            th, sh = weighted_spectrum(sys_lr, Z)
            spectra.append((th, sh))
            prs.append(lr.participation_ratio(Z))
            taus.append(t)
        sys_lr.set_theta(np.angle(Z))
        sig = sys_lr.sigma_spectrum()[0]
        Z = sys_lr.cayley_step(Z, sig)
    return np.array(taus), spectra, np.array(prs), sys_lr.m


def main():
    n = 16
    taus, spectra, prs, m = expansion_with_pr(n)

    # ---- 図A: 機構（熱化と凝集の同時進行） ----
    fig = plt.figure(figsize=(12.5, 4.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0])

    # (a) 重み付き振動数スペクトルのスナップショット
    ax = fig.add_subplot(gs[0])
    snap_idx = [0, 20, len(taus) - 1]  # τ=0, 500, 4000
    colors = ["tab:blue", "tab:orange", "tab:red"]
    floor = 1e-7
    for k, (si, col) in enumerate(zip(snap_idx, colors)):
        th, sh = spectra[si]
        ax.vlines(th / OMEGA, floor, np.maximum(sh, floor), color=col,
                  alpha=0.75, lw=2, label=rf"$\tau={taus[si]}$")
    ax.axvline(0.5, color="k", ls=":", alpha=0.5)
    ax.axvline(1.0, color="k", ls=":", alpha=0.5)
    ax.set_yscale("log")
    ax.set_ylim(floor, 2.0)
    ax.set_xlabel(r"$\theta/\Omega$ (normalized frequency)")
    ax.set_ylabel("share (log)")
    ax.set_title(rf"$N={n}$: weighted frequency spectrum "
                 "(clusters toward band edges)", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    # (b) 熱化（PR/M上昇）と凝集（周波数クラス数減少）の同時進行
    ax = fig.add_subplot(gs[1])
    class_counts = [n_read(th, sh, 64, 1e-4) for th, sh in spectra]
    ax.plot(taus, np.array(prs) / m, "tab:green", lw=1.5,
            label=r"amplitude spread PR/$M$ (thermalization)")
    ax.set_ylabel(r"PR/$M$", color="tab:green")
    ax.set_ylim(0, 1.05)
    ax.tick_params(axis="y", labelcolor="tab:green")
    ax2 = ax.twinx()
    ax2.plot(taus, class_counts, "tab:purple", lw=1.5,
             label=r"occupied frequency classes ($B=64$)")
    ax2.set_ylabel("frequency classes", color="tab:purple")
    ax2.tick_params(axis="y", labelcolor="tab:purple")
    ax.set_xlabel(r"$\tau$")
    ax.set_title("Amplitude thermalizes while frequency condenses", fontsize=10)
    ax.grid(True, alpha=0.3)
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "mechanism_band_clustering_v1.png"), dpi=160)
    plt.close(fig)

    # ---- 図B: U2改良版（中間 vs 終時刻の規約格子） ----
    Bs = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]
    eps2s = np.logspace(-4, math.log10(0.5), 16)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6), sharey=True)
    viol = 0
    for ax, si, tag in [(axes[0], 20, "mid-expansion"),
                        (axes[1], len(taus) - 1, "late (condensed)")]:
        th, sh = spectra[si]
        grid = np.zeros((len(Bs), len(eps2s)), dtype=int)
        for i, B in enumerate(Bs):
            for j, e2 in enumerate(eps2s):
                c = n_read(th, sh, B, e2)
                grid[i, j] = c
                if c > min(B, int(1.0 / e2), n):
                    viol += 1
        im = ax.imshow(grid, aspect="auto", origin="lower", cmap="viridis",
                       vmin=0, vmax=10)
        ax.set_xticks(range(0, len(eps2s), 3))
        ax.set_xticklabels([f"{e:.0e}" for e in eps2s[::3]], fontsize=7)
        ax.set_yticks(range(len(Bs)))
        ax.set_yticklabels(Bs, fontsize=7)
        ax.set_xlabel(r"$\varepsilon^2$")
        ax.set_title(rf"$\tau={taus[si]}$ ({tag}): max count = {grid.max()}",
                     fontsize=10)
    axes[0].set_ylabel(r"$B=\Omega/\omega_0$")
    fig.colorbar(im, ax=axes, shrink=0.9, label=r"$n_{\mathrm{read}}$")
    fig.suptitle(rf"$N={n}$: same run, two instants — count varies with convention "
                 f"(mid) and condenses (late); hard-bound violations = {viol}",
                 fontsize=10)
    fig.savefig(os.path.join(RESULT_DIR, "u2_floor_grid_v2.png"), dpi=160)
    plt.close(fig)

    print(f"図2点を生成。硬い上界の違反 = {viol}")
    print(f"中間時刻 τ={taus[20]} の最大読出し数 = "
          f"{max(n_read(*spectra[20], B, e) for B in Bs for e in eps2s)}")


if __name__ == "__main__":
    main()
