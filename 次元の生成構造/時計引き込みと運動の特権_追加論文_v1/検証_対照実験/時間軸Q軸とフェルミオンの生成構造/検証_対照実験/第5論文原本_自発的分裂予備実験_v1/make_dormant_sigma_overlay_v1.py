#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""休眠フラクション f(τ) の元図（τ=0 からの semilogy・指数拡大が見える）に、
σ から作った実効散乱係数 R_eff を右軸（0=ボゾン, 0.697=灰色α, 1=フェルミオン）で重ねる、
単一の図化プログラム。

データは run_sigma_spectrum_largeN_v1.py が生成した
  sigma_spectrum_result_v1/sigmaspec_N*.csv   （列: tau, f, sigma_1..4, n_active）
を読むだけ（このプログラムは計算をしない・加工しない）。

  R_eff = cos^2(pi * sigma_2/sigma_1)   ← 白黒猫ブリッジ（そのまま、狙い値なし）

時間軸は元図と同じく τ=0 から全区間。指数拡大→交差→準振動の全過程を残す。
出力: sigma_spectrum_result_v1/dormant_sigma_overlay_v1.png
"""

import csv
import glob
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(BASE, "sigma_spectrum_result_v1")
GRAY = math.cos(math.pi * 23.0 / 124.0) ** 2  # α137根 R_eff = 0.697...


def load(path):
    tau_f, f_all = [], []
    tau_s, reff = [], []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            t = float(row["tau"])
            tau_f.append(t)
            f_all.append(float(row["f"]))
            s1, s2 = row.get("sigma_1", ""), row.get("sigma_2", "")
            if s1 not in ("", None) and s2 not in ("", None):
                s1 = float(s1); s2 = float(s2)
                if s1 > 0.0:
                    tau_s.append(t)
                    reff.append(math.cos(math.pi * s2 / s1) ** 2)
    return (np.array(tau_f), np.array(f_all), np.array(tau_s), np.array(reff))


def main():
    paths = sorted(
        glob.glob(os.path.join(DIR, "sigmaspec_N*.csv")),
        key=lambda p: int(os.path.basename(p).split("_")[1][1:]),
    )
    if not paths:
        raise SystemExit(f"no sigmaspec_N*.csv in {DIR}")
    cmap = plt.get_cmap("turbo")
    colors = [cmap(x) for x in np.linspace(0.05, 0.95, len(paths))]

    fig, axL = plt.subplots(figsize=(10.5, 6.2))
    axR = axL.twinx()

    for path, col in zip(paths, colors):
        n = int(os.path.basename(path).split("_")[1][1:])
        tau_f, f_all, tau_s, reff = load(path)
        # 左軸: f(τ) 元図と同じ semilogy（τ=0 から・指数拡大が見える）
        axL.semilogy(tau_f, f_all, color=col, lw=1.5, label=f"N={n}")
        # 右軸: R_eff を各N自分の最大振幅で正規化し、中央(0)を原点に（振動の"形"だけ比較）
        rc = reff - reff.mean()
        amp = np.max(np.abs(rc))
        if amp > 0.0:
            rc = rc / amp
        axR.plot(tau_s, rc, color=col, lw=1.4, ls=":", alpha=0.9)

    axR.axhline(0.0, color="gray", lw=0.7, alpha=0.6)  # 中央=原点

    axL.set_xlabel(r"$\tau$ (step)")
    axL.set_ylabel(r"dormant fraction $f(\tau)$  (solid, left, log)")
    axR.set_ylabel(r"$-R_{\rm eff}$ (flipped) normalized per-$N$  $\sim$ expansion velocity (dotted)")
    axL.set_ylim(1e-31, 2.0)
    axR.set_ylim(1.1, -1.1)  # 上下反転：速度として読む（拡大中=高, 飽和=低）
    axL.set_title(r"Exponential onset $\to$ quasi-oscillation:"
                  r" $f(\tau)$ (left, solid) and effective scattering $R_{\rm eff}$ (right, dotted)")
    axL.grid(alpha=0.3, which="both")
    axL.legend(loc="lower right", fontsize=9, title="solid=f, dotted=R_eff")
    fig.tight_layout()
    out = os.path.join(DIR, "dormant_sigma_overlay_v1.png")
    fig.savefig(out, dpi=160)
    print("wrote", out)


if __name__ == "__main__":
    main()
