#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""準安定系列 N=5,40,300 の図：指数増幅→飽和→準安定プラトー。

run_metastable_series_v1.py から呼ばれ、RESULT_DIR / N_VALUES / DELTA / SEED を
外から差し替えて使う。原本 make_largeN_figure_v1.py とは別物（O6は交差後1500で
プラトーが写らないため、本図は準安定域まで描き、交差時刻と後期中央値を明示する）。
"""

import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULT_DIR = None
N_VALUES = (5, 40, 300)
DELTA = 1e-15
SEED = 0


def tag_for(n):
    return f"N{n:05d}_delta{DELTA:.0e}_seed{SEED}"


def load(n):
    tag = tag_for(n)
    taus, fs = [], []
    with open(os.path.join(RESULT_DIR, f"fcurve_{tag}.csv")) as fh:
        for row in csv.DictReader(fh):
            taus.append(float(row["tau"]))
            fs.append(float(row["f"]))
    with open(os.path.join(RESULT_DIR, f"summary_{tag}.json")) as fh:
        s = json.load(fh)
    meta_path = os.path.join(RESULT_DIR, f"metastable_{tag}.json")
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
    return np.array(taus), np.array(fs), s, meta


def main():
    colors = {5: "#9467bd", 40: "#1f77b4", 300: "#d62728"}
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    for n in N_VALUES:
        taus, f, s, meta = load(n)
        col = colors.get(n, "#333333")
        cross = s.get("crossing_tau")
        plateau = meta.get("plateau", {})
        med = plateau.get("tail_median_f")
        label = f"N={n} (M={s['m']:,}), rate={s['onset_rate_per_step']:.3f}/step"
        if med is not None:
            label += f", plateau median f≈{med:.2f}"
        ax.semilogy(taus, f, color=col, lw=1.4, label=label)
        if cross is not None:
            ax.axvline(cross, color=col, ls=":", lw=0.9, alpha=0.7)
        if med is not None:
            ax.hlines(med, plateau["tail_start_tau"], taus[-1],
                      color=col, ls="--", lw=0.9, alpha=0.8)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"dormant fraction $f(\tau)$ (complement projection)")
    ax.set_title(r"Onset to metastable plateau ($\delta=10^{-15}$, "
                 r"$f(0)\approx10^{-30}$): dotted=threshold crossing, "
                 r"dashed=plateau median")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="lower right", fontsize=8.5)
    ax.set_ylim(1e-31, 2.0)
    fig.tight_layout()
    out = os.path.join(RESULT_DIR, "onset_to_metastable_N5_40_300_v1.png")
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
