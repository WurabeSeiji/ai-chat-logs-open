#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""実験O6（大N開始走行）の論文用図を生成する。"""

import csv
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(BASE, "largeN_splitting_result_v1")


def load(tag):
    taus, fs = [], []
    with open(os.path.join(DIR, f"fcurve_{tag}.csv")) as fh:
        for row in csv.DictReader(fh):
            taus.append(float(row["tau"]))
            fs.append(float(row["f"]))
    with open(os.path.join(DIR, f"summary_{tag}.json")) as fh:
        s = json.load(fh)
    return np.array(taus), np.array(fs), s


def main():
    tags = sorted(
        (os.path.basename(p)[7:-4] for p in glob.glob(os.path.join(DIR, "fcurve_*.csv"))),
        key=lambda t: int(t.split("_")[0][1:]))
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    for tag, col in zip(tags, colors):
        taus, f, s = load(tag)
        ax.semilogy(taus, f, color=col, lw=1.6,
                    label=(f"N={s['n']} (M={s['m']:,}), rate={s['onset_rate_per_step']:.3f}/step,"
                           f" {s['steps_per_decade']:.0f} steps/decade"))
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"dormant fraction $f(\tau)$ (complement projection)")
    ax.set_title(r"Large-$N$ onset under sequential reconstruction ($\delta=10^{-15}$,"
                 r" $f(0)\approx\delta^2=10^{-30}$)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="lower right", fontsize=9)
    ax.set_ylim(1e-31, 2.0)
    fig.tight_layout()
    out = os.path.join(DIR, "dormant_growth_large_n_v1.png")
    fig.savefig(out, dpi=160)
    print("wrote", out)


if __name__ == "__main__":
    main()
