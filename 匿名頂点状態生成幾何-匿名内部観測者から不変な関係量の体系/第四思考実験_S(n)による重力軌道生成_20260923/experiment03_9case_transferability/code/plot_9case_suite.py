#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot Experiment 03 from already-saved row data only.

No dynamics are recomputed here. This separation guarantees that every plotted
curve can be traced back to a CSV row file in raw/.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def load_metrics(path):
    rows = []
    with open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=None)
    args = ap.parse_args()
    base = Path(args.base) if args.base else Path(__file__).resolve().parent.parent
    raw, summary, figs = base / "raw", base / "summary", base / "figures"
    figs.mkdir(parents=True, exist_ok=True)

    with open(summary / "cases.json", encoding="utf-8") as f:
        cfg = json.load(f)
    cases = cfg["cases"]

    # 3x3 orbit comparison grid: rows p, columns e.
    fig, axs = plt.subplots(3, 3, figsize=(13.5, 13.0))
    for c in cases:
        i, j = c["p_index"], c["e_index"]
        ax = axs[i, j]
        cid = c["case_id"]
        pn = read_csv(raw / f"pn_plot_{cid}.csv")
        sn = read_csv(raw / f"sn_{cid}.csv")
        ax.plot(pn["x"], pn["y"], lw=1.15, label="Paper 3 PN")
        ax.plot(sn["x"], sn["y"], lw=1.0, ls="--", label="same S(n)")
        ax.plot([0.0], [0.0], marker="+", ms=6)
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_title(f"p={c['p0']:g}, e={c['e0']:.2f}", fontsize=10)
        ax.set_xlabel("x [M]", fontsize=8)
        ax.set_ylabel("y [M]", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.2)
        if i == 0 and j == 0:
            ax.legend(fontsize=7)
    fig.suptitle("Experiment 03: one S(n) rule across 9 orbital conditions (q=1)", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(figs / "experiment03_9case_orbit_grid.png", dpi=180)
    fig.savefig(figs / "experiment03_9case_orbit_grid.svg")
    plt.close(fig)

    rows = load_metrics(summary / "case_metrics.csv")
    pvals = [24.3, 60.0, 120.0]
    evals = [0.20, 0.45, 0.55]
    rms = np.full((3, 3), np.nan)
    dprec = np.full((3, 3), np.nan)
    shrink_err = np.full((3, 3), np.nan)
    for r in rows:
        i = pvals.index(float(r["p0"]))
        j = evals.index(float(r["e0"]))
        rms[i, j] = float(r["radial_rms_over_a"])
        dprec[i, j] = float(r["precession_error_rad"])
        shrink_err[i, j] = float(r["sn_peri_shrink_fraction"]) - float(r["pn_peri_shrink_fraction"])

    def heat(data, title, stem, fmt):
        fig, ax = plt.subplots(figsize=(7.4, 5.8))
        im = ax.imshow(data, aspect="auto")
        ax.set_xticks(range(3), [f"{e:.2f}" for e in evals])
        ax.set_yticks(range(3), [f"{p:g}" for p in pvals])
        ax.set_xlabel("initial eccentricity e")
        ax.set_ylabel("initial semilatus rectum p [M]")
        ax.set_title(title)
        for i in range(3):
            for j in range(3):
                ax.text(j, i, format(data[i, j], fmt), ha="center", va="center", fontsize=9)
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(figs / f"{stem}.png", dpi=180)
        fig.savefig(figs / f"{stem}.svg")
        plt.close(fig)

    heat(rms, "Radial RMS error / a: same S(n) rule vs Paper 3 PN",
         "experiment03_radial_rms_over_a", ".3e")
    heat(dprec, "Periapsis-precession error [rad/orbit]",
         "experiment03_precession_error", ".3e")
    heat(shrink_err, "Periapsis shrink-fraction error over common passages",
         "experiment03_peri_shrink_error", ".3e")

    print("wrote orbit grid and three metric figures")


if __name__ == "__main__":
    main()
