#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot the first 10 continuous orbits with logarithmic radial scale.

Input:
  ../02_RAW_DATA/raw_macro_trajectory.csv

Outputs:
  ../04_FIGURES_SVG/figure07_first_10_orbits_log_radial_scale.svg
  figure07_first_10_orbits_log_radial_scale.png (local/review output)

The raw trajectory is never recomputed here. This script is readout/plotting only.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RAW = ROOT / "02_RAW_DATA" / "raw_macro_trajectory.csv"
SVG_OUT = ROOT / "04_FIGURES_SVG" / "figure07_first_10_orbits_log_radial_scale.svg"
PNG_OUT = HERE / "figure07_first_10_orbits_log_radial_scale.png"

STEPS_PER_ORBIT = 4000
N_ORBITS = 10


def main():
    df = pd.read_csv(RAW)
    end_step = N_ORBITS * STEPS_PER_ORBIT
    seg = df.iloc[: end_step + 1].copy()

    # Use only stored readout coordinates and stored P; no trajectory regeneration.
    theta = np.unwrap(
        np.arctan2(
            seg["y_readout"].to_numpy(),
            seg["x_readout"].to_numpy(),
        )
    )
    r = seg["P"].to_numpy()

    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="polar")

    ax.plot(theta, r, linewidth=1.2)
    ax.set_yscale("log")

    # Show only the actual radial band traversed by the first 10 orbits.
    ax.set_ylim(r[-1] * 0.9995, r[0] * 1.0005)

    ax.set_title("First 10 orbits — logarithmic radial scale", pad=22)
    ax.set_ylabel("r / M  (log scale)", labelpad=35)
    ax.grid(True)

    ax.scatter([theta[0]], [r[0]], s=40, label="start: orbit 0")
    ax.scatter([theta[-1]], [r[-1]], s=40, marker="x", label="end: orbit 10")
    ax.legend(loc="lower left", bbox_to_anchor=(1.02, 0.0))

    SVG_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(SVG_OUT, bbox_inches="tight")
    fig.savefig(PNG_OUT, dpi=220, bbox_inches="tight")

    print(f"r(start) = {r[0]:.15g}")
    print(f"r(after 10 orbits) = {r[-1]:.15g}")
    print(f"SVG: {SVG_OUT}")
    print(f"PNG: {PNG_OUT}")


if __name__ == "__main__":
    main()
