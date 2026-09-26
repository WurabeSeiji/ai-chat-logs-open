#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproduce the Paper-5 C1 orbit figure from saved Paper-4 data and the
six-persistent-state internalized-nu generator.

Inputs (kept in the same directory):
  - run_nu_internal_state_control_v1.py
  - reference_audit_two_paths_rk4_snapshot.py
  - experiment02_sn_2p5pn_rr_C1_N12.csv

Outputs:
  - figures/figure_orbit_c1_paper4_vs_6state.svg
  - figures/figure_orbit_c1_paper4_vs_6state.png
  - figures/figure_orbit_c1_paper4_vs_6state_summary.txt

No fitted parameters and no reference-data feedback are used during generation.
The saved Paper-4 C1 CSV is read only after the six-state trajectory has been
generated, for plotting and residual reporting.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CONTROL = HERE / "run_nu_internal_state_control_v1.py"
SAVED_C1 = HERE / "experiment02_sn_2p5pn_rr_C1_N12.csv"
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)


def load_control():
    spec = importlib.util.spec_from_file_location("nu_control", CONTROL)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def main():
    ctrl = load_control()
    ref = np.genfromtxt(SAVED_C1, delimiter=",", names=True)

    p0 = float(ref["p"][0])
    e0 = float(ref["e"][0])
    nu0 = 0.25
    nsteps = len(ref) - 1

    x0 = np.array([1 + 0j, p0, e0, 1 + 0j, 1 + 0j], dtype=np.complex128)
    state = ctrl.init_B6(*x0, nu0)

    x6 = np.empty(nsteps + 1, dtype=float)
    y6 = np.empty(nsteps + 1, dtype=float)
    r6 = np.empty(nsteps + 1, dtype=float)
    t6 = np.empty(nsteps + 1, dtype=float)

    # Generate the full six-state trajectory first, without reading reference rows.
    for k in range(nsteps + 1):
        z = ctrl.core6(state)[:5]
        r, x, y, t = ctrl.readout_core(z)
        x6[k], y6[k], r6[k], t6[k] = x, y, r, t
        if k < nsteps:
            state = ctrl.macro_fast6(state)

    dx = np.max(np.abs(x6 - ref["x"]))
    dy = np.max(np.abs(y6 - ref["y"]))
    dr = np.max(np.abs(r6 - ref["r"]))
    dt = np.max(np.abs(t6 - ref["t_osc"]))

    fig, ax = plt.subplots(figsize=(7.5, 7.5))
    ax.plot(ref["x"], ref["y"], linewidth=2.2, label="Paper 4 saved C1 orbit")
    stride = max(1, len(ref) // 240)
    ax.plot(x6[::stride], y6[::stride], linestyle="none", marker="o",
            markersize=2.2, label="six-state internalized-nu generator")
    ax.plot([x6[0]], [y6[0]], marker="s", linestyle="none", markersize=6,
            label="start")
    ax.plot([x6[-1]], [y6[-1]], marker="x", linestyle="none", markersize=7,
            label="end")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x / M")
    ax.set_ylabel("y / M")
    ax.set_title("C1 orbit: Paper 4 saved trajectory vs six-state reconstruction\n"
                 "12000 macrosteps = 3 radial periods")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()

    svg = OUT / "figure_orbit_c1_paper4_vs_6state.svg"
    png = OUT / "figure_orbit_c1_paper4_vs_6state.png"
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(png, dpi=220, bbox_inches="tight")
    plt.close(fig)

    summary = OUT / "figure_orbit_c1_paper4_vs_6state_summary.txt"
    summary.write_text(
        "Paper-5 C1 orbit figure reproducibility summary\n"
        f"rows={len(ref)}\n"
        f"macrosteps={nsteps}\n"
        f"p0={p0:.17g}\n"
        f"e0={e0:.17g}\n"
        f"nu0={nu0:.17g}\n"
        f"max_abs_dx={dx:.17e}\n"
        f"max_abs_dy={dy:.17e}\n"
        f"max_abs_dr={dr:.17e}\n"
        f"max_abs_dt={dt:.17e}\n",
        encoding="utf-8",
    )
    print(summary.read_text(encoding="utf-8"))
    print(svg)
    print(png)


if __name__ == "__main__":
    main()
