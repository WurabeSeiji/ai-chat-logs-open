#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figures P1-P7 for the precision/seed isolation experiment.

Deep values (down to 1e-200) are plotted from the log10 column directly;
figs P1/P4 use a log y-axis of 10**log10 clipped at 1e-300 (float64-safe),
figs P2/P3/P5/P6 plot log10(Hperp/H) on a linear axis.
"""
from __future__ import annotations
import csv, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
FIG = HERE / "figures"
CONDS = [("A_IC64_DYN64", "tab:blue", "A: IC64+Dyn64"),
         ("B_IC64_DYN100", "tab:orange", "B: IC64+Dyn100"),
         ("C_IC100_DYN100", "tab:green", "C: IC100+Dyn100")]
LOG_ONSET = math.log10(0.05)


def load(N, cond):
    p = HERE / "data" / f"N{N}_D{N}" / cond / "timeseries.csv"
    if not p.exists():
        return None
    with open(p, newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        hdr = next(r)
        rows = list(r)
    i_t, i_l = hdr.index("tau"), hdr.index("log10_Hperp_frac")
    tau = np.array([float(x[i_t]) for x in rows])
    lg = np.array([(-320.0 if x[i_l] == "-inf" else float(x[i_l]))
                   for x in rows])
    return tau, lg


def onset_tau(tau, lg):
    i = np.flatnonzero(lg > LOG_ONSET)
    return float(tau[i[0]]) if i.size else None


def fig_abc_tau(N, fname):
    fig, ax = plt.subplots(figsize=(9, 6))
    for cond, c, lab in CONDS:
        d = load(N, cond)
        if d is None:
            continue
        tau, lg = d
        ax.semilogy(tau, 10.0 ** np.clip(lg, -300, 1), lw=1.0, color=c,
                    label=lab)
    ax.set_xlabel(r"$\tau=(2\pi/D)\,$step")
    ax.set_ylabel(r"$H_\perp/H$")
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title(f"N={N}, D={N}: A/B/C overlay")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def fig_abc_log(N, fname):
    fig, ax = plt.subplots(figsize=(9, 6))
    for cond, c, lab in CONDS:
        d = load(N, cond)
        if d is None:
            continue
        tau, lg = d
        ax.plot(tau, lg, lw=1.0, color=c, label=lab)
    ax.axhline(LOG_ONSET, color="gray", lw=0.7, ls=":")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$\log_{10}(H_\perp/H)$")
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title(f"N={N}, D={N}: early growth region (linear log10 axis)")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def fig_shifted(N, fname):
    fig, ax = plt.subplots(figsize=(9, 6))
    for cond, c, lab in (CONDS[1], CONDS[2]):
        d = load(N, cond)
        if d is None:
            continue
        tau, lg = d
        ot = onset_tau(tau, lg)
        if ot is None:
            ax.plot(tau, lg, lw=1.0, color=c,
                    label=lab + " (no onset; unshifted)")
        else:
            ax.plot(tau - ot, lg, lw=1.0, color=c,
                    label=lab + f" (shifted by onset tau={ot:.1f})")
    ax.axhline(LOG_ONSET, color="gray", lw=0.7, ls=":")
    ax.set_xlabel(r"$\tau-\tau_{\rm onset}$ (auxiliary shifted view)")
    ax.set_ylabel(r"$\log_{10}(H_\perp/H)$")
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title(f"N={N}: B vs C slope comparison after onset alignment "
                 f"(auxiliary figure)")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def fig_seed_floor(fname):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    marks = {7: "o", 8: "s"}
    for N in (7, 8):
        for cond, c, lab in CONDS:
            d = load(N, cond)
            if d is None:
                continue
            tau, lg = d
            if len(lg) < 2:
                continue
            x = -lg[1]  # -log10 of step-1 Hperp/H (effective seed floor)
            y = onset_tau(tau, lg)
            if y is None:
                continue
            ax.scatter([x], [y], marker=marks[N], color=c, s=60,
                       label=f"N={N} {lab.split(':')[0]}")
    ax.set_xlabel(r"$-\log_{10}$(step-1 $H_\perp/H$)  (effective seed floor)")
    ax.set_ylabel(r"onset $\tau$ ($H_\perp/H>0.05$)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title("Onset vs initial effective seed floor (few points; "
                 "no law claimed)")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    fig_abc_tau(8, "figP1_N8_ABC_tau.png")
    fig_abc_log(8, "figP2_N8_ABC_log_growth.png")
    fig_shifted(8, "figP3_N8_shifted_growth.png")
    fig_abc_tau(7, "figP4_N7_ABC_tau.png")
    fig_abc_log(7, "figP5_N7_ABC_log_growth.png")
    fig_shifted(7, "figP6_N7_shifted_growth.png")
    fig_seed_floor("figP7_onset_vs_seed_floor.png")
    print("figures done")


if __name__ == "__main__":
    main()
