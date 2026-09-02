#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figures for the D-sweep (instruction section 9 naming).

Main figures use tau, D/N, chi axes only. Step-axis figures go to
figures/audit_step_axis/. Heatmaps: per-column zero-order hold of the
actual samples on a common tau grid (no smoothing, no interpolation
between D columns); cells beyond a run's coverage are masked.
"""
from __future__ import annotations
import csv, math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
SELECTED_D = [4, 5, 6, 7, 8, 9, 10, 12, 16, 24, 32, 64, 124, 256]


def load_ts(N, D, stage):
    p = ROOT / "data" / f"N{N}" / f"D{D:04d}" / f"timeseries_stage{stage}.csv"
    if not p.exists():
        return None
    raw = list(csv.reader(open(p, newline="", encoding="utf-8")))
    hdr, data = raw[0], raw[1:]
    g = lambda k: np.array([r[hdr.index(k)] for r in data], dtype=float)
    return {"step": g("step"), "tau": g("tau"), "chi": g("chi"),
            "f": g("Hperp_frac")}


def best_ts(N, D):
    return load_ts(N, D, "B") or load_ts(N, D, "A")


def load_summary(stage):
    p = ROOT / "results" / f"stage{stage}_summary.csv"
    raw = list(csv.DictReader(open(p, newline="", encoding="utf-8")))
    return raw


def fig_tau_curves(N, fname):
    fig, ax = plt.subplots(figsize=(9, 6))
    cmap = plt.get_cmap("viridis")
    for i, D in enumerate(SELECTED_D):
        ts = best_ts(N, D)
        if ts is None:
            continue
        f = np.clip(ts["f"], 1e-36, None)
        ax.semilogy(ts["tau"], f, lw=1.0,
                    color=cmap(i / (len(SELECTED_D) - 1)),
                    label=f"D={D}" + (" (=N)" if D == N else ""))
    ax.set_xlabel(r"$\tau = (2\pi/D)\,\mathrm{step}$")
    ax.set_ylabel(r"$H_\perp/H$")
    ax.set_xlim(0, 500 * 2 * math.pi / N)
    ax.set_ylim(1e-34, 3)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower right")
    ax.set_title(f"N={N}: Hperp/H vs tau, selected D (stage B window)")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def _num(rows, key):
    return np.array([float(r[key]) if r[key] != "NA" else np.nan
                     for r in rows])


def fig_onset(fname):
    sb = load_summary("B")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for N, mk, c in ((7, "o", "tab:blue"), (8, "s", "tab:orange")):
        rows = sorted((r for r in sb if int(r["N"]) == N),
                      key=lambda r: int(r["D"]))
        x = np.array([int(r["D"]) / N for r in rows])
        y = _num(rows, "onset_tau_0p05")
        ax.plot(x, y, mk + "-", ms=4, lw=0.8, color=c, label=f"N={N}")
    ax.axvline(1.0, color="gray", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("D/N")
    ax.set_ylabel(r"onset $\tau$  ($H_\perp/H>0.05$)")
    ax.grid(alpha=0.3)
    ax.legend()
    ax.set_title("Onset tau vs D/N (stage B, same tau window)")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def fig_gamma(fname):
    sb = load_summary("B")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for N, c in ((7, "tab:blue"), (8, "tab:orange")):
        rows = sorted((r for r in sb if int(r["N"]) == N),
                      key=lambda r: int(r["D"]))
        x = np.array([int(r["D"]) / N for r in rows])
        y = _num(rows, "gamma_tau")
        r2 = _num(rows, "fit_R2")
        good = r2 >= 0.98
        ax.plot(x[good], y[good], "o-", ms=4, lw=0.8, color=c,
                label=f"N={N} (R2>=0.98)")
        if np.any(~good & np.isfinite(y)):
            ax.plot(x[~good], y[~good], "x", ms=6, color=c,
                    label=f"N={N} (R2<0.98)")
    ax.axvline(1.0, color="gray", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("D/N")
    ax.set_ylabel(r"$\gamma_\tau$ (fit window $10^{-12}\leq f\leq 10^{-4}$)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title("Early exponential growth rate vs D/N (stage B)")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def fig_saturation(fname):
    sb = load_summary("B")
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for N, c in ((7, "tab:blue"), (8, "tab:orange")):
        rows = sorted((r for r in sb if int(r["N"]) == N),
                      key=lambda r: int(r["D"]))
        x = np.array([int(r["D"]) / N for r in rows])
        m = _num(rows, "sat_mean_Hperp_frac")
        q05 = _num(rows, "sat_q05_Hperp_frac")
        q95 = _num(rows, "sat_q95_Hperp_frac")
        ax.plot(x, m, "o-", ms=4, lw=0.8, color=c, label=f"N={N} mean")
        ax.fill_between(x, q05, q95, color=c, alpha=0.2,
                        label=f"N={N} q05-q95")
    ax.axvline(1.0, color="gray", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("D/N")
    ax.set_ylabel(r"tail-window $H_\perp/H$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_title("Saturation statistics vs D/N (tail 20%, stage B)")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def fig_chi(fname):
    fig, ax = plt.subplots(figsize=(9, 6))
    cmap = plt.get_cmap("viridis")
    for i, D in enumerate(SELECTED_D):
        for N, ls in ((7, "-"), (8, "--")):
            ts = best_ts(N, D)
            if ts is None:
                continue
            f = np.clip(ts["f"], 1e-36, None)
            ax.semilogy(ts["chi"], f, ls, lw=0.9,
                        color=cmap(i / (len(SELECTED_D) - 1)),
                        label=f"N={N} D={D}" if D in (4, 7, 8, 124) else None)
    ax.set_xlabel(r"$\chi = 2 r_N^2 (N-2)\,\tau$")
    ax.set_ylabel(r"$H_\perp/H$")
    ax.set_ylim(1e-34, 3)
    ax.set_xlim(0, 2 * 500 * 2 * math.pi / 7 * 2 * (1 / 15) * 5 / 10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2, loc="lower right")
    ax.set_title("chi-axis collapse test (N=7 solid, N=8 dashed)")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def fig_heatmap(N, fname):
    Ds = sorted(int(p.name[1:]) for p in (ROOT / "data" / f"N{N}").iterdir()
                if p.is_dir() and p.name.startswith("D"))
    tau_max = 500 * 2 * math.pi / N
    nbins = 400
    edges = np.linspace(0, tau_max, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    grid = np.full((nbins, len(Ds)), np.nan)
    for j, D in enumerate(Ds):
        ts = best_ts(N, D)
        if ts is None:
            continue
        tau, f = ts["tau"], np.clip(ts["f"], 1e-36, None)
        # zero-order hold of actual samples within covered range
        idx = np.searchsorted(tau, centers, side="right") - 1
        ok = (idx >= 0) & (centers <= tau[-1])
        grid[ok, j] = np.log10(f[idx[ok]])
    x = np.array(Ds) / N
    fig, ax = plt.subplots(figsize=(10, 6))
    # non-uniform x grid: build edges between successive D/N values
    xe = np.concatenate([[x[0] - (x[1] - x[0]) / 2],
                         (x[:-1] + x[1:]) / 2,
                         [x[-1] + (x[-1] - x[-2]) / 2]])
    pcm = ax.pcolormesh(xe, edges, np.ma.masked_invalid(grid),
                        cmap="magma", vmin=-32, vmax=0, shading="flat")
    ax.axvline(1.0, color="cyan", lw=0.7, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("D/N")
    ax.set_ylabel(r"$\tau$")
    ax.set_title(f"N={N}: log10(Hperp/H) on (D/N, tau) grid "
                 f"(zero-order hold of samples, masked beyond coverage)")
    fig.colorbar(pcm, ax=ax, label=r"$\log_{10}(H_\perp/H)$")
    fig.tight_layout()
    fig.savefig(FIG / fname, dpi=170)
    plt.close(fig)


def fig_audit_step(N, fname):
    fig, ax = plt.subplots(figsize=(9, 6))
    cmap = plt.get_cmap("viridis")
    for i, D in enumerate(SELECTED_D):
        ts = best_ts(N, D)
        if ts is None:
            continue
        f = np.clip(ts["f"], 1e-36, None)
        ax.semilogy(ts["step"], f, lw=1.0,
                    color=cmap(i / (len(SELECTED_D) - 1)), label=f"D={D}")
    ax.set_xlabel("step (audit axis)")
    ax.set_ylabel(r"$H_\perp/H$")
    ax.set_ylim(1e-34, 3)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2, loc="lower right")
    ax.set_title(f"N={N}: audit figure, step axis")
    fig.tight_layout()
    fig.savefig(FIG / "audit_step_axis" / fname, dpi=170)
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    (FIG / "audit_step_axis").mkdir(exist_ok=True)
    fig_tau_curves(7, "fig01_N7_tau_curves.png")
    fig_tau_curves(8, "fig02_N8_tau_curves.png")
    fig_onset("fig03_onset_tau_vs_D_over_N.png")
    fig_gamma("fig04_growth_gamma_tau_vs_D_over_N.png")
    fig_saturation("fig05_saturation_vs_D_over_N.png")
    fig_chi("fig06_chi_collapse.png")
    fig_heatmap(7, "fig07_N7_heatmap.png")
    fig_heatmap(8, "fig08_N8_heatmap.png")
    fig_audit_step(7, "audit_N7_step_curves.png")
    fig_audit_step(8, "audit_N8_step_curves.png")
    print("figures done")


if __name__ == "__main__":
    main()
