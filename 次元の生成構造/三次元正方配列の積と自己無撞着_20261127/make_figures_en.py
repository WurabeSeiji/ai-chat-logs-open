#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
English figure generator (plotting only; reads saved CSV/NPZ data):
Generates all 5 English figures into results/figures_en/:
  1. fig_n64_512melt_evolution_en.png
  2. fig_n64_512melt_snapshots_en.png
  3. fig_n64_512melt_spectrum_en.png
  4. fig_n64_viewangle_vs_condensation_en.png
  5. fig_n64_viewangle_spectrum_en.png
"""
import csv
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["axes.unicode_minus"] = False


def load_csv_typed(path, int_cols, str_cols=()):
    data = {}
    with open(path, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            for k, v in row.items():
                data.setdefault(k, []).append(v)
    out = {}
    for k, vals in data.items():
        if k in str_cols:
            out[k] = vals
        elif k in int_cols:
            out[k] = np.array([int(x) for x in vals], dtype=np.int64)
        else:
            out[k] = np.array([float(x) for x in vals], dtype=np.float64)
    return out


def plot_evolution_en(d, out_path):
    steps = d["step"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True)
    fig.suptitle(r"N = 64 ($O_8 \otimes O_8$), Real Single-Cell Defect (0,1,0)+1e-4: Crystal Melting & Recombination (Steps 0–500)", fontsize=12.5, fontweight="bold")

    ax = axes[0, 0]
    ax.plot(steps, d["nz_1e5"], label="Nonzero count (|U| > 1e-5)", color="#1f77b4", lw=1.6)
    ax.plot(steps, d["n_eff"], label=r"Effective count $N_{\mathrm{eff}} = 1/\Sigma|U|^4$", color="#d62728", lw=1.6, ls="--")
    ax.set_yscale("log")
    ax.set_ylabel("Component count (log)")
    ax.set_title("1. Crystal Melting and Cluster Size")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)

    ax = axes[0, 1]
    ax.plot(steps, d["norm_P"], label="Growth factor ||P||", color="#2ca02c", lw=1.6)
    ax.set_ylabel("Growth factor ||P||", color="#2ca02c")
    ax.tick_params(axis="y", labelcolor="#2ca02c")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(steps, d["dist_crys"], label="Distance from crystal dist(T0)", color="#9467bd", lw=1.4, alpha=0.85)
    ax2.set_ylabel("Distance from crystal dist(T0)", color="#9467bd")
    ax2.tick_params(axis="y", labelcolor="#9467bd")
    ax2.set_ylim(-0.05, 1.08)
    ax.set_title("2. Stepwise Growth Factor ||P|| and Departure from Crystal")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="center right", fontsize=8.5)

    ax = axes[1, 0]
    ax.plot(steps, d["dim_i"], label="Active axes dim_i (= dim_k)", color="#ff7f0e", lw=1.6)
    ax.plot(steps, d["dim_j"], label="Active axes dim_j", color="#17becf", lw=1.3, ls="--")
    ax.plot(steps, d["max_comp"], label="Largest connected component", color="#8c564b", lw=1.3, ls=":")
    ax.set_ylabel("Axes / Vertices")
    ax.set_xlabel("Step")
    ax.set_title("3. Active Axes and Largest Connected Component")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8.5)

    ax = axes[1, 1]
    ax.plot(steps[1:], np.maximum(d["d1"][1:], 1e-16), label="d1 (1-step projective distance)", color="#1f77b4", lw=1.1, alpha=0.85)
    ax.plot(steps[2:], np.maximum(d["d2"][2:], 1e-16), label="d2 (2-step projective distance)", color="#e377c2", lw=1.1, alpha=0.85)
    ax.set_yscale("log")
    ax.set_ylim(1e-7, 2.0)
    ax.set_ylabel("Projective distance (log)")
    ax.set_xlabel("Step")
    ax.set_title("4. Step-to-Step Evolution d1, d2")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right", fontsize=8.5)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_snapshots_en(npz, d, out_path):
    chosen = [0, 30, 40, 60, 80, 500]
    step_to_idx = {int(s): idx for idx, s in enumerate(d["step"])}

    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5))
    fig.suptitle(r"N = 64 ($O_8 \otimes O_8$) 64×64 Projected Intensity Maps  $\sqrt{\Sigma_j |U[i,j,k]|^2}$", fontsize=13, fontweight="bold")

    for ax, st in zip(axes.flat, chosen):
        M = npz[f"map_ik_{st}"]
        M_log = np.log10(np.maximum(M, 1e-6))
        im = ax.imshow(M_log, origin="upper", cmap="inferno", vmin=-5.5, vmax=0.0, aspect="equal")
        idx = step_to_idx[st]
        nz = d["nz_1e5"][idx]
        neff = d["n_eff"][idx]
        di, dj, dk = d["dim_i"][idx], d["dim_j"][idx], d["dim_k"][idx]
        ax.set_title(f"step {st}: nz={nz}, N_eff={neff:.1f}, axes=({di},{dj},{dk})", fontsize=9.5)
        ax.set_xlabel("k (0..63)", fontsize=8.5)
        ax.set_ylabel("i (0..63)", fontsize=8.5)
        ax.set_xticks([0, 16, 32, 48, 63])
        ax.set_yticks([0, 16, 32, 48, 63])

    cbar_ax = fig.add_axes([0.92, 0.12, 0.018, 0.76])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(r"$\log_{10}\sqrt{\Sigma_j |U[i,j,k]|^2}$", fontsize=9.5)
    fig.subplots_adjust(left=0.06, right=0.89, top=0.90, bottom=0.07, wspace=0.25, hspace=0.28)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_spectrum_en(npz, out_path):
    chosen = [0, 30, 40, 80, 500]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    fig.suptitle(r"N = 64 ($O_8 \otimes O_8$) Amplitude Spectrum and 64-Element Energy Concentration", fontsize=13, fontweight="bold")

    ax = axes[0]
    bins = np.linspace(-5.5, 0.0, 56)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    for st, col in zip(chosen, colors):
        vals = npz[f"abs_{st}"]
        if len(vals) > 0:
            ax.hist(np.log10(np.maximum(vals, 1e-6)), bins=bins, histtype="step", lw=1.8, label=f"step {st} (nz={len(vals)})", color=col)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\log_{10} |U_{ijk}|$")
    ax.set_ylabel("Component count (log)")
    ax.set_title("1. Evolution of Component Amplitude Spectrum")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5)

    ax = axes[1]
    x = np.arange(64)
    ax.plot(x, npz["node_e_0"], label="step 0 (uniform)", color="#7f7f7f", lw=1.2, ls=":")
    ax.bar(x - 0.18, npz["node_e_80"], width=0.36, label="step 80", color="#d62728", alpha=0.7)
    ax.bar(x + 0.18, npz["node_e_500"], width=0.36, label="step 500", color="#1f77b4", alpha=0.7)
    ax.set_xlabel("Element index (0..63)")
    ax.set_ylabel("Energy per element (3-axis avg)")
    ax.set_title("2. Energy Concentration across 64 Elements (step 80 vs 500)")
    ax.set_xticks(np.arange(0, 65, 8))
    ax.set_ylim(0, 0.052)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", ncol=3, fontsize=8.5)

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_viewangle_vs_condensation_en(d, out_path):
    st = d["step"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True)
    fig.suptitle(r"N = 64 ($O_8 \otimes O_8$) Particle Condensation vs. Triangle View-Angle ($\rho = L_3/L_1$) & Exponent Stretching (Steps 0–100)", fontsize=12.2, fontweight="bold")

    ax = axes[0, 0]
    ax.plot(st, d["e_512"], label="Energy fraction in 512 cells E_512", color="#d62728", lw=2.0)
    ax.set_ylabel("Energy fraction E_512", color="#d62728")
    ax.tick_params(axis="y", labelcolor="#d62728")
    ax.set_ylim(-0.05, 1.08)
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(st, d["n_eff"], label="Effective count N_eff", color="#1f77b4", lw=1.5, ls="--")
    ax2.set_yscale("log")
    ax2.set_ylabel("Effective count N_eff (log)", color="#1f77b4")
    ax2.tick_params(axis="y", labelcolor="#1f77b4")
    ax.set_title("1. Particle Condensation (E_512 and N_eff)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="center left", fontsize=8.5)

    ax = axes[0, 1]
    ax.plot(st, d["dln_contrast"], label="ln(mean |U_512| / mean |U_melt|)", color="#2ca02c", lw=2.0)
    ax.set_ylabel("Exponent contrast Δln (particle vs background)")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_title("2. Exponent Contrast Between Particle and Melted Background")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)

    ax = axes[1, 0]
    ax.plot(st, d["tri_part_rho_med"], label="Particle-side triangles (median)", color="#d62728", lw=1.8)
    ax.fill_between(st, np.maximum(d["tri_part_rho_p10"], 1e-30), d["tri_part_rho_p90"], color="#d62728", alpha=0.15, label="Particle 10%–90% band")
    ax.plot(st, d["tri_melt_rho_med"], label="Background-side triangles (median)", color="#1f77b4", lw=1.8, ls="--")
    ax.fill_between(st, np.maximum(d["tri_melt_rho_p10"], 1e-30), d["tri_melt_rho_p90"], color="#1f77b4", alpha=0.15, label="Background 10%–90% band")
    ax.set_yscale("log")
    ax.set_ylim(1e-25, 2.0)
    ax.set_xlabel("Step")
    ax.set_ylabel(r"View-angle ratio $\rho = L_3 / L_1$ (log)")
    ax.set_title(r"3. Triangle View-Angle $\rho = L_3/L_1$ (1 = equilateral)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", fontsize=8)

    ax = axes[1, 1]
    ax.plot(st, d["tri_part_sig_dln_med"], label="Particle (3 sides > 1e-12, median)", color="#d62728", lw=1.8)
    ax.plot(st, d["tri_melt_sig_dln_med"], label="Background (3 sides > 1e-12, median)", color="#1f77b4", lw=1.8, ls="--")
    ax.plot(st, d["tri_part_dln_med"], label="Particle (all nonzero, median)", color="#ff7f0e", lw=1.2, ls=":")
    ax.plot(st, d["tri_melt_dln_med"], label="Background (all nonzero, median)", color="#17becf", lw=1.2, ls=":")
    ax.set_yscale("symlog", linthresh=0.5)
    ax.set_xlabel("Step")
    ax.set_ylabel(r"Exponent stretching $\Delta\ln = \ln(L_1/L_3)$")
    ax.set_title(r"4. Exponent Stretching $\Delta\ln$ (equilateral = 0)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_viewangle_spectrum_en(npz, out_path):
    chosen = [0, 30, 40, 60, 70, 80]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), sharey=True)
    fig.suptitle(r"N = 64 ($O_8 \otimes O_8$) Evolution of Triangle View-Angle Spectrum $\log_{10}(L_3 / L_1)$ (0 = equilateral)", fontsize=12.5, fontweight="bold")

    bins = np.linspace(-25.0, 0.0, 60)
    colors = ["#7f7f7f", "#ff7f0e", "#2ca02c", "#1f77b4", "#9467bd", "#d62728"]

    for ax, prefix, title in [
        (axes[0], "rho_part", "1. Particle Representative Position"),
        (axes[1], "rho_melt", "2. Melted Background Representative Position"),
    ]:
        for st, col in zip(chosen, colors):
            vals = npz[f"{prefix}_{st}"]
            if len(vals) > 0:
                logv = np.log10(np.maximum(vals, 1e-25))
                ax.hist(logv, bins=bins, histtype="step", lw=1.7, label=f"step {st} (n={len(vals)})", color=col)
        ax.set_yscale("log")
        ax.set_xlabel(r"$\log_{10}(\rho) = \log_{10}(L_3 / L_1)$   (0: equilateral, left: elongated / distant)")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=8.5)
    axes[0].set_ylabel("Triangle count (log)")

    plt.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def run():
    outdir = os.path.join("results", "figures_en")
    os.makedirs(outdir, exist_ok=True)

    d_melt = load_csv_typed(
        os.path.join("results", "tri_loop_n64_512melt_timeseries.csv"),
        int_cols=("step", "nz_1e5", "nz_1e3", "nz_1e1", "dim_i", "dim_j", "dim_k", "n_comp", "max_comp"),
        str_cols=("comp_sizes",),
    )
    npz_melt = np.load(os.path.join("results", "tri_loop_n64_512melt_snapshots.npz"))

    d_va = load_csv_typed(
        os.path.join("results", "tri_loop_n64_viewangle_timeseries.csv"),
        int_cols=("step", "tri_part_cnt", "tri_melt_cnt"),
    )
    npz_va = np.load(os.path.join("results", "tri_loop_n64_viewangle_snapshots.npz"))

    plot_evolution_en(d_melt, os.path.join(outdir, "fig_n64_512melt_evolution_en.png"))
    plot_snapshots_en(npz_melt, d_melt, os.path.join(outdir, "fig_n64_512melt_snapshots_en.png"))
    plot_spectrum_en(npz_melt, os.path.join(outdir, "fig_n64_512melt_spectrum_en.png"))
    plot_viewangle_vs_condensation_en(d_va, os.path.join(outdir, "fig_n64_viewangle_vs_condensation_en.png"))
    plot_viewangle_spectrum_en(npz_va, os.path.join(outdir, "fig_n64_viewangle_spectrum_en.png"))
    print("Generated 5 English figures in", outdir)


if __name__ == "__main__":
    run()
