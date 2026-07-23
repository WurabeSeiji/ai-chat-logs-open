"""平面分解読出し予備実験の論文用図生成 v1（第4論文用）"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
RESULT_DIR = BASE_DIR / "nbody_plane_decomposition_readout_result_v1"
MPL_DIR = RESULT_DIR / ".matplotlib"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_nbody_plane_decomposition_readout_preliminary_v1 as exp

BLUE = "#2a78d6"
GREEN = "#008300"
MAGENTA = "#e87ba4"
GRAY = "#52514e"


def representative(body_count: int, seed_offset: int = 0):
    """代表試行を再計算して系列を得る。"""
    rng = np.random.default_rng(exp.SEED + 1000 + seed_offset)
    pairs = exp.relation_pairs(body_count)
    adjacency = exp.relation_adjacency(pairs)
    state0 = exp.initial_state(len(pairs), rng)
    generator = exp.build_generator(state0, adjacency)
    update, gamma = exp.cayley(generator)
    states = np.empty((exp.STEP_COUNT + 1, len(pairs)), dtype=complex)
    states[0] = state0
    for s in range(exp.STEP_COUNT):
        states[s + 1] = update @ states[s]
    planes, kernel_basis = exp.plane_decomposition(generator, exp.RANK_TOL)
    return states, planes, kernel_basis, gamma


def figure_plane_isomorphism() -> None:
    states, planes, kernel_basis, gamma = representative(6)
    show = 240
    fig, ax = plt.subplots(figsize=(8, 5))
    shades = plt.cm.Blues(np.linspace(0.45, 0.95, len(planes)))
    for idx, pl in enumerate(planes):
        _, increments = exp.phase_series(states, pl["p"], pl["q"])
        label = "plane phase increments (6 planes)" if idx == 0 else None
        ax.plot(np.arange(show), np.abs(increments[:show]), color=shades[idx],
                linewidth=1.6, label=label)
    if kernel_basis.shape[1] > 0:
        kernel_proj = states @ kernel_basis
        kernel_change = np.max(np.abs(np.diff(kernel_proj, axis=0)), axis=1)
        ax.plot(np.arange(show), kernel_change[:show], color=MAGENTA,
                linewidth=2.0, label="kernel component change (static)")
    ax.set_yscale("symlog", linthresh=1e-16)
    ax.set_xlabel(r"step $\tau$")
    ax.set_ylabel(r"per-step phase increment $|\Delta\varphi_j|$ / kernel change")
    ax.set_title("Each rotation plane advances at its own constant frequency;\n"
                 "the kernel does not advance at all (N=6)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="center right")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "plane_ab_isomorphism_v1.png", dpi=180)
    plt.close(fig)


def figure_conservation_drifts() -> None:
    data = json.load(open(RESULT_DIR / "plane_decomposition_readout_result_v1.json",
                          encoding="utf-8"))
    rows = data["summaries"]
    n_values = [r["body_count"] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(n_values, [r["F_max_amp_drift"] for r in rows], marker="o",
                color=BLUE, label="plane amplitude drift")
    ax.semilogy(n_values, [r["H_max_planeH_drift"] for r in rows], marker="s",
                color=GREEN, label=r"plane energy $H_j$ drift")
    ker_rows = [r for r in rows if r["kernel_dim_values"][0] > 0]
    ax.semilogy([r["body_count"] for r in ker_rows],
                [r["H_max_kerH_drift"] for r in ker_rows],
                marker="^", color=MAGENTA, linestyle="none", markersize=9,
                label=r"kernel energy $H_{\ker}$ drift (N with kernel)")
    ax.semilogy(n_values, [r["H_max_identity_err"] for r in rows], marker="d",
                color=GRAY, label=r"identity $\sum H_j+H_{\ker}-H$")
    ax.axhline(1e-10, color=GRAY, linestyle="--", linewidth=1.2,
               label=r"tolerance $10^{-10}$")
    ax.set_xticks(n_values)
    ax.set_xlabel("body count N")
    ax.set_ylabel("maximum over 32 trials x 720 steps")
    ax.set_title("Per-plane and kernel energies are individually conserved")
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="upper right", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "plane_energy_decomposition_v1.png", dpi=180)
    plt.close(fig)


def figure_inverse_square_across_planes() -> None:
    states, planes, kernel_basis, gamma = representative(9)
    thetas = [2.0 * math.atan(gamma * pl["sigma"]) for pl in planes]
    theta_min = min(thetas)
    dtheta = [2.0 * math.pi / (t / theta_min) for t in thetas]
    alpha = [4.0 * math.sin(t / 2.0) ** 2 for t in thetas]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(dtheta, alpha, marker="o", linestyle="none", color=BLUE,
              markersize=8, label="9 rotation planes of one N=9 trial")
    ref_x = np.array(sorted(dtheta))
    ref_y = alpha[0] * (ref_x / dtheta[0]) ** -2
    ax.loglog(ref_x, ref_y, color=GRAY, linestyle="--", linewidth=1.2,
              label=r"reference slope $-2$")
    ax.set_xlabel(r"spectral cell width  $\Delta\theta_j=2\pi\,\theta_{\min}/\theta_j$  (ratios generally non-integer)")
    ax.set_ylabel(r"acceleration coefficient  $4\sin^2(\theta_j/2)$")
    ax.set_title("Continuous-spectrum re-expression:\n"
                 "plane coefficients follow slope $-2$ by definition of the cell width")
    ax.grid(True, alpha=0.25, which="both")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "inverse_square_across_planes_v1.png", dpi=180)
    plt.close(fig)


def figure_kernel_asymmetry() -> None:
    data = json.load(open(RESULT_DIR / "plane_decomposition_readout_result_v1.json",
                          encoding="utf-8"))
    rows = [r for r in data["summaries"] if r["kernel_dim_values"][0] >= 2]
    n_values = [r["body_count"] for r in rows]
    angles = [r["I_min_kernel_direction_angle_deg"] for r in rows]
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar([str(n) for n in n_values], angles, width=0.5, color=BLUE)
    for bar, angle in zip(bars, angles):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.2,
                f"{angle:.1f}" + "°", ha="center", va="bottom", color=GRAY)
    max_proj = max(r["I_max_kernel_projector_gauge_err"] for r in rows)
    ax.set_xlabel("body count N (multi-dimensional kernel)")
    ax.set_ylabel("min angle between gauge-equivalent kernel directions (deg)")
    ax.set_title("Kernel directions are not unique (angles shown),\n"
                 f"while the kernel projector is unique (gauge error <= {max_proj:.1e})")
    ax.set_ylim(0, max(angles) * 1.3)
    ax.grid(True, alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "kernel_readout_asymmetry_v1.png", dpi=180)
    plt.close(fig)


def main() -> None:
    figure_plane_isomorphism()
    figure_conservation_drifts()
    figure_inverse_square_across_planes()
    figure_kernel_asymmetry()
    for name in ("plane_ab_isomorphism_v1.png", "plane_energy_decomposition_v1.png",
                 "inverse_square_across_planes_v1.png", "kernel_readout_asymmetry_v1.png"):
        print(RESULT_DIR / name)


if __name__ == "__main__":
    main()
