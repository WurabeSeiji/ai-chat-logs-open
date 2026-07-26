#!/usr/bin/env python3
"""Stage A2c必須図1〜15をPNG/SVGで生成する。"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
PROCESSED = HERE / "processed"
FIGURES = HERE / "figures"
CFG = json.loads((HERE / "config_locked.json").read_text(encoding="utf-8"))


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def arr(data: list[dict[str, str]], field: str) -> np.ndarray:
    return np.asarray([float(r[field]) for r in data], dtype=float)


def save(fig, stem: str) -> None:
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FIGURES / f"{stem}.{ext}", dpi=180 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def crossing_line(ax) -> None:
    ax.axvline(1167, color="black", ls=":", lw=0.9, label="existing crossing=1167")


def main() -> None:
    summary_path = PROCESSED / "lineage_analysis_summary.json"
    if not summary_path.is_file():
        raise SystemExit("EXECUTION_FAILED: analyze_direction_lineage.py gate missing")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "ANALYSIS_COMPLETE":
        raise SystemExit("EXECUTION_FAILED: lineage analysis incomplete")
    if FIGURES.exists() and any(FIGURES.iterdir()):
        raise SystemExit("EXECUTION_FAILED: figures/が空ではないため上書きを拒否")
    FIGURES.mkdir(parents=True, exist_ok=True)

    metrics = rows(RAW / "trajectory_direction_metrics.csv")
    snapshots = rows(PROCESSED / "direction_basis_snapshots.csv")
    late = rows(PROCESSED / "early_vs_late_direction_overlap.csv")
    transverse = rows(PROCESSED / "early_vs_transverse_overlap.csv")
    late_transverse = rows(PROCESSED / "late_direction_vs_transverse_overlap.csv")
    continuity = rows(PROCESSED / "consecutive_subspace_continuity.csv")
    flevels = [r for r in rows(PROCESSED / "lineage_by_f_level.csv") if r["status"] == "found"]
    qbands = rows(PROCESSED / "lineage_by_q_resolution_band.csv")

    step = arr(metrics, "step")
    f = arr(metrics, "f")
    q3 = arr(metrics, "q3_over_q1")
    q4 = arr(metrics, "q4_over_q1")
    d3 = arr(metrics, "direction_3_occupation")
    d4 = arr(metrics, "direction_4_occupation")
    late_overlap = arr(late, "overlap")
    late_theta1 = arr(late, "theta_1_rad")
    late_theta2 = arr(late, "theta_2_rad")
    trans_by_seed = {
        seed: [r for r in transverse if int(r["seed"]) == seed]
        for seed in range(CFG["transverse"]["seeds"])
    }

    # 1. f and all comparison coordinates
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.semilogy(step, f, lw=1.0, color="#1f77b4", label="Stage A0 trajectory f")
    fixed = sorted(set(int(x) for x in CFG["fixed_steps"]))
    for i, s in enumerate(fixed):
        ax.axvline(s, color="0.72", lw=0.45, alpha=0.55,
                   label="fixed steps" if i == 0 else None)
    level_steps = sorted(set(int(r["step"]) for r in flevels))
    ax.scatter(level_steps, f[level_steps], s=17, color="#d62728", zorder=4, label="fixed f first passages")
    ax.axvspan(CFG["late_reference"]["start_step"], CFG["late_reference"]["end_step"],
               color="#2ca02c", alpha=0.10, label="D34_late averaging window")
    crossing_line(ax)
    ax.set(xlabel="absolute step", ylabel="f",
           title="Figure 1: f and fixed direction-lineage comparison coordinates")
    ax.grid(alpha=0.22, which="both")
    ax.legend(fontsize=8, ncol=2)
    save(fig, "figure01_f_and_lineage_sampling_coordinates")

    # 2. q ratios and numerical-resolution bands
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.semilogy(step, np.maximum(q3, 1e-18), lw=0.9, label="q3/q1")
    ax.semilogy(step, np.maximum(q4, 1e-18), lw=0.9, label="q4/q1")
    for threshold in (1e-8, 1e-6, 1e-4, 1e-2):
        ax.axhline(threshold, color="0.4", lw=0.6, ls="--")
        ax.text(5010, threshold, f"{threshold:.0e}", va="center", fontsize=7)
    crossing_line(ax)
    ax.set(xlabel="absolute step", ylabel="relative q (display floor 1e-18)",
           title="Figure 2: q3/q1, q4/q1 and fixed numerical-resolution bands")
    ax.set_xlim(0, 5050)
    ax.grid(alpha=0.22, which="both")
    ax.legend(fontsize=8)
    save(fig, "figure02_q_ratios_and_resolution_bands")

    # 3. direction 3/4 occupation
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.semilogy(step, np.maximum(d3, 1e-34), lw=0.9, label="direction 3 occupation")
    ax.semilogy(step, np.maximum(d4, 1e-34), lw=0.9, label="direction 4 occupation")
    crossing_line(ax)
    ax.set(xlabel="absolute step", ylabel="occupation (display floor 1e-34)",
           title="Figure 3: existing aligned direction 3/4 occupations")
    ax.grid(alpha=0.22, which="both")
    ax.legend(fontsize=8)
    save(fig, "figure03_direction3_direction4_occupation")

    # 4. D34(t) vs D34_late principal angles
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.plot(step, late_theta1, lw=0.8, label="theta1")
    ax.plot(step, late_theta2, lw=0.8, label="theta2=max")
    ax.axhline(0.1, color="0.35", ls="--", lw=0.8, label="classification coordinate 0.1 rad")
    crossing_line(ax)
    ax.set(xlabel="absolute step", ylabel="principal angle (rad)",
           title="Figure 4: D34(t) versus D34_late principal angles")
    ax.grid(alpha=0.22)
    ax.legend(fontsize=8)
    save(fig, "figure04_D34_vs_late_principal_angles")

    # 5. D34(t) vs D34_late overlap
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.plot(step, late_overlap, lw=0.9)
    ax.axhline(0.95, color="0.35", ls="--", lw=0.8, label="classification coordinate 0.95")
    crossing_line(ax)
    ax.set(xlabel="absolute step", ylabel="trace(P34 Plate)/2", ylim=(-0.02, 1.02),
           title="Figure 5: D34(t) versus D34_late projector overlap")
    ax.grid(alpha=0.22)
    ax.legend(fontsize=8)
    save(fig, "figure05_D34_vs_late_overlap")

    # 6. D34(t) vs Tperp principal angles
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    for seed, data in trans_by_seed.items():
        ax.plot(arr(data, "step"), arr(data, "theta_2_rad"), lw=0.75, label=f"seed {seed}: max angle")
    crossing_line(ax)
    ax.set(xlabel="absolute step", ylabel="maximum principal angle (rad)",
           title="Figure 6: D34(t) versus existing Tperp spans")
    ax.grid(alpha=0.22)
    ax.legend(fontsize=8)
    save(fig, "figure06_D34_vs_Tperp_principal_angles")

    # 7. D34(t) vs Tperp overlap
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    for seed, data in trans_by_seed.items():
        ax.plot(arr(data, "step"), arr(data, "overlap"), lw=0.75, label=f"seed {seed}")
    ax.axhline(0.1, color="0.35", ls="--", lw=0.8, label="classification coordinate 0.1")
    crossing_line(ax)
    ax.set(xlabel="absolute step", ylabel="projector overlap", ylim=(-0.02, 1.02),
           title="Figure 7: D34(t) versus existing Tperp projector overlap")
    ax.grid(alpha=0.22)
    ax.legend(fontsize=8)
    save(fig, "figure07_D34_vs_Tperp_overlap")

    # 8. D34_late vs each Tperp
    fig, ax = plt.subplots(figsize=(7.4, 5.1))
    seeds = [int(r["seed"]) for r in late_transverse]
    values = [float(r["overlap"]) for r in late_transverse]
    bars = ax.bar([str(x) for x in seeds], values, color="#9467bd")
    ax.bar_label(bars, fmt="%.4f", fontsize=8)
    ax.set(xlabel="existing transverse seed", ylabel="projector overlap", ylim=(0, 1),
           title="Figure 8: D34_late versus each Tperp")
    ax.grid(alpha=0.22, axis="y")
    save(fig, "figure08_late_D34_vs_each_Tperp_overlap")

    cstep = arr(continuity, "step_to")
    cangle = arr(continuity, "maximum_principal_angle_rad")
    cdist = arr(continuity, "projector_distance")

    # 9. consecutive maximum angle
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.plot(cstep, cangle, lw=0.75)
    crossing_line(ax)
    ax.set(xlabel="ending absolute step", ylabel="maximum principal angle (rad)",
           title="Figure 9: consecutive-step D34 ambient-subspace rotation")
    ax.grid(alpha=0.22)
    save(fig, "figure09_consecutive_D34_max_angle")

    # 10. consecutive projector distance
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.plot(cstep, cdist, lw=0.75)
    crossing_line(ax)
    ax.set(xlabel="ending absolute step", ylabel="projector distance",
           title="Figure 10: consecutive-step D34 projector distance")
    ax.grid(alpha=0.22)
    save(fig, "figure10_consecutive_D34_projector_distance")

    # 11. f-level lineage
    xlevel = arr(flevels, "level")
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    ax.plot(xlevel, arr(flevels, "D34_vs_late_overlap"), "o-", ms=3, label="D34 vs late")
    ax.plot(xlevel, arr(flevels, "D34_vs_Tperp_max_overlap"), "s-", ms=3,
            label="D34 vs best of 3 Tperp seeds")
    ax.set_xscale("log")
    ax.set(xlabel="fixed f first-passage level", ylabel="projector overlap", ylim=(-0.02, 1.02),
           title="Figure 11: direction lineage at fixed f coordinates")
    ax.grid(alpha=0.22, which="both")
    ax.legend(fontsize=8)
    save(fig, "figure11_lineage_by_f_level")

    # 12. q-resolution band lineage
    labels = [r["q_resolution_band"] for r in qbands]
    x = np.arange(len(labels))
    width = 0.36
    late_band = np.asarray([
        float(r["median_D34_vs_late_overlap"]) if r.get("median_D34_vs_late_overlap") else np.nan
        for r in qbands
    ])
    trans_band = np.asarray([
        float(r["median_max_D34_vs_Tperp_overlap"]) if r.get("median_max_D34_vs_Tperp_overlap") else np.nan
        for r in qbands
    ])
    fig, ax = plt.subplots(figsize=(10.0, 5.5))
    ax.bar(x - width / 2, late_band, width, label="median D34 vs late")
    ax.bar(x + width / 2, trans_band, width, label="median max D34 vs Tperp")
    ax.set_xticks(x, labels, rotation=18, ha="right")
    ax.set(ylabel="median projector overlap", ylim=(0, 1),
           title="Figure 12: lineage stratified by numerical q-resolution band")
    ax.grid(alpha=0.22, axis="y")
    ax.legend(fontsize=8)
    save(fig, "figure12_lineage_by_q_resolution_band")

    # 13. crossing neighborhood
    lo, hi = CFG["display_windows"]["crossing_zoom"]
    mask = (step >= lo) & (step <= hi)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.0, 7.1), sharex=True)
    ax1.semilogy(step[mask], f[mask], lw=0.9, label="f")
    crossing_line(ax1)
    ax1.set_ylabel("f")
    ax1.grid(alpha=0.22, which="both")
    ax2.plot(step[mask], late_overlap[mask], lw=0.9, label="D34 vs late")
    for seed, data in trans_by_seed.items():
        ts = arr(data, "step")
        mm = (ts >= lo) & (ts <= hi)
        ax2.plot(ts[mm], arr(data, "overlap")[mm], lw=0.65, label=f"D34 vs Tperp seed {seed}")
    crossing_line(ax2)
    ax2.set(xlabel="absolute step", ylabel="projector overlap", ylim=(-0.02, 1.02))
    ax2.grid(alpha=0.22)
    ax2.legend(fontsize=7, ncol=2)
    fig.suptitle("Figure 13: crossing=1167 neighborhood, fixed step 900–1400")
    save(fig, "figure13_crossing_neighborhood_zoom")

    # 14. late/reference region
    lo, hi = CFG["display_windows"]["metastable_zoom"]
    mask = (step >= lo) & (step <= hi)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.0, 7.1), sharex=True)
    ax1.plot(step[mask], late_theta2[mask], lw=0.8, label="D34 vs late max angle")
    ax1.set_ylabel("angle (rad)")
    ax1.grid(alpha=0.22)
    ax1.axvspan(CFG["late_reference"]["start_step"], CFG["late_reference"]["end_step"],
                color="#2ca02c", alpha=0.10)
    ax2.plot(step[mask], late_overlap[mask], lw=0.8, label="D34 vs late")
    for seed, data in trans_by_seed.items():
        ts = arr(data, "step")
        mm = (ts >= lo) & (ts <= hi)
        ax2.plot(ts[mm], arr(data, "overlap")[mm], lw=0.65, label=f"Tperp seed {seed}")
    ax2.set(xlabel="absolute step", ylabel="overlap", ylim=(-0.02, 1.02))
    ax2.grid(alpha=0.22)
    ax2.legend(fontsize=7, ncol=2)
    fig.suptitle("Figure 14: fixed step 1400–2500 direction-lineage view")
    save(fig, "figure14_metastable_region_zoom")

    # 15. column bookkeeping vs 2D subspace
    swaps = np.asarray([1.0 if r["basis_column_swap"] == "True" else 0.0 for r in continuity])
    flips = arr(continuity, "sign_flip_count")
    col1 = arr(continuity, "aligned_column_1_abs_inner")
    col2 = arr(continuity, "aligned_column_2_abs_inner")
    fig, axes = plt.subplots(3, 1, figsize=(10.0, 8.2), sharex=True)
    axes[0].plot(cstep, swaps, lw=0.55, label="column swap flag")
    axes[0].plot(cstep, flips / 2.0, lw=0.55, label="sign-flip count / 2")
    axes[0].set_ylabel("column bookkeeping")
    axes[0].legend(fontsize=7)
    axes[1].plot(cstep, col1, lw=0.6, label="aligned column 1 |inner|")
    axes[1].plot(cstep, col2, lw=0.6, label="aligned column 2 |inner|")
    axes[1].set_ylabel("aligned column continuity")
    axes[1].legend(fontsize=7)
    axes[2].plot(cstep, cangle, lw=0.65, label="2D subspace max angle")
    axes[2].set(xlabel="ending absolute step", ylabel="ambient rotation (rad)")
    axes[2].legend(fontsize=7)
    for ax in axes:
        crossing_line(ax)
        ax.grid(alpha=0.22)
    fig.suptitle("Figure 15: direction-column exchange versus 2D subspace continuity")
    save(fig, "figure15_column_exchange_vs_subspace_continuity")

    manifest = {
        "stage": "A2c",
        "status": "FIGURES_COMPLETE",
        "figure_count": 15,
        "formats": ["png", "svg"],
        "files": sorted(p.name for p in FIGURES.iterdir()),
    }
    (FIGURES / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("FIGURES_COMPLETE: 15 PNG + 15 SVG")


if __name__ == "__main__":
    main()
