#!/usr/bin/env python3
"""Stage A2a必須図1〜14をPNG/SVGで生成する。イベント判定は行わない。"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
# ★変更1: 原本パッケージを明示（コピー先では HERE が別物になるため）
PKG = (HERE / ".." / ".." / "第8論文_二段階seed除去による準安定相の因果分離"
       / "paper8_stage_A2a_seedless_N5").resolve()
RAW = PKG / "raw"
PROCESSED = PKG / "processed"
# ★変更2: 出力先を本フォルダへ（公開図を上書きしないため）
FIGURES = HERE / "figures_control"; FIGURES.mkdir(parents=True, exist_ok=True)
CONFIG = json.loads((PKG / "config_locked.json").read_text(encoding="utf-8"))
RUN1, RUN2 = CONFIG["run_ids"]
A0_F = (
    PKG.parent / "paper7_N5_reproduction" / "reproduced" /
    "metastable_series_result_v1" / "fcurve_N00005_delta1e-15_seed0.csv"
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def arrays(path: Path, columns: list[str]) -> dict[str, np.ndarray]:
    data = rows(path)
    return {c: np.asarray([float(r[c]) for r in data], dtype=float) for c in columns}


def save(fig, stem: str) -> None:
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(FIGURES / f"{stem}.{ext}", dpi=180 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


def plot_level_brackets(table: list[dict], field: str, ylabel: str, stem: str) -> None:
    x = np.asarray([float(r["level"]) for r in table])
    before = np.asarray([float(r[f"seeded_{field}_before"]) for r in table])
    after = np.asarray([float(r[f"seeded_{field}_after"]) for r in table])
    seedless = np.asarray([float(r[f"seedless_{field}_exact"]) for r in table])
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    ax.plot(x, before, "o-", ms=3, lw=0.9, label="seeded actual record before/at")
    ax.plot(x, after, "o-", ms=3, lw=0.9, label="seeded actual record after/at")
    ax.fill_between(x, np.minimum(before, after), np.maximum(before, after), alpha=0.16,
                    label="seeded actual-record bracket (not interpolation)")
    ax.plot(x, seedless, "s-", ms=3, lw=1.1, label="seedless exact first-passage step")
    ax.set_xscale("log")
    if np.all(before >= 0) and np.all(after >= 0) and np.all(seedless >= 0):
        positive = np.concatenate([before[before > 0], after[after > 0], seedless[seedless > 0]])
        if positive.size and positive.max() / positive.min() > 1e4:
            ax.set_yscale("symlog", linthresh=max(positive.min(), 1e-32))
    ax.set_xlabel("f level")
    ax.set_ylabel(ylabel)
    ax.set_title(ylabel + " by fixed f level")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8)
    save(fig, stem)


def main() -> None:
    gate = PROCESSED / "seeded_comparison_summary.json"
    if not gate.is_file() or json.loads(gate.read_text(encoding="utf-8")).get("status") != "COMPARISON_TABLES_COMPLETE":
        raise SystemExit("EXECUTION_FAILED: compare_with_seeded_reference.pyの完了記録がない")
    FIGURES.mkdir(parents=True, exist_ok=True)
    if any(FIGURES.iterdir()):
        raise SystemExit("EXECUTION_FAILED: figures/が空でないため上書きを拒否")

    f1 = arrays(RAW / RUN1 / "f_timeseries.csv", ["step", "f", "log10_f"])
    f2 = arrays(RAW / RUN2 / "f_timeseries.csv", ["step", "f", "log10_f"])
    q1 = arrays(RAW / RUN1 / "q_timeseries.csv", ["step", "q1", "q2", "q3", "q4"])
    q2 = arrays(RAW / RUN2 / "q_timeseries.csv", ["step", "q1", "q2", "q3", "q4"])
    o1 = arrays(RAW / RUN1 / "occupation_timeseries.csv", [
        "step", "direction_1_occupation", "direction_2_occupation",
        "direction_3_occupation", "direction_4_occupation",
        "other_rotating_occupation", "kernel_occupation",
    ])
    o2 = arrays(RAW / RUN2 / "occupation_timeseries.csv", [
        "step", "direction_1_occupation", "direction_2_occupation",
        "direction_3_occupation", "direction_4_occupation",
        "other_rotating_occupation", "kernel_occupation",
    ])
    seeded = arrays(A0_F, ["tau", "f"])

    # 1
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    ax.semilogy(f1["step"], f1["f"], lw=1.0)
    ax.set(xlabel="absolute step", ylabel="f", title="Figure 1: seedless f, full step 0–5000")
    ax.grid(alpha=0.25, which="both")
    save(fig, "figure01_seedless_f_full")

    # 2
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    ax.plot(f1["step"], f1["log10_f"], lw=1.0)
    ax.set(xlabel="absolute step", ylabel="log10(f)", title="Figure 2: seedless log10(f), step 0–5000")
    ax.grid(alpha=0.25)
    save(fig, "figure02_seedless_log10_f")

    # 3
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    ax.semilogy(seeded["tau"], seeded["f"], lw=1.0, label="seeded Stage A1b reference")
    ax.semilogy(f1["step"], f1["f"], lw=1.0, label="seedless A2a")
    ax.set_xlim(0, 5000)
    ax.set(xlabel="absolute step", ylabel="f", title="Figure 3: seeded vs seedless, absolute step")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    save(fig, "figure03_seeded_vs_seedless_absolute_step")

    # 4
    aligned = arrays(PROCESSED / "time_aligned_f_comparison.csv", [
        "relative_step_from_f_ge_1e-12", "seeded_f", "seedless_f"
    ])
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    ax.semilogy(aligned["relative_step_from_f_ge_1e-12"], aligned["seeded_f"], label="seeded")
    ax.semilogy(aligned["relative_step_from_f_ge_1e-12"], aligned["seedless_f"], label="seedless")
    ax.axvline(0, color="k", lw=0.8, ls=":")
    ax.set(xlabel="step shifted so first f≥1e-12 is 0", ylabel="f",
           title="Figure 4: parallel-shift comparison (display rule only)")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    save(fig, "figure04_first_passage_aligned_f")

    # 5
    growth = rows(PROCESSED / "seeded_vs_seedless_growth_rate.csv")
    x = np.sqrt(np.asarray([float(r["lower_level"]) for r in growth]) *
                np.asarray([float(r["upper_level"]) for r in growth]))
    seeded_rate = np.asarray([float(r["seeded_mean_exponential_rate_per_step"]) for r in growth])
    seedless_rate = np.asarray([
        float(r["seedless_mean_exponential_rate_per_step"]) if r["seedless_mean_exponential_rate_per_step"] else np.nan
        for r in growth
    ])
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    ax.plot(x, seeded_rate, "o-", ms=3, label="seeded")
    ax.plot(x, seedless_rate, "s-", ms=3, label="seedless")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set(xlabel="geometric midpoint of adjacent f levels", ylabel="mean exponential rate / step",
           title="Figure 5: amplification rate between fixed f coordinates")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    save(fig, "figure05_decade_growth_rate_comparison")

    by_level = rows(PROCESSED / "seeded_vs_seedless_by_f_level.csv")
    plot_level_brackets(by_level, "q3_over_q1", "q3/q1", "figure06_q3_over_q1_by_f_level")
    plot_level_brackets(by_level, "q4_over_q1", "q4/q1", "figure07_q4_over_q1_by_f_level")
    plot_level_brackets(
        by_level, "direction_3_occupation", "direction 3 occupation",
        "figure08_direction3_occupation_by_f_level"
    )
    plot_level_brackets(
        by_level, "direction_4_occupation", "direction 4 occupation",
        "figure09_direction4_occupation_by_f_level"
    )
    plot_level_brackets(
        by_level, "kernel_occupation", "kernel occupation",
        "figure10_kernel_occupation_by_f_level"
    )

    # 11
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    for field, label in (
        ("direction_1_occupation", "direction 1"),
        ("direction_2_occupation", "direction 2"),
        ("direction_3_occupation", "direction 3"),
        ("direction_4_occupation", "direction 4"),
        ("other_rotating_occupation", "other rotating"),
        ("kernel_occupation", "kernel"),
    ):
        ax.plot(o1["step"], o1[field], lw=1.0, label=label)
    ax.set(xlabel="absolute step", ylabel="occupation",
           title="Figure 11: seedless direction/rotation/kernel occupations")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    save(fig, "figure11_seedless_occupation_timeseries")

    # 12
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    for field in ("q1", "q2", "q3", "q4"):
        ax.plot(q1["step"], q1[field], lw=1.0, label=field)
    ax.set(xlabel="absolute step", ylabel="q", title="Figure 12: seedless q1–q4")
    ax.grid(alpha=0.25)
    ax.legend()
    save(fig, "figure12_seedless_q1_q4_timeseries")

    # 13: 固定観察窓で、採用イベント線を置かない。
    mask = f1["step"] <= 3000
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.4, 7.0), sharex=True)
    ax1.semilogy(f1["step"][mask], f1["f"][mask], lw=1.0)
    ax1.set_ylabel("f (log)")
    ax1.grid(alpha=0.25, which="both")
    ax2.plot(f1["step"][mask], f1["f"][mask], lw=1.0)
    ax2.set(xlabel="absolute step", ylabel="f (linear)")
    ax2.grid(alpha=0.25)
    fig.suptitle("Figure 13: seedless amplification-to-oscillation observation, fixed step 0–3000")
    save(fig, "figure13_seedless_amplification_to_oscillation_zoom")

    # 14
    fig, axes = plt.subplots(3, 1, figsize=(9.4, 8.0), sharex=False)
    axes[0].plot(f1["step"], f2["f"] - f1["f"], lw=0.9)
    axes[0].set(ylabel="Δf", title="f: exec2 − exec1")
    for field in ("q1", "q2", "q3", "q4"):
        axes[1].plot(q1["step"], q2[field] - q1[field], lw=0.8, label=f"Δ{field}")
    axes[1].set(ylabel="Δq")
    axes[1].legend(ncol=4, fontsize=8)
    for field in ("direction_3_occupation", "direction_4_occupation", "kernel_occupation"):
        axes[2].plot(o1["step"], o2[field] - o1[field], lw=0.8, label=f"Δ{field}")
    axes[2].set(xlabel="absolute step", ylabel="Δoccupation")
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    fig.suptitle("Figure 14: deterministic exec 1/2 differences")
    save(fig, "figure14_exec1_exec2_difference")

    manifest = {
        "stage": "A2a",
        "status": "FIGURES_COMPLETE",
        "figure_count": 14,
        "formats": ["png", "svg"],
        "files": sorted(p.name for p in FIGURES.iterdir()),
    }
    (FIGURES / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print("FIGURES_COMPLETE: 14 PNG + 14 SVG")


if __name__ == "__main__":
    main()
