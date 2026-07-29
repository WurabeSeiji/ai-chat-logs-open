#!/usr/bin/env python3
"""AB波からRを読む散乱本体を変えず、R=0.0,...,1.0を1500衝突まで比較する。"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
INITIAL_SEARCH_PATH = HERE.parent / "search_initial_conditions_and_plot_v1.py"
TARGETS = tuple(index / 10.0 for index in range(11))
MAX_COLLISION = 1500
DISPLAY_COLLISIONS = (
    0,
    1,
    2,
    3,
    5,
    10,
    20,
    42,
    100,
    200,
    300,
    500,
    1000,
    1500,
)
SEARCH_TOLERANCE = 1.0e-14
MAX_SEARCH_ITERATIONS = 200
INVARIANT_TOLERANCE = 1.0e-10
EMPTY_CHANNEL_NORM_SQUARED = 1.0e-20


def load_initial_search_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "initial_state_search_for_R_0_to_1_longrun_sweep_v1",
        INITIAL_SEARCH_PATH,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load initial-state search: {INITIAL_SEARCH_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


initial_search = load_initial_search_module()
toy = initial_search.toy
plt = initial_search.plt


@dataclass(frozen=True)
class SectorPair:
    a_bosonic: np.ndarray
    b_bosonic: np.ndarray
    a_fermionic: np.ndarray
    b_fermionic: np.ndarray
    original_bosonic_power: float
    original_fermionic_power: float


@dataclass(frozen=True)
class MixSearchResult:
    target_reflection_rate: float
    sector_mix_parameter: float
    achieved_reflection_rate: float
    achieved_theta: float
    absolute_error: float
    iterations: int
    method: str


def pair_norm(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.vdot(a, a).real + np.vdot(b, b).real)


def target_tag(target_r: float) -> str:
    return f"R{int(round(10.0 * target_r)):02d}"


def split_pair_into_readout_sectors(
    a: np.ndarray,
    b: np.ndarray,
    source_params: Any,
) -> SectorPair:
    """現行theta読出しと同じFFTマスクで、元のAB波を二保存セクターへ分ける。"""

    shape = (source_params.chi_grid_n, source_params.eta_grid_n)
    frequencies = np.rint(
        np.fft.fftfreq(source_params.chi_grid_n, d=1.0 / source_params.chi_grid_n)
    ).astype(int)
    abs_frequency = np.abs(frequencies)
    fermionic_mask = (abs_frequency >= 4) & ((abs_frequency % 2) == 0)

    def split(vector: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        transformed = np.fft.fft(
            vector.reshape(shape),
            axis=0,
            norm="ortho",
        )
        fermionic = np.fft.ifft(
            transformed * fermionic_mask[:, None],
            axis=0,
            norm="ortho",
        ).reshape(-1)
        bosonic = np.fft.ifft(
            transformed * (~fermionic_mask)[:, None],
            axis=0,
            norm="ortho",
        ).reshape(-1)
        return bosonic, fermionic

    a_bosonic, a_fermionic = split(a)
    b_bosonic, b_fermionic = split(b)
    bosonic_power = pair_norm(a_bosonic, b_bosonic)
    fermionic_power = pair_norm(a_fermionic, b_fermionic)
    if bosonic_power <= 0.0 or fermionic_power <= 0.0:
        raise ValueError("both readout sectors must have positive pair power")

    bosonic_scale = math.sqrt(bosonic_power)
    fermionic_scale = math.sqrt(fermionic_power)
    return SectorPair(
        a_bosonic=a_bosonic / bosonic_scale,
        b_bosonic=b_bosonic / bosonic_scale,
        a_fermionic=a_fermionic / fermionic_scale,
        b_fermionic=b_fermionic / fermionic_scale,
        original_bosonic_power=bosonic_power,
        original_fermionic_power=fermionic_power,
    )


def candidate_pair(
    sector_pair: SectorPair,
    mix_parameter: float,
) -> tuple[np.ndarray, np.ndarray]:
    if not 0.0 <= mix_parameter <= 1.0:
        raise ValueError("mix_parameter must be in [0, 1]")
    bosonic_weight = math.sqrt(max(1.0 - mix_parameter, 0.0))
    fermionic_weight = math.sqrt(max(mix_parameter, 0.0))
    a = (
        bosonic_weight * sector_pair.a_bosonic
        + fermionic_weight * sector_pair.a_fermionic
    )
    b = (
        bosonic_weight * sector_pair.b_bosonic
        + fermionic_weight * sector_pair.b_fermionic
    )
    return a, b


def search_initial_sector_mix(
    target_r: float,
    sector_pair: SectorPair,
    source_params: Any,
) -> tuple[MixSearchResult, np.ndarray, np.ndarray]:
    """候補AB配列を作り、無変更theta_from_abの返値だけで初期混合量を探索する。"""

    if not 0.0 <= target_r <= 1.0:
        raise ValueError("target_r must be in [0, 1]")
    if target_r == 0.0 or target_r == 1.0:
        mix_parameter = target_r
        a, b = candidate_pair(sector_pair, mix_parameter)
        readout = toy.theta_from_ab(a, b, source_params)
        return (
            MixSearchResult(
                target_reflection_rate=target_r,
                sector_mix_parameter=mix_parameter,
                achieved_reflection_rate=readout.reflection_rate,
                achieved_theta=readout.theta,
                absolute_error=abs(readout.reflection_rate - target_r),
                iterations=0,
                method="readout-sector endpoint initial state",
            ),
            a,
            b,
        )

    low = 0.0
    high = 1.0
    best_mix = 0.5
    best_a, best_b = candidate_pair(sector_pair, best_mix)
    best_readout = toy.theta_from_ab(best_a, best_b, source_params)
    best_error = abs(best_readout.reflection_rate - target_r)
    iterations = 0
    for iterations in range(1, MAX_SEARCH_ITERATIONS + 1):
        midpoint = 0.5 * (low + high)
        a, b = candidate_pair(sector_pair, midpoint)
        readout = toy.theta_from_ab(a, b, source_params)
        error = abs(readout.reflection_rate - target_r)
        if error < best_error:
            best_mix = midpoint
            best_a = a
            best_b = b
            best_readout = readout
            best_error = error
        if error <= SEARCH_TOLERANCE:
            break
        if readout.reflection_rate < target_r:
            low = midpoint
        else:
            high = midpoint

    return (
        MixSearchResult(
            target_reflection_rate=target_r,
            sector_mix_parameter=best_mix,
            achieved_reflection_rate=best_readout.reflection_rate,
            achieved_theta=best_readout.theta,
            absolute_error=best_error,
            iterations=iterations,
            method="black-box bisection over initial readout-sector mixture",
        ),
        best_a,
        best_b,
    )


def normalized_pair_return_residual(
    a: np.ndarray,
    b: np.ndarray,
    initial_a: np.ndarray,
    initial_b: np.ndarray,
    initial_pair_norm: float,
) -> float:
    difference_norm_squared = float(
        np.vdot(a - initial_a, a - initial_a).real
        + np.vdot(b - initial_b, b - initial_b).real
    )
    return math.sqrt(max(difference_norm_squared, 0.0) / initial_pair_norm)


def safe_state_metrics(vector: np.ndarray, metric_context: Any) -> dict[str, float]:
    if float(np.vdot(vector, vector).real) <= EMPTY_CHANNEL_NORM_SQUARED:
        return {"L": 0.0, "N_eff": 0.0, "N_eff_2": 0.0}
    return toy.state_metrics(vector, metric_context)


def normalized_chi_density(source_params: Any, vector: np.ndarray) -> np.ndarray:
    if float(np.vdot(vector, vector).real) <= EMPTY_CHANNEL_NORM_SQUARED:
        return np.zeros(source_params.chi_grid_n, dtype=float)
    density = toy.src.chi_density(source_params, vector)
    maximum = float(np.max(density))
    return density / maximum if maximum > 0.0 else np.zeros_like(density)


def run_forward(
    initial_a: np.ndarray,
    initial_b: np.ndarray,
    target_tag_value: str,
    source_params: Any,
    metric_context: Any,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """目標Rを受け取らず、完成した初期AB配列だけから1500衝突を進める。"""

    a = initial_a.copy()
    b = initial_b.copy()
    initial_pair_norm = toy.pair_hermitian_norm(a, b)
    initial_closure = toy.pair_zero_closure(a, b)
    initial_readout = toy.theta_from_ab(a, b, source_params)
    _, initial_spectrum = toy.combined_chi_power(a, b, source_params)

    rows: list[dict[str, Any]] = []
    snapshots: dict[int, dict[str, Any]] = {}
    for collision in range(MAX_COLLISION + 1):
        readout = toy.theta_from_ab(a, b, source_params)
        current_pair_norm = toy.pair_hermitian_norm(a, b)
        current_closure = toy.pair_zero_closure(a, b)
        metrics_a = safe_state_metrics(a, metric_context)
        metrics_b = safe_state_metrics(b, metric_context)
        a_norm_fraction = float(np.vdot(a, a).real / current_pair_norm)
        return_residual = normalized_pair_return_residual(
            a,
            b,
            initial_a,
            initial_b,
            initial_pair_norm,
        )
        row = {
            "target_tag": target_tag_value,
            "collision": collision,
            "R_derived_from_current_AB": readout.reflection_rate,
            "theta_derived_from_current_AB": readout.theta,
            "R_drift_from_initial": abs(
                readout.reflection_rate - initial_readout.reflection_rate
            ),
            "theta_drift_from_initial": abs(readout.theta - initial_readout.theta),
            "pair_norm_drift": abs(current_pair_norm - initial_pair_norm),
            "closure_drift_abs": abs(current_closure - initial_closure),
            "pair_return_residual": return_residual,
            "A_pair_norm_fraction": a_norm_fraction,
            "B_pair_norm_fraction": 1.0 - a_norm_fraction,
            "L_A": metrics_a["L"],
            "N_eff_A": metrics_a["N_eff"],
            "L_B": metrics_b["L"],
            "N_eff_B": metrics_b["N_eff"],
        }
        rows.append(row)
        if collision in DISPLAY_COLLISIONS:
            _, current_spectrum = toy.combined_chi_power(a, b, source_params)
            snapshots[collision] = {
                "rho_A": normalized_chi_density(source_params, a),
                "rho_B": normalized_chi_density(source_params, b),
                "combined_spectrum_max_drift": float(
                    np.max(np.abs(current_spectrum - initial_spectrum))
                ),
                "metrics": row,
            }
        if collision < MAX_COLLISION:
            a, b = toy.rotate_ab(a, b, readout.theta)

    positive_rows = rows[1:]
    nearest_return = min(
        positive_rows,
        key=lambda row: float(row["pair_return_residual"]),
    )
    max_snapshot_spectrum_drift = max(
        float(snapshot["combined_spectrum_max_drift"])
        for snapshot in snapshots.values()
    )
    max_r_drift = max(float(row["R_drift_from_initial"]) for row in rows)
    max_theta_drift = max(float(row["theta_drift_from_initial"]) for row in rows)
    max_norm_drift = max(float(row["pair_norm_drift"]) for row in rows)
    max_closure_drift = max(float(row["closure_drift_abs"]) for row in rows)
    invariant_pass = (
        max_r_drift <= INVARIANT_TOLERANCE
        and max_theta_drift <= INVARIANT_TOLERANCE
        and max_norm_drift <= INVARIANT_TOLERANCE
        and max_closure_drift <= INVARIANT_TOLERANCE
        and max_snapshot_spectrum_drift <= INVARIANT_TOLERANCE
    )
    summary = {
        "initial_AB_readout": asdict(initial_readout),
        "initial_pair_norm": initial_pair_norm,
        "initial_closure": {
            "real": initial_closure.real,
            "imag": initial_closure.imag,
            "abs": abs(initial_closure),
        },
        "nearest_return_collision_through_1500": int(
            nearest_return["collision"]
        ),
        "nearest_return_residual_through_1500": float(
            nearest_return["pair_return_residual"]
        ),
        "max_R_drift": max_r_drift,
        "max_theta_drift": max_theta_drift,
        "max_pair_norm_drift": max_norm_drift,
        "max_closure_drift_abs": max_closure_drift,
        "max_snapshot_combined_spectrum_drift": max_snapshot_spectrum_drift,
        "invariant_verdict": "PASS" if invariant_pass else "CHECK",
        "selected_collision_metrics": [
            snapshots[collision]["metrics"] for collision in DISPLAY_COLLISIONS
        ],
        "_snapshots_for_plot": snapshots,
    }
    return rows, summary


def make_waveform_figure(
    target_r: float,
    summary: dict[str, Any],
    source_params: Any,
) -> tuple[Path, Path]:
    snapshots = summary.pop("_snapshots_for_plot")
    chi, _ = toy.src.make_grids(source_params)
    x = chi / math.pi
    fig, axes = plt.subplots(
        7,
        2,
        figsize=(12, 21),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    for ax, collision in zip(axes.flatten(), DISPLAY_COLLISIONS):
        snapshot = snapshots[collision]
        metrics = snapshot["metrics"]
        ax.plot(
            x,
            snapshot["rho_A"],
            label=f"A L={metrics['L_A']:.3g}, N={metrics['N_eff_A']:.3g}",
        )
        ax.plot(
            x,
            snapshot["rho_B"],
            label=f"B L={metrics['L_B']:.3g}, N={metrics['N_eff_B']:.3g}",
        )
        if collision > 0 and float(metrics["pair_return_residual"]) <= 1.0e-10:
            ax.set_facecolor("#eef8ee")
            suffix = " (return)"
        else:
            suffix = ""
        ax.set_title(f"R={target_r:.1f}, collision={collision}{suffix}")
        ax.set_ylabel("rho_chi / channel max")
        ax.legend(fontsize=6)
    for ax in axes[-1]:
        ax.set_xlabel("chi / pi")
    fig.suptitle(
        f"AB-derived R={summary['initial_AB_readout']['reflection_rate']:.15g}; "
        f"theta={summary['initial_AB_readout']['theta']:.12g}",
        fontsize=14,
    )
    stem = f"{target_tag(target_r)}_longrun_waveforms_v1"
    png_path = HERE / f"{stem}.png"
    svg_path = HERE / f"{stem}.svg"
    fig.savefig(png_path, dpi=160)
    fig.savefig(svg_path, dpi=160)
    plt.close(fig)
    return png_path, svg_path


def make_early_comparison_figure(
    summaries: list[dict[str, Any]],
    source_params: Any,
) -> tuple[Path, Path]:
    chi, _ = toy.src.make_grids(source_params)
    x = chi / math.pi
    fig, axes = plt.subplots(
        len(TARGETS),
        3,
        figsize=(15, 29),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    for row_index, (target_r, summary) in enumerate(zip(TARGETS, summaries)):
        snapshots = summary["_early_snapshots_for_plot"]
        for column_index, collision in enumerate((0, 1, 2)):
            ax = axes[row_index, column_index]
            ax.plot(x, snapshots[collision]["rho_A"], label="A")
            ax.plot(x, snapshots[collision]["rho_B"], label="B")
            if row_index == 0:
                ax.set_title(f"collision={collision}")
            if column_index == 0:
                ax.set_ylabel(f"R={target_r:.1f}\nrho/max")
            if row_index == len(TARGETS) - 1:
                ax.set_xlabel("chi / pi")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "R=0.0 to 1.0: waveform comparison before and after the first two collisions",
        fontsize=15,
    )
    png_path = HERE / "R_0_to_1_early_waveform_comparison_v1.png"
    svg_path = HERE / "R_0_to_1_early_waveform_comparison_v1.svg"
    fig.savefig(png_path, dpi=160)
    fig.savefig(svg_path, dpi=160)
    plt.close(fig)
    return png_path, svg_path


def make_longrun_heatmap(
    rows_by_target: dict[str, list[dict[str, Any]]],
) -> tuple[Path, Path]:
    matrices: dict[str, np.ndarray] = {}
    for key in (
        "A_pair_norm_fraction",
        "B_pair_norm_fraction",
        "L_A",
        "L_B",
        "N_eff_A",
        "N_eff_B",
    ):
        matrices[key] = np.asarray(
            [
                [float(row[key]) for row in rows_by_target[target_tag(target)]]
                for target in TARGETS
            ]
        )

    fig, axes = plt.subplots(
        3,
        2,
        figsize=(15, 11),
        sharex=True,
        constrained_layout=True,
    )
    specifications = (
        ("A_pair_norm_fraction", "A pair-norm fraction", 0.0, 1.0),
        ("B_pair_norm_fraction", "B pair-norm fraction", 0.0, 1.0),
        (
            "L_A",
            "A localization L",
            min(matrices["L_A"].min(), matrices["L_B"].min()),
            max(matrices["L_A"].max(), matrices["L_B"].max()),
        ),
        (
            "L_B",
            "B localization L",
            min(matrices["L_A"].min(), matrices["L_B"].min()),
            max(matrices["L_A"].max(), matrices["L_B"].max()),
        ),
        (
            "N_eff_A",
            "A effective harmonic number",
            min(matrices["N_eff_A"].min(), matrices["N_eff_B"].min()),
            max(matrices["N_eff_A"].max(), matrices["N_eff_B"].max()),
        ),
        (
            "N_eff_B",
            "B effective harmonic number",
            min(matrices["N_eff_A"].min(), matrices["N_eff_B"].min()),
            max(matrices["N_eff_A"].max(), matrices["N_eff_B"].max()),
        ),
    )
    for ax, (key, title, minimum, maximum) in zip(axes.flatten(), specifications):
        image = ax.imshow(
            matrices[key],
            origin="lower",
            aspect="auto",
            extent=(0, MAX_COLLISION, -0.05, 1.05),
            interpolation="nearest",
            vmin=minimum,
            vmax=maximum,
            cmap="viridis",
        )
        ax.set_title(title)
        ax.set_ylabel("AB-derived R")
        ax.set_yticks(TARGETS)
        ax.grid(False)
        fig.colorbar(image, ax=ax, shrink=0.85)
    for ax in axes[-1]:
        ax.set_xlabel("collision")
    fig.suptitle(
        "R dependence of AB exchange dynamics through collision 1500",
        fontsize=15,
    )
    png_path = HERE / "R_0_to_1_longrun_heatmap_v1.png"
    svg_path = HERE / "R_0_to_1_longrun_heatmap_v1.svg"
    fig.savefig(png_path, dpi=160)
    fig.savefig(svg_path, dpi=160)
    plt.close(fig)
    return png_path, svg_path


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# R=0.0–1.0 初期状態系列・1500衝突実験 v1",
        "",
        "散乱本体は変更していない。目標Rは独立した初期状態探索だけに使用し、",
        "前進処理へ渡したのは探索済みの初期A/B複素配列だけである。",
        "",
        "固定A基本波・B63全体振幅だけの族は R<61/64 なので、R=1を含む全域には",
        "元のAB波を現行読出しのボゾン／フェルミオン保存セクターへ分解した初期状態族を使用した。",
        "",
        "| target R | sector mix | AB-derived R | theta | nearest return collision <=1500 | residual | max R drift | invariant |",
        "|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for result, summary in zip(
        payload["initial_state_search_results"],
        payload["target_summaries"],
    ):
        lines.append(
            "| {target:.1f} | {mix:.15g} | {derived:.15g} | {theta:.15g} | "
            "{collision} | {residual:.3e} | {drift:.3e} | {verdict} |".format(
                target=result["target_reflection_rate"],
                mix=result["sector_mix_parameter"],
                derived=summary["initial_AB_readout"]["reflection_rate"],
                theta=summary["initial_AB_readout"]["theta"],
                collision=summary["nearest_return_collision_through_1500"],
                residual=summary["nearest_return_residual_through_1500"],
                drift=summary["max_R_drift"],
                verdict=summary["invariant_verdict"],
            )
        )
    lines.extend(
        [
            "",
            "表示衝突回:",
            "",
            "$$",
            ",\\ ".join(str(value) for value in DISPLAY_COLLISIONS),
            "$$",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    params = toy.base.Params(
        high_n=63,
        recursive_collision_count=MAX_COLLISION,
    )
    source_params = toy.base.build_source_params(params)
    metric_context = toy.base.MetricContext(source_params)
    original_a, original_b, case = initial_search.make_unit_templates(source_params)
    sector_pair = split_pair_into_readout_sectors(
        original_a,
        original_b,
        source_params,
    )

    all_rows: list[dict[str, Any]] = []
    rows_by_target: dict[str, list[dict[str, Any]]] = {}
    search_results: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    figure_paths: list[str] = []

    for target_r in TARGETS:
        search_result, initial_a, initial_b = search_initial_sector_mix(
            target_r,
            sector_pair,
            source_params,
        )
        tag = target_tag(target_r)
        rows, summary = run_forward(
            initial_a,
            initial_b,
            tag,
            source_params,
            metric_context,
        )
        early_snapshots = {
            collision: summary["_snapshots_for_plot"][collision]
            for collision in (0, 1, 2)
        }
        summary["_early_snapshots_for_plot"] = early_snapshots
        png_path, svg_path = make_waveform_figure(
            target_r,
            summary,
            source_params,
        )
        summary["waveform_figure_png"] = png_path.name
        summary["waveform_figure_svg"] = svg_path.name
        figure_paths.extend((png_path.name, svg_path.name))
        rows_by_target[tag] = rows
        all_rows.extend(rows)
        search_results.append(asdict(search_result))
        summaries.append(summary)
        print(
            f"target_R={target_r:.1f}",
            f"mix={search_result.sector_mix_parameter:.17g}",
            f"AB_derived_R={summary['initial_AB_readout']['reflection_rate']:.17g}",
            f"nearest_return={summary['nearest_return_collision_through_1500']}",
            f"residual={summary['nearest_return_residual_through_1500']:.3e}",
            f"invariant={summary['invariant_verdict']}",
        )

    early_png, early_svg = make_early_comparison_figure(
        summaries,
        source_params,
    )
    heatmap_png, heatmap_svg = make_longrun_heatmap(rows_by_target)
    figure_paths.extend(
        (early_png.name, early_svg.name, heatmap_png.name, heatmap_svg.name)
    )

    for summary in summaries:
        summary.pop("_early_snapshots_for_plot")

    payload = {
        "experiment": "R_0_to_1_longrun_sweep_v1",
        "design_boundary": {
            "target_R_used_only_in": "standalone initial-state error minimization",
            "core_theta_readout_modified": False,
            "forward_scattering_external_R_or_theta": False,
            "forward_function_accepts_target_R": False,
            "fixed_A_plus_scaled_full_B63_reachability": "[0, 61/64)",
            "full_range_initial_state_extension": (
                "normalized bosonic/fermionic readout-sector components "
                "derived from the original A fundamental plus B63 pair"
            ),
        },
        "case": asdict(case),
        "conditions": {
            "targets": list(TARGETS),
            "max_collision": MAX_COLLISION,
            "display_collisions": list(DISPLAY_COLLISIONS),
            "search_tolerance": SEARCH_TOLERANCE,
            "invariant_tolerance": INVARIANT_TOLERANCE,
            "forward_scattering_rule": (
                "unchanged real orthogonal AB rotation with theta re-derived "
                "from current AB at every collision"
            ),
        },
        "source_sector_powers_before_normalization": {
            "bosonic": sector_pair.original_bosonic_power,
            "fermionic": sector_pair.original_fermionic_power,
        },
        "core_runner": {
            "path": str(
                initial_search.TOY_RUNNER_PATH.relative_to(initial_search.TOY_DIR)
            ),
            "sha256": toy.sha256(initial_search.TOY_RUNNER_PATH),
        },
        "initial_state_search_results": search_results,
        "target_summaries": summaries,
        "figures": figure_paths,
    }
    rows_path = HERE / "R_0_to_1_longrun_rows_v1.csv"
    result_path = HERE / "R_0_to_1_longrun_result_v1.json"
    report_path = HERE / "R_0_to_1_longrun_report_v1.md"
    write_rows(rows_path, all_rows)
    result_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    write_report(report_path, payload)


if __name__ == "__main__":
    main()
