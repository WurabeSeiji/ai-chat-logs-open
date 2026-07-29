#!/usr/bin/env python3
"""2つの厳密有限位数根を初期AB状態へ逆探索し、1500衝突まで追跡する。"""

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
    248,
    300,
    500,
    1000,
    1240,
    1500,
)
SEARCH_TOLERANCE = 1.0e-15
INVARIANT_TOLERANCE = 1.0e-10
RECURRENCE_TOLERANCE = 1.0e-8


def load_initial_search_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "initial_state_search_for_exact_roots_longrun_v1",
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
class RootSpec:
    root_id: str
    n: int
    m: int
    interpretation: str

    @property
    def target_r(self) -> float:
        return math.cos(math.pi * self.m / self.n) ** 2

    @property
    def theta_for_current_real_rotation(self) -> float:
        return math.asin(math.sqrt(self.target_r))

    @property
    def current_real_rotation_period(self) -> int:
        numerator = self.n - 2 * self.m
        return 4 * self.n // math.gcd(numerator, 4 * self.n)


ROOTS = (
    RootSpec(
        root_id="R124_23",
        n=124,
        m=23,
        interpretation="low-energy alpha inverse near 137",
    ),
    RootSpec(
        root_id="R620_117",
        n=620,
        m=117,
        interpretation="high-energy alpha inverse near 128.946",
    ),
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


def run_root(
    root: RootSpec,
    initial_b_amplitude: float,
    source_params: Any,
    metric_context: Any,
    a_template: np.ndarray,
    b_template: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """完成した初期AB配列だけを受け取り、無変更の前進処理を実行する。"""

    a = a_template.copy()
    b = initial_b_amplitude * b_template
    initial_a = a.copy()
    initial_b = b.copy()
    initial_readout = toy.theta_from_ab(a, b, source_params)
    initial_pair_norm = toy.pair_hermitian_norm(a, b)
    initial_closure = toy.pair_zero_closure(a, b)
    _, initial_spectrum = toy.combined_chi_power(a, b, source_params)

    rows: list[dict[str, Any]] = []
    snapshots: dict[int, dict[str, Any]] = {}
    for collision in range(MAX_COLLISION + 1):
        readout = toy.theta_from_ab(a, b, source_params)
        pair_norm = toy.pair_hermitian_norm(a, b)
        closure = toy.pair_zero_closure(a, b)
        metrics_a = toy.state_metrics(a, metric_context)
        metrics_b = toy.state_metrics(b, metric_context)
        return_residual = normalized_pair_return_residual(
            a,
            b,
            initial_a,
            initial_b,
            initial_pair_norm,
        )
        row = {
            "root_id": root.root_id,
            "n_old_exchange_operator": root.n,
            "m": root.m,
            "current_real_rotation_period": root.current_real_rotation_period,
            "collision": collision,
            "R_derived_from_current_AB": readout.reflection_rate,
            "theta_derived_from_current_AB": readout.theta,
            "R_drift_from_initial": abs(
                readout.reflection_rate - initial_readout.reflection_rate
            ),
            "theta_drift_from_initial": abs(readout.theta - initial_readout.theta),
            "pair_norm_drift": abs(pair_norm - initial_pair_norm),
            "closure_drift_abs": abs(closure - initial_closure),
            "pair_return_residual": return_residual,
            "L_A": metrics_a["L"],
            "N_eff_A": metrics_a["N_eff"],
            "L_B": metrics_b["L"],
            "N_eff_B": metrics_b["N_eff"],
        }
        rows.append(row)
        if collision in DISPLAY_COLLISIONS:
            _, current_spectrum = toy.combined_chi_power(a, b, source_params)
            snapshots[collision] = {
                "rho_A": initial_search.chi_density(source_params, a),
                "rho_B": initial_search.chi_density(source_params, b),
                "combined_spectrum_max_drift": float(
                    np.max(np.abs(current_spectrum - initial_spectrum))
                ),
                "metrics": row,
            }
        if collision < MAX_COLLISION:
            a, b = toy.rotate_ab(a, b, readout.theta)

    period = root.current_real_rotation_period
    period_rows = [
        row
        for row in rows
        if int(row["collision"]) > 0 and int(row["collision"]) % period == 0
    ]
    first_period_row = rows[period]
    max_r_drift = max(float(row["R_drift_from_initial"]) for row in rows)
    max_theta_drift = max(float(row["theta_drift_from_initial"]) for row in rows)
    max_norm_drift = max(float(row["pair_norm_drift"]) for row in rows)
    max_closure_drift = max(float(row["closure_drift_abs"]) for row in rows)
    max_snapshot_spectrum_drift = max(
        float(snapshot["combined_spectrum_max_drift"])
        for snapshot in snapshots.values()
    )
    invariant_pass = (
        max_r_drift <= INVARIANT_TOLERANCE
        and max_theta_drift <= INVARIANT_TOLERANCE
        and max_norm_drift <= INVARIANT_TOLERANCE
        and max_closure_drift <= INVARIANT_TOLERANCE
        and max_snapshot_spectrum_drift <= INVARIANT_TOLERANCE
    )
    recurrence_pass = (
        float(first_period_row["pair_return_residual"])
        <= RECURRENCE_TOLERANCE
    )
    summary = {
        "root": asdict(root),
        "target_R_exact_double": root.target_r,
        "target_theta_for_current_real_rotation": (
            root.theta_for_current_real_rotation
        ),
        "old_complex_exchange_operator_order": root.n,
        "current_real_rotation_period": period,
        "initial_B_amplitude": initial_b_amplitude,
        "initial_AB_readout": asdict(initial_readout),
        "initial_pair_norm": initial_pair_norm,
        "initial_closure": {
            "real": initial_closure.real,
            "imag": initial_closure.imag,
            "abs": abs(initial_closure),
        },
        "first_period_collision": period,
        "first_period_pair_return_residual": float(
            first_period_row["pair_return_residual"]
        ),
        "periodic_return_residuals_through_1500": [
            {
                "collision": int(row["collision"]),
                "pair_return_residual": float(row["pair_return_residual"]),
            }
            for row in period_rows
        ],
        "max_R_drift": max_r_drift,
        "max_theta_drift": max_theta_drift,
        "max_pair_norm_drift": max_norm_drift,
        "max_closure_drift_abs": max_closure_drift,
        "max_snapshot_combined_spectrum_drift": max_snapshot_spectrum_drift,
        "invariant_verdict": "PASS" if invariant_pass else "CHECK",
        "recurrence_verdict": "PASS" if recurrence_pass else "CHECK",
        "selected_collision_metrics": [
            snapshots[collision]["metrics"] for collision in DISPLAY_COLLISIONS
        ],
    }
    summary["_snapshots_for_plot"] = snapshots
    return rows, summary


def make_waveform_figure(
    root: RootSpec,
    summary: dict[str, Any],
    source_params: Any,
) -> tuple[Path, Path]:
    snapshots = summary.pop("_snapshots_for_plot")
    chi, _ = toy.src.make_grids(source_params)
    x = chi / math.pi
    fig, axes = plt.subplots(
        8,
        2,
        figsize=(12, 23),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    period = root.current_real_rotation_period
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
        if collision > 0 and collision % period == 0:
            ax.set_facecolor("#eef8ee")
            suffix = " (exact return)"
        else:
            suffix = ""
        ax.set_title(f"{root.root_id}, collision={collision}{suffix}")
        ax.set_ylabel("rho_chi / max")
        ax.legend(fontsize=6)
    for ax in axes[-1]:
        ax.set_xlabel("chi / pi")
    fig.suptitle(
        f"{root.root_id}: AB-derived R={summary['initial_AB_readout']['reflection_rate']:.15g}; "
        f"real-rotation period={period}",
        fontsize=14,
    )
    png_path = HERE / f"{root.root_id}_longrun_waveforms_v1.png"
    svg_path = HERE / f"{root.root_id}_longrun_waveforms_v1.svg"
    fig.savefig(png_path, dpi=160)
    fig.savefig(svg_path, dpi=160)
    plt.close(fig)
    return png_path, svg_path


def mark_periods(ax: Any, period: int) -> None:
    for collision in range(period, MAX_COLLISION + 1, period):
        ax.axvline(collision, color="0.75", linewidth=0.7, linestyle=":")


def make_diagnostics_figure(
    rows_by_root: dict[str, list[dict[str, Any]]],
    summaries: list[dict[str, Any]],
) -> tuple[Path, Path]:
    fig, axes = plt.subplots(
        3,
        2,
        figsize=(14, 10),
        sharex="col",
        constrained_layout=True,
    )
    for column, summary in enumerate(summaries):
        root_id = str(summary["root"]["root_id"])
        rows = rows_by_root[root_id]
        collisions = np.asarray([int(row["collision"]) for row in rows])
        residuals = np.asarray(
            [float(row["pair_return_residual"]) for row in rows]
        )
        l_a = np.asarray([float(row["L_A"]) for row in rows])
        l_b = np.asarray([float(row["L_B"]) for row in rows])
        n_a = np.asarray([float(row["N_eff_A"]) for row in rows])
        n_b = np.asarray([float(row["N_eff_B"]) for row in rows])
        period = int(summary["current_real_rotation_period"])

        nonzero_collisions = collisions[1:]
        nonzero_residuals = residuals[1:]
        axes[0, column].semilogy(
            nonzero_collisions,
            nonzero_residuals,
            color="tab:purple",
            label="return residual",
        )
        period_mask = nonzero_collisions % period == 0
        axes[0, column].scatter(
            nonzero_collisions[period_mask],
            nonzero_residuals[period_mask],
            color="tab:red",
            s=22,
            zorder=3,
            label="period return",
        )
        axes[0, column].set_title(f"{root_id}: pair return residual")
        axes[0, column].set_ylabel("normalized residual")
        axes[0, column].set_ylim(1e-15, 3.0)
        axes[0, column].legend()

        axes[1, column].plot(collisions, l_a, label="A")
        axes[1, column].plot(collisions, l_b, label="B")
        axes[1, column].set_title("localization")
        axes[1, column].set_ylabel("L")
        axes[1, column].legend()

        axes[2, column].plot(collisions, n_a, label="A")
        axes[2, column].plot(collisions, n_b, label="B")
        axes[2, column].set_title("effective harmonic number")
        axes[2, column].set_ylabel("N_eff")
        axes[2, column].set_xlabel("collision")
        axes[2, column].legend()

        for row in range(3):
            mark_periods(axes[row, column], period)
            axes[row, column].grid(alpha=0.25)

    fig.suptitle(
        "Exact finite-order-root initial states: dynamics through collision 1500",
        fontsize=15,
    )
    png_path = HERE / "exact_finite_order_roots_longrun_diagnostics_v1.png"
    svg_path = HERE / "exact_finite_order_roots_longrun_diagnostics_v1.svg"
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
        "# 厳密有限位数根の初期状態逆探索・1500衝突実験 v1",
        "",
        "散乱本体は変更していない。目標Rは独立した初期B振幅探索だけに使用し、",
        "前進処理には探索済みの初期A/B配列だけを渡した。",
        "",
        "| root | target R | initial B amplitude | AB-derived R | real-rotation period | residual at first period | invariant | recurrence |",
        "|---|---:|---:|---:|---:|---:|---|---|",
    ]
    for summary in payload["root_summaries"]:
        lines.append(
            "| {root} | {target:.15g} | {amplitude:.15g} | {derived:.15g} | "
            "{period} | {residual:.3e} | {invariant} | {recurrence} |".format(
                root=summary["root"]["root_id"],
                target=summary["target_R_exact_double"],
                amplitude=summary["initial_B_amplitude"],
                derived=summary["initial_AB_readout"]["reflection_rate"],
                period=summary["current_real_rotation_period"],
                residual=summary["first_period_pair_return_residual"],
                invariant=summary["invariant_verdict"],
                recurrence=summary["recurrence_verdict"],
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
            "248と1240は、現在の実直交回転における厳密閉鎖点として追加した。",
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
    a_template, b_template, case = initial_search.make_unit_templates(
        source_params
    )

    all_rows: list[dict[str, Any]] = []
    rows_by_root: dict[str, list[dict[str, Any]]] = {}
    root_summaries: list[dict[str, Any]] = []
    search_results: list[dict[str, Any]] = []
    figure_paths: list[str] = []

    for root in ROOTS:
        search_result = initial_search.search_initial_b_amplitude(
            root.target_r,
            a_template,
            b_template,
            source_params,
            tolerance=SEARCH_TOLERANCE,
        )
        rows, summary = run_root(
            root,
            search_result.initial_b_amplitude,
            source_params,
            metric_context,
            a_template,
            b_template,
        )
        png_path, svg_path = make_waveform_figure(
            root,
            summary,
            source_params,
        )
        summary["waveform_figure_png"] = png_path.name
        summary["waveform_figure_svg"] = svg_path.name
        figure_paths.extend((png_path.name, svg_path.name))
        all_rows.extend(rows)
        rows_by_root[root.root_id] = rows
        root_summaries.append(summary)
        search_results.append(asdict(search_result))

    diagnostics_png, diagnostics_svg = make_diagnostics_figure(
        rows_by_root,
        root_summaries,
    )
    figure_paths.extend((diagnostics_png.name, diagnostics_svg.name))

    payload = {
        "experiment": "exact_finite_order_roots_longrun_v1",
        "design_boundary": {
            "target_R_used_only_in": "standalone initial-state error minimization",
            "core_theta_readout_modified": False,
            "forward_scattering_external_R_or_theta": False,
            "forward_function_accepts_target_R": False,
        },
        "case": asdict(case),
        "conditions": {
            "A_template": "unit broad fundamental N=1",
            "B_template": "unit equal-amplitude odd harmonics 1,3,...,63",
            "searched_initial_parameter": "B amplitude only",
            "max_collision": MAX_COLLISION,
            "display_collisions": list(DISPLAY_COLLISIONS),
            "search_tolerance": SEARCH_TOLERANCE,
            "invariant_tolerance": INVARIANT_TOLERANCE,
            "recurrence_tolerance": RECURRENCE_TOLERANCE,
            "forward_scattering_rule": (
                "unchanged real orthogonal AB rotation with theta re-derived from current AB"
            ),
        },
        "core_runner": {
            "path": str(initial_search.TOY_RUNNER_PATH.relative_to(initial_search.TOY_DIR)),
            "sha256": toy.sha256(initial_search.TOY_RUNNER_PATH),
        },
        "initial_state_search_results": search_results,
        "root_summaries": root_summaries,
        "figures": figure_paths,
    }
    rows_path = HERE / "exact_finite_order_roots_longrun_rows_v1.csv"
    result_path = HERE / "exact_finite_order_roots_longrun_result_v1.json"
    report_path = HERE / "exact_finite_order_roots_longrun_report_v1.md"
    write_rows(rows_path, all_rows)
    result_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_report(report_path, payload)

    for search_result, summary in zip(search_results, root_summaries):
        print(
            f"root={summary['root']['root_id']}",
            f"target_R={summary['target_R_exact_double']:.17g}",
            f"initial_B_amplitude={search_result['initial_b_amplitude']:.17g}",
            f"AB_derived_R={summary['initial_AB_readout']['reflection_rate']:.17g}",
            f"period={summary['current_real_rotation_period']}",
            f"return_residual={summary['first_period_pair_return_residual']:.3e}",
            f"invariant={summary['invariant_verdict']}",
            f"recurrence={summary['recurrence_verdict']}",
        )


if __name__ == "__main__":
    main()
