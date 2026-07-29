#!/usr/bin/env python3
"""初期AB状態だけを数値探索し、無変更のAB→theta処理で前進検証する。"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_DIR = HERE.parent
TOY_RUNNER_PATH = TOY_DIR / "run_ab_invariant_theta_toy_v1.py"

ALPHA_DENOMINATOR = 137.035999084
R_ALPHA = 1.0 - math.sqrt(4.0 * math.pi / ALPHA_DENOMINATOR)
DEFAULT_TARGETS = (0.50, R_ALPHA)
EVOLUTION_COLLISIONS = (0, 1, 2, 3, 5, 10, 20, 42)
SEARCH_TOLERANCE = 1.0e-12
INVARIANT_TOLERANCE = 1.0e-10
MAX_SEARCH_ITERATIONS = 200
MAX_INITIAL_B_AMPLITUDE = 1.0e6


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "unchanged_ab_theta_toy_for_initial_state_search_v1",
        TOY_RUNNER_PATH,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load toy runner: {TOY_RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


toy = load_toy_module()
plt = toy.base.plt


@dataclass(frozen=True)
class InitialStateSearchResult:
    target_reflection_rate: float
    initial_a_amplitude: float
    initial_b_amplitude: float
    achieved_reflection_rate: float
    achieved_theta: float
    absolute_error: float
    iterations: int
    method: str


def target_name(target_r: float) -> str:
    if math.isclose(target_r, 0.50, rel_tol=0.0, abs_tol=1.0e-14):
        return "R050"
    if math.isclose(target_r, R_ALPHA, rel_tol=0.0, abs_tol=1.0e-14):
        return "Ralpha"
    return "R" + f"{target_r:.12f}".replace(".", "p")


def chi_density(source_params: Any, vector: np.ndarray) -> np.ndarray:
    density = toy.src.chi_density(source_params, vector)
    maximum = float(np.max(density))
    return density / maximum if maximum > 0.0 else density


def make_unit_templates(source_params: Any) -> tuple[np.ndarray, np.ndarray, Any]:
    case = next(
        item
        for item in toy.build_cases(63)
        if item.name == "odd_fermion_candidate_B63"
    )
    base_case = toy.to_base_case(case)
    a_template = toy.base.make_case_state(
        source_params,
        base_case,
        "A",
        hair_enabled=True,
    )
    b_template = toy.base.make_case_state(
        source_params,
        base_case,
        "B",
        hair_enabled=True,
    )
    return a_template, b_template, case


def read_initial_r(
    a_template: np.ndarray,
    b_template: np.ndarray,
    source_params: Any,
    initial_b_amplitude: float,
) -> Any:
    """候補初期状態を作り、既存のtheta_from_abをそのまま呼ぶ。"""

    a = a_template.copy()
    b = initial_b_amplitude * b_template
    return toy.theta_from_ab(a, b, source_params)


def search_initial_b_amplitude(
    target_r: float,
    a_template: np.ndarray,
    b_template: np.ndarray,
    source_params: Any,
) -> InitialStateSearchResult:
    """目標値との誤差だけを見て、初期B振幅をブラックボックス探索する。

    theta_from_ab の式は複製・変更しない。目標Rはこの前処理の比較値であり、
    前進散乱へ渡さない。
    """

    if not 0.0 < target_r < 1.0:
        raise ValueError("target_r must be strictly between 0 and 1")

    low = 0.0
    high = 1.0
    low_readout = read_initial_r(
        a_template,
        b_template,
        source_params,
        low,
    )
    high_readout = read_initial_r(
        a_template,
        b_template,
        source_params,
        high,
    )
    while high_readout.reflection_rate < target_r:
        high *= 2.0
        if high > MAX_INITIAL_B_AMPLITUDE:
            raise ValueError(
                "target R was not bracketed by the allowed initial B-amplitude search"
            )
        high_readout = read_initial_r(
            a_template,
            b_template,
            source_params,
            high,
        )

    best_amplitude = high
    best_readout = high_readout
    best_error = abs(high_readout.reflection_rate - target_r)
    low_error = abs(low_readout.reflection_rate - target_r)
    if low_error < best_error:
        best_amplitude = low
        best_readout = low_readout
        best_error = low_error

    iterations = 0
    for iterations in range(1, MAX_SEARCH_ITERATIONS + 1):
        midpoint = 0.5 * (low + high)
        readout = read_initial_r(
            a_template,
            b_template,
            source_params,
            midpoint,
        )
        error = abs(readout.reflection_rate - target_r)
        if error < best_error:
            best_amplitude = midpoint
            best_readout = readout
            best_error = error
        if error <= SEARCH_TOLERANCE:
            break
        if readout.reflection_rate < target_r:
            low = midpoint
        else:
            high = midpoint

    return InitialStateSearchResult(
        target_reflection_rate=target_r,
        initial_a_amplitude=1.0,
        initial_b_amplitude=best_amplitude,
        achieved_reflection_rate=best_readout.reflection_rate,
        achieved_theta=best_readout.theta,
        absolute_error=best_error,
        iterations=iterations,
        method="black-box bisection over initial B amplitude",
    )


def run_forward_from_initial_state(
    initial_b_amplitude: float,
    source_params: Any,
    metric_context: Any,
    a_template: np.ndarray,
    b_template: np.ndarray,
    figure_tag: str,
) -> dict[str, Any]:
    """目標Rを引数に取らず、完成した初期AB状態だけで前進計算する。"""

    a = a_template.copy()
    b = initial_b_amplitude * b_template
    initial_readout = toy.theta_from_ab(a, b, source_params)
    initial_norm = toy.pair_hermitian_norm(a, b)
    initial_closure = toy.pair_zero_closure(a, b)
    _, initial_spectrum = toy.combined_chi_power(a, b, source_params)
    chi, _ = toy.src.make_grids(source_params)
    x = chi / math.pi

    snapshots: dict[int, dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []
    for collision in range(max(EVOLUTION_COLLISIONS) + 1):
        readout = toy.theta_from_ab(a, b, source_params)
        current_norm = toy.pair_hermitian_norm(a, b)
        current_closure = toy.pair_zero_closure(a, b)
        _, current_spectrum = toy.combined_chi_power(a, b, source_params)
        metrics_a = toy.state_metrics(a, metric_context)
        metrics_b = toy.state_metrics(b, metric_context)
        row = {
            "collision": collision,
            "R_derived_from_current_AB": readout.reflection_rate,
            "theta_derived_from_current_AB": readout.theta,
            "R_drift_from_initial": abs(
                readout.reflection_rate - initial_readout.reflection_rate
            ),
            "theta_drift_from_initial": abs(readout.theta - initial_readout.theta),
            "pair_norm_drift": abs(current_norm - initial_norm),
            "closure_drift_abs": abs(current_closure - initial_closure),
            "combined_spectrum_max_drift": float(
                np.max(np.abs(current_spectrum - initial_spectrum))
            ),
            "L_A": metrics_a["L"],
            "N_eff_A": metrics_a["N_eff"],
            "L_B": metrics_b["L"],
            "N_eff_B": metrics_b["N_eff"],
        }
        rows.append(row)
        if collision in EVOLUTION_COLLISIONS:
            snapshots[collision] = {
                "rho_A": chi_density(source_params, a),
                "rho_B": chi_density(source_params, b),
                "metrics": row,
            }
        if collision < max(EVOLUTION_COLLISIONS):
            a, b = toy.rotate_ab(a, b, readout.theta)

    fig, axes = plt.subplots(
        4,
        2,
        figsize=(12, 12),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    for ax, collision in zip(axes.flatten(), EVOLUTION_COLLISIONS):
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
        ax.set_title(
            "AB-derived R={rate:.8f}, collision={collision}".format(
                rate=initial_readout.reflection_rate,
                collision=collision,
            )
        )
        ax.set_ylabel("rho_chi / max")
        ax.legend(fontsize=7)
    for ax in axes[-1]:
        ax.set_xlabel("chi / pi")

    stem = f"ab_initial_state_{figure_tag}_waveform_evolution_v1"
    png_path = HERE / f"{stem}.png"
    svg_path = HERE / f"{stem}.svg"
    fig.savefig(png_path, dpi=160)
    fig.savefig(svg_path, dpi=160)
    plt.close(fig)

    max_r_drift = max(float(row["R_drift_from_initial"]) for row in rows)
    max_theta_drift = max(float(row["theta_drift_from_initial"]) for row in rows)
    max_norm_drift = max(float(row["pair_norm_drift"]) for row in rows)
    max_closure_drift = max(float(row["closure_drift_abs"]) for row in rows)
    max_spectrum_drift = max(
        float(row["combined_spectrum_max_drift"]) for row in rows
    )
    invariants_pass = (
        max_r_drift <= INVARIANT_TOLERANCE
        and max_theta_drift <= INVARIANT_TOLERANCE
        and max_norm_drift <= INVARIANT_TOLERANCE
        and max_closure_drift <= INVARIANT_TOLERANCE
        and max_spectrum_drift <= INVARIANT_TOLERANCE
    )
    return {
        "forward_input": {
            "initial_A_array": "unit broad fundamental template",
            "initial_B_array": "searched amplitude times unit odd-harmonic B63 template",
            "initial_B_amplitude": initial_b_amplitude,
            "external_R_or_theta": False,
        },
        "initial_AB_readout": asdict(initial_readout),
        "initial_pair_norm": initial_norm,
        "initial_closure": {
            "real": initial_closure.real,
            "imag": initial_closure.imag,
            "abs": abs(initial_closure),
        },
        "max_R_drift": max_r_drift,
        "max_theta_drift": max_theta_drift,
        "max_pair_norm_drift": max_norm_drift,
        "max_closure_drift_abs": max_closure_drift,
        "max_combined_spectrum_drift": max_spectrum_drift,
        "invariant_verdict": "PASS" if invariants_pass else "CHECK",
        "figure_png": png_path.name,
        "figure_svg": svg_path.name,
        "selected_collision_metrics": [
            rows[collision] for collision in EVOLUTION_COLLISIONS
        ],
    }


def write_report(path: Path, result: dict[str, Any]) -> None:
    lines = [
        "# 初期状態の数値探索と無変更AB読出しによる前進検証 v1",
        "",
        "目標Rは独立した初期状態探索の誤差評価にだけ使用した。",
        "`run_ab_invariant_theta_toy_v1.py` は変更せず、前進計算には",
        "探索済みの初期A/B配列だけを渡した。",
        "",
        "| target R | initial B amplitude | AB-derived R | error | search iterations | invariant |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for item in result["targets"]:
        search = item["initial_state_search"]
        forward = item["forward_verification"]
        lines.append(
            "| {target:.15g} | {amplitude:.15g} | {derived:.15g} | "
            "{error:.3e} | {iterations} | {verdict} |".format(
                target=search["target_reflection_rate"],
                amplitude=search["initial_b_amplitude"],
                derived=forward["initial_AB_readout"]["reflection_rate"],
                error=search["absolute_error"],
                iterations=search["iterations"],
                verdict=forward["invariant_verdict"],
            )
        )
    lines.extend(
        [
            "",
            "探索は既存の `theta_from_ab(A,B)` の返り値を観測するブラックボックス二分探索であり、",
            "その式を探索側へ複製していない。前進処理は目標Rを引数に取らない。",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Search only the initial AB state, then run the unchanged AB-derived theta model"
        )
    )
    parser.add_argument(
        "--target-r",
        type=float,
        action="append",
        help="comparison target used only by the standalone initial-state search",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    targets = tuple(args.target_r) if args.target_r else DEFAULT_TARGETS
    params = toy.base.Params(
        high_n=63,
        recursive_collision_count=max(EVOLUTION_COLLISIONS),
    )
    source_params = toy.base.build_source_params(params)
    metric_context = toy.base.MetricContext(source_params)
    a_template, b_template, case = make_unit_templates(source_params)

    target_results: list[dict[str, Any]] = []
    for target_r in targets:
        search = search_initial_b_amplitude(
            target_r,
            a_template,
            b_template,
            source_params,
        )
        forward = run_forward_from_initial_state(
            search.initial_b_amplitude,
            source_params,
            metric_context,
            a_template,
            b_template,
            target_name(target_r),
        )
        target_results.append(
            {
                "initial_state_search": asdict(search),
                "forward_verification": forward,
            }
        )

    payload = {
        "experiment": "initial_state_search_v1",
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
            "collisions": list(EVOLUTION_COLLISIONS),
            "display_normalization": "each rho_chi divided by its own maximum",
            "forward_scattering_rule": (
                "unchanged real orthogonal AB rotation with theta re-derived from current AB"
            ),
            "alpha_denominator": ALPHA_DENOMINATOR,
            "R_alpha": R_ALPHA,
        },
        "search": {
            "method": "black-box bisection",
            "tolerance": SEARCH_TOLERANCE,
            "maximum_iterations": MAX_SEARCH_ITERATIONS,
        },
        "core_runner": {
            "path": str(TOY_RUNNER_PATH.relative_to(TOY_DIR)),
            "sha256": toy.sha256(TOY_RUNNER_PATH),
        },
        "targets": target_results,
    }
    result_path = HERE / "initial_state_search_result_v1.json"
    report_path = HERE / "initial_state_search_report_v1.md"
    result_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_report(report_path, payload)

    for item in target_results:
        search = item["initial_state_search"]
        forward = item["forward_verification"]
        print(
            f"target_R={search['target_reflection_rate']:.15g}",
            f"initial_B_amplitude={search['initial_b_amplitude']:.15g}",
            f"AB_derived_R={forward['initial_AB_readout']['reflection_rate']:.15g}",
            f"search_error={search['absolute_error']:.3e}",
            f"invariant={forward['invariant_verdict']}",
        )


if __name__ == "__main__":
    main()
