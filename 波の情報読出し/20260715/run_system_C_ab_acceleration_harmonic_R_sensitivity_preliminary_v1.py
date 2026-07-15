from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
SOURCE_PATH = BASE_DIR.parent / "20260711" / "run_ab_two_body_fermionic_reflection_harmonic_readout_v4.py"
DEFAULT_OUT_DIR = BASE_DIR / "system_C_ab_acceleration_harmonic_R_sensitivity_result_v1"
OUT_DIR = DEFAULT_OUT_DIR
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_source_module() -> Any:
    spec = importlib.util.spec_from_file_location("ab_acceleration_harmonic_v4", SOURCE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load source module: {SOURCE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


src = load_source_module()


def configure_output_dir(path: Path) -> None:
    global OUT_DIR
    OUT_DIR = path
    OUT_DIR.mkdir(parents=True, exist_ok=True)


@dataclass
class Params:
    coarse_r_min: float = 0.600
    coarse_r_max: float = 0.900
    coarse_r_step: float = 0.010
    fine_r_min: float = 0.680
    fine_r_max: float = 0.710
    fine_r_step: float = 0.001
    adaptive_refine_top_k: int = 8
    adaptive_refine_half_width: float = 0.012
    adaptive_refine_step: float = 0.001
    r137: float = 1.0 - math.sqrt(4.0 * math.pi / 137.035999084)
    r128_nominal: float = 1.0 - math.sqrt(4.0 * math.pi / 128.0)
    epsilon: float = 1.0e-300
    peak_top_k: int = 12
    depth_plot_cap: float = 6.0


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, float(value)))


def delta_from_reflection_rate(r_value: float) -> float:
    return 2.0 * math.asin(math.sqrt(clamp01(r_value)))


def plot_depth(value: float, params: Params) -> float:
    return max(0.0, min(float(value), params.depth_plot_cap))


def relative_depth(value: float, max_value: float, params: Params) -> float:
    denominator = max(float(max_value), params.epsilon)
    ratio = max(float(value), params.epsilon) / denominator
    return plot_depth(-math.log10(max(ratio, params.epsilon)), params)


def r_values(params: Params) -> List[float]:
    values = [0.0, 0.5, 1.0, params.r137, params.r128_nominal]
    coarse_steps = int(round((params.coarse_r_max - params.coarse_r_min) / params.coarse_r_step))
    values.extend(params.coarse_r_min + params.coarse_r_step * i for i in range(coarse_steps + 1))
    fine_steps = int(round((params.fine_r_max - params.fine_r_min) / params.fine_r_step))
    values.extend(params.fine_r_min + params.fine_r_step * i for i in range(fine_steps + 1))
    return sorted({round(float(value), 12) for value in values})


def parse_float_list(text: str) -> List[float]:
    if not text.strip():
        return []
    return sorted({round(float(token.strip()), 12) for token in text.split(",") if token.strip()})


def selected_r_values(params: Params, args: argparse.Namespace) -> List[float]:
    values = parse_float_list(args.r_values) if args.r_values else r_values(params)
    if args.r_min is not None:
        values = [value for value in values if value >= args.r_min]
    if args.r_max is not None:
        values = [value for value in values if value <= args.r_max]
    if not values:
        raise ValueError("no R values selected")
    return values


def adaptive_refinement_values(peaks: List[Dict[str, Any]], params: Params) -> List[float]:
    centers = [float(row["R_peak"]) for row in peaks[: params.adaptive_refine_top_k]]
    centers.extend([params.r137, params.r128_nominal])
    values: List[float] = []
    steps = int(round((2.0 * params.adaptive_refine_half_width) / params.adaptive_refine_step))
    for center in centers:
        lo = max(params.coarse_r_min, center - params.adaptive_refine_half_width)
        hi = min(params.coarse_r_max, center + params.adaptive_refine_half_width)
        for idx in range(steps + 1):
            value = lo + params.adaptive_refine_step * idx
            if value <= hi + 1.0e-15:
                values.append(value)
    return sorted({round(float(value), 12) for value in values})


def output_dir_from_args(args: argparse.Namespace) -> Path:
    if args.output_dir:
        return Path(args.output_dir).expanduser().resolve()
    if args.run_id:
        return DEFAULT_OUT_DIR / args.run_id
    return DEFAULT_OUT_DIR


def safe_slug(text: str, max_len: int = 120) -> str:
    allowed: List[str] = []
    for ch in text:
        if ch.isalnum() or ch in {"-", "_"}:
            allowed.append(ch)
        else:
            allowed.append("-")
    slug = "".join(allowed).strip("-_")
    while "--" in slug:
        slug = slug.replace("--", "-")
    if len(slug) <= max_len:
        return slug
    digest = hashlib.sha1(slug.encode("utf-8")).hexdigest()[:8]
    return f"{slug[: max_len - 9]}-{digest}"


def compact_float(value: float) -> str:
    return f"{value:.12g}".replace(".", "p").replace("-", "m")


def r_values_slug(values: List[float], params: Params) -> str:
    if values == r_values(params):
        return "Rdefault"
    if len(values) <= 4:
        return safe_slug("R" + "-".join(compact_float(value) for value in values), max_len=90)
    return safe_slug(f"R{compact_float(min(values))}-{compact_float(max(values))}_n{len(values)}", max_len=90)


def build_file_stem(values: List[float], params: Params, run_id: Optional[str], explicit_stem: Optional[str]) -> str:
    if explicit_stem:
        return safe_slug(explicit_stem, max_len=150)
    parts = ["system_C", r_values_slug(values, params)]
    if run_id:
        parts.insert(1, safe_slug(run_id, max_len=40))
    return safe_slug("_".join(parts), max_len=150)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="System C AB acceleration-like harmonic R sensitivity sweep")
    parser.add_argument("--r-min", type=float, help="minimum R to include")
    parser.add_argument("--r-max", type=float, help="maximum R to include")
    parser.add_argument("--r-values", help="comma-separated explicit R values")
    parser.add_argument("--run-id", help="write outputs under the default result directory with this subdirectory name")
    parser.add_argument("--file-stem", help="explicit output file stem")
    parser.add_argument("--output-dir", help="write outputs to this exact directory")
    parser.add_argument("--no-plots", action="store_true", help="skip png generation")
    return parser.parse_args()


def sweep_region(r_value: float, params: Params) -> str:
    if params.fine_r_min <= r_value <= params.fine_r_max:
        return "fine"
    if params.coarse_r_min <= r_value <= params.coarse_r_max:
        return "coarse"
    return "control"


def scattering_state_for_r(chi_pass: np.ndarray, complement: np.ndarray, r_value: float) -> Dict[str, Any]:
    delta_f = delta_from_reflection_rate(r_value)
    protocol = src.ScatteringProtocol("fermionic_partial_R", "fermionic_reflection", delta_f, True)
    return src.two_channel_scattering_state(protocol, chi_pass, complement)


def harmonic_rows_r(case: Any, readout_mode: Any, r_value: float, source_params: Any) -> List[Dict[str, Any]]:
    initial_deviation_rad = math.radians(case.deviation_deg)
    chi_pass, complement = src.base_rotation(initial_deviation_rad, readout_mode, source_params)
    scatter = scattering_state_for_r(chi_pass, complement, r_value)
    chi_read = scatter["chi_read"]
    complement_read = scatter["closure_complement"]
    envelope_values = scatter["envelope_V_AB"]
    event_mask = src.reflection_event_mask(chi_pass)
    delta_f = delta_from_reflection_rate(r_value)
    t_delta = scatter["t_delta"]
    r_delta = scatter["r_delta"]
    transmission_rate = scatter["transmission_rate"]
    reflection_rate = scatter["reflection_rate"]
    omega_discrete_sq = 4.0 * math.sin(source_params.omega_step / 2.0) ** 2
    f_circle = src.second_difference_abs(chi_read)
    rows: List[Dict[str, Any]] = []
    for step in range(source_params.step_count + 1):
        d_near, d_far, v_ab = src.label_free_readout(float(chi_read[step]))
        envelope_v = float(envelope_values[step])
        f_center = omega_discrete_sq * abs(float(chi_read[step]))
        q_raw = readout_mode.per_step_leak * envelope_v
        rows.append(
            {
                "case_id": case.case_id,
                "initial_deviation_deg": case.deviation_deg,
                "scattering_protocol": "fermionic_partial_R",
                "scattering_kind": "fermionic_partial_reflection",
                "scattering_matrix_used": True,
                "delta_f": delta_f,
                "t_delta_real": float(t_delta.real),
                "t_delta_imag": float(t_delta.imag),
                "r_delta_real": float(r_delta.real),
                "r_delta_imag": float(r_delta.imag),
                "transmission_rate": transmission_rate,
                "reflection_rate": reflection_rate,
                "q_out_factor": scatter["q_out_factor_diagnostic"],
                "q_out_factor_applied": False,
                "full_two_channel_scattering_used": scatter["full_two_channel_scattering_used"],
                "display_only_abs_reflection": scatter["display_only_abs_reflection"],
                "R_input": r_value,
                "readout_mode": readout_mode.name,
                "active_readout": readout_mode.active_readout,
                "step": step,
                "chi_pass": float(chi_pass[step]),
                "chi_read": float(chi_read[step]),
                "closure_complement": float(complement_read[step]),
                "psi_A_in_real": float(scatter["psi_a_in"][step].real),
                "psi_A_in_imag": float(scatter["psi_a_in"][step].imag),
                "psi_B_in_real": float(scatter["psi_b_in"][step].real),
                "psi_B_in_imag": float(scatter["psi_b_in"][step].imag),
                "psi_A_out_real": float(scatter["psi_a_out"][step].real),
                "psi_A_out_imag": float(scatter["psi_a_out"][step].imag),
                "psi_B_out_real": float(scatter["psi_b_out"][step].real),
                "psi_B_out_imag": float(scatter["psi_b_out"][step].imag),
                "channel_norm_in": float(scatter["channel_norm_in"][step]),
                "channel_norm_out": float(scatter["channel_norm_out"][step]),
                "scattering_unitarity_error": float(scatter["scattering_unitarity_error"][step]),
                "reflection_event_cell": bool(event_mask[step]),
                "D_AB_near_rad": d_near,
                "D_AB_far_rad": d_far,
                "D_AB_near_deg": math.degrees(d_near),
                "D_AB_far_deg": math.degrees(d_far),
                "V_AB": v_ab,
                "envelope_AB_abs": math.sqrt(envelope_v),
                "envelope_V_AB": envelope_v,
                "f_AB_center": f_center,
                "f_AB_circle": float(f_circle[step]),
                "f_AB_projection_consistency_error": abs(f_center - float(f_circle[step])),
                "Q_raw": q_raw,
                "Q_closed": 0.0,
                "closure_relaxation": q_raw,
                "absolute_background_axis_used": False,
                "f_A_or_f_B_used": False,
            }
        )
    return rows


def c1_rows_r(case: Any, tau_mode: Any, readout_mode: Any, r_value: float, source_params: Any) -> List[Dict[str, Any]]:
    initial_deviation_rad = math.radians(case.deviation_deg)
    chi_pass, complement = src.base_rotation(initial_deviation_rad, readout_mode, source_params)
    scatter = scattering_state_for_r(chi_pass, complement, r_value)
    chi_read = scatter["chi_read"]
    complement_read = scatter["closure_complement"]
    envelope_values = scatter["envelope_V_AB"]
    event_mask = src.reflection_event_mask(chi_pass)
    tau = src.tau_series(initial_deviation_rad, chi_read, readout_mode, tau_mode, source_params)
    tau_for_geometry = np.zeros_like(chi_read) if not tau_mode.tau_available else tau
    area = src.rolling_closed_area(chi_read, tau_for_geometry, source_params.period_steps)
    ranks = src.rolling_rank(chi_read, tau_for_geometry, source_params.period_steps, tau_mode.tau_available, source_params.rank_tol)
    c_errors = src.rolling_c_error(chi_read, tau_for_geometry, source_params.period_steps, tau_mode.tau_available)
    delta_f = delta_from_reflection_rate(r_value)
    transmission_rate = scatter["transmission_rate"]
    reflection_rate = scatter["reflection_rate"]
    rows: List[Dict[str, Any]] = []
    for step in range(source_params.step_count + 1):
        d_near, d_far, v_ab = src.label_free_readout(float(chi_read[step]))
        envelope_v = float(envelope_values[step])
        q_raw = readout_mode.per_step_leak * envelope_v * (1.0 + (0.0 if not tau_mode.tau_available else 0.25))
        rows.append(
            {
                "case_id": case.case_id,
                "initial_deviation_deg": case.deviation_deg,
                "scattering_protocol": "fermionic_partial_R",
                "scattering_kind": "fermionic_partial_reflection",
                "scattering_matrix_used": True,
                "delta_f": delta_f,
                "transmission_rate": transmission_rate,
                "reflection_rate": reflection_rate,
                "q_out_factor": scatter["q_out_factor_diagnostic"],
                "q_out_factor_applied": False,
                "full_two_channel_scattering_used": scatter["full_two_channel_scattering_used"],
                "display_only_abs_reflection": scatter["display_only_abs_reflection"],
                "R_input": r_value,
                "readout_mode": readout_mode.name,
                "tau_mode": tau_mode.name,
                "tau_kind": tau_mode.kind,
                "active_readout": readout_mode.active_readout,
                "step": step,
                "chi_pass": float(chi_pass[step]),
                "chi_read": float(chi_read[step]),
                "tau_read": "" if not tau_mode.tau_available else float(tau[step]),
                "tau_available": tau_mode.tau_available,
                "closure_complement": float(complement_read[step]),
                "psi_A_in_real": float(scatter["psi_a_in"][step].real),
                "psi_A_in_imag": float(scatter["psi_a_in"][step].imag),
                "psi_B_in_real": float(scatter["psi_b_in"][step].real),
                "psi_B_in_imag": float(scatter["psi_b_in"][step].imag),
                "psi_A_out_real": float(scatter["psi_a_out"][step].real),
                "psi_A_out_imag": float(scatter["psi_a_out"][step].imag),
                "psi_B_out_real": float(scatter["psi_b_out"][step].real),
                "psi_B_out_imag": float(scatter["psi_b_out"][step].imag),
                "channel_norm_in": float(scatter["channel_norm_in"][step]),
                "channel_norm_out": float(scatter["channel_norm_out"][step]),
                "scattering_unitarity_error": float(scatter["scattering_unitarity_error"][step]),
                "reflection_event_cell": bool(event_mask[step]),
                "D_AB_near_rad": d_near,
                "D_AB_far_rad": d_far,
                "V_AB": v_ab,
                "envelope_V_AB": envelope_v,
                "A_chi_tau": float(area[step]),
                "abs_A_chi_tau": abs(float(area[step])),
                "d_eff_chi_tau": int(ranks[step]),
                "epsilon_c": "" if not np.isfinite(c_errors[step]) else float(c_errors[step]),
                "Q_raw": q_raw,
                "Q_closed": 0.0,
                "tau_is_step_used": False,
                "external_c_used": False,
                "absolute_background_axis_used": False,
                "f_A_or_f_B_used": False,
            }
        )
    return rows


def bool_rate(values: Iterable[bool]) -> float:
    items = [bool(value) for value in values]
    if not items:
        return 0.0
    return sum(1 for value in items if value) / float(len(items))


def score_from_metrics(row: Dict[str, Any]) -> float:
    harmonic_penalty = 1.0 - float(row["harmonic_valid_rate_nonstrong"])
    area_penalty = 1.0 - float(row["c1_area_valid_rate"])
    c1_penalty = 1.0 - float(row["c1_calibrated_rate"])
    projection_penalty = math.log10(
        1.0
        + float(row["fermionic_max_f_AB_projection_error_regular_nonstrong"]) / float(row["f_consistency_tol"])
    )
    return harmonic_penalty + area_penalty + c1_penalty + projection_penalty


def classify_by_terrain(rows: List[Dict[str, Any]], best_row: Dict[str, Any], params: Params) -> str:
    sweep_rows = [row for row in rows if str(row["sweep_region"]) in {"coarse", "fine"}]
    if not sweep_rows:
        return "broken"
    scores = [float(row["score_C"]) for row in sweep_rows]
    min_score = min(scores)
    max_score = max(scores)
    if min_score > 3.0:
        return "broken"
    if max_score - min_score <= max(1.0e-9, abs(min_score) * 1.0e-3):
        return "flat"
    r_best = float(best_row["R"])
    if not (params.fine_r_min <= r_best <= params.fine_r_max):
        return "outside"
    local_min_count = len([row for row in peak_rows(rows, params) if "local_min" in str(row["peak_kind"])])
    if local_min_count >= 2:
        return "multi_peak"
    band = r_band_width(sweep_rows, min_score, 1.10)
    if band <= 0.005:
        return "sharp"
    return "broad"


def r_band_width(rows: List[Dict[str, Any]], min_score: float, multiplier: float) -> float:
    selected = [row for row in rows if float(row["score_C"]) <= min_score * multiplier + 1.0e-15]
    if not selected:
        return float("nan")
    r_list = [float(row["R"]) for row in selected]
    return max(r_list) - min(r_list)


def run_for_r(r_value: float, params: Params, source_params: Any) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    harmonic_summaries: List[Dict[str, Any]] = []
    condition_depth_rows: List[Dict[str, Any]] = []

    for case in src.INITIAL_CASES:
        for readout_mode in src.READOUT_MODES:
            rows = harmonic_rows_r(case, readout_mode, r_value, source_params)
            summary = src.summarize_harmonic_case(rows, source_params)
            harmonic_summaries.append(summary)
            error = (
                float(summary["max_f_AB_projection_consistency_error_regular_cells"]) / source_params.f_consistency_tol
                + (0.0 if bool(summary["regular_cell_harmonic_consistent"]) else 1.0)
            )
            condition_depth_rows.append(
                {
                    "sweep_region": sweep_region(r_value, params),
                    "R": r_value,
                    "case_id": case.case_id,
                    "case_group": f"harmonic|{case.case_id}",
                    "group_label": f"harmonic {case.case_id}",
                    "condition_kind": "harmonic",
                    "condition_id": f"{case.case_id}|{readout_mode.name}",
                    "readout_mode": readout_mode.name,
                    "tau_mode": "",
                    "condition_error": error,
                    "stability_depth": -math.log10(max(error, params.epsilon)),
                }
            )

    c1_summaries: List[Dict[str, Any]] = []
    for case in src.INITIAL_CASES:
        for tau_mode in src.TAU_MODES:
            for readout_mode in src.READOUT_MODES:
                rows = c1_rows_r(case, tau_mode, readout_mode, r_value, source_params)
                summary = src.summarize_c1_case(rows, tau_mode, source_params)
                c1_summaries.append(summary)
                c_error = float(summary["max_epsilon_c_abs"])
                error = (
                    c_error / source_params.c_tol
                    + (0.0 if bool(summary["area_sweep_detected"]) == bool(summary["area_expected"]) else 1.0)
                    + (0.0 if bool(summary["rank_matches_expectation"]) else 1.0)
                )
                condition_depth_rows.append(
                    {
                        "sweep_region": sweep_region(r_value, params),
                        "R": r_value,
                        "case_id": case.case_id,
                        "case_group": f"c1|{case.case_id}",
                        "group_label": f"c1 {case.case_id}",
                        "condition_kind": "c1",
                        "condition_id": f"{case.case_id}|{tau_mode.name}|{readout_mode.name}",
                        "readout_mode": readout_mode.name,
                        "tau_mode": tau_mode.name,
                        "condition_error": error,
                        "stability_depth": -math.log10(max(error, params.epsilon)),
                    }
                )

    harmonic_nonstrong = [row for row in harmonic_summaries if row["readout_mode"] != "readout_strong"]
    harmonic_strong = [row for row in harmonic_summaries if row["readout_mode"] == "readout_strong"]
    independent_c1 = [row for row in c1_summaries if row["tau_mode"] == "tau_independent_c1"]
    independent_off = [
        row for row in independent_c1 if row["readout_mode"] == "readout_off"
    ]
    delta_f = delta_from_reflection_rate(r_value)
    _t, _r, transmission_rate, reflection_rate = src.scattering_coefficients(delta_f)
    aggregate = {
        "sweep_region": sweep_region(r_value, params),
        "R": r_value,
        "T": float(transmission_rate),
        "Delta_F": delta_f,
        "q_out_factor": float(transmission_rate - reflection_rate),
        "q_out_factor_diagnostic": float(transmission_rate - reflection_rate),
        "q_out_factor_applied": False,
        "full_two_channel_scattering_used_rate": bool_rate(
            bool(row["full_two_channel_scattering_used"]) for row in harmonic_summaries + c1_summaries
        ),
        "max_scattering_unitarity_error": max(
            float(row["max_scattering_unitarity_error"]) for row in harmonic_summaries + c1_summaries
        ),
        "max_Q_closed_abs": max(float(row["max_Q_closed_abs"]) for row in harmonic_summaries + c1_summaries),
        "harmonic_case_count": len(harmonic_summaries),
        "c1_case_count": len(c1_summaries),
        "harmonic_valid_rate_all": bool_rate(row["regular_cell_harmonic_consistent"] for row in harmonic_summaries),
        "harmonic_valid_rate_nonstrong": bool_rate(row["regular_cell_harmonic_consistent"] for row in harmonic_nonstrong),
        "harmonic_valid_rate_strong": bool_rate(row["regular_cell_harmonic_consistent"] for row in harmonic_strong),
        "fermionic_regular_cell_harmonic_consistent_all_cases": src.bool_all(
            bool(row["regular_cell_harmonic_consistent"]) for row in harmonic_summaries
        ),
        "fermionic_regular_cell_harmonic_consistent_nonstrong_modes": src.bool_all(
            bool(row["regular_cell_harmonic_consistent"]) for row in harmonic_nonstrong
        ),
        "fermionic_strong_readout_perturbs_harmonic_projection": any(
            not bool(row["regular_cell_harmonic_consistent"]) for row in harmonic_strong
        ),
        "fermionic_max_f_AB_projection_error_regular_nonstrong": max(
            float(row["max_f_AB_projection_consistency_error_regular_cells"]) for row in harmonic_nonstrong
        ),
        "fermionic_max_f_AB_projection_error_regular_strong": max(
            float(row["max_f_AB_projection_consistency_error_regular_cells"]) for row in harmonic_strong
        ),
        "readout_off_decay_max_abs": max(
            abs(float(row["decay_rate_V_AB"])) for row in harmonic_summaries if row["readout_mode"] == "readout_off"
        ),
        "readout_strong_decay_min_abs": min(
            abs(float(row["decay_rate_V_AB"])) for row in harmonic_summaries if row["readout_mode"] == "readout_strong"
        ),
        "c1_area_valid_rate": bool_rate(
            bool(row["area_sweep_detected"]) for row in independent_c1
        ),
        "c1_calibrated_rate": bool_rate(
            bool(row["c1_calibrated"]) for row in independent_c1
        ),
        "c1_readout_off_max_epsilon_c_abs": max(
            float(row["max_epsilon_c_abs"]) for row in independent_off
        ),
        "fermionic_c1_area_sweep_detected_all_cases": src.bool_all(
            bool(row["area_sweep_detected"]) for row in independent_c1
        ),
        "tau_is_step_used_any": any(bool(row["tau_is_step_used"]) for row in c1_summaries),
        "external_c_used_any": any(bool(row["external_c_used"]) for row in c1_summaries),
        "f_A_or_f_B_used_any": any(bool(row["f_A_or_f_B_used"]) for row in c1_summaries),
        "f_consistency_tol": source_params.f_consistency_tol,
        "c_tol": source_params.c_tol,
        "distance_to_R_137": abs(r_value - params.r137),
        "distance_to_R_128_nominal": abs(r_value - params.r128_nominal),
    }
    aggregate["score_C"] = score_from_metrics(aggregate)
    aggregate["stability_depth"] = -math.log10(max(float(aggregate["score_C"]), params.epsilon))
    return aggregate, condition_depth_rows


def peak_kind_for_r(row: Dict[str, Any], params: Params, is_local: bool) -> str:
    r_value = float(row["R"])
    kinds: List[str] = []
    if is_local:
        kinds.append("local_min")
    if abs(r_value - params.r137) <= 5.0e-12:
        kinds.append("probe_R137")
    if abs(r_value - params.r128_nominal) <= 5.0e-12:
        kinds.append("probe_alpha128_nominal")
    return "+".join(kinds) if kinds else "selected"


def peak_rows(rows: List[Dict[str, Any]], params: Params) -> List[Dict[str, Any]]:
    sweep_rows = sorted(
        [row for row in rows if str(row["sweep_region"]) in {"coarse", "fine"}],
        key=lambda row: float(row["R"]),
    )
    candidates: Dict[float, Tuple[Dict[str, Any], bool]] = {}
    for i, row in enumerate(sweep_rows):
        score = float(row["score_C"])
        left = float(sweep_rows[i - 1]["score_C"]) if i > 0 else float("inf")
        right = float(sweep_rows[i + 1]["score_C"]) if i + 1 < len(sweep_rows) else float("inf")
        is_local = score <= left and score <= right and (score < left or score < right)
        is_probe = abs(float(row["R"]) - params.r137) <= 5.0e-12 or abs(float(row["R"]) - params.r128_nominal) <= 5.0e-12
        if is_local or is_probe:
            candidates[float(row["R"])] = (row, is_local)
    ranked = sorted(candidates.values(), key=lambda item: float(item[0]["score_C"]))
    selected = ranked[: params.peak_top_k]
    selected_r = {float(row["R"]) for row, _is_local in selected}
    for row, is_local in ranked:
        r_value = float(row["R"])
        is_probe = abs(r_value - params.r137) <= 5.0e-12 or abs(r_value - params.r128_nominal) <= 5.0e-12
        if is_probe and r_value not in selected_r:
            selected.append((row, is_local))
            selected_r.add(r_value)
    out: List[Dict[str, Any]] = []
    for rank, (row, is_local) in enumerate(selected, start=1):
        out.append(
            {
                "peak_rank": rank,
                "peak_kind": peak_kind_for_r(row, params, is_local),
                "R_peak": float(row["R"]),
                "score_C": float(row["score_C"]),
                "stability_depth": float(row["stability_depth"]),
                "q_out_factor": float(row["q_out_factor"]),
                "q_out_factor_diagnostic": float(row["q_out_factor_diagnostic"]),
                "q_out_factor_applied": bool(row["q_out_factor_applied"]),
                "full_two_channel_scattering_used_rate": float(row["full_two_channel_scattering_used_rate"]),
                "max_scattering_unitarity_error": float(row["max_scattering_unitarity_error"]),
                "harmonic_valid_rate_nonstrong": float(row["harmonic_valid_rate_nonstrong"]),
                "c1_area_valid_rate": float(row["c1_area_valid_rate"]),
                "c1_calibrated_rate": float(row["c1_calibrated_rate"]),
                "distance_to_R_137": abs(float(row["R"]) - params.r137),
                "distance_to_R_128_nominal": abs(float(row["R"]) - params.r128_nominal),
            }
        )
    return out


def best_rows(rows: List[Dict[str, Any]], params: Params) -> List[Dict[str, Any]]:
    sweep_rows = [row for row in rows if str(row["sweep_region"]) in {"coarse", "fine"}]
    best_pool = sweep_rows or rows
    best = min(best_pool, key=lambda row: float(row["score_C"]))
    min_score = float(best["score_C"])
    within_5 = r_band_width(best_pool, min_score, 1.05)
    within_10 = r_band_width(best_pool, min_score, 1.10)
    out = {
        "R_star_C": float(best["R"]),
        "score_C_min": min_score,
        "classification_C": classify_by_terrain(rows, best, params),
        "R_band_width_5": within_5,
        "R_band_width_10": within_10,
        "distance_to_R_137": abs(float(best["R"]) - params.r137),
        "distance_to_R_128_nominal": abs(float(best["R"]) - params.r128_nominal),
        "q_out_factor_at_star": float(best["q_out_factor"]),
        "q_out_factor_diagnostic_at_star": float(best["q_out_factor_diagnostic"]),
        "q_out_factor_applied_at_star": bool(best["q_out_factor_applied"]),
        "full_two_channel_scattering_used_rate_at_star": float(best["full_two_channel_scattering_used_rate"]),
        "max_scattering_unitarity_error_at_star": float(best["max_scattering_unitarity_error"]),
        "harmonic_valid_rate_nonstrong_at_star": float(best["harmonic_valid_rate_nonstrong"]),
        "c1_area_valid_rate_at_star": float(best["c1_area_valid_rate"]),
        "c1_calibrated_rate_at_star": float(best["c1_calibrated_rate"]),
    }
    return [out]


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def prepare_condition_terrain_rows(condition_rows: List[Dict[str, Any]], params: Params) -> List[Dict[str, Any]]:
    rows = [dict(row) for row in condition_rows]
    by_group: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        row["R"] = float(row["R"])
        row["condition_error"] = float(row["condition_error"])
        by_group.setdefault(str(row["case_group"]), []).append(row)
    for group_rows in by_group.values():
        max_error = max([float(row["condition_error"]) for row in group_rows] or [1.0])
        max_error = max(max_error, params.epsilon)
        for row in group_rows:
            terrain_score = float(row["condition_error"]) / max_error
            depth = -math.log10(max(terrain_score, params.epsilon))
            row["terrain_score"] = terrain_score
            row["log10_terrain_score"] = math.log10(max(terrain_score, params.epsilon))
            row["depth"] = depth
            row["depth_plot"] = plot_depth(depth, params)
    return rows


def score_components(row: Dict[str, Any]) -> Dict[str, float]:
    projection_penalty = math.log10(
        1.0
        + float(row["fermionic_max_f_AB_projection_error_regular_nonstrong"]) / float(row["f_consistency_tol"])
    )
    return {
        "harmonic_penalty": 1.0 - float(row["harmonic_valid_rate_nonstrong"]),
        "area_penalty": 1.0 - float(row["c1_area_valid_rate"]),
        "c1_penalty": 1.0 - float(row["c1_calibrated_rate"]),
        "projection_penalty": projection_penalty,
    }


def selected_plot_rows(
    rows: List[Dict[str, Any]],
    params: Params,
    deep: bool,
    focus_values: List[float],
    half_width: float = 0.020,
) -> List[Dict[str, Any]]:
    sweep_rows = sorted(
        [row for row in rows if str(row["sweep_region"]) in {"coarse", "fine"}],
        key=lambda row: float(row["R"]),
    )
    if not deep:
        return sweep_rows
    available_min = min(float(row["R"]) for row in sweep_rows)
    available_max = max(float(row["R"]) for row in sweep_rows)
    focus = list(focus_values)
    for probe in [params.r137, params.r128_nominal]:
        if available_min <= probe <= available_max:
            focus.append(probe)
    if not focus:
        focus = [(available_min + available_max) / 2.0]
    x_min = max(available_min, min(focus) - half_width)
    x_max = min(available_max, max(focus) + half_width)
    return [row for row in sweep_rows if x_min <= float(row["R"]) <= x_max]


def apply_r_axis_marks(ax: Any, x_values: List[float], params: Params) -> None:
    x_low = min(x_values) if x_values else params.coarse_r_min
    x_high = max(x_values) if x_values else params.coarse_r_max
    if x_low <= params.r137 <= x_high:
        ax.axvline(params.r137, color="black", linestyle=":", linewidth=1.1, label=None)
    if x_low <= params.r128_nominal <= x_high:
        ax.axvline(params.r128_nominal, color="0.35", linestyle="--", linewidth=1.1, label=None)
    ax.set_xlim(x_low, x_high)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.grid(alpha=0.22)


def annotate_projection_local_minima(
    ax: Any,
    plot_rows: List[Dict[str, Any]],
    projection_values: List[float],
    peaks: List[Dict[str, Any]],
) -> None:
    if not plot_rows or not projection_values:
        return
    local_minima = [peak for peak in peaks if "local_min" in str(peak["peak_kind"])][:3]
    if not local_minima:
        return
    x_values = [float(row["R"]) for row in plot_rows]
    x_min = min(x_values)
    x_max = max(x_values)
    y_min, y_max = ax.get_ylim()
    y_span = max(y_max - y_min, 1.0e-12)
    label_levels = [0.13, 0.28, 0.43]
    for i, peak in enumerate(local_minima):
        r_peak = float(peak["R_peak"])
        if not (x_min <= r_peak <= x_max):
            continue
        nearest_index = min(range(len(x_values)), key=lambda idx: abs(x_values[idx] - r_peak))
        y_peak = float(projection_values[nearest_index])
        y_label = y_min + label_levels[i % len(label_levels)] * y_span
        ax.axvline(r_peak, color="tab:orange", linestyle=":", linewidth=1.0, alpha=0.8)
        ax.scatter([r_peak], [y_peak], marker="D", color="tab:orange", s=34, zorder=5)
        ax.annotate(
            f"R={r_peak:.12f}",
            xy=(r_peak, y_peak),
            xytext=(r_peak + 0.0015, y_label),
            textcoords="data",
            fontsize=8,
            color="tab:orange",
            arrowprops={"arrowstyle": "-", "color": "tab:orange", "lw": 0.8},
            bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": "tab:orange", "alpha": 0.82},
        )


def make_projection_penalty_landscape_plot(rows: List[Dict[str, Any]], peaks: List[Dict[str, Any]], params: Params, file_stem: str, deep: bool) -> str:
    sweep_rows = [row for row in rows if str(row["sweep_region"]) in {"coarse", "fine"}]
    best_projection_row_all = min(sweep_rows, key=lambda row: score_components(row)["projection_penalty"])
    plot_rows = selected_plot_rows(rows, params, deep, [float(best_projection_row_all["R"])])
    component_rows = [(row, score_components(row)) for row in plot_rows]
    best_projection_row = min(plot_rows, key=lambda row: score_components(row)["projection_penalty"])

    x_values = [float(row["R"]) for row in plot_rows]
    projection_values = [components["projection_penalty"] for _row, components in component_rows]

    fig, ax = plt.subplots(figsize=(12, 4.8), constrained_layout=True)
    ax.plot(x_values, projection_values, color="black", linewidth=1.6, label="projection penalty")
    ax.scatter(x_values, projection_values, color="0.65", s=18, alpha=0.55, linewidths=0)
    ax.scatter(
        [float(best_projection_row["R"])],
        [score_components(best_projection_row)["projection_penalty"]],
        marker="x",
        s=80,
        color="black",
        linewidths=1.6,
    )
    if projection_values:
        y_min = min(projection_values)
        y_max = max(projection_values)
        padding = max((y_max - y_min) * 0.08, 1.0e-7)
        ax.set_ylim(y_min - padding, y_max + padding)
    apply_r_axis_marks(ax, x_values, params)
    annotate_projection_local_minima(ax, plot_rows, projection_values, peaks)
    ax.set_title("System C projection penalty landscape")
    ax.set_xlabel("R")
    ax.set_ylabel("projection penalty\nlower is better")
    ax.legend(loc="best", fontsize=8)
    suffix = "deep" if deep else "overview"
    filename = f"{file_stem}_stability_depth_distribution_{suffix}_v1.png"
    fig.savefig(OUT_DIR / filename, dpi=170)
    plt.close(fig)
    return filename


def make_score_plot(condition_rows: List[Dict[str, Any]], rows: List[Dict[str, Any]], peaks: List[Dict[str, Any]], params: Params, file_stem: str) -> str:
    plot_conditions = [
        row
        for row in condition_rows
        if params.coarse_r_min <= float(row["R"]) <= params.coarse_r_max
    ]
    plot_rows = sorted(
        [row for row in rows if str(row["sweep_region"]) in {"coarse", "fine"}],
        key=lambda row: float(row["R"]),
    )
    condition_max_error = max([float(row["condition_error"]) for row in plot_conditions] or [1.0])
    aggregate_max_score = max([float(row["score_C"]) for row in plot_rows] or [1.0])
    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.scatter(
        [float(row["R"]) for row in plot_conditions],
        [relative_depth(float(row["condition_error"]), condition_max_error, params) for row in plot_conditions],
        color="0.68",
        alpha=0.22,
        s=14,
        linewidths=0,
        label="conditions",
    )
    ax.plot(
        [float(row["R"]) for row in plot_rows],
        [relative_depth(float(row["score_C"]), aggregate_max_score, params) for row in plot_rows],
        color="black",
        linewidth=1.4,
        label="aggregate score depth",
    )
    if peaks:
        ax.scatter(
            [row["R_peak"] for row in peaks if int(row["peak_rank"]) <= 8],
            [relative_depth(float(row["score_C"]), aggregate_max_score, params) for row in peaks if int(row["peak_rank"]) <= 8],
            color="tab:orange",
            s=42,
            zorder=4,
            label="local/probe minima",
        )
    ax.axvline(params.r137, color="tab:red", linestyle=":", linewidth=1.2, label="R_137")
    ax.axvline(params.r128_nominal, color="tab:purple", linestyle="--", linewidth=1.2, label="alpha~128 nominal")
    ax.set_title("System C AB acceleration-like harmonic readout R-depth distribution")
    ax.set_xlabel("R")
    ax.set_ylabel("depth = -log10(relative error), capped at 6")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    filename = f"{file_stem}_scores_v1.png"
    fig.savefig(OUT_DIR / filename, dpi=170)
    plt.close(fig)
    return filename


def make_peak_zoom_plot(condition_rows: List[Dict[str, Any]], rows: List[Dict[str, Any]], peaks: List[Dict[str, Any]], params: Params, file_stem: str) -> str:
    focus = [float(row["R_peak"]) for row in peaks[:8]]
    focus.extend([params.r137, params.r128_nominal])
    x_min = max(params.coarse_r_min, min(focus) - 0.015)
    x_max = min(params.coarse_r_max, max(focus) + 0.015)
    plot_conditions = [
        row
        for row in condition_rows
        if x_min <= float(row["R"]) <= x_max
    ]
    sweep_conditions = [row for row in condition_rows if params.coarse_r_min <= float(row["R"]) <= params.coarse_r_max]
    sweep_rows = [row for row in rows if str(row["sweep_region"]) in {"coarse", "fine"}]
    condition_max_error = max([float(row["condition_error"]) for row in sweep_conditions] or [1.0])
    aggregate_max_score = max([float(row["score_C"]) for row in sweep_rows] or [1.0])
    plot_rows = sorted(
        [row for row in rows if str(row["sweep_region"]) in {"coarse", "fine"} and x_min <= float(row["R"]) <= x_max],
        key=lambda row: float(row["R"]),
    )
    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.scatter(
        [float(row["R"]) for row in plot_conditions],
        [relative_depth(float(row["condition_error"]), condition_max_error, params) for row in plot_conditions],
        color="0.68",
        alpha=0.22,
        s=14,
        linewidths=0,
        label="conditions",
    )
    ax.plot(
        [float(row["R"]) for row in plot_rows],
        [relative_depth(float(row["score_C"]), aggregate_max_score, params) for row in plot_rows],
        color="black",
        linewidth=1.5,
        label="aggregate score depth",
    )
    ax.axvline(params.r137, color="tab:red", linestyle=":", linewidth=1.2)
    ax.axvline(params.r128_nominal, color="tab:purple", linestyle="--", linewidth=1.2)
    for i, peak in enumerate(peaks[:6]):
        r_peak = float(peak["R_peak"])
        if not (x_min <= r_peak <= x_max):
            continue
        y_peak = relative_depth(float(peak["score_C"]), aggregate_max_score, params)
        color = "tab:purple" if "alpha128" in str(peak["peak_kind"]) else "tab:orange"
        ax.scatter([r_peak], [y_peak], marker="D", color=color, s=46, zorder=5)
        ax.annotate(
            f"#{int(peak['peak_rank'])} {peak['peak_kind']}\nR={r_peak:.12g}\nscore={float(peak['score_C']):.6g}",
            xy=(r_peak, y_peak),
            xytext=(8, 10 + 18 * (i % 3)),
            textcoords="offset points",
            fontsize=7,
            color=color,
            arrowprops={"arrowstyle": "-", "color": color, "lw": 0.8},
        )
    ax.set_title("System C peak distribution around top local/probe minima")
    ax.set_xlabel("R")
    ax.set_ylabel("depth = -log10(relative error), capped at 6")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    filename = f"{file_stem}_peak_zoom_v1.png"
    fig.savefig(OUT_DIR / filename, dpi=170)
    plt.close(fig)
    return filename


def write_report(result: Dict[str, Any], file_stem: str) -> str:
    best = result["best"][0]
    peaks = result["peaks"]
    lines = [
        "# 系統C AB加速度様調和読出しR感度スイープ 予備実験レポート",
        "",
        "## 代表結果",
        "",
        "| 量 | 値 |",
        "|---|---:|",
        f"| `R_star_C` | `{best['R_star_C']}` |",
        f"| `score_C_min` | `{best['score_C_min']}` |",
        f"| `classification_C` | `{best['classification_C']}` |",
        f"| `q_out_factor_diagnostic_at_star` | `{best['q_out_factor_diagnostic_at_star']}` |",
        f"| `q_out_factor_applied_at_star` | `{best['q_out_factor_applied_at_star']}` |",
        f"| `full_two_channel_scattering_used_rate_at_star` | `{best['full_two_channel_scattering_used_rate_at_star']}` |",
        f"| `max_scattering_unitarity_error_at_star` | `{best['max_scattering_unitarity_error_at_star']}` |",
        f"| `distance_to_R_137` | `{best['distance_to_R_137']}` |",
        f"| `distance_to_R_128_nominal` | `{best['distance_to_R_128_nominal']}` |",
        f"| `harmonic_valid_rate_nonstrong_at_star` | `{best['harmonic_valid_rate_nonstrong_at_star']}` |",
        f"| `c1_area_valid_rate_at_star` | `{best['c1_area_valid_rate_at_star']}` |",
        f"| `c1_calibrated_rate_at_star` | `{best['c1_calibrated_rate_at_star']}` |",
        "",
        "## 上位ピーク",
        "",
        "| rank | kind | R | score_C | q_out_factor_diagnostic |",
        "|---:|---|---:|---:|---:|",
    ]
    for row in peaks[:12]:
        lines.append(
            f"| {row['peak_rank']} | `{row['peak_kind']}` | `{row['R_peak']}` | `{row['score_C']}` | `{row['q_out_factor_diagnostic']}` |"
        )
    lines.extend(
        [
            "",
            "## 読み方",
            "",
            "本実験は、20260711 のAB二体加速度様調和読出しに部分反射率 `R` を導入し、調和読出しと `c=1` 面積読出しがどの程度 `R` に敏感かを確認する。",
            "",
            "`q_out_factor_diagnostic` は透過率と反射率から読まれる診断量であり、演算子として `chi_read` へ掛けていない。実装ではA/B二チャネル散乱行列を入射チャネルへ作用させる。",
            "",
            "`score_C` は、調和読出し、c=1面積読出し、射影誤差を合わせた管理用集約値である。物理的なR地形は `projection_penalty` 図で読む。",
            "",
            "`R_128_nominal` は高エネルギー側で `1/alpha` が128近傍へ走ることを読むための名目プローブであり、精密測定値そのものではない。",
            "",
            "## 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
            f"| rows | `{file_stem}_rows_v1.csv` |",
            f"| condition rows | `{file_stem}_condition_rows_v1.csv` |",
            f"| summary | `{file_stem}_summary_v1.csv` |",
            f"| best | `{file_stem}_best_v1.csv` |",
            f"| peaks | `{file_stem}_peaks_v1.csv` |",
            f"| result | `{file_stem}_result_v1.json` |",
            f"| scores | `{file_stem}_scores_v1.png` |",
            f"| projection landscape overview | `{file_stem}_stability_depth_distribution_overview_v1.png` |",
            f"| projection landscape deep | `{file_stem}_stability_depth_distribution_deep_v1.png` |",
            f"| peak zoom | `{file_stem}_peak_zoom_v1.png` |",
        ]
    )
    filename = f"{file_stem}_report_v1.md"
    (OUT_DIR / filename).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return filename


def run(args: argparse.Namespace) -> Dict[str, Any]:
    params = Params()
    source_params = src.Params()
    values = selected_r_values(params, args)
    configure_output_dir(output_dir_from_args(args))
    file_stem = build_file_stem(values, params, args.run_id, args.file_stem)
    rows_by_r: Dict[float, Dict[str, Any]] = {}
    condition_rows_by_r: Dict[float, List[Dict[str, Any]]] = {}

    def run_values(target_values: List[float]) -> None:
        for r_value in target_values:
            key = round(float(r_value), 12)
            if key in rows_by_r:
                continue
            aggregate, condition_depth = run_for_r(key, params, source_params)
            rows_by_r[key] = aggregate
            condition_rows_by_r[key] = condition_depth

    run_values(values)
    initial_rows = sorted(rows_by_r.values(), key=lambda row: float(row["R"]))
    initial_peaks = peak_rows(initial_rows, params)
    if not args.r_values:
        run_values(adaptive_refinement_values(initial_peaks, params))

    rows = sorted(rows_by_r.values(), key=lambda row: float(row["R"]))
    condition_rows = [
        row
        for r_value in sorted(condition_rows_by_r.keys())
        for row in condition_rows_by_r[r_value]
    ]
    condition_rows = prepare_condition_terrain_rows(condition_rows, params)
    peaks = peak_rows(rows, params)
    best = best_rows(rows, params)
    for row in rows:
        row["classification_C"] = best[0]["classification_C"]
    if not args.no_plots:
        make_score_plot(condition_rows, rows, peaks, params, file_stem)
        make_projection_penalty_landscape_plot(rows, peaks, params, file_stem, deep=False)
        make_projection_penalty_landscape_plot(rows, peaks, params, file_stem, deep=True)
        make_peak_zoom_plot(condition_rows, rows, peaks, params, file_stem)
    result = {
        "experiment": "system_C_ab_acceleration_harmonic_R_sensitivity_v1",
        "params": asdict(params),
        "source": str(SOURCE_PATH.relative_to(BASE_DIR.parent)),
        "file_stem": file_stem,
        "initial_R_count": len(values),
        "final_R_count": len(rows),
        "adaptive_refinement_used": not bool(args.r_values),
        "rows": rows,
        "condition_rows": condition_rows,
        "best": best,
        "peaks": peaks,
    }
    write_csv(OUT_DIR / f"{file_stem}_rows_v1.csv", rows)
    write_csv(OUT_DIR / f"{file_stem}_condition_rows_v1.csv", condition_rows)
    write_csv(OUT_DIR / f"{file_stem}_summary_v1.csv", rows)
    write_csv(OUT_DIR / f"{file_stem}_best_v1.csv", best)
    write_csv(OUT_DIR / f"{file_stem}_peaks_v1.csv", peaks)
    (OUT_DIR / f"{file_stem}_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_report(result, file_stem)
    return result


def main() -> None:
    args = parse_args()
    result = run(args)
    print(json.dumps({"best": result["best"], "top_peaks": result["peaks"][:5]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
