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
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
SYSTEM_A_PATH = BASE_DIR / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py"
DEFAULT_OUT_DIR = BASE_DIR / "system_B_gray_cat_metastable_R_sweep_result_v1"
OUT_DIR = DEFAULT_OUT_DIR
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_system_a_module() -> Any:
    spec = importlib.util.spec_from_file_location("system_a_r_sweep_v1", SYSTEM_A_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load system A module: {SYSTEM_A_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


system_a = load_system_a_module()
src = system_a.src


def configure_output_dir(path: Path) -> None:
    global OUT_DIR
    OUT_DIR = path
    OUT_DIR.mkdir(parents=True, exist_ok=True)


@dataclass
class Params:
    steps: int = 4096
    tail_fraction: float = 0.5
    s_gray_limit: float = 0.05
    eigen_amp_tol: float = 1.0e-4
    mean_tol: float = 5.0e-3
    drift_tol: float = 1.0e-2
    selection_limit: float = 0.95
    q_tol: float = 1.0e-10
    metastable_amp_target: float = 0.02
    coarse_r_min: float = 0.600
    coarse_r_max: float = 0.900
    coarse_r_step: float = 0.010
    fine_r_min: float = 0.680
    fine_r_max: float = 0.710
    fine_r_step: float = 0.001
    r137: float = 1.0 - math.sqrt(4.0 * math.pi / 137.035999084)
    r128: float = 1.0 - math.sqrt(4.0 * math.pi / 128.0)
    pairs: Tuple[Tuple[int, int], ...] = (
        (1, 1),
        (1, 2),
        (1, 3),
        (1, 5),
        (1, 15),
        (1, 63),
    )
    phi_values: Tuple[float, ...] = (0.0, math.pi / 12.0, math.pi / 2.0, math.pi)
    s0_values: Tuple[float, ...] = (0.0, 1.0e-3, 1.0e-2)
    stability_gain_values: Tuple[float, ...] = (-2.0e-3, 0.0, 2.0e-3)
    noise_values: Tuple[float, ...] = (0.0,)
    seed_values: Tuple[int, ...] = (0,)
    peak_top_k: int = 12


def normalize_pair(a: complex, b: complex) -> Tuple[complex, complex]:
    q = abs(a) ** 2 + abs(b) ** 2
    if q <= 0.0:
        raise ValueError("zero AB norm")
    scale = 1.0 / math.sqrt(q)
    return a * scale, b * scale


def state_from_s_phi(s: float, phi: float) -> Tuple[complex, complex]:
    if abs(s) >= 0.5:
        raise ValueError(f"|s| must be < 0.5: {s}")
    a = math.sqrt(0.5 + s)
    b = math.sqrt(0.5 - s) * complex(math.cos(phi), math.sin(phi))
    return normalize_pair(a, b)


def apply_scattering(a: complex, b: complex, r_value: float) -> Tuple[complex, complex]:
    delta_f = src.delta_from_reflection_rate(r_value)
    t, r, _t_power, _r_power = src.scattering_coefficients(delta_f)
    return normalize_pair(r * a + t * b, t * a + r * b)


def apply_stability_gain(a: complex, b: complex, gain: float) -> Tuple[complex, complex]:
    if gain == 0.0:
        return a, b
    p_a = abs(a) ** 2
    p_b = abs(b) ** 2
    q = p_a + p_b
    if q <= 0.0:
        raise ValueError("zero AB norm")
    s = (p_a - p_b) / q
    s_next = s + gain * s * (1.0 - s * s)
    s_next = max(min(s_next, 1.0 - 1.0e-12), -1.0 + 1.0e-12)
    phase_a = a / abs(a) if abs(a) > 0.0 else 1.0 + 0.0j
    phase_b = b / abs(b) if abs(b) > 0.0 else 1.0 + 0.0j
    a_next = math.sqrt(0.5 * (1.0 + s_next)) * phase_a
    b_next = math.sqrt(0.5 * (1.0 - s_next)) * phase_b
    return normalize_pair(a_next, b_next)


def apply_noise(a: complex, b: complex, noise_amp: float, rng: np.random.Generator) -> Tuple[complex, complex]:
    if noise_amp <= 0.0:
        return a, b
    da = noise_amp * (rng.normal() + 1j * rng.normal())
    db = noise_amp * (rng.normal() + 1j * rng.normal())
    return normalize_pair(a + da, b + db)


def r_values(params: Params) -> List[float]:
    values = [0.0, 0.5, 1.0, params.r137, params.r128]
    coarse_steps = int(round((params.coarse_r_max - params.coarse_r_min) / params.coarse_r_step))
    values.extend(params.coarse_r_min + params.coarse_r_step * i for i in range(coarse_steps + 1))
    fine_steps = int(round((params.fine_r_max - params.fine_r_min) / params.fine_r_step))
    values.extend(params.fine_r_min + params.fine_r_step * i for i in range(fine_steps + 1))
    return sorted({round(float(v), 12) for v in values})


def selected_cases_from_args(params: Params, args: argparse.Namespace) -> Tuple[Any, ...]:
    proxy_params = system_a.Params(
        recursive_collision_count=params.steps,
        coarse_r_min=params.coarse_r_min,
        coarse_r_max=params.coarse_r_max,
        coarse_r_step=params.coarse_r_step,
        fine_r_min=params.fine_r_min,
        fine_r_max=params.fine_r_max,
        fine_r_step=params.fine_r_step,
        r137=params.r137,
        r128=params.r128,
        pairs=params.pairs,
    )
    return system_a.selected_cases_from_args(proxy_params, args)


def selected_r_values_from_args(params: Params, args: argparse.Namespace) -> List[float]:
    values = system_a.parse_float_list(args.r_values) if args.r_values else r_values(params)
    if args.r_min is not None:
        values = [value for value in values if value >= args.r_min]
    if args.r_max is not None:
        values = [value for value in values if value <= args.r_max]
    if not values:
        raise ValueError("no R values selected")
    return values


def output_dir_from_args(args: argparse.Namespace) -> Path:
    if args.output_dir:
        return Path(args.output_dir).expanduser().resolve()
    if args.run_id:
        return DEFAULT_OUT_DIR / args.run_id
    return DEFAULT_OUT_DIR


def safe_slug(text: str, max_len: int = 80) -> str:
    return system_a.safe_slug(text, max_len=max_len)


def build_file_stem(cases: Tuple[Any, ...], values: List[float], params: Params, run_id: Optional[str], explicit_stem: Optional[str]) -> str:
    if explicit_stem:
        return safe_slug(explicit_stem, max_len=120)
    proxy_params = system_a.Params(
        recursive_collision_count=params.steps,
        coarse_r_min=params.coarse_r_min,
        coarse_r_max=params.coarse_r_max,
        coarse_r_step=params.coarse_r_step,
        fine_r_min=params.fine_r_min,
        fine_r_max=params.fine_r_max,
        fine_r_step=params.fine_r_step,
        r137=params.r137,
        r128=params.r128,
        pairs=params.pairs,
    )
    parts = ["system_B", system_a.cases_slug(cases), system_a.r_values_slug(values, proxy_params), f"S{params.steps}"]
    if run_id:
        parts.insert(1, safe_slug(run_id, max_len=40))
    return safe_slug("_".join(parts), max_len=150)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="System B gray-cat metastable R sweep")
    parser.add_argument("--pairs", help="comma-separated pair list, e.g. 1:2,1:3,1:63")
    parser.add_argument("--pair-start", type=int, help="zero-based start index in default pair list")
    parser.add_argument("--pair-stop", type=int, help="zero-based stop index in default pair list")
    parser.add_argument(
        "--harmonic-mode",
        choices=("single", "even-packet", "alternating-packet-1", "alternating-packet-2", "all"),
        default="single",
        help="harmonic construction mode",
    )
    parser.add_argument("--packet-a", help="custom explicit harmonic packet for A, e.g. 1 or 1,3,5")
    parser.add_argument("--packet-b", help="custom explicit harmonic packet for B, e.g. 1,2,4,6")
    parser.add_argument("--r-min", type=float, help="minimum R to include")
    parser.add_argument("--r-max", type=float, help="maximum R to include")
    parser.add_argument("--r-values", help="comma-separated explicit R values")
    parser.add_argument("--max-steps", type=int, help="override time step count")
    parser.add_argument("--run-id", help="write outputs under the default result directory with this subdirectory name")
    parser.add_argument("--file-stem", help="explicit output file stem")
    parser.add_argument("--output-dir", help="write outputs to this exact directory")
    parser.add_argument("--no-plots", action="store_true", help="skip png generation")
    return parser.parse_args()


def classify_series(params: Params, s_values: np.ndarray, q_values: np.ndarray) -> Dict[str, Any]:
    tail_start = int(round(params.steps * (1.0 - params.tail_fraction)))
    tail = s_values[tail_start:]
    q_tail = q_values[tail_start:]
    first_half_tail = tail[: len(tail) // 2]
    second_half_tail = tail[len(tail) // 2 :]
    s_mean = float(np.mean(tail))
    s_amp = float(0.5 * (np.max(tail) - np.min(tail)))
    s_max_abs = float(np.max(np.abs(tail)))
    s_final = float(s_values[-1])
    s_drift = float(abs(np.mean(second_half_tail) - np.mean(first_half_tail))) if len(first_half_tail) else 0.0
    q_max_error = float(np.max(np.abs(q_values - 1.0)))
    q_tail_max_error = float(np.max(np.abs(q_tail - 1.0)))

    if q_max_error > params.q_tol:
        phase = "norm_error"
    elif abs(s_mean) >= params.selection_limit and s_amp < params.s_gray_limit and abs(s_final) >= params.selection_limit:
        phase = "natural_selection"
    elif abs(s_mean) <= params.mean_tol and s_amp <= params.eigen_amp_tol and s_drift <= params.drift_tol:
        phase = "gray_eigen"
    elif abs(s_mean) <= params.s_gray_limit and params.eigen_amp_tol < s_amp < params.s_gray_limit and s_drift <= params.drift_tol:
        phase = "gray_metastable"
    elif s_amp >= params.s_gray_limit and s_drift <= params.drift_tol:
        phase = "large_oscillation"
    else:
        phase = "unstable_or_drifting"

    if phase == "gray_metastable":
        phase_penalty = 0.0
    elif phase == "gray_eigen":
        phase_penalty = 0.25
    elif phase == "large_oscillation":
        phase_penalty = 0.5
    else:
        phase_penalty = 1.0
    gray_error = abs(s_mean) + abs(s_amp - params.metastable_amp_target) + s_drift + phase_penalty
    gray_depth = -math.log10(max(gray_error, 1.0e-300))
    return {
        "phase": phase,
        "S_mean": s_mean,
        "S_amp": s_amp,
        "S_max_abs": s_max_abs,
        "S_final": s_final,
        "S_drift": s_drift,
        "Q_max_error": q_max_error,
        "Q_tail_max_error": q_tail_max_error,
        "gray_error": gray_error,
        "gray_depth": gray_depth,
        "gray_depth_plot": min(gray_depth, 6.0),
    }


def run_condition(
    params: Params,
    case: Any,
    r_value: float,
    phi: float,
    s0: float,
    stability_gain: float,
    noise_amp: float,
    seed: int,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    rng = np.random.default_rng(seed)
    a, b = state_from_s_phi(s0, phi)
    s_values = np.empty(params.steps + 1)
    q_values = np.empty(params.steps + 1)
    samples: List[Dict[str, Any]] = []
    sample_steps = {0, 1, 2, 3, 5, 10, 20, 42, 64, 128, 256, params.steps}
    delta_f = src.delta_from_reflection_rate(r_value)
    _t, _r, t_power, r_power = src.scattering_coefficients(delta_f)
    condition_id = f"phi{phi / math.pi:.6g}_s{s0:.6g}_g{stability_gain:.6g}_noise{noise_amp:.6g}_seed{seed}"

    for step in range(params.steps + 1):
        p_a = abs(a) ** 2
        p_b = abs(b) ** 2
        q = p_a + p_b
        s = (p_a - p_b) / q
        s_values[step] = s
        q_values[step] = q
        if step in sample_steps:
            samples.append(
                {
                    "case_id": case.case_id,
                    "mode": case.mode,
                    "harmonic_packet_A": case.packet_a_text,
                    "harmonic_packet_B": case.packet_b_text,
                    "N_A": case.n_a,
                    "N_B": case.n_b,
                    "R": r_value,
                    "T": float(t_power),
                    "Delta_F": delta_f,
                    "condition_id": condition_id,
                    "step": step,
                    "p_A": p_a / q,
                    "p_B": p_b / q,
                    "S": s,
                    "Q": q,
                }
            )
        if step >= params.steps:
            break
        a, b = apply_scattering(a, b, r_value)
        a, b = apply_stability_gain(a, b, stability_gain)
        a, b = apply_noise(a, b, noise_amp, rng)

    metrics = classify_series(params, s_values, q_values)
    if params.fine_r_min <= r_value <= params.fine_r_max:
        sweep_region = "fine"
    elif params.coarse_r_min <= r_value <= params.coarse_r_max:
        sweep_region = "coarse"
    else:
        sweep_region = "control"
    row = {
        "case_id": case.case_id,
        "mode": case.mode,
        "state_family": case.state_family,
        "harmonic_packet_A": case.packet_a_text,
        "harmonic_packet_B": case.packet_b_text,
        "N_A": case.n_a,
        "N_B": case.n_b,
        "pair_kind": "one_side_high_harmonic" if case.n_a != case.n_b else "same_harmonic_control",
        "R": r_value,
        "T": float(t_power),
        "Delta_F": delta_f,
        "sweep_region": sweep_region,
        "condition_id": condition_id,
        "phi": phi,
        "phi_over_pi": phi / math.pi,
        "s0": s0,
        "stability_gain": stability_gain,
        "noise_amp": noise_amp,
        "seed": seed,
        "p_A_final": 0.5 * (1.0 + metrics["S_final"]),
        "p_B_final": 0.5 * (1.0 - metrics["S_final"]),
        **metrics,
    }
    return row, samples


def run_case(params: Params, case: Any, r_value: float) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    rows: List[Dict[str, Any]] = []
    samples: List[Dict[str, Any]] = []
    for phi in params.phi_values:
        for s0 in params.s0_values:
            for stability_gain in params.stability_gain_values:
                for noise_amp in params.noise_values:
                    seeds = params.seed_values if noise_amp > 0.0 else (0,)
                    for seed in seeds:
                        row, sample = run_condition(params, case, r_value, phi, s0, stability_gain, noise_amp, seed)
                        rows.append(row)
                        samples.extend(sample)
    return rows, samples


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def md_cell(value: Any) -> str:
    return str(value).replace("|", "\\|")


def aggregate_by_case_r(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, float], List[Dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault((str(row["case_id"]), float(row["R"])), []).append(row)
    out: List[Dict[str, Any]] = []
    for (case_id, r_value), case_rows in sorted(grouped.items()):
        first = case_rows[0]
        phase_counts: Dict[str, int] = {}
        for row in case_rows:
            phase_counts[str(row["phase"])] = phase_counts.get(str(row["phase"]), 0) + 1
        best = min(case_rows, key=lambda item: float(item["gray_error"]))
        gray_count = phase_counts.get("gray_eigen", 0) + phase_counts.get("gray_metastable", 0)
        selection_count = phase_counts.get("natural_selection", 0)
        joint_gray_score = float(best["gray_depth"]) + 0.25 * gray_count - 0.25 * selection_count
        out.append(
            {
                "case_id": case_id,
                "mode": str(first["mode"]),
                "state_family": str(first["state_family"]),
                "harmonic_packet_A": str(first["harmonic_packet_A"]),
                "harmonic_packet_B": str(first["harmonic_packet_B"]),
                "N_A": int(first["N_A"]),
                "N_B": int(first["N_B"]),
                "pair_kind": str(first["pair_kind"]),
                "R": r_value,
                "T": float(first["T"]),
                "Delta_F": float(first["Delta_F"]),
                "sweep_region": str(first["sweep_region"]),
                "gray_eigen_count": phase_counts.get("gray_eigen", 0),
                "gray_metastable_count": phase_counts.get("gray_metastable", 0),
                "large_oscillation_count": phase_counts.get("large_oscillation", 0),
                "natural_selection_count": phase_counts.get("natural_selection", 0),
                "unstable_or_drifting_count": phase_counts.get("unstable_or_drifting", 0),
                "norm_error_count": phase_counts.get("norm_error", 0),
                "best_condition_id": str(best["condition_id"]),
                "best_phase": str(best["phase"]),
                "best_gray_depth": float(best["gray_depth"]),
                "best_gray_error": float(best["gray_error"]),
                "best_S_mean": float(best["S_mean"]),
                "best_S_amp": float(best["S_amp"]),
                "best_S_drift": float(best["S_drift"]),
                "joint_gray_score": joint_gray_score,
            }
        )
    return out


def best_rows_for_case(rows: List[Dict[str, Any]], params: Params) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for case_id in sorted({str(row["case_id"]) for row in rows}):
        all_case_rows = [row for row in rows if str(row["case_id"]) == case_id]
        sweep_rows = [row for row in all_case_rows if str(row["sweep_region"]) in {"coarse", "fine"}]
        if not sweep_rows:
            sweep_rows = all_case_rows
        best_joint = max(sweep_rows, key=lambda row: float(row["joint_gray_score"]))
        best_metastable = max(sweep_rows, key=lambda row: (int(row["gray_metastable_count"]), float(row["best_gray_depth"])))
        best_eigen = max(sweep_rows, key=lambda row: (int(row["gray_eigen_count"]), float(row["best_gray_depth"])))
        best_selection = max(sweep_rows, key=lambda row: int(row["natural_selection_count"]))
        max_score = float(best_joint["joint_gray_score"])
        within_5 = [row for row in sweep_rows if float(row["joint_gray_score"]) >= max_score * 0.95 - 1.0e-15]
        within_10 = [row for row in sweep_rows if float(row["joint_gray_score"]) >= max_score * 0.90 - 1.0e-15]

        def width(selected: List[Dict[str, Any]]) -> float:
            if not selected:
                return float("nan")
            r_list = [float(row["R"]) for row in selected]
            return max(r_list) - min(r_list)

        coarse_rows = [row for row in all_case_rows if str(row["sweep_region"]) == "coarse"]
        fine_rows = [row for row in all_case_rows if str(row["sweep_region"]) == "fine"]
        control_rows = [row for row in all_case_rows if str(row["sweep_region"]) == "control"]
        best_coarse = max(coarse_rows, key=lambda row: float(row["joint_gray_score"])) if coarse_rows else None
        best_fine = max(fine_rows, key=lambda row: float(row["joint_gray_score"])) if fine_rows else None
        best_control = max(control_rows, key=lambda row: float(row["joint_gray_score"])) if control_rows else None
        out.append(
            {
                "case_id": case_id,
                "mode": str(best_joint["mode"]),
                "state_family": str(best_joint["state_family"]),
                "harmonic_packet_A": str(best_joint["harmonic_packet_A"]),
                "harmonic_packet_B": str(best_joint["harmonic_packet_B"]),
                "N_A": int(best_joint["N_A"]),
                "N_B": int(best_joint["N_B"]),
                "pair_kind": str(best_joint["pair_kind"]),
                "R_star_gray_metastable": float(best_metastable["R"]),
                "gray_metastable_count_at_star": int(best_metastable["gray_metastable_count"]),
                "R_star_gray_eigen": float(best_eigen["R"]),
                "gray_eigen_count_at_star": int(best_eigen["gray_eigen_count"]),
                "R_star_selection": float(best_selection["R"]),
                "natural_selection_count_at_star": int(best_selection["natural_selection_count"]),
                "R_star_joint": float(best_joint["R"]),
                "joint_gray_score_max": float(best_joint["joint_gray_score"]),
                "best_phase_at_joint": str(best_joint["best_phase"]),
                "best_condition_at_joint": str(best_joint["best_condition_id"]),
                "best_gray_depth_at_joint": float(best_joint["best_gray_depth"]),
                "best_S_mean_at_joint": float(best_joint["best_S_mean"]),
                "best_S_amp_at_joint": float(best_joint["best_S_amp"]),
                "best_S_drift_at_joint": float(best_joint["best_S_drift"]),
                "R_band_width_5": width(within_5),
                "R_band_width_10": width(within_10),
                "distance_to_R_137": abs(float(best_joint["R"]) - params.r137),
                "distance_to_R_128": abs(float(best_joint["R"]) - params.r128),
                "best_control_R": float(best_control["R"]) if best_control is not None else float("nan"),
                "best_coarse_R": float(best_coarse["R"]) if best_coarse is not None else float("nan"),
                "best_fine_R": float(best_fine["R"]) if best_fine is not None else float("nan"),
            }
        )
    return out


def peak_kind_for_r(r_value: float, params: Params, is_local: bool) -> str:
    kinds: List[str] = []
    if is_local:
        kinds.append("local_max")
    if abs(r_value - params.r137) <= 5.0e-12:
        kinds.append("probe_R137")
    if abs(r_value - params.r128) <= 5.0e-12:
        kinds.append("probe_alpha128_nominal")
    return "+".join(kinds) if kinds else "selected"


def peak_rows_for_case(rows: List[Dict[str, Any]], params: Params) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for case_id in sorted({str(row["case_id"]) for row in rows}):
        case_rows = sorted(
            [row for row in rows if str(row["case_id"]) == case_id and str(row["sweep_region"]) in {"coarse", "fine"}],
            key=lambda row: float(row["R"]),
        )
        if not case_rows:
            continue
        candidates: Dict[float, Tuple[Dict[str, Any], bool]] = {}
        for i, row in enumerate(case_rows):
            score = float(row["joint_gray_score"])
            left = float(case_rows[i - 1]["joint_gray_score"]) if i > 0 else float("-inf")
            right = float(case_rows[i + 1]["joint_gray_score"]) if i + 1 < len(case_rows) else float("-inf")
            is_local = score >= left and score >= right and (score > left or score > right)
            is_probe = abs(float(row["R"]) - params.r137) <= 5.0e-12 or abs(float(row["R"]) - params.r128) <= 5.0e-12
            if is_local or is_probe:
                candidates[float(row["R"])] = (row, is_local)
        ranked = sorted(candidates.values(), key=lambda item: float(item[0]["joint_gray_score"]), reverse=True)
        selected: List[Tuple[Dict[str, Any], bool]] = ranked[: params.peak_top_k]
        selected_r = {float(row["R"]) for row, _is_local in selected}
        for row, is_local in ranked:
            r_value = float(row["R"])
            is_probe = abs(r_value - params.r137) <= 5.0e-12 or abs(r_value - params.r128) <= 5.0e-12
            if is_probe and r_value not in selected_r:
                selected.append((row, is_local))
                selected_r.add(r_value)
        for rank, (row, is_local) in enumerate(selected, start=1):
            r_value = float(row["R"])
            out.append(
                {
                    "case_id": case_id,
                    "mode": str(row["mode"]),
                    "state_family": str(row["state_family"]),
                    "harmonic_packet_A": str(row["harmonic_packet_A"]),
                    "harmonic_packet_B": str(row["harmonic_packet_B"]),
                    "N_A": int(row["N_A"]),
                    "N_B": int(row["N_B"]),
                    "pair_kind": str(row["pair_kind"]),
                    "peak_rank": rank,
                    "peak_kind": peak_kind_for_r(r_value, params, is_local),
                    "R_peak": r_value,
                    "joint_gray_score": float(row["joint_gray_score"]),
                    "best_phase": str(row["best_phase"]),
                    "best_condition_id": str(row["best_condition_id"]),
                    "best_gray_depth": float(row["best_gray_depth"]),
                    "best_gray_error": float(row["best_gray_error"]),
                    "best_S_mean": float(row["best_S_mean"]),
                    "best_S_amp": float(row["best_S_amp"]),
                    "best_S_drift": float(row["best_S_drift"]),
                    "gray_eigen_count": int(row["gray_eigen_count"]),
                    "gray_metastable_count": int(row["gray_metastable_count"]),
                    "natural_selection_count": int(row["natural_selection_count"]),
                    "distance_to_R_137": abs(r_value - params.r137),
                    "distance_to_R_128": abs(r_value - params.r128),
                }
            )
    return out


def primary_case_ids(rows: List[Dict[str, Any]]) -> List[str]:
    return sorted({str(row["case_id"]) for row in rows if str(row["sweep_region"]) in {"coarse", "fine"}})


def case_label_from_rows(rows: List[Dict[str, Any]], case_id: str) -> str:
    row = next(item for item in rows if str(item["case_id"]) == case_id)
    return f"{row['mode']} A=[{row['harmonic_packet_A']}] B=[{row['harmonic_packet_B']}]"


def best_depth_rows_by_r(rows: List[Dict[str, Any]], case_id: str) -> List[Dict[str, Any]]:
    by_r: Dict[float, List[Dict[str, Any]]] = {}
    for row in rows:
        if str(row["case_id"]) != case_id:
            continue
        if str(row["sweep_region"]) not in {"coarse", "fine"}:
            continue
        by_r.setdefault(float(row["R"]), []).append(row)
    return [min(r_rows, key=lambda row: float(row["gray_error"])) for _, r_rows in sorted(by_r.items())]


def summary_row_by_r(summaries: List[Dict[str, Any]], case_id: str, r_value: float, tol: float = 5.0e-12) -> Optional[Dict[str, Any]]:
    for row in summaries:
        if str(row["case_id"]) == case_id and abs(float(row["R"]) - r_value) <= tol:
            return row
    return None


def nearest_peak_to_r(peaks: List[Dict[str, Any]], case_id: str, r_value: float) -> Optional[Dict[str, Any]]:
    case_peaks = [row for row in peaks if str(row["case_id"]) == case_id and "local_max" in str(row["peak_kind"])]
    if not case_peaks:
        return None
    return min(case_peaks, key=lambda row: abs(float(row["R_peak"]) - r_value))


def make_score_plot(best: List[Dict[str, Any]], summaries: List[Dict[str, Any]], peaks: List[Dict[str, Any]], params: Params, file_stem: str) -> str:
    cols = 2
    rows_n = math.ceil(len(best) / cols)
    fig, axes = plt.subplots(rows_n, cols, figsize=(12, 3.2 * rows_n), sharex=True)
    axes_flat = np.atleast_1d(axes).reshape(-1)
    for ax, best_row in zip(axes_flat, best):
        case_id = str(best_row["case_id"])
        case_rows = sorted([row for row in summaries if str(row["case_id"]) == case_id], key=lambda row: float(row["R"]))
        ax.plot([row["R"] for row in case_rows], [row["joint_gray_score"] for row in case_rows], marker="o", linewidth=1.2)
        case_peaks = [row for row in peaks if str(row["case_id"]) == case_id and int(row["peak_rank"]) <= 5]
        if case_peaks:
            ax.scatter(
                [float(row["R_peak"]) for row in case_peaks],
                [float(row["joint_gray_score"]) for row in case_peaks],
                color="tab:orange",
                s=42,
                zorder=4,
                label="top local peaks",
            )
        ax.axvline(best_row["R_star_joint"], color="black", linestyle="--", linewidth=1)
        ax.axvline(params.r137, color="tab:red", linestyle=":", linewidth=1)
        ax.axvline(params.r128, color="tab:purple", linestyle=":", linewidth=1)
        ax.set_title(f"{best_row['mode']} B=[{best_row['harmonic_packet_B']}] R*={best_row['R_star_joint']:.6f}")
        ax.set_ylabel("joint gray score")
        ax.grid(alpha=0.25)
        if case_peaks:
            ax.legend(loc="best", fontsize=8)
    for ax in axes_flat[len(best):]:
        ax.axis("off")
    for ax in axes_flat[-cols:]:
        ax.set_xlabel("R")
    fig.tight_layout()
    filename = f"{file_stem}_scores_v1.png"
    fig.savefig(OUT_DIR / filename, dpi=160)
    plt.close(fig)
    return filename


def make_gray_depth_distribution_plot(rows: List[Dict[str, Any]], summaries: List[Dict[str, Any]], peaks: List[Dict[str, Any]], params: Params, file_stem: str, deep: bool) -> str:
    case_ids = primary_case_ids(rows)
    fig, axes = plt.subplots(len(case_ids), 1, figsize=(12, 2.6 * len(case_ids)), sharex=not deep, sharey=True, constrained_layout=True)
    axes_list = list(np.atleast_1d(axes))
    for ax, case_id in zip(axes_list, case_ids):
        case_rows = [row for row in rows if str(row["case_id"]) == case_id and str(row["sweep_region"]) in {"coarse", "fine"}]
        best_row = max(case_rows, key=lambda row: float(row["gray_depth"]))
        if deep:
            r_center = float(best_row["R"])
            plot_rows = [row for row in case_rows if abs(float(row["R"]) - r_center) <= 0.020]
        else:
            plot_rows = case_rows
        ax.scatter(
            [float(row["R"]) for row in plot_rows],
            [float(row["gray_depth_plot"]) for row in plot_rows],
            s=12 if not deep else 16,
            color="0.70" if not deep else "0.68",
            alpha=0.22 if not deep else 0.25,
            linewidths=0,
            label="conditions",
        )
        envelope = best_depth_rows_by_r(rows, case_id)
        if deep:
            envelope = [row for row in envelope if abs(float(row["R"]) - float(best_row["R"])) <= 0.020]
        ax.plot([float(row["R"]) for row in envelope], [float(row["gray_depth_plot"]) for row in envelope], color="black", linewidth=1.5)
        ax.scatter([float(best_row["R"])], [float(best_row["gray_depth_plot"])], marker="x", s=90, color="black", linewidths=1.8)
        ax.axvline(params.r137, color="black", linestyle=":", linewidth=1.1)
        ax.axvline(params.r128, color="0.35", linestyle="--", linewidth=1.1)
        if deep:
            case_peaks = [row for row in peaks if str(row["case_id"]) == case_id and int(row["peak_rank"]) <= 3]
            r128_nearest = nearest_peak_to_r(peaks, case_id, params.r128)
            if r128_nearest is not None and all(abs(float(row["R_peak"]) - float(r128_nearest["R_peak"])) > 5.0e-12 for row in case_peaks):
                case_peaks.append(r128_nearest)
            for i, peak in enumerate(case_peaks):
                r_peak = float(peak["R_peak"])
                if not plot_rows or not (min(float(row["R"]) for row in plot_rows) <= r_peak <= max(float(row["R"]) for row in plot_rows)):
                    continue
                y_peak = min(float(peak["best_gray_depth"]), 6.0)
                is_r128_near = abs(r_peak - params.r128) == min(abs(float(row["R_peak"]) - params.r128) for row in peaks if str(row["case_id"]) == case_id and "local_max" in str(row["peak_kind"]))
                marker_color = "tab:purple" if is_r128_near else "tab:orange"
                ax.scatter([r_peak], [y_peak], marker="D", s=42, color=marker_color, zorder=5)
                label_prefix = "near alpha~128" if is_r128_near and int(peak["peak_rank"]) > 3 else f"#{int(peak['peak_rank'])}"
                ax.annotate(
                    f"{label_prefix} R={r_peak:.12g}\nscore={float(peak['joint_gray_score']):.6g}",
                    xy=(r_peak, y_peak),
                    xytext=(8, 10 + 16 * (i % 3)),
                    textcoords="offset points",
                    fontsize=7,
                    color=marker_color,
                    arrowprops={"arrowstyle": "-", "color": marker_color, "lw": 0.8},
                )
            r128_row = summary_row_by_r(summaries, case_id, params.r128)
            if r128_row is not None:
                y_r128 = min(float(r128_row["best_gray_depth"]), 6.0)
                ax.scatter([params.r128], [y_r128], marker="s", s=34, facecolors="none", edgecolors="tab:purple", zorder=5)
                ax.annotate(
                    f"alpha~128 nominal\nR={params.r128:.12g}\nscore={float(r128_row['joint_gray_score']):.6g}",
                    xy=(params.r128, y_r128),
                    xytext=(-86, -34),
                    textcoords="offset points",
                    fontsize=7,
                    color="tab:purple",
                    arrowprops={"arrowstyle": "-", "color": "tab:purple", "lw": 0.8},
                )
        title = f"{case_label_from_rows(rows, case_id)} gray-depth distribution"
        if deep:
            title = f"{case_label_from_rows(rows, case_id)} deepest gray R={float(best_row['R']):.12g}"
        ax.set_title(title)
        ax.set_ylabel("depth = -log10(gray error), capped at 6")
        ax.grid(alpha=0.20)
    axes_list[-1].set_xlabel("R")
    axes_list[0].legend(loc="upper right", fontsize=8)
    suffix = "deep" if deep else "overview"
    filename = f"{file_stem}_gray_depth_distribution_{suffix}_v1.png"
    fig.savefig(OUT_DIR / filename, dpi=170)
    plt.close(fig)
    return filename


def make_peak_zoom_plot(best: List[Dict[str, Any]], summaries: List[Dict[str, Any]], peaks: List[Dict[str, Any]], params: Params, file_stem: str) -> str:
    cols = 2
    rows_n = math.ceil(len(best) / cols)
    fig, axes = plt.subplots(rows_n, cols, figsize=(13, 3.6 * rows_n), sharex=True)
    axes_flat = np.atleast_1d(axes).reshape(-1)
    x_min = min(0.680, params.r128 - 0.003)
    x_max = max(0.703, params.r137 + 0.004)
    for ax, best_row in zip(axes_flat, best):
        case_id = str(best_row["case_id"])
        case_rows = sorted(
            [
                row
                for row in summaries
                if str(row["case_id"]) == case_id and x_min <= float(row["R"]) <= x_max
            ],
            key=lambda row: float(row["R"]),
        )
        ax.plot(
            [float(row["R"]) for row in case_rows],
            [float(row["joint_gray_score"]) for row in case_rows],
            color="black",
            linewidth=1.4,
            marker="o",
            markersize=3,
        )
        ax.axvline(params.r137, color="tab:red", linestyle=":", linewidth=1.1)
        ax.axvline(params.r128, color="tab:purple", linestyle="--", linewidth=1.1)
        case_peaks = [row for row in peaks if str(row["case_id"]) == case_id and int(row["peak_rank"]) <= 3]
        r128_nearest = nearest_peak_to_r(peaks, case_id, params.r128)
        if r128_nearest is not None and all(abs(float(row["R_peak"]) - float(r128_nearest["R_peak"])) > 5.0e-12 for row in case_peaks):
            case_peaks.append(r128_nearest)
        for i, peak in enumerate(case_peaks):
            r_peak = float(peak["R_peak"])
            if not (x_min <= r_peak <= x_max):
                continue
            score = float(peak["joint_gray_score"])
            near_r128 = r128_nearest is not None and abs(r_peak - float(r128_nearest["R_peak"])) <= 5.0e-12
            color = "tab:purple" if near_r128 else "tab:orange"
            marker = "D" if near_r128 else "o"
            ax.scatter([r_peak], [score], s=54, color=color, marker=marker, zorder=5)
            label = "near alpha~128" if near_r128 and int(peak["peak_rank"]) > 3 else f"#{int(peak['peak_rank'])}"
            ax.annotate(
                f"{label}\nR={r_peak:.12g}\nscore={score:.6g}",
                xy=(r_peak, score),
                xytext=(8, 12 + 18 * (i % 2)),
                textcoords="offset points",
                fontsize=7,
                color=color,
                arrowprops={"arrowstyle": "-", "color": color, "lw": 0.8},
            )
        r128_row = summary_row_by_r(summaries, case_id, params.r128)
        if r128_row is not None:
            r128_score = float(r128_row["joint_gray_score"])
            ax.scatter([params.r128], [r128_score], marker="s", s=42, facecolors="white", edgecolors="tab:purple", zorder=6)
            ax.annotate(
                f"alpha~128 nominal\nR={params.r128:.12g}\nscore={r128_score:.6g}",
                xy=(params.r128, r128_score),
                xytext=(-88, -48),
                textcoords="offset points",
                fontsize=7,
                color="tab:purple",
                arrowprops={"arrowstyle": "-", "color": "tab:purple", "lw": 0.8},
            )
        y_values = [float(row["joint_gray_score"]) for row in case_rows]
        if y_values:
            y_min = min(y_values)
            y_max = max(y_values)
            pad = max(0.25, (y_max - y_min) * 0.18)
            ax.set_ylim(y_min - pad * 0.45, y_max + pad * 1.75)
        ax.set_title(f"{best_row['mode']} B=[{best_row['harmonic_packet_B']}] peak zoom", fontsize=10)
        ax.set_ylabel("joint gray score")
        ax.grid(alpha=0.22)
    for ax in axes_flat[len(best):]:
        ax.axis("off")
    for ax in axes_flat[-cols:]:
        ax.set_xlabel("R")
    fig.tight_layout()
    filename = f"{file_stem}_peak_zoom_v1.png"
    fig.savefig(OUT_DIR / filename, dpi=180)
    plt.close(fig)
    return filename


def build_report(params: Params, best: List[Dict[str, Any]], peaks: List[Dict[str, Any]], outputs: Dict[str, str]) -> str:
    lines: List[str] = [
        "# 系統B 灰色猫準安定界面R近傍スイープ 予備実験結果 v1",
        "",
        "## 1. 実験条件",
        "",
        "```text",
        f"coarse R range = {params.coarse_r_min} .. {params.coarse_r_max}",
        f"coarse Delta R = {params.coarse_r_step}",
        f"fine R range = {params.fine_r_min} .. {params.fine_r_max}",
        f"fine Delta R = {params.fine_r_step}",
        f"R_137 = {params.r137}",
        f"R_alpha128_nominal = {params.r128}",
        f"steps = {params.steps}",
        "C = off",
        "D = off",
        "```",
        "",
        "## 2. 系統B 判定サマリー",
        "",
        "| mode | A packet | B packet | N_A | N_B | kind | R*_metastable | metastable count | R*_eigen | eigen count | R*_selection | selection count | R*_joint | phase@joint | depth | S_amp | band5 | band10 | d137 | d128 | coarse R | fine R | control R |",
        "|---|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in best:
        lines.append(
            "| {mode} | {harmonic_packet_A} | {harmonic_packet_B} | {N_A} | {N_B} | {pair_kind} | "
            "{R_star_gray_metastable:.12g} | {gray_metastable_count_at_star} | {R_star_gray_eigen:.12g} | {gray_eigen_count_at_star} | "
            "{R_star_selection:.12g} | {natural_selection_count_at_star} | {R_star_joint:.12g} | {best_phase_at_joint} | "
            "{best_gray_depth_at_joint:.12g} | {best_S_amp_at_joint:.12g} | {R_band_width_5:.12g} | {R_band_width_10:.12g} | "
            "{distance_to_R_137:.12g} | {distance_to_R_128:.12g} | {best_coarse_R:.12g} | {best_fine_R:.12g} | {best_control_R:.12g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## 3. R地形ピーク候補",
            "",
            "`R*_joint` は単一の代表点である。ただし系統BではR地形が多峰的になるため、局所ピークを別表として保存する。",
            "",
            "ここでは、各ケースで `joint_gray_score` の局所最大点を順位付きで記録する。",
            "",
            "`R_137` と `alpha~128 nominal` は、局所最大でない場合でも固定プローブ点として表に残す。",
            "",
            "| case | rank | kind | R_peak | score | phase | depth | S_amp | d137 | d128 |",
            "|---|---:|---|---:|---:|---|---:|---:|---:|---:|",
        ]
    )
    for row in peaks:
        lines.append(
            "| {case_id} | {peak_rank} | {peak_kind} | {R_peak:.12g} | {joint_gray_score:.12g} | {best_phase} | "
            "{best_gray_depth:.12g} | {best_S_amp:.12g} | {distance_to_R_137:.12g} | {distance_to_R_128:.12g} |".format(
                **{**row, "case_id": md_cell(row["case_id"])}
            )
        )
    lines.extend(
        [
            "",
            "## 4. 読み",
            "",
            "本予備実験では、20260714 の `S=p_A-p_B` 準安定力学を母体とし、AB交換を `R` 指定の散乱写像に置き換えた。",
            "",
            "Nと倍音パケットは、系統Aと同じ比較軸として記録した。",
            "",
            "この二成分準安定モデルでは、Nと奇偶は `S` 写像に直接入らないため、N依存性が出る場合は、R写像または条件選別を通じた間接効果として読む。",
            "",
            "`R_137` は低エネルギー側の精密プローブ、`alpha~128 nominal` は高エネルギー側の名目プローブである。",
            "",
            "`alpha~128 nominal` は精密値ではないため、固定点そのものだけでなく、その近傍局所ピークをあわせて読む。",
            "",
            "したがって、最深ピークが別のRにある場合でも、alpha近傍に独立したピークが立つなら、それ自体をR地形上の重要な候補として扱う。",
            "",
            "分布図では、`gray_error=|S_mean|+|S_amp-S_amp_target|+S_drift+penalty` を `gray_depth=-log10(gray_error)` として縦軸に置く。",
            "",
            "灰色点は各 `R` における各初期条件、黒線は各 `R` で最も深い条件の包絡線である。",
            "",
            "## 5. 出力",
            "",
            "| 種類 | ファイル |",
            "|---|---|",
        ]
    )
    for label, filename in outputs.items():
        lines.append(f"| {label} | `{filename}` |")
    lines.append("")
    return "\n".join(lines)


def run(
    selected_cases: Optional[Tuple[Any, ...]] = None,
    selected_r_values: Optional[List[float]] = None,
    output_dir: Optional[Path] = None,
    max_steps: Optional[int] = None,
    make_plots: bool = True,
    run_id: Optional[str] = None,
    file_stem: Optional[str] = None,
) -> Dict[str, Any]:
    params = Params()
    if max_steps is not None:
        params.steps = max_steps
    cases = selected_cases if selected_cases is not None else tuple(system_a.odd_kernel_case(pair) for pair in params.pairs)
    r_list = selected_r_values if selected_r_values is not None else r_values(params)
    if output_dir is not None:
        configure_output_dir(output_dir)
    resolved_file_stem = build_file_stem(cases, r_list, params, run_id, file_stem)

    rows: List[Dict[str, Any]] = []
    samples: List[Dict[str, Any]] = []
    for case in cases:
        for r_value in r_list:
            case_rows, case_samples = run_case(params, case, r_value)
            rows.extend(case_rows)
            samples.extend(case_samples)

    summaries = aggregate_by_case_r(rows)
    best = best_rows_for_case(summaries, params)
    peaks = peak_rows_for_case(summaries, params)

    outputs = {
        "rows": f"{resolved_file_stem}_rows_v1.csv",
        "samples": f"{resolved_file_stem}_samples_v1.csv",
        "summary": f"{resolved_file_stem}_summary_v1.csv",
        "best": f"{resolved_file_stem}_best_v1.csv",
        "peaks": f"{resolved_file_stem}_peaks_v1.csv",
        "result": f"{resolved_file_stem}_result_v1.json",
        "report": f"{resolved_file_stem}_report_v1.md",
    }
    if make_plots:
        outputs["scores"] = make_score_plot(best, summaries, peaks, params, resolved_file_stem)
        outputs["gray_depth_overview"] = make_gray_depth_distribution_plot(rows, summaries, peaks, params, resolved_file_stem, deep=False)
        outputs["gray_depth_deep"] = make_gray_depth_distribution_plot(rows, summaries, peaks, params, resolved_file_stem, deep=True)
        outputs["peak_zoom"] = make_peak_zoom_plot(best, summaries, peaks, params, resolved_file_stem)

    write_csv(OUT_DIR / outputs["rows"], rows)
    write_csv(OUT_DIR / outputs["samples"], samples)
    write_csv(OUT_DIR / outputs["summary"], summaries)
    write_csv(OUT_DIR / outputs["best"], best)
    write_csv(OUT_DIR / outputs["peaks"], peaks)

    result = {
        "experiment": "system_B_gray_cat_metastable_R_sweep_preliminary_v1",
        "params": asdict(params),
        "cases": [
            {
                "case_id": case.case_id,
                "mode": case.mode,
                "state_family": case.state_family,
                "N_A": case.n_a,
                "N_B": case.n_b,
                "harmonic_packet_A": case.packet_a_text,
                "harmonic_packet_B": case.packet_b_text,
            }
            for case in cases
        ],
        "r_values": r_list,
        "best": best,
        "peaks": peaks,
        "outputs": outputs,
    }
    (OUT_DIR / outputs["result"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT_DIR / outputs["report"]).write_text(build_report(params, best, peaks, outputs), encoding="utf-8")
    return result


if __name__ == "__main__":
    args = parse_args()
    params = Params()
    if args.max_steps is not None:
        params.steps = args.max_steps
    cases = selected_cases_from_args(params, args)
    r_list = selected_r_values_from_args(params, args)
    out = run(
        selected_cases=cases,
        selected_r_values=r_list,
        output_dir=output_dir_from_args(args),
        max_steps=args.max_steps,
        make_plots=not args.no_plots,
        run_id=args.run_id,
        file_stem=args.file_stem,
    )
    print(json.dumps({"best": out["best"], "outputs": out["outputs"]}, ensure_ascii=False, indent=2))
