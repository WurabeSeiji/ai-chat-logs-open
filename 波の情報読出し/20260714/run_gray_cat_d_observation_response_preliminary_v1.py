from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np

from run_gray_cat_ab_metastable_interface_preliminary_v1 import (
    Params as ABParams,
    apply_exchange,
    apply_noise,
    apply_stability_gain,
    classify,
    normalize_pair,
    state_from_s_phi,
)
from run_gray_cat_c_readout_window_preliminary_v1 import (
    Params as CParams,
    apply_c_backaction,
    c_visibility,
    candidate_rows,
    read_json,
)


BASE_DIR = Path(__file__).resolve().parent
AB_RESULT_PATH = (
    BASE_DIR
    / "gray_cat_ab_metastable_interface_preliminary_result_v1"
    / "gray_cat_ab_metastable_interface_preliminary_result_v1.json"
)
OUT_DIR = BASE_DIR / "gray_cat_d_observation_response_preliminary_result_v1"


@dataclass(frozen=True)
class Params:
    pre_steps_values: Tuple[int, ...] = (0, 1, 2, 5, 10, 20, 50, 100, 250, 500, 1000, 2000)
    d_steps: int = 2048
    c_readout_kappa: float = 0.02
    d_readout_kappa: float = 0.02
    c_g: float = 1.0
    c_backaction_scale: float = 1.0e-5
    c_modes: Tuple[str, ...] = ("record_only", "weak_C_window")
    g_d_values: Tuple[float, ...] = (0.0, 0.05, 0.10, 0.20, 0.50, 1.00)
    d_backaction_scale_values: Tuple[float, ...] = (0.0, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00)
    d_selection_limit: float = 0.95
    d_selection_amp_limit: float = 0.05
    small_s_limit: float = 0.01
    large_s_limit: float = 0.05
    candidate_limit: int = 5


def measure_s_q(a: complex, b: complex) -> Tuple[float, float]:
    p_a = abs(a) ** 2
    p_b = abs(b) ** 2
    q = p_a + p_b
    if q <= 0.0:
        raise ValueError("zero AB norm")
    return float((p_a - p_b) / q), float(q)


def d_visibility(g_d: float, kappa: float) -> float:
    if g_d <= 0.0:
        return 0.0
    return float(g_d / (g_d + kappa))


def apply_d_backaction(a: complex, b: complex, d_gain: float, s_d: float) -> Tuple[complex, complex]:
    if d_gain == 0.0 or s_d == 0.0:
        return a, b
    s, _ = measure_s_q(a, b)
    s_next = s + d_gain * s_d * (1.0 - s * s)
    s_next = max(min(s_next, 1.0 - 1.0e-12), -1.0 + 1.0e-12)
    phase_a = a / abs(a) if abs(a) > 0.0 else 1.0 + 0.0j
    phase_b = b / abs(b) if abs(b) > 0.0 else 1.0 + 0.0j
    a_next = math.sqrt(0.5 * (1.0 + s_next)) * phase_a
    b_next = math.sqrt(0.5 * (1.0 - s_next)) * phase_b
    return normalize_pair(a_next, b_next)


def selected_candidates(ab_result: Dict[str, Any], params: Params) -> List[Dict[str, Any]]:
    c_params = CParams(candidate_limit=10)
    base_candidates = candidate_rows(ab_result, c_params)
    metastable = [row for row in base_candidates if row["candidate_kind"] == "gray_metastable"]
    eigen = [row for row in base_candidates if row["candidate_kind"] == "gray_eigen_control"]
    large = [row for row in base_candidates if row["candidate_kind"] == "large_oscillation_control"]
    selected = metastable[: params.candidate_limit] + eigen[:1] + large[:1]
    for i, row in enumerate(selected):
        row["case_id"] = (
            f'{row["candidate_kind"]}_{i}_eps{row["epsilon"]:.6g}'
            f'_phi{row["phi_over_pi"]:.6g}_s{row["s0"]:.6g}_g{row["stability_gain"]:.6g}'
        )
    return selected


def pre_evolve(
    candidate: Dict[str, Any],
    params: Params,
    pre_steps: int,
    c_mode: str,
) -> Tuple[complex, complex]:
    a, b = state_from_s_phi(float(candidate["s0"]), float(candidate["phi"]))
    rng = np.random.default_rng(int(candidate["seed"]))
    c_gain = params.c_g * params.c_backaction_scale if c_mode == "weak_C_window" else 0.0
    c_vis = c_visibility(params.c_g, params.c_readout_kappa)
    for _ in range(pre_steps):
        s, _ = measure_s_q(a, b)
        s_c = c_vis * s
        a, b = apply_exchange(a, b, float(candidate["epsilon"]))
        a, b = apply_stability_gain(a, b, float(candidate["stability_gain"]))
        a, b = apply_c_backaction(a, b, c_gain, s_c)
        a, b = apply_noise(a, b, float(candidate["noise_amp"]), rng)
    return a, b


def evolve_from_state(
    ab_params: ABParams,
    params: Params,
    candidate: Dict[str, Any],
    a0: complex,
    b0: complex,
    g_d: float,
    d_backaction_scale: float,
) -> Dict[str, Any]:
    visibility = d_visibility(g_d, params.d_readout_kappa)
    d_gain = g_d * d_backaction_scale
    a, b = a0, b0
    s_values = np.empty(params.d_steps + 1)
    q_values = np.empty(params.d_steps + 1)
    d_a_values = np.empty(params.d_steps + 1)
    d_b_values = np.empty(params.d_steps + 1)
    for k in range(params.d_steps + 1):
        s, q = measure_s_q(a, b)
        s_d = visibility * s
        s_values[k] = s
        q_values[k] = q
        d_a_values[k] = 0.5 * (1.0 + s_d)
        d_b_values[k] = 0.5 * (1.0 - s_d)
        if k >= params.d_steps:
            break
        a, b = apply_exchange(a, b, float(candidate["epsilon"]))
        a, b = apply_stability_gain(a, b, float(candidate["stability_gain"]))
        a, b = apply_d_backaction(a, b, d_gain, s_d)
    metrics = classify(ab_params, s_values, q_values)
    tail_start = int(round(params.d_steps * (1.0 - ab_params.tail_fraction)))
    tail = s_values[tail_start:]
    tail_d_a = d_a_values[tail_start:]
    tail_d_b = d_b_values[tail_start:]
    tail_mean = float(np.mean(tail))
    tail_amp = float(0.5 * (np.max(tail) - np.min(tail)))
    if tail_mean >= params.d_selection_limit and tail_amp <= params.d_selection_amp_limit:
        d_outcome = "white_selected"
    elif tail_mean <= -params.d_selection_limit and tail_amp <= params.d_selection_amp_limit:
        d_outcome = "black_selected"
    elif metrics["phase"] == "gray_eigen":
        d_outcome = "gray_kept_eigen"
    elif metrics["phase"] == "gray_metastable":
        d_outcome = "gray_kept_metastable"
    elif metrics["phase"] == "large_oscillation":
        d_outcome = "large_oscillation"
    else:
        d_outcome = "unresolved"
    return {
        "phase_after": metrics["phase"],
        "D_outcome": d_outcome,
        "S_mean_after": metrics["S_mean"],
        "S_amp_after": metrics["S_amp"],
        "S_drift_after": metrics["S_drift"],
        "S_final_after": metrics["S_final"],
        "Q_max_error": metrics["Q_max_error"],
        "D_A_mean_tail": float(np.mean(tail_d_a)),
        "D_B_mean_tail": float(np.mean(tail_d_b)),
        "D_tail_mean": tail_mean,
        "D_tail_amp": tail_amp,
    }


def sign_label(x: float, limit: float) -> str:
    if x > limit:
        return "A"
    if x < -limit:
        return "B"
    return "gray"


def run_case(
    ab_params: ABParams,
    params: Params,
    candidate: Dict[str, Any],
    pre_steps: int,
    c_mode: str,
    g_d: float,
    d_backaction_scale: float,
) -> Dict[str, Any]:
    a0, b0 = pre_evolve(candidate, params, pre_steps, c_mode)
    s_start, q_start = measure_s_q(a0, b0)
    c_vis = c_visibility(params.c_g, params.c_readout_kappa)
    c_s = c_vis * s_start
    c_a = 0.5 * (1.0 + c_s)
    c_b = 0.5 * (1.0 - c_s)
    baseline = evolve_from_state(ab_params, params, candidate, a0, b0, 0.0, 0.0)
    observed = evolve_from_state(ab_params, params, candidate, a0, b0, g_d, d_backaction_scale)
    baseline_selected = baseline["D_outcome"] in ("white_selected", "black_selected")
    observed_selected = observed["D_outcome"] in ("white_selected", "black_selected")
    induced_selection = bool(observed_selected and not baseline_selected)
    c_sign = sign_label(c_s, params.large_s_limit)
    d_sign = sign_label(observed["S_mean_after"], params.large_s_limit)
    if c_sign == "gray" and d_sign == "gray":
        agreement = "gray_kept"
    elif c_sign == "gray" and d_sign != "gray":
        agreement = "D_selected_from_small_S"
    elif c_sign == d_sign:
        agreement = "same_sign"
    else:
        agreement = "opposite_sign"
    return {
        "case_id": candidate["case_id"],
        "candidate_kind": candidate["candidate_kind"],
        "base_phase": candidate["phase"],
        "epsilon": candidate["epsilon"],
        "phi": candidate["phi"],
        "phi_over_pi": candidate["phi_over_pi"],
        "s0": candidate["s0"],
        "stability_gain": candidate["stability_gain"],
        "pre_steps": pre_steps,
        "C_mode": c_mode,
        "C_g": params.c_g,
        "C_gain": params.c_g * params.c_backaction_scale if c_mode == "weak_C_window" else 0.0,
        "S_start": s_start,
        "Q_start": q_start,
        "C_A_start": c_a,
        "C_B_start": c_b,
        "C_sign_start": c_sign,
        "g_D": g_d,
        "D_visibility": d_visibility(g_d, params.d_readout_kappa),
        "D_backaction_scale": d_backaction_scale,
        "D_gain": g_d * d_backaction_scale,
        "baseline_outcome": baseline["D_outcome"],
        "baseline_phase_after": baseline["phase_after"],
        "baseline_S_mean_after": baseline["S_mean_after"],
        "baseline_S_amp_after": baseline["S_amp_after"],
        "D_outcome": observed["D_outcome"],
        "phase_after_D": observed["phase_after"],
        "S_mean_after_D": observed["S_mean_after"],
        "S_amp_after_D": observed["S_amp_after"],
        "S_drift_after_D": observed["S_drift_after"],
        "S_final_after_D": observed["S_final_after"],
        "D_A_mean_tail": observed["D_A_mean_tail"],
        "D_B_mean_tail": observed["D_B_mean_tail"],
        "Q_max_error_after_D": observed["Q_max_error"],
        "D_induced_selection": induced_selection,
        "D_vs_C_agreement": agreement,
    }


def count_by(rows: Iterable[Dict[str, Any]], key: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        counts[value] = counts.get(value, 0) + 1
    return counts


def summarise(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    induced_rows = [row for row in rows if row["D_induced_selection"]]
    high_gain_rows = [row for row in rows if row["g_D"] == 1.0 and row["D_backaction_scale"] == 1.0]
    gray_eigen_rows = [row for row in high_gain_rows if row["candidate_kind"] == "gray_eigen_control"]
    metastable_rows = [row for row in induced_rows if row["candidate_kind"] == "gray_metastable"]
    large_rows = [row for row in high_gain_rows if row["candidate_kind"] == "large_oscillation_control"]
    return {
        "total_cases": len(rows),
        "D_outcome_counts": count_by(rows, "D_outcome"),
        "baseline_outcome_counts": count_by(rows, "baseline_outcome"),
        "D_induced_selection_count": len(induced_rows),
        "D_vs_C_agreement_counts": count_by(rows, "D_vs_C_agreement"),
        "phase_after_D_counts": count_by(rows, "phase_after_D"),
        "high_gain_D_outcome_counts": count_by(high_gain_rows, "D_outcome"),
        "top_induced_selection_cases": sorted(
            induced_rows,
            key=lambda r: (r["candidate_kind"] != "gray_metastable", abs(r["S_start"]), -abs(r["S_mean_after_D"])),
        )[:20],
        "gray_eigen_high_gain_cases": gray_eigen_rows[:12],
        "metastable_induced_selection_cases": metastable_rows[:20],
        "large_separation_high_gain_cases": large_rows[:12],
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_report(result: Dict[str, Any]) -> str:
    params = result["params"]
    summary = result["summary"]
    lines = [
        "# 白猫・黒猫・灰色猫 D観測応答 予備実験結果 v1",
        "",
        "## 1. 実験条件",
        "",
        "```text",
        "Stage 1: AB準安定界面探索済み",
        "Stage 2: C読出し窓確認済み",
        "Stage 3: D強観測応答",
        f"pre_steps_values = {params['pre_steps_values']}",
        f"d_steps = {params['d_steps']}",
        f"c_modes = {params['c_modes']}",
        f"g_D_values = {params['g_d_values']}",
        f"d_backaction_scale_values = {params['d_backaction_scale_values']}",
        "```",
        "",
        "Dあり条件は、同じD開始状態からDなし対照を並走させて判定した。",
        "Dなし対照でも同じ白猫または黒猫選択が起こる場合は、D起因選択とは数えない。",
        "",
        "## 2. 全体集計",
        "",
        f"total_cases = {summary['total_cases']}",
        f"D_induced_selection_count = {summary['D_induced_selection_count']}",
        "",
        "### D結果分類",
        "",
        "| D_outcome | count |",
        "|---|---:|",
    ]
    for key, count in sorted(summary["D_outcome_counts"].items()):
        lines.append(f"| {key} | {count} |")
    lines += [
        "",
        "### C読出し符号との対応",
        "",
        "| D_vs_C_agreement | count |",
        "|---|---:|",
    ]
    for key, count in sorted(summary["D_vs_C_agreement_counts"].items()):
        lines.append(f"| {key} | {count} |")
    lines += [
        "",
        "## 3. D起因選択の代表例",
        "",
        "| kind | pre | C_mode | eps | phi/pi | s0 | gain | S_start | g_D | D_gain | outcome | baseline | S_mean_after_D | agreement |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---|",
    ]
    for row in summary["top_induced_selection_cases"]:
        lines.append(
            "| {candidate_kind} | {pre_steps} | {C_mode} | {epsilon:.6g} | {phi_over_pi:.6g} | {s0:.6g} | {stability_gain:.6g} | {S_start:.6g} | {g_D:.6g} | {D_gain:.6g} | {D_outcome} | {baseline_outcome} | {S_mean_after_D:.6g} | {D_vs_C_agreement} |".format(
                **row
            )
        )
    lines += [
        "",
        "## 4. 灰色猫固有相の強D応答",
        "",
        "| pre | C_mode | S_start | outcome | S_mean_after_D | S_amp_after_D | Q_err |",
        "|---:|---|---:|---|---:|---:|---:|",
    ]
    for row in summary["gray_eigen_high_gain_cases"]:
        lines.append(
            "| {pre_steps} | {C_mode} | {S_start:.6g} | {D_outcome} | {S_mean_after_D:.6g} | {S_amp_after_D:.6g} | {Q_max_error_after_D:.3g} |".format(
                **row
            )
        )
    lines += [
        "",
        "## 5. 大振幅分離領域の強D応答",
        "",
        "| pre | C_mode | S_start | C_sign | outcome | S_mean_after_D | agreement |",
        "|---:|---|---:|---|---|---:|---|",
    ]
    for row in summary["large_separation_high_gain_cases"]:
        lines.append(
            "| {pre_steps} | {C_mode} | {S_start:.6g} | {C_sign_start} | {D_outcome} | {S_mean_after_D:.6g} | {D_vs_C_agreement} |".format(
                **row
            )
        )
    lines += [
        "",
        "## 6. 判定",
        "",
        "灰色猫固有相、灰色猫準安定相、大振幅分離領域を同じD写像で比較できるデータが得られた。",
        "次段階では、D起因選択が現れた準安定候補を中心に、D結合強度と観測開始位相の境界を細かく調べる。",
    ]
    return "\n".join(lines) + "\n"


def run() -> Dict[str, Any]:
    params = Params()
    ab_params = ABParams()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ab_result = read_json(AB_RESULT_PATH)
    candidates = selected_candidates(ab_result, params)
    rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        for pre_steps in params.pre_steps_values:
            for c_mode in params.c_modes:
                for g_d in params.g_d_values:
                    for d_backaction_scale in params.d_backaction_scale_values:
                        rows.append(run_case(ab_params, params, candidate, pre_steps, c_mode, g_d, d_backaction_scale))
    result = {
        "experiment": "gray_cat_d_observation_response_preliminary_v1",
        "params": asdict(params),
        "ab_source": str(AB_RESULT_PATH.relative_to(BASE_DIR)),
        "candidate_count": len(candidates),
        "candidate_ids": [candidate["case_id"] for candidate in candidates],
        "summary": summarise(rows),
        "rows": rows,
        "outputs": {
            "json": "gray_cat_d_observation_response_preliminary_result_v1.json",
            "csv": "gray_cat_d_observation_response_rows_v1.csv",
            "report": "gray_cat_d_observation_response_report_v1.md",
        },
    }
    write_csv(OUT_DIR / result["outputs"]["csv"], rows)
    (OUT_DIR / result["outputs"]["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT_DIR / result["outputs"]["report"]).write_text(build_report(result), encoding="utf-8")
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["summary"], ensure_ascii=False, indent=2))
