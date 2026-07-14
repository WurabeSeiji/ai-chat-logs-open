from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

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


BASE_DIR = Path(__file__).resolve().parent
AB_RESULT_PATH = (
    BASE_DIR
    / "gray_cat_ab_metastable_interface_preliminary_result_v1"
    / "gray_cat_ab_metastable_interface_preliminary_result_v1.json"
)
OUT_DIR = BASE_DIR / "gray_cat_c_readout_window_preliminary_result_v1"


@dataclass(frozen=True)
class Params:
    steps: int = 4096
    readout_kappa: float = 0.02
    c_read_rel_tol: float = 0.10
    c_read_abs_tol: float = 1.0e-4
    c_bias_delta_tol: float = 1.0e-2
    g_c_values: Tuple[float, ...] = (0.0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00)
    backaction_scale_values: Tuple[float, ...] = (0.0, 1.0e-5, 5.0e-5, 1.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 5.0e-3, 1.0e-2)
    candidate_limit: int = 10


def read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def c_visibility(g_c: float, kappa: float) -> float:
    if g_c <= 0.0:
        return 0.0
    return float(g_c / (g_c + kappa))


def apply_c_backaction(a: complex, b: complex, c_gain: float, s_c: float) -> Tuple[complex, complex]:
    if c_gain == 0.0 or s_c == 0.0:
        return a, b
    p_a = abs(a) ** 2
    p_b = abs(b) ** 2
    q = p_a + p_b
    if q <= 0.0:
        raise ValueError("zero AB norm")
    s = (p_a - p_b) / q
    s_next = s + c_gain * s_c * (1.0 - s * s)
    s_next = max(min(s_next, 1.0 - 1.0e-12), -1.0 + 1.0e-12)
    phase_a = a / abs(a) if abs(a) > 0.0 else 1.0 + 0.0j
    phase_b = b / abs(b) if abs(b) > 0.0 else 1.0 + 0.0j
    a_next = math.sqrt(0.5 * (1.0 + s_next)) * phase_a
    b_next = math.sqrt(0.5 * (1.0 - s_next)) * phase_b
    return normalize_pair(a_next, b_next)


def candidate_rows(ab_result: Dict[str, Any], params: Params) -> List[Dict[str, Any]]:
    candidates = list(ab_result["summary"]["top_metastable_candidates"][: params.candidate_limit])
    rows = ab_result["rows"]
    eigen = [row for row in rows if row["phase"] == "gray_eigen" and row["noise_amp"] == 0.0]
    large = [row for row in rows if row["phase"] == "large_oscillation" and row["noise_amp"] == 0.0]
    if eigen:
        candidates.append(sorted(eigen, key=lambda r: abs(r["S_amp"]))[0] | {"candidate_kind": "gray_eigen_control"})
    if large:
        candidates.append(sorted(large, key=lambda r: abs(r["S_amp"] - 0.08))[0] | {"candidate_kind": "large_oscillation_control"})
    for row in candidates:
        row.setdefault("candidate_kind", row["phase"])
    return candidates


def run_case(
    ab_params: ABParams,
    params: Params,
    candidate: Dict[str, Any],
    g_c: float,
    backaction_scale: float,
) -> Dict[str, Any]:
    visibility = c_visibility(g_c, params.readout_kappa)
    c_gain = g_c * backaction_scale
    a, b = state_from_s_phi(float(candidate["s0"]), float(candidate["phi"]))
    rng = np.random.default_rng(int(candidate["seed"]))
    s_values = np.empty(params.steps + 1)
    q_values = np.empty(params.steps + 1)
    s_c_values = np.empty(params.steps + 1)
    c_error_values = np.empty(params.steps + 1)
    c_a_values = np.empty(params.steps + 1)
    c_b_values = np.empty(params.steps + 1)
    for k in range(params.steps + 1):
        p_a = abs(a) ** 2
        p_b = abs(b) ** 2
        q = p_a + p_b
        s = (p_a - p_b) / q
        s_c = visibility * s
        c_a = 0.5 * (1.0 + s_c)
        c_b = 0.5 * (1.0 - s_c)
        s_values[k] = s
        q_values[k] = q
        s_c_values[k] = s_c
        c_error_values[k] = abs(s_c - s)
        c_a_values[k] = c_a
        c_b_values[k] = c_b
        if k >= params.steps:
            break
        a, b = apply_exchange(a, b, float(candidate["epsilon"]))
        a, b = apply_stability_gain(a, b, float(candidate["stability_gain"]))
        a, b = apply_c_backaction(a, b, c_gain, s_c)
        a, b = apply_noise(a, b, float(candidate["noise_amp"]), rng)

    metrics = classify(ab_params, s_values, q_values)
    tail_start = int(round(params.steps * (1.0 - ab_params.tail_fraction)))
    tail_s = s_values[tail_start:]
    tail_c_error = c_error_values[tail_start:]
    tail_c_a = c_a_values[tail_start:]
    tail_c_b = c_b_values[tail_start:]
    mean_abs_s = float(np.mean(np.abs(tail_s)))
    c_read_abs_error = float(np.mean(tail_c_error))
    c_read_rel_error = float(c_read_abs_error / max(mean_abs_s, 1.0e-12))
    c_bias_delta = float(abs(metrics["S_mean"] - float(candidate["S_mean"])) + abs(metrics["S_amp"] - float(candidate["S_amp"])))
    c_readout_ok = bool(
        (mean_abs_s <= ab_params.eigen_amp_tol and c_read_abs_error <= params.c_read_abs_tol)
        or (mean_abs_s > ab_params.eigen_amp_tol and c_read_rel_error <= params.c_read_rel_tol)
    )
    c_nonselective_ok = bool(metrics["phase"] != "natural_selection" and c_bias_delta <= params.c_bias_delta_tol)
    c_window_ok = bool(c_readout_ok and c_nonselective_ok)
    c_informative_window_ok = bool(c_window_ok and candidate["candidate_kind"] == "gray_metastable" and g_c > 0.0)
    return {
        "candidate_kind": candidate["candidate_kind"],
        "base_phase": candidate["phase"],
        "epsilon": candidate["epsilon"],
        "phi": candidate["phi"],
        "phi_over_pi": candidate["phi_over_pi"],
        "s0": candidate["s0"],
        "stability_gain": candidate["stability_gain"],
        "base_S_mean": candidate["S_mean"],
        "base_S_amp": candidate["S_amp"],
        "g_C": g_c,
        "visibility_C": visibility,
        "backaction_scale": backaction_scale,
        "c_gain": c_gain,
        "phase_after_C": metrics["phase"],
        "S_mean_after_C": metrics["S_mean"],
        "S_amp_after_C": metrics["S_amp"],
        "S_drift_after_C": metrics["S_drift"],
        "S_final_after_C": metrics["S_final"],
        "Q_max_error": metrics["Q_max_error"],
        "C_A_mean_tail": float(np.mean(tail_c_a)),
        "C_B_mean_tail": float(np.mean(tail_c_b)),
        "C_read_abs_error": c_read_abs_error,
        "C_read_rel_error": c_read_rel_error,
        "C_induced_bias_delta": c_bias_delta,
        "C_readout_ok": c_readout_ok,
        "C_nonselective_ok": c_nonselective_ok,
        "C_window_ok": c_window_ok,
        "C_informative_window_ok": c_informative_window_ok,
    }


def summarise(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    window_rows = [row for row in rows if row["C_window_ok"]]
    informative_window_rows = [row for row in rows if row["C_informative_window_ok"]]
    by_candidate: Dict[str, int] = {}
    for row in window_rows:
        key = f'{row["candidate_kind"]}|eps={row["epsilon"]:.6g}|phi={row["phi_over_pi"]:.6g}|s0={row["s0"]:.6g}|gain={row["stability_gain"]:.6g}'
        by_candidate[key] = by_candidate.get(key, 0) + 1
    phase_after_counts: Dict[str, int] = {}
    for row in rows:
        phase_after_counts[row["phase_after_C"]] = phase_after_counts.get(row["phase_after_C"], 0) + 1
    top_windows = sorted(
        informative_window_rows,
        key=lambda r: (r["C_read_rel_error"], r["C_induced_bias_delta"], r["g_C"]),
    )[:15]
    nonzero_backaction_windows = [row for row in informative_window_rows if row["c_gain"] > 0.0]
    top_nonzero_backaction_windows = sorted(
        nonzero_backaction_windows,
        key=lambda r: (r["C_read_rel_error"], r["C_induced_bias_delta"], r["c_gain"]),
    )[:15]
    return {
        "total_cases": len(rows),
        "C_window_count": len(window_rows),
        "C_informative_window_count": len(informative_window_rows),
        "C_nonzero_backaction_window_count": len(nonzero_backaction_windows),
        "phase_after_C_counts": phase_after_counts,
        "C_window_count_by_candidate": by_candidate,
        "top_C_readout_windows": top_windows,
        "top_C_nonzero_backaction_windows": top_nonzero_backaction_windows,
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_report(result: Dict[str, Any]) -> str:
    summary = result["summary"]
    params = result["params"]
    lines = [
        "# 白猫・黒猫・灰色猫 C読出し窓 予備実験結果 v1",
        "",
        "## 1. 実験条件",
        "",
        "```text",
        "C = on",
        "D = off",
        f"steps = {params['steps']}",
        f"g_C_values = {params['g_c_values']}",
        f"backaction_scale_values = {params['backaction_scale_values']}",
        f"readout_kappa = {params['readout_kappa']}",
        "```",
        "",
        "## 2. 判定",
        "",
        "AB二体準安定候補にCを加え、A/B配分を読めるが一方選択を起こさない結合窓を確認した。",
        "",
        f"C_window_count = {summary['C_window_count']} / {summary['total_cases']}",
        f"C_informative_window_count = {summary['C_informative_window_count']} / {summary['total_cases']}",
        f"C_nonzero_backaction_window_count = {summary['C_nonzero_backaction_window_count']} / {summary['total_cases']}",
        "",
        "## 3. C後の相分類",
        "",
        "| phase_after_C | count |",
        "|---|---:|",
    ]
    for phase, count in sorted(summary["phase_after_C_counts"].items()):
        lines.append(f"| {phase} | {count} |")
    lines += [
        "",
        "## 4. 上位C読出し窓",
        "",
        "| kind | epsilon | phi/pi | s0 | base_gain | g_C | c_gain | C_rel_err | C_bias_delta | phase_after_C |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary["top_C_readout_windows"]:
        lines.append(
            "| {candidate_kind} | {epsilon:.6g} | {phi_over_pi:.6g} | {s0:.6g} | {stability_gain:.6g} | {g_C:.6g} | {c_gain:.6g} | {C_read_rel_error:.6g} | {C_induced_bias_delta:.6g} | {phase_after_C} |".format(
                **row
            )
        )
    lines += [
        "",
        "## 5. 非ゼロCバックアクションを持つ読出し窓",
        "",
        "| kind | epsilon | phi/pi | s0 | base_gain | g_C | c_gain | C_rel_err | C_bias_delta | phase_after_C |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary["top_C_nonzero_backaction_windows"]:
        lines.append(
            "| {candidate_kind} | {epsilon:.6g} | {phi_over_pi:.6g} | {s0:.6g} | {stability_gain:.6g} | {g_C:.6g} | {c_gain:.6g} | {C_read_rel_error:.6g} | {C_induced_bias_delta:.6g} | {phase_after_C} |".format(
                **row
            )
        )
    lines += [
        "",
        "## 6. 次段階",
        "",
        "C読出し窓を固定し、Dを加えた観測実験で、灰色猫固有相・灰色猫準安定相・大振幅分離領域の応答を比較する。",
    ]
    return "\n".join(lines) + "\n"


def run() -> Dict[str, Any]:
    params = Params()
    ab_params = ABParams()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ab_result = read_json(AB_RESULT_PATH)
    candidates = candidate_rows(ab_result, params)
    rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        for g_c in params.g_c_values:
            for backaction_scale in params.backaction_scale_values:
                rows.append(run_case(ab_params, params, candidate, g_c, backaction_scale))
    result = {
        "experiment": "gray_cat_c_readout_window_preliminary_v1",
        "params": asdict(params),
        "ab_source": str(AB_RESULT_PATH.relative_to(BASE_DIR)),
        "candidate_count": len(candidates),
        "summary": summarise(rows),
        "rows": rows,
        "outputs": {
            "json": "gray_cat_c_readout_window_preliminary_result_v1.json",
            "csv": "gray_cat_c_readout_window_rows_v1.csv",
            "report": "gray_cat_c_readout_window_report_v1.md",
        },
    }
    write_csv(OUT_DIR / result["outputs"]["csv"], rows)
    (OUT_DIR / result["outputs"]["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT_DIR / result["outputs"]["report"]).write_text(build_report(result), encoding="utf-8")
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["summary"], ensure_ascii=False, indent=2))
