from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "gray_cat_ab_metastable_interface_preliminary_result_v1"


@dataclass(frozen=True)
class Params:
    steps: int = 4096
    tail_fraction: float = 0.5
    s_gray_limit: float = 0.05
    eigen_amp_tol: float = 1.0e-4
    mean_tol: float = 5.0e-3
    drift_tol: float = 1.0e-2
    selection_limit: float = 0.95
    q_tol: float = 1.0e-10
    epsilon_values: Tuple[float, ...] = (0.0, 1.0e-5, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2)
    phi_values: Tuple[float, ...] = (
        0.0,
        math.pi / 12.0,
        math.pi / 6.0,
        math.pi / 4.0,
        math.pi / 2.0,
        3.0 * math.pi / 4.0,
        math.pi,
    )
    s0_values: Tuple[float, ...] = (0.0, 1.0e-4, 1.0e-3, 1.0e-2, 3.0e-2)
    stability_gain_values: Tuple[float, ...] = (-1.0e-2, -2.0e-3, 0.0, 2.0e-3, 1.0e-2)
    noise_values: Tuple[float, ...] = (0.0, 1.0e-8)
    seed_values: Tuple[int, ...] = (0, 1, 2)


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


def apply_exchange(a: complex, b: complex, epsilon: float) -> Tuple[complex, complex]:
    c = math.cos(epsilon)
    s = math.sin(epsilon)
    a_next = c * a + 1j * s * b
    b_next = 1j * s * a + c * b
    return normalize_pair(a_next, b_next)


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


def classify(params: Params, s_values: np.ndarray, q_values: np.ndarray) -> Dict[str, Any]:
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

    return {
        "phase": phase,
        "S_mean": s_mean,
        "S_amp": s_amp,
        "S_max_abs": s_max_abs,
        "S_final": s_final,
        "S_drift": s_drift,
        "Q_max_error": q_max_error,
        "Q_tail_max_error": q_tail_max_error,
    }


def run_case(
    params: Params,
    epsilon: float,
    phi: float,
    s0: float,
    stability_gain: float,
    noise_amp: float,
    seed: int,
) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    a, b = state_from_s_phi(s0, phi)
    s_values = np.empty(params.steps + 1)
    q_values = np.empty(params.steps + 1)
    for k in range(params.steps + 1):
        p_a = abs(a) ** 2
        p_b = abs(b) ** 2
        q = p_a + p_b
        s_values[k] = (p_a - p_b) / q
        q_values[k] = q
        if k >= params.steps:
            break
        a, b = apply_exchange(a, b, epsilon)
        a, b = apply_stability_gain(a, b, stability_gain)
        a, b = apply_noise(a, b, noise_amp, rng)
    metrics = classify(params, s_values, q_values)
    return {
        "epsilon": epsilon,
        "phi": phi,
        "phi_over_pi": phi / math.pi,
        "s0": s0,
        "stability_gain": stability_gain,
        "noise_amp": noise_amp,
        "seed": seed,
        **metrics,
    }


def representative_series(params: Params) -> Dict[str, List[Dict[str, float]]]:
    configs = {
        "gray_eigen": {"epsilon": 1.0e-3, "phi": 0.0, "s0": 0.0, "stability_gain": 0.0, "noise_amp": 0.0, "seed": 0},
        "gray_metastable": {"epsilon": 1.0e-5, "phi": math.pi / 2.0, "s0": 0.0, "stability_gain": 0.0, "noise_amp": 0.0, "seed": 0},
        "gray_restoring": {"epsilon": 1.0e-3, "phi": math.pi / 2.0, "s0": 0.0, "stability_gain": -1.0e-2, "noise_amp": 0.0, "seed": 0},
        "natural_selection": {"epsilon": 1.0e-4, "phi": math.pi / 12.0, "s0": 1.0e-3, "stability_gain": 1.0e-2, "noise_amp": 0.0, "seed": 0},
    }
    out: Dict[str, List[Dict[str, float]]] = {}
    sample_steps = sorted(set([0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, params.steps]))
    for name, cfg in configs.items():
        rng = np.random.default_rng(int(cfg["seed"]))
        a, b = state_from_s_phi(float(cfg["s0"]), float(cfg["phi"]))
        rows: List[Dict[str, float]] = []
        for k in range(params.steps + 1):
            if k in sample_steps:
                p_a = abs(a) ** 2
                p_b = abs(b) ** 2
                q = p_a + p_b
                rows.append({"step": k, "p_A": p_a / q, "p_B": p_b / q, "S": (p_a - p_b) / q, "Q": q})
            if k >= params.steps:
                break
            a, b = apply_exchange(a, b, float(cfg["epsilon"]))
            a, b = apply_stability_gain(a, b, float(cfg["stability_gain"]))
            a, b = apply_noise(a, b, float(cfg["noise_amp"]), rng)
        out[name] = rows
    return out


def summarise(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    phase_counts: Dict[str, int] = {}
    for row in rows:
        phase_counts[row["phase"]] = phase_counts.get(row["phase"], 0) + 1
    by_gain: Dict[str, Dict[str, int]] = {}
    for row in rows:
        key = f'{row["stability_gain"]:.6g}'
        by_gain.setdefault(key, {})
        by_gain[key][row["phase"]] = by_gain[key].get(row["phase"], 0) + 1
    candidates = [
        row
        for row in rows
        if row["phase"] == "gray_metastable"
        and row["noise_amp"] == 0.0
        and abs(row["S_mean"]) < 0.01
        and 0.005 <= row["S_amp"] <= 0.03
    ]
    candidates = sorted(candidates, key=lambda r: (abs(r["S_amp"] - 0.02), abs(r["S_drift"])))
    eigen = [row for row in rows if row["phase"] == "gray_eigen"]
    selection = [row for row in rows if row["phase"] == "natural_selection"]
    return {
        "total_cases": len(rows),
        "phase_counts": phase_counts,
        "phase_counts_by_stability_gain": by_gain,
        "metastable_candidate_count": len(candidates),
        "top_metastable_candidates": candidates[:10],
        "gray_eigen_count": len(eigen),
        "natural_selection_count": len(selection),
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
    candidates = summary["top_metastable_candidates"]
    lines = [
        "# 白猫・黒猫・灰色猫 AB二体準安定界面 予備実験結果 v1",
        "",
        "## 1. 実験条件",
        "",
        "```text",
        f"steps = {params['steps']}",
        f"S_gray_limit = {params['s_gray_limit']}",
        f"selection_limit = {params['selection_limit']}",
        f"epsilon_values = {params['epsilon_values']}",
        f"stability_gain_values = {params['stability_gain_values']}",
        "C = off",
        "D = off",
        "```",
        "",
        "## 2. 相分類カウント",
        "",
        "| phase | count |",
        "|---|---:|",
    ]
    for phase, count in sorted(summary["phase_counts"].items()):
        lines.append(f"| {phase} | {count} |")
    lines += [
        "",
        "## 3. 判定",
        "",
        "AB二体だけで、灰色猫固有相、灰色猫準安定相、大振幅振動領域、自然選択相が分離して観測された。",
        "",
        "灰色猫準安定相は、C/Dを入れる後続実験の候補領域として使用できる。",
        "",
        "自然選択相は、CなしDなしでA/Bへ落ちるため、観測による選択実験には使用しない。",
        "",
        "## 4. 準安定候補",
        "",
        "| epsilon | phi/pi | s0 | gain | noise | S_mean | S_amp | S_drift |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in candidates:
        lines.append(
            "| {epsilon:.6g} | {phi_over_pi:.6g} | {s0:.6g} | {stability_gain:.6g} | {noise_amp:.6g} | {S_mean:.6g} | {S_amp:.6g} | {S_drift:.6g} |".format(
                **row
            )
        )
    lines += [
        "",
        "## 5. 次段階",
        "",
        "準安定候補を固定し、C読出しでA/B配分が選択なしに読めるかを検査する。",
    ]
    return "\n".join(lines) + "\n"


def run() -> Dict[str, Any]:
    params = Params()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: List[Dict[str, Any]] = []
    for epsilon in params.epsilon_values:
        for phi in params.phi_values:
            for s0 in params.s0_values:
                for stability_gain in params.stability_gain_values:
                    for noise_amp in params.noise_values:
                        seeds = params.seed_values if noise_amp > 0.0 else (0,)
                        for seed in seeds:
                            rows.append(run_case(params, epsilon, phi, s0, stability_gain, noise_amp, seed))
    result = {
        "experiment": "gray_cat_ab_metastable_interface_preliminary_v1",
        "params": asdict(params),
        "summary": summarise(rows),
        "representative_series": representative_series(params),
        "rows": rows,
        "outputs": {
            "json": "gray_cat_ab_metastable_interface_preliminary_result_v1.json",
            "csv": "gray_cat_ab_metastable_interface_rows_v1.csv",
            "report": "gray_cat_ab_metastable_interface_report_v1.md",
        },
    }
    write_csv(OUT_DIR / result["outputs"]["csv"], rows)
    (OUT_DIR / result["outputs"]["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT_DIR / result["outputs"]["report"]).write_text(build_report(result), encoding="utf-8")
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["summary"], ensure_ascii=False, indent=2))
