from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np

from run_abc_multigauge_interference_readout_v1 import Gauge, default_gauges, summarize_stage_readouts
from run_abc_multigauge_generalized_elastic_collision_velocity_sweep_v1 import (
    Params,
    readout_all_float,
    simulate,
    stage_quantities,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_multigauge_generalized_elastic_collision_noise_robustness_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


NOISE_CASES: List[Dict[str, float]] = [
    {"A_A": 1.0, "A_B": 1.0, "q_A0": 1.0, "q_B0": -1.0},
    {"A_A": 1.0, "A_B": 2.0, "q_A0": 1.4, "q_B0": -0.6},
    {"A_A": 1.0, "A_B": 2.0, "q_A0": 1.2, "q_B0": 0.2},
    {"A_A": 1.5, "A_B": 1.0, "q_A0": 1.8, "q_B0": -0.2},
]

NOISE_LEVELS: List[float] = [0.0, 1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4]
NOISE_MODES = ["zero_mean_gauge_noise", "common_bias_control"]


def dense_gauges(params: Params) -> List[Gauge]:
    gauges = list(default_gauges(params))
    h_values = [3.5e-4, 5.0e-4, 6.5e-4]
    phi_values = [0.0, 0.17, 0.41, 0.73]
    offsets = [-7.5e-4, 0.0, 7.5e-4]
    for h in h_values:
        for phi in phi_values:
            for delta_chi in offsets:
                for delta_tau in offsets:
                    gauges.append(
                        Gauge(
                            name=f"dense_h{h:.4f}_phi{phi:.2f}_chi{delta_chi:+.4f}_tau{delta_tau:+.4f}",
                            delta_chi=delta_chi,
                            delta_tau=delta_tau,
                            delta_phi=phi,
                            h_chi=h,
                            h_tau=h,
                            nh_chi_c=params.Nh_chi_C,
                            nh_tau_c=params.Nh_tau_C,
                            c_gain=params.A_C,
                        )
                    )
    return gauges


def stable_noise_value(*parts: Any) -> float:
    text = "|".join(str(part) for part in parts)
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    value = int.from_bytes(digest[:8], "big") / float(2**64 - 1)
    return 2.0 * value - 1.0


def group_rows(rows: Iterable[Dict[str, Any]], *keys: str) -> Dict[tuple[Any, ...], List[Dict[str, Any]]]:
    grouped: Dict[tuple[Any, ...], List[Dict[str, Any]]] = {}
    for row in rows:
        key = tuple(row[name] for name in keys)
        grouped.setdefault(key, []).append(row)
    return grouped


def perturb_rows(
    rows: List[Dict[str, Any]],
    case_name: str,
    noise_mode: str,
    noise_level: float,
) -> List[Dict[str, Any]]:
    perturbed = [dict(row) for row in rows]
    grouped = group_rows(perturbed, "stage", "particle")
    for (stage, particle), selected in grouped.items():
        raw_p = np.array(
            [stable_noise_value(case_name, noise_mode, noise_level, stage, particle, row["gauge"], "p") for row in selected]
        )
        raw_e = np.array(
            [stable_noise_value(case_name, noise_mode, noise_level, stage, particle, row["gauge"], "e") for row in selected]
        )
        raw_r = np.array(
            [stable_noise_value(case_name, noise_mode, noise_level, stage, particle, row["gauge"], "r") for row in selected]
        )
        if noise_mode == "zero_mean_gauge_noise":
            raw_p = raw_p - float(np.mean(raw_p))
            raw_e = raw_e - float(np.mean(raw_e))
            raw_r = raw_r - float(np.mean(raw_r))
        elif noise_mode == "common_bias_control":
            raw_p = np.full_like(raw_p, stable_noise_value(case_name, noise_mode, noise_level, stage, particle, "p"))
            raw_e = np.full_like(raw_e, stable_noise_value(case_name, noise_mode, noise_level, stage, particle, "e"))
            raw_r = np.full_like(raw_r, stable_noise_value(case_name, noise_mode, noise_level, stage, particle, "r"))
        else:
            raise ValueError(noise_mode)
        for index, row in enumerate(selected):
            p_scale = max(abs(float(row["p_expected"])), 1.0)
            e_scale = max(abs(float(row["E_expected"])), 1.0)
            r_scale = max(abs(float(row["R_expected"])), 1.0)
            row["noise_mode"] = noise_mode
            row["noise_level"] = noise_level
            row["p_clean"] = float(row["p_read"])
            row["E_clean"] = float(row["E_read"])
            row["R_clean"] = float(row["R_read"])
            row["p_read"] = float(row["p_read"]) + noise_level * p_scale * float(raw_p[index])
            row["E_read"] = float(row["E_read"]) + noise_level * e_scale * float(raw_e[index])
            row["R_read"] = float(row["R_read"]) + noise_level * r_scale * float(raw_r[index])
            row["p_abs_error"] = abs(float(row["p_read"]) - float(row["p_expected"]))
            row["E_abs_error"] = abs(float(row["E_read"]) - float(row["E_expected"]))
            row["R_abs_error"] = abs(float(row["R_read"]) - float(row["R_expected"]))
    return perturbed


def max_abs(rows: Iterable[Dict[str, Any]], key: str) -> float:
    return float(max(abs(float(row[key])) for row in rows))


def summary_lookup(summaries: List[Dict[str, Any]]) -> Dict[tuple[str, str], Dict[str, Any]]:
    return {(str(row["stage"]), str(row["particle"])): row for row in summaries}


def mean_readout_errors(gauge_rows: List[Dict[str, Any]], summaries: List[Dict[str, Any]]) -> Dict[str, float]:
    grouped = group_rows(gauge_rows, "stage", "particle")
    lookup = summary_lookup(summaries)
    p_errors: List[float] = []
    e_errors: List[float] = []
    r_errors: List[float] = []
    for (stage, particle), selected in grouped.items():
        expected_p = float(selected[0]["p_expected"])
        expected_e = float(selected[0]["E_expected"])
        expected_r = float(selected[0]["R_expected"])
        summary = lookup[(str(stage), str(particle))]
        p_errors.append(abs(float(summary["p_mean"]) - expected_p))
        e_errors.append(abs(float(summary["E_mean"]) - expected_e))
        r_errors.append(abs(float(summary["R_mean"]) - expected_r))
    return {
        "p_mean_abs_error": float(max(p_errors)),
        "E_mean_abs_error": float(max(e_errors)),
        "R_mean_abs_error": float(max(r_errors)),
    }


def within_particle_tr_separation(summaries: List[Dict[str, Any]]) -> float:
    ratios: List[float] = []
    for particle in ["A", "B"]:
        selected = [row for row in summaries if row["particle"] == particle]
        r_values = np.array([float(row["R_mean"]) for row in selected])
        t_values = np.array([float(row["t_mean"]) for row in selected])
        var_r = float(np.var(r_values))
        var_t = float(np.var(t_values))
        ratios.append(float(var_r / var_t) if var_t > 0.0 else float("inf"))
    return float(max(ratios))


def conservation_errors(quantities: List[Dict[str, Any]]) -> Dict[str, float]:
    by_stage = {str(row["stage"]): row for row in quantities}
    initial = by_stage["initial"]
    after = by_stage["collision_map"]
    final = by_stage["final"]
    relative_initial = float(initial["relative_p"])
    relative_after = float(after["relative_p"])
    return {
        "P_R_error": max(
            abs(float(after["P_R_total"]) - float(initial["P_R_total"])),
            abs(float(final["P_R_total"]) - float(initial["P_R_total"])),
        ),
        "K_R_error": max(
            abs(float(after["K_R_phase_total"]) - float(initial["K_R_phase_total"])),
            abs(float(final["K_R_phase_total"]) - float(initial["K_R_phase_total"])),
        ),
        "relative_flip_error": abs(relative_after + relative_initial),
        "E_tau_R_error": abs(float(final["E_tau_R_total"]) - float(initial["E_tau_R_total"])),
        "R_total_error": abs(float(final["R_total"]) - float(initial["R_total"])),
    }


def case_noise_summary(
    case_name: str,
    params: Params,
    noise_mode: str,
    noise_level: float,
    gauge_rows: List[Dict[str, Any]],
    summaries: List[Dict[str, Any]],
    quantities: List[Dict[str, Any]],
) -> Dict[str, Any]:
    mean_errors = mean_readout_errors(gauge_rows, summaries)
    conservation = conservation_errors(quantities)
    gauge_count = int(max(int(row["gauge_count"]) for row in summaries))
    max_r_std = float(max(float(row["R_std"]) for row in summaries))
    tr_ratio = within_particle_tr_separation(summaries)
    max_scale = max(1.0, params.A_A**2, params.A_B**2)
    unbiased_tol = max(1.0e-11, 1.0e-8 * max_scale)
    bias_detect_tol = max(1.0e-12, noise_level * 0.05)
    zero_mean = noise_mode == "zero_mean_gauge_noise"
    common_bias = noise_mode == "common_bias_control"
    zero_mean_valid = bool(
        zero_mean
        and mean_errors["p_mean_abs_error"] <= unbiased_tol
        and mean_errors["E_mean_abs_error"] <= unbiased_tol
        and mean_errors["R_mean_abs_error"] <= unbiased_tol
        and conservation["P_R_error"] <= 5.0e-8 * max_scale
        and conservation["K_R_error"] <= 5.0e-8 * max_scale
        and conservation["relative_flip_error"] <= 5.0e-8
    )
    bias_detected = bool(
        common_bias
        and noise_level > 0.0
        and max(
            mean_errors["p_mean_abs_error"],
            mean_errors["E_mean_abs_error"],
            mean_errors["R_mean_abs_error"] / max_scale,
        )
        >= bias_detect_tol
    )
    return {
        "case": case_name,
        "noise_mode": noise_mode,
        "noise_level": noise_level,
        "A_A": params.A_A,
        "A_B": params.A_B,
        "q_A0": params.q_A0,
        "q_B0": params.q_B0,
        "R_A": params.A_A**2,
        "R_B": params.A_B**2,
        "gauge_count": gauge_count,
        "p_max_abs_error": max_abs(gauge_rows, "p_abs_error"),
        "E_max_abs_error": max_abs(gauge_rows, "E_abs_error"),
        "R_max_abs_error": max_abs(gauge_rows, "R_abs_error"),
        "R_max_gauge_std": max_r_std,
        "within_particle_separation_ratio_time": tr_ratio,
        **mean_errors,
        **conservation,
        "zero_mean_multigauge_valid": zero_mean_valid,
        "common_bias_detected": bias_detected,
    }


def run_case(index: int, settings: Dict[str, float]) -> Dict[str, Any]:
    params = Params(**settings)
    gauges = dense_gauges(params)
    stages, events = simulate(params)
    clean_rows = readout_all_float(stages, gauges, params)
    case_name = (
        f"c{index:02d}_A{params.A_A:.2f}_B{params.A_B:.2f}_"
        f"u{params.q_A0:.2f}_v{params.q_B0:.2f}"
    )
    case_summaries: List[Dict[str, Any]] = []
    all_gauge_rows: List[Dict[str, Any]] = []
    all_quantity_rows: List[Dict[str, Any]] = []
    for noise_mode in NOISE_MODES:
        for noise_level in NOISE_LEVELS:
            rows = perturb_rows(clean_rows, case_name, noise_mode, noise_level)
            for row in rows:
                row["case"] = case_name
                row["A_A"] = params.A_A
                row["A_B"] = params.A_B
                row["q_A0"] = params.q_A0
                row["q_B0"] = params.q_B0
            summaries = summarize_stage_readouts(rows)
            quantities = stage_quantities(case_name, summaries)
            for row in quantities:
                row["noise_mode"] = noise_mode
                row["noise_level"] = noise_level
            case_summaries.append(case_noise_summary(case_name, params, noise_mode, noise_level, rows, summaries, quantities))
            all_gauge_rows.extend(rows)
            all_quantity_rows.extend(quantities)
    return {
        "case": case_name,
        "parameters": asdict(params),
        "gauges": [asdict(gauge) for gauge in gauges],
        "events": events,
        "case_summaries": case_summaries,
        "gauge_rows": all_gauge_rows,
        "stage_quantities": all_quantity_rows,
    }


def aggregate(case_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    rows = [row for entry in case_results for row in entry["case_summaries"]]
    zero_rows = [row for row in rows if row["noise_mode"] == "zero_mean_gauge_noise"]
    bias_rows = [
        row
        for row in rows
        if row["noise_mode"] == "common_bias_control" and float(row["noise_level"]) >= 1.0e-10
    ]
    return {
        "case_count": len(case_results),
        "noise_mode_count": len(NOISE_MODES),
        "noise_level_count": len(NOISE_LEVELS),
        "total_summary_rows": len(rows),
        "max_gauge_count": int(max(int(row["gauge_count"]) for row in rows)),
        "zero_mean_multigauge_valid_all": all(bool(row["zero_mean_multigauge_valid"]) for row in zero_rows),
        "common_bias_detection_floor": 1.0e-10,
        "common_bias_detected_all_above_floor": all(bool(row["common_bias_detected"]) for row in bias_rows),
        "zero_mean_max_p_mean_abs_error": float(max(float(row["p_mean_abs_error"]) for row in zero_rows)),
        "zero_mean_max_E_mean_abs_error": float(max(float(row["E_mean_abs_error"]) for row in zero_rows)),
        "zero_mean_max_R_mean_abs_error": float(max(float(row["R_mean_abs_error"]) for row in zero_rows)),
        "zero_mean_max_P_R_error": float(max(float(row["P_R_error"]) for row in zero_rows)),
        "zero_mean_max_K_R_error": float(max(float(row["K_R_error"]) for row in zero_rows)),
        "zero_mean_max_relative_flip_error": float(max(float(row["relative_flip_error"]) for row in zero_rows)),
        "biased_control_max_p_mean_abs_error": float(max(float(row["p_mean_abs_error"]) for row in bias_rows)),
        "biased_control_max_E_mean_abs_error": float(max(float(row["E_mean_abs_error"]) for row in bias_rows)),
        "biased_control_max_R_mean_abs_error": float(max(float(row["R_mean_abs_error"]) for row in bias_rows)),
        "single_gauge_only_used": False,
        "noise_robustness_valid": bool(
            all(bool(row["zero_mean_multigauge_valid"]) for row in zero_rows)
            and all(bool(row["common_bias_detected"]) for row in bias_rows)
        ),
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    keys = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def make_plots(rows: List[Dict[str, Any]]) -> None:
    zero_rows = [row for row in rows if row["noise_mode"] == "zero_mean_gauge_noise"]
    grouped = group_rows(zero_rows, "noise_level")
    levels = sorted(float(key[0]) for key in grouped.keys())
    p_errors = [max(float(row["p_mean_abs_error"]) for row in grouped[(level,)]) for level in levels]
    k_errors = [max(float(row["K_R_error"]) for row in grouped[(level,)]) for level in levels]
    r_errors = [max(float(row["R_mean_abs_error"]) for row in grouped[(level,)]) for level in levels]
    floor = 1.0e-18
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(levels, np.maximum(p_errors, floor), marker="o", label="p mean error")
    ax.plot(levels, np.maximum(k_errors, floor), marker="o", label="R*p^2 conservation error")
    ax.plot(levels, np.maximum(r_errors, floor), marker="o", label="R mean error")
    ax.set_xscale("symlog", linthresh=1.0e-13)
    ax.set_yscale("log")
    ax.set_xlabel("zero-mean gauge noise level")
    ax.set_ylabel("max error")
    ax.set_title("zero-mean multigauge noise cancellation")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_noise_robustness_zero_mean_v1.png", dpi=160)
    plt.close(fig)

    bias_rows = [row for row in rows if row["noise_mode"] == "common_bias_control" and float(row["noise_level"]) > 0.0]
    grouped_bias = group_rows(bias_rows, "noise_level")
    bias_levels = sorted(float(key[0]) for key in grouped_bias.keys())
    bias_p_errors = [max(float(row["p_mean_abs_error"]) for row in grouped_bias[(level,)]) for level in bias_levels]
    bias_r_errors = [max(float(row["R_mean_abs_error"]) for row in grouped_bias[(level,)]) for level in bias_levels]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(bias_levels, bias_p_errors, marker="o", label="p mean error")
    ax.plot(bias_levels, bias_r_errors, marker="o", label="R mean error")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("common bias level")
    ax.set_ylabel("max mean error")
    ax.set_title("common readout bias remains detectable")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_multigauge_generalized_elastic_collision_noise_robustness_bias_control_v1.png", dpi=160)
    plt.close(fig)


def write_report(result: Dict[str, Any]) -> None:
    lines = [
        "# ABC Multigauge Generalized Elastic Collision Noise Robustness v1",
        "",
        "## Purpose",
        "",
        "This experiment injects deterministic readout-side noise into p/E/R gauge rows after the physical state has been simulated.",
        "Zero-mean gauge noise is expected to cancel by multigauge averaging; common bias is expected to remain detectable.",
        "",
        "## Aggregate Verdict",
        "",
    ]
    for key, value in result["aggregate_verdict"].items():
        lines.append(f"- {key}: `{value}`")
    lines.extend(
        [
            "",
            "## Zero-Mean Summary",
            "",
            "| noise level | max p mean err | max R mean err | max R*p^2 err | valid all cases |",
            "|---:|---:|---:|---:|---|",
        ]
    )
    zero_rows = [row for row in result["summary_rows"] if row["noise_mode"] == "zero_mean_gauge_noise"]
    for level in NOISE_LEVELS:
        selected = [row for row in zero_rows if math.isclose(float(row["noise_level"]), level, rel_tol=0.0, abs_tol=0.0)]
        lines.append(
            f"| {level:.1e} | "
            f"{max(float(row['p_mean_abs_error']) for row in selected):.16e} | "
            f"{max(float(row['R_mean_abs_error']) for row in selected):.16e} | "
            f"{max(float(row['K_R_error']) for row in selected):.16e} | "
            f"`{all(bool(row['zero_mean_multigauge_valid']) for row in selected)}` |"
        )
    lines.extend(
        [
            "",
            "## Common Bias Control",
            "",
            "| noise level | max p mean err | max R mean err | detected all cases |",
            "|---:|---:|---:|---|",
        ]
    )
    bias_rows = [row for row in result["summary_rows"] if row["noise_mode"] == "common_bias_control" and float(row["noise_level"]) > 0.0]
    for level in [level for level in NOISE_LEVELS if level > 0.0]:
        selected = [row for row in bias_rows if math.isclose(float(row["noise_level"]), level, rel_tol=0.0, abs_tol=0.0)]
        lines.append(
            f"| {level:.1e} | "
            f"{max(float(row['p_mean_abs_error']) for row in selected):.16e} | "
            f"{max(float(row['R_mean_abs_error']) for row in selected):.16e} | "
            f"`{all(bool(row['common_bias_detected']) for row in selected)}` |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_generalized_elastic_collision_noise_robustness_result_v1.json` |",
            "| summary CSV | `abc_multigauge_generalized_elastic_collision_noise_robustness_summary_v1.csv` |",
            "| quantity CSV | `abc_multigauge_generalized_elastic_collision_noise_robustness_stage_quantities_v1.csv` |",
            "| gauge CSV | `abc_multigauge_generalized_elastic_collision_noise_robustness_gauge_rows_v1.csv` |",
            "| zero mean plot | `abc_multigauge_generalized_elastic_collision_noise_robustness_zero_mean_v1.png` |",
            "| bias plot | `abc_multigauge_generalized_elastic_collision_noise_robustness_bias_control_v1.png` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_noise_robustness_report_v1.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    case_results = [run_case(index + 1, settings) for index, settings in enumerate(NOISE_CASES)]
    summary_rows = [row for entry in case_results for row in entry["case_summaries"]]
    return {
        "experiment": "abc_multigauge_generalized_elastic_collision_noise_robustness_v1",
        "noise_cases": NOISE_CASES,
        "noise_modes": NOISE_MODES,
        "noise_levels": NOISE_LEVELS,
        "summary_rows": summary_rows,
        "aggregate_verdict": aggregate(case_results),
        "case_results": [
            {
                "case": entry["case"],
                "parameters": entry["parameters"],
                "events": entry["events"],
                "case_summaries": entry["case_summaries"],
            }
            for entry in case_results
        ],
        "note": (
            "Noise is injected after state simulation into readout rows only. "
            "Zero-mean gauge noise tests multigauge cancellation; common bias is a detection control."
        ),
    }


def write_outputs(result: Dict[str, Any]) -> None:
    case_results = [run_case(index + 1, settings) for index, settings in enumerate(NOISE_CASES)]
    gauge_rows: List[Dict[str, Any]] = []
    quantity_rows: List[Dict[str, Any]] = []
    for entry in case_results:
        gauge_rows.extend(entry["gauge_rows"])
        quantity_rows.extend(entry["stage_quantities"])
    (OUT_DIR / "abc_multigauge_generalized_elastic_collision_noise_robustness_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(
        OUT_DIR / "abc_multigauge_generalized_elastic_collision_noise_robustness_summary_v1.csv",
        result["summary_rows"],
    )
    write_csv(
        OUT_DIR / "abc_multigauge_generalized_elastic_collision_noise_robustness_stage_quantities_v1.csv",
        quantity_rows,
    )
    write_csv(
        OUT_DIR / "abc_multigauge_generalized_elastic_collision_noise_robustness_gauge_rows_v1.csv",
        gauge_rows,
    )
    make_plots(result["summary_rows"])
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
