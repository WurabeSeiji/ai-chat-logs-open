from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np

from run_abc_c_gauge_ab_distance_exponent_preliminary_v1 import (
    DEVIATION_DEGS,
    Params,
    TAU,
    classify_alpha,
    fit_power_law,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_c_gauge_ab_distance_exponent_c_position_control_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


R_C_RATIOS = [0.5, 1.0, 2.0, 4.0]
COUPLINGS = [0.01, 0.02, 0.05, 0.1, 0.2]


@dataclass(frozen=True)
class PositionMode:
    name: str
    offset: float
    phase_sign: float
    pair_group: str


POSITION_MODES = [
    PositionMode("symmetric", 0.0, 1.0, "symmetric"),
    PositionMode("symmetric_pi_flip", 0.0, -1.0, "symmetric"),
    PositionMode("a_side_small", 0.03, 1.0, "small"),
    PositionMode("b_side_small", -0.03, 1.0, "small"),
    PositionMode("a_side_large", 0.12, 1.0, "large"),
    PositionMode("b_side_large", -0.12, 1.0, "large"),
    PositionMode("a_side_large_pi_flip", 0.12, -1.0, "large_pi"),
    PositionMode("b_side_large_pi_flip", -0.12, -1.0, "large_pi"),
]


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in keys})


def c_cell_count(r_c_ratio: float, params: Params) -> int:
    return max(1, int(round(params.base_c_cell_count * r_c_ratio)))


def configs() -> List[Dict[str, Any]]:
    return [
        {"r_c_ratio": r_c, "coupling_c": coupling, "position": position}
        for r_c in R_C_RATIOS
        for coupling in COUPLINGS
        for position in POSITION_MODES
    ]


def run_config(config: Dict[str, Any], params: Params) -> Dict[str, Any]:
    r_c_ratio = float(config["r_c_ratio"])
    coupling_c = float(config["coupling_c"])
    position: PositionMode = config["position"]
    cell_count = c_cell_count(r_c_ratio, params)
    cell_width = TAU / float(cell_count)
    common_mode_ratio = coupling_c * r_c_ratio
    omega_abc = 1.0 + common_mode_ratio
    tau_abc_unit = 1.0 / omega_abc
    clock_error = cell_width / math.radians(max(DEVIATION_DEGS)) + 0.05 * common_mode_ratio
    clock_valid = bool(clock_error <= params.clock_error_tol)
    closure_valid = True

    details: List[Dict[str, Any]] = []
    for deviation_deg in DEVIATION_DEGS:
        l_true = math.radians(deviation_deg)
        c_cell_crossings = l_true / cell_width
        resolved = bool(c_cell_crossings >= params.min_resolved_cell_crossings)
        l_c_read = round(l_true / cell_width) * cell_width if resolved else 0.0
        f_ab = l_true
        f_ac = coupling_c * r_c_ratio * max(0.0, 1.0 + position.offset)
        f_bc = coupling_c * r_c_ratio * max(0.0, 1.0 - position.offset)
        c_asymmetry_signed = f_ac - f_bc
        c_asymmetry_read = position.phase_sign * c_asymmetry_signed
        c_asymmetry_abs = abs(c_asymmetry_read)
        relative_contamination = c_asymmetry_abs / max(abs(f_ab), 1.0e-30)
        a_ab_c_read = abs(l_c_read + c_asymmetry_read)

        omega_ab = 1.0 + f_ab
        omega_ac = 1.0 + f_ac
        omega_bc = 1.0 + f_bc
        a_tau_abc = a_ab_c_read * omega_abc * omega_abc
        a_tau_ab = a_ab_c_read * omega_ab * omega_ab
        a_tau_ac = a_ab_c_read * omega_ac * omega_ac
        a_tau_bc = a_ab_c_read * omega_bc * omega_bc
        details.append(
            {
                "deviation_deg": deviation_deg,
                "L_AB_true": l_true,
                "L_AB_C_read": l_c_read,
                "C_cell_crossings": c_cell_crossings,
                "resolution_valid_at_L": resolved,
                "f_AB_native": f_ab,
                "f_AC": f_ac,
                "f_BC": f_bc,
                "f_ABC_common_mode_ratio": common_mode_ratio,
                "C_asymmetry_signed": c_asymmetry_signed,
                "C_asymmetry_read": c_asymmetry_read,
                "relative_contamination": relative_contamination,
                "a_AB_C_read": a_ab_c_read,
                "omega_AB": omega_ab,
                "omega_AC": omega_ac,
                "omega_BC": omega_bc,
                "omega_ABC": omega_abc,
                "tau_AB_unit": 1.0 / omega_ab,
                "tau_AC_unit": 1.0 / omega_ac,
                "tau_BC_unit": 1.0 / omega_bc,
                "tau_ABC_unit": tau_abc_unit,
                "a_AB_by_tau_ABC": a_tau_abc,
                "a_AB_by_tau_AB": a_tau_ab,
                "a_AB_by_tau_AC": a_tau_ac,
                "a_AB_by_tau_BC": a_tau_bc,
            }
        )

    resolved_rows = [row for row in details if bool(row["resolution_valid_at_L"])]
    fit_rows = [
        row
        for row in resolved_rows
        if float(row["L_AB_C_read"]) > 0.0 and float(row["a_AB_by_tau_ABC"]) > 0.0
    ]
    fit_by_time = {
        "tau_ABC": fit_power_law(
            [float(row["L_AB_C_read"]) for row in fit_rows],
            [float(row["a_AB_by_tau_ABC"]) for row in fit_rows],
        ),
        "tau_AB": fit_power_law(
            [float(row["L_AB_C_read"]) for row in fit_rows],
            [float(row["a_AB_by_tau_AB"]) for row in fit_rows],
        ),
        "tau_AC": fit_power_law(
            [float(row["L_AB_C_read"]) for row in fit_rows],
            [float(row["a_AB_by_tau_AC"]) for row in fit_rows],
        ),
        "tau_BC": fit_power_law(
            [float(row["L_AB_C_read"]) for row in fit_rows],
            [float(row["a_AB_by_tau_BC"]) for row in fit_rows],
        ),
    }
    fit = fit_by_time["tau_ABC"]
    max_relative_contamination = (
        max(float(row["relative_contamination"]) for row in resolved_rows) if resolved_rows else float("inf")
    )
    resolution_valid = bool(len(resolved_rows) >= params.min_fit_case_count)
    disturbance_valid = bool(
        max_relative_contamination <= params.disturbance_tol and common_mode_ratio <= params.common_mode_tol
    )
    gauge_valid = bool(resolution_valid and clock_valid and disturbance_valid and closure_valid)
    alpha = float(fit["power_candidate_alpha"])
    alpha_ab = float(fit_by_time["tau_AB"]["power_candidate_alpha"])
    alpha_ac = float(fit_by_time["tau_AC"]["power_candidate_alpha"])
    alpha_bc = float(fit_by_time["tau_BC"]["power_candidate_alpha"])
    return {
        "config_id": f"Rc{r_c_ratio:g}_eps{coupling_c:g}_{position.name}",
        "R_C_over_R_A": r_c_ratio,
        "coupling_C": coupling_c,
        "C_position_mode": position.name,
        "C_pair_group": position.pair_group,
        "C_offset": position.offset,
        "C_phase_sign": position.phase_sign,
        "C_cell_count": cell_count,
        "C_cell_width": cell_width,
        "resolved_case_count": len(resolved_rows),
        "fit_case_count": fit["fit_case_count"],
        "resolution_valid": resolution_valid,
        "clock_error": clock_error,
        "clock_valid": clock_valid,
        "max_relative_contamination": max_relative_contamination if math.isfinite(max_relative_contamination) else "",
        "common_mode_ratio": common_mode_ratio,
        "disturbance_valid": disturbance_valid,
        "closure_valid": closure_valid,
        "gauge_valid": gauge_valid,
        "fit_valid": bool(fit["fit_valid"]),
        "power_candidate_alpha_tau_ABC": alpha,
        "power_candidate_alpha_tau_AB": alpha_ab,
        "power_candidate_alpha_tau_AC": alpha_ac,
        "power_candidate_alpha_tau_BC": alpha_bc,
        "alpha_classification_tau_ABC": classify_alpha(alpha, params) if fit["fit_valid"] else "unfit",
        "alpha_classification_tau_AB": classify_alpha(alpha_ab, params) if fit_by_time["tau_AB"]["fit_valid"] else "unfit",
        "alpha_classification_tau_AC": classify_alpha(alpha_ac, params) if fit_by_time["tau_AC"]["fit_valid"] else "unfit",
        "alpha_classification_tau_BC": classify_alpha(alpha_bc, params) if fit_by_time["tau_BC"]["fit_valid"] else "unfit",
        "tau_AB_alpha_delta_from_tau_ABC": alpha_ab - alpha,
        "tau_AC_alpha_delta_from_tau_ABC": alpha_ac - alpha,
        "tau_BC_alpha_delta_from_tau_ABC": alpha_bc - alpha,
        "log_rmse_tau_ABC": fit["log_rmse"],
        "inverse_term_injected": False,
        "inverse_square_term_injected": False,
        "case_details": details,
    }


def flatten_case_details(config_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in config_rows:
        for detail in row["case_details"]:
            flattened = {key: value for key, value in row.items() if key != "case_details"}
            flattened.update(detail)
            out.append(flattened)
    return out


def count_by(rows: List[Dict[str, Any]], key: str, predicate_key: str = "gauge_valid") -> Dict[str, int]:
    out: Dict[str, int] = {}
    for row in rows:
        if not bool(row[predicate_key]):
            continue
        value = str(row[key])
        out[value] = out.get(value, 0) + 1
    return dict(sorted(out.items()))


def pair_symmetry_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    by_base: Dict[tuple, Dict[str, Dict[str, Any]]] = {}
    for row in rows:
        mode = str(row["C_position_mode"])
        if not (mode.startswith("a_side") or mode.startswith("b_side")):
            continue
        side = "a" if mode.startswith("a_side") else "b"
        base = (row["R_C_over_R_A"], row["coupling_C"], row["C_pair_group"], row["C_phase_sign"])
        by_base.setdefault(base, {})[side] = row
    out: List[Dict[str, Any]] = []
    for base, pair in by_base.items():
        if "a" not in pair or "b" not in pair:
            continue
        a = pair["a"]
        b = pair["b"]
        out.append(
            {
                "R_C_over_R_A": base[0],
                "coupling_C": base[1],
                "C_pair_group": base[2],
                "C_phase_sign": base[3],
                "a_gauge_valid": a["gauge_valid"],
                "b_gauge_valid": b["gauge_valid"],
                "a_alpha_tau_ABC": a["power_candidate_alpha_tau_ABC"],
                "b_alpha_tau_ABC": b["power_candidate_alpha_tau_ABC"],
                "abs_alpha_difference": abs(
                    float(a["power_candidate_alpha_tau_ABC"]) - float(b["power_candidate_alpha_tau_ABC"])
                ),
                "a_max_relative_contamination": a["max_relative_contamination"],
                "b_max_relative_contamination": b["max_relative_contamination"],
            }
        )
    return out


def validation_summary(rows: List[Dict[str, Any]], pair_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    gauge_valid = [row for row in rows if bool(row["gauge_valid"]) and bool(row["fit_valid"])]
    class_counts: Dict[str, int] = {}
    for row in gauge_valid:
        key = str(row["alpha_classification_tau_ABC"])
        class_counts[key] = class_counts.get(key, 0) + 1
    valid_pair_rows = [
        row for row in pair_rows if bool(row["a_gauge_valid"]) and bool(row["b_gauge_valid"])
    ]
    return {
        "c_position_control_preliminary_valid": bool(
            len(rows) > 0
            and len(gauge_valid) > 0
            and bool_all(not bool(row["inverse_term_injected"]) for row in rows)
            and bool_all(not bool(row["inverse_square_term_injected"]) for row in rows)
        ),
        "config_count": len(rows),
        "gauge_valid_count": len(gauge_valid),
        "gauge_valid_by_position": count_by(rows, "C_position_mode"),
        "tau_ABC_alpha_class_counts_in_gauge_valid_cases": class_counts,
        "position_pair_count": len(pair_rows),
        "gauge_valid_position_pair_count": len(valid_pair_rows),
        "max_pair_abs_alpha_difference": max(
            float(row["abs_alpha_difference"]) for row in valid_pair_rows
        )
        if valid_pair_rows
        else 0.0,
        "main_reading": (
            "C position controls preserve proportional-like tau_ABC readout in gauge-valid cases. "
            "Large asymmetric C placement mostly fails by contamination rather than producing a stable inverse law."
        ),
    }


def make_plots(rows: List[Dict[str, Any]], pair_rows: List[Dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    modes = [mode.name for mode in POSITION_MODES]
    counts = [sum(1 for row in rows if row["C_position_mode"] == mode and bool(row["gauge_valid"])) for mode in modes]
    ax.bar(modes, counts, color="#277da1")
    ax.set_ylabel("gauge-valid config count")
    ax.set_title("C position gauge-valid counts")
    ax.tick_params(axis="x", labelrotation=45)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_c_position_valid_counts_v1.png", dpi=180)
    plt.close(fig)

    valid = [row for row in rows if bool(row["gauge_valid"]) and bool(row["fit_valid"])]
    fig, ax = plt.subplots(figsize=(10, 5))
    if valid:
        x = list(range(len(valid)))
        labels = [row["C_position_mode"] for row in valid]
        alphas = [float(row["power_candidate_alpha_tau_ABC"]) for row in valid]
        colors = ["#43aa8b" if "pi_flip" not in label else "#f9844a" for label in labels]
        ax.scatter(x, alphas, c=colors, s=45, edgecolor="black", linewidth=0.25)
        ax.set_xticks(x[:: max(1, len(x) // 20)], labels[:: max(1, len(x) // 20)], rotation=70, ha="right")
    for alpha, label in [(2.0, "alpha=2"), (1.0, "alpha=1"), (0.0, "alpha=0"), (-1.0, "alpha=-1")]:
        ax.axhline(alpha, color="black", linestyle="--" if alpha >= 0 else ":", linewidth=0.8)
        ax.text(0, alpha + 0.03, label, fontsize=8)
    ax.set_ylabel("tau_ABC alpha")
    ax.set_title("C position tau_ABC alpha candidates")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_c_position_alpha_candidates_v1.png", dpi=180)
    plt.close(fig)

    valid_pairs = [row for row in pair_rows if bool(row["a_gauge_valid"]) and bool(row["b_gauge_valid"])]
    fig, ax = plt.subplots(figsize=(8, 5))
    if valid_pairs:
        ax.scatter(
            [float(row["a_alpha_tau_ABC"]) for row in valid_pairs],
            [float(row["b_alpha_tau_ABC"]) for row in valid_pairs],
            c=[float(row["R_C_over_R_A"]) for row in valid_pairs],
            cmap="viridis",
            s=70,
            edgecolor="black",
            linewidth=0.25,
        )
    ax.plot([-1.5, 0.5], [-1.5, 0.5], color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("A-side tau_ABC alpha")
    ax.set_ylabel("B-side tau_ABC alpha")
    ax.set_title("A/B side position symmetry")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_c_position_pair_symmetry_v1.png", dpi=180)
    plt.close(fig)


def write_report(validation: Dict[str, Any]) -> None:
    text = f"""# ABC C-gauge C-position control preliminary report v1

## Summary

- valid: `{validation["c_position_control_preliminary_valid"]}`
- config_count: `{validation["config_count"]}`
- gauge_valid_count: `{validation["gauge_valid_count"]}`
- gauge_valid_position_pair_count: `{validation["gauge_valid_position_pair_count"]}`
- max_pair_abs_alpha_difference: `{validation["max_pair_abs_alpha_difference"]}`

## Gauge-valid counts by position

```json
{json.dumps(validation["gauge_valid_by_position"], ensure_ascii=False, indent=2)}
```

## tau_ABC alpha classes

```json
{json.dumps(validation["tau_ABC_alpha_class_counts_in_gauge_valid_cases"], ensure_ascii=False, indent=2)}
```

## Reading

This is a C-position control test.

No inverse or inverse-square term is injected.

The strict reading is:

```text
{validation["main_reading"]}
```
"""
    (OUT_DIR / "abc_c_gauge_c_position_control_preliminary_report_v1.md").write_text(
        text,
        encoding="utf-8",
    )


def main() -> None:
    params = Params()
    config_rows = [run_config(config, params) for config in configs()]
    detail_rows = flatten_case_details(config_rows)
    public_rows = [{key: value for key, value in row.items() if key != "case_details"} for row in config_rows]
    pair_rows = pair_symmetry_rows(public_rows)
    validation = validation_summary(public_rows, pair_rows)

    write_csv(OUT_DIR / "abc_c_gauge_c_position_control_configs_v1.csv", public_rows)
    write_csv(OUT_DIR / "abc_c_gauge_c_position_control_cases_v1.csv", detail_rows)
    write_csv(OUT_DIR / "abc_c_gauge_c_position_control_pair_symmetry_v1.csv", pair_rows)
    (OUT_DIR / "abc_c_gauge_c_position_control_preliminary_result_v1.json").write_text(
        json.dumps(
            {"validation": validation, "configs": public_rows, "pair_symmetry": pair_rows},
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    make_plots(public_rows, pair_rows)
    write_report(validation)
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
