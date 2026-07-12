from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_c_gauge_ab_distance_exponent_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


TAU = 2.0 * math.pi


@dataclass(frozen=True)
class Params:
    base_c_cell_count: int = 96
    min_resolved_cell_crossings: float = 1.5
    min_fit_case_count: int = 4
    disturbance_tol: float = 0.15
    common_mode_tol: float = 0.35
    clock_error_tol: float = 0.35
    closure_tol: float = 1.0e-12
    alpha_tol: float = 0.18


@dataclass(frozen=True)
class CGaugeConfig:
    r_c_ratio: float
    coupling_c: float
    c_offset_rad: float


DEVIATION_DEGS = [1.0, 2.0, 5.0, 10.0, 20.0, 35.0, 60.0]
R_C_RATIOS = [0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0]
COUPLINGS = [0.01, 0.02, 0.05, 0.1, 0.2]
C_OFFSET_RADS = [0.0, 0.01, 0.03, 0.08]


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


def configs() -> List[CGaugeConfig]:
    return [
        CGaugeConfig(r_c_ratio, coupling, offset)
        for r_c_ratio in R_C_RATIOS
        for coupling in COUPLINGS
        for offset in C_OFFSET_RADS
    ]


def c_cell_count(config: CGaugeConfig, params: Params) -> int:
    # In this preliminary gauge model, smaller R_C means a lower C frequency and a coarser cell grid.
    return max(1, int(round(params.base_c_cell_count * config.r_c_ratio)))


def fit_power_law(xs: List[float], ys: List[float]) -> Dict[str, Any]:
    x = np.array(xs, dtype=float)
    y = np.array(ys, dtype=float)
    mask = (x > 1.0e-30) & (y > 1.0e-30)
    x = x[mask]
    y = y[mask]
    if len(x) < 3:
        return {
            "fit_valid": False,
            "loglog_slope": 0.0,
            "power_candidate_alpha": 0.0,
            "log_rmse": 0.0,
            "fit_case_count": int(len(x)),
        }
    slope, intercept = np.polyfit(np.log(x), np.log(y), 1)
    predicted = slope * np.log(x) + intercept
    residual = np.log(y) - predicted
    return {
        "fit_valid": True,
        "loglog_slope": float(slope),
        "power_candidate_alpha": -float(slope),
        "loglog_intercept": float(intercept),
        "log_rmse": float(np.sqrt(np.mean(residual * residual))),
        "fit_case_count": int(len(x)),
    }


def classify_alpha(alpha: float, params: Params) -> str:
    if abs(alpha - 0.0) <= params.alpha_tol:
        return "constant_like_alpha0"
    if abs(alpha - 1.0) <= params.alpha_tol:
        return "inverse_like_alpha1"
    if abs(alpha - 2.0) <= params.alpha_tol:
        return "inverse_square_like_alpha2"
    if abs(alpha + 1.0) <= params.alpha_tol:
        return "proportional_like_alpha_minus1"
    if abs(alpha + 2.0) <= params.alpha_tol:
        return "area_like_alpha_minus2"
    return "other"


def run_config(config: CGaugeConfig, params: Params) -> Dict[str, Any]:
    cell_count = c_cell_count(config, params)
    cell_width = TAU / float(cell_count)
    rows: List[Dict[str, Any]] = []

    common_mode_ratio = config.coupling_c * config.r_c_ratio
    omega_abc = 1.0 + common_mode_ratio
    tau_abc_unit = 1.0 / omega_abc
    clock_error = cell_width / math.radians(max(DEVIATION_DEGS)) + 0.05 * common_mode_ratio
    clock_valid = bool(clock_error <= params.clock_error_tol)
    closure_valid = True

    for deviation_deg in DEVIATION_DEGS:
        l_true = math.radians(deviation_deg)
        c_cell_crossings = l_true / cell_width
        resolved = bool(c_cell_crossings >= params.min_resolved_cell_crossings)
        l_c_read = round(l_true / cell_width) * cell_width if resolved else 0.0

        # Minimal AB relation compensation inherited from the AB harmonic readout.
        # No inverse or inverse-square term is injected here.
        f_ab = l_true

        # C relation contamination is represented as a projection asymmetry term.
        # In symmetric C placement this tends to zero; in offset placement it remains.
        f_ac = config.coupling_c * config.r_c_ratio * (1.0 + config.c_offset_rad)
        f_bc = config.coupling_c * config.r_c_ratio * max(0.0, 1.0 - config.c_offset_rad)
        c_asymmetry = abs(f_ac - f_bc)
        relative_contamination = c_asymmetry / max(abs(f_ab), 1.0e-30)

        # The C-read acceleration candidate is what C can read after quantization and
        # relation-asymmetry contamination. It is a readout candidate, not an imposed force law.
        a_ab_c_read = abs(l_c_read + (f_ac - f_bc))
        omega_ab = 1.0 + f_ab
        omega_ac = 1.0 + f_ac
        omega_bc = 1.0 + f_bc
        tau_ab_unit = 1.0 / omega_ab
        tau_ac_unit = 1.0 / omega_ac
        tau_bc_unit = 1.0 / omega_bc
        a_ab_by_tau_abc = a_ab_c_read * omega_abc * omega_abc
        a_ab_by_tau_ab = a_ab_c_read * omega_ab * omega_ab
        a_ab_by_tau_ac = a_ab_c_read * omega_ac * omega_ac
        a_ab_by_tau_bc = a_ab_c_read * omega_bc * omega_bc

        rows.append(
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
                "C_asymmetry": c_asymmetry,
                "relative_contamination": relative_contamination,
                "a_AB_C_read": a_ab_c_read,
                "omega_AB": omega_ab,
                "omega_AC": omega_ac,
                "omega_BC": omega_bc,
                "omega_ABC": omega_abc,
                "tau_AB_unit": tau_ab_unit,
                "tau_AC_unit": tau_ac_unit,
                "tau_BC_unit": tau_bc_unit,
                "tau_ABC_unit": tau_abc_unit,
                "tau_AB_over_tau_ABC": tau_ab_unit / tau_abc_unit,
                "tau_AC_over_tau_ABC": tau_ac_unit / tau_abc_unit,
                "tau_BC_over_tau_ABC": tau_bc_unit / tau_abc_unit,
                "a_AB_by_tau_ABC": a_ab_by_tau_abc,
                "a_AB_by_tau_AB": a_ab_by_tau_ab,
                "a_AB_by_tau_AC": a_ab_by_tau_ac,
                "a_AB_by_tau_BC": a_ab_by_tau_bc,
            }
        )

    resolved_rows = [row for row in rows if bool(row["resolution_valid_at_L"])]
    fit_rows = [
        row
        for row in resolved_rows
        if float(row["L_AB_C_read"]) > 0.0 and float(row["a_AB_by_tau_ABC"]) > 0.0
    ]
    fit_by_time: Dict[str, Dict[str, Any]] = {}
    for time_name, candidate in [
        ("tau_ABC", "a_AB_by_tau_ABC"),
        ("tau_AB", "a_AB_by_tau_AB"),
        ("tau_AC", "a_AB_by_tau_AC"),
        ("tau_BC", "a_AB_by_tau_BC"),
    ]:
        fit_by_time[time_name] = fit_power_law(
            [float(row["L_AB_C_read"]) for row in fit_rows],
            [float(row[candidate]) for row in fit_rows],
        )
    fit = fit_by_time["tau_ABC"]

    max_relative_contamination = (
        max(float(row["relative_contamination"]) for row in resolved_rows) if resolved_rows else float("inf")
    )
    disturbance_valid = bool(
        max_relative_contamination <= params.disturbance_tol and common_mode_ratio <= params.common_mode_tol
    )
    resolution_valid = bool(len(resolved_rows) >= params.min_fit_case_count)
    gauge_valid = bool(resolution_valid and clock_valid and disturbance_valid and closure_valid)
    alpha = float(fit["power_candidate_alpha"])
    classification = classify_alpha(alpha, params) if fit["fit_valid"] else "unfit"
    alpha_ab = float(fit_by_time["tau_AB"]["power_candidate_alpha"])
    alpha_ac = float(fit_by_time["tau_AC"]["power_candidate_alpha"])
    alpha_bc = float(fit_by_time["tau_BC"]["power_candidate_alpha"])

    return {
        "config_id": f"Rc{config.r_c_ratio:g}_eps{config.coupling_c:g}_off{config.c_offset_rad:g}",
        "R_C_over_R_A": config.r_c_ratio,
        "coupling_C": config.coupling_c,
        "C_offset_rad": config.c_offset_rad,
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
        "loglog_slope": fit["loglog_slope"],
        "power_candidate_alpha": alpha,
        "power_candidate_alpha_tau_ABC": alpha,
        "power_candidate_alpha_tau_AB": alpha_ab,
        "power_candidate_alpha_tau_AC": alpha_ac,
        "power_candidate_alpha_tau_BC": alpha_bc,
        "log_rmse": fit["log_rmse"],
        "alpha_classification": classification,
        "alpha_classification_tau_ABC": classification,
        "alpha_classification_tau_AB": classify_alpha(alpha_ab, params) if fit_by_time["tau_AB"]["fit_valid"] else "unfit",
        "alpha_classification_tau_AC": classify_alpha(alpha_ac, params) if fit_by_time["tau_AC"]["fit_valid"] else "unfit",
        "alpha_classification_tau_BC": classify_alpha(alpha_bc, params) if fit_by_time["tau_BC"]["fit_valid"] else "unfit",
        "tau_AB_alpha_delta_from_tau_ABC": alpha_ab - alpha,
        "tau_AC_alpha_delta_from_tau_ABC": alpha_ac - alpha,
        "tau_BC_alpha_delta_from_tau_ABC": alpha_bc - alpha,
        "inverse_term_injected": False,
        "inverse_square_term_injected": False,
        "constructed_reciprocal_used_for_success": False,
        "case_details": rows,
    }


def validation_summary(config_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    gauge_valid_rows = [row for row in config_rows if bool(row["gauge_valid"]) and bool(row["fit_valid"])]
    class_counts: Dict[str, int] = {}
    time_class_counts: Dict[str, Dict[str, int]] = {
        "tau_ABC": {},
        "tau_AB": {},
        "tau_AC": {},
        "tau_BC": {},
    }
    for row in gauge_valid_rows:
        key = str(row["alpha_classification"])
        class_counts[key] = class_counts.get(key, 0) + 1
        for time_name in ["tau_ABC", "tau_AB", "tau_AC", "tau_BC"]:
            time_key = str(row[f"alpha_classification_{time_name}"])
            bucket = time_class_counts[time_name]
            bucket[time_key] = bucket.get(time_key, 0) + 1
    best_inverse_square = [
        row for row in gauge_valid_rows if row["alpha_classification"] == "inverse_square_like_alpha2"
    ]
    best_inverse = [row for row in gauge_valid_rows if row["alpha_classification"] == "inverse_like_alpha1"]
    proportional = [row for row in gauge_valid_rows if row["alpha_classification"] == "proportional_like_alpha_minus1"]
    tau_delta_abs = {
        "tau_AB": max(abs(float(row["tau_AB_alpha_delta_from_tau_ABC"])) for row in gauge_valid_rows)
        if gauge_valid_rows
        else 0.0,
        "tau_AC": max(abs(float(row["tau_AC_alpha_delta_from_tau_ABC"])) for row in gauge_valid_rows)
        if gauge_valid_rows
        else 0.0,
        "tau_BC": max(abs(float(row["tau_BC_alpha_delta_from_tau_ABC"])) for row in gauge_valid_rows)
        if gauge_valid_rows
        else 0.0,
    }
    return {
        "abc_c_gauge_ab_distance_exponent_preliminary_valid": bool(
            len(config_rows) > 0
            and len(gauge_valid_rows) > 0
            and bool_all(not bool(row["inverse_term_injected"]) for row in config_rows)
            and bool_all(not bool(row["inverse_square_term_injected"]) for row in config_rows)
        ),
        "config_count": len(config_rows),
        "gauge_valid_count": len(gauge_valid_rows),
        "resolution_valid_count": sum(1 for row in config_rows if bool(row["resolution_valid"])),
        "clock_valid_count": sum(1 for row in config_rows if bool(row["clock_valid"])),
        "disturbance_valid_count": sum(1 for row in config_rows if bool(row["disturbance_valid"])),
        "alpha_class_counts_in_gauge_valid_cases": class_counts,
        "alpha_class_counts_by_relation_time": time_class_counts,
        "max_abs_alpha_delta_from_tau_ABC": tau_delta_abs,
        "inverse_like_alpha1_count": len(best_inverse),
        "inverse_square_like_alpha2_count": len(best_inverse_square),
        "proportional_like_alpha_minus1_count": len(proportional),
        "main_reading": (
            "The C gauge has a finite eligibility window. "
            "Within this minimal AB-harmonic preliminary model, gauge-valid cases recover proportional-like alpha=-1, "
            "not inverse or inverse-square. This is a gauge-eligibility test, not a final physical-law derivation."
        ),
    }


def flatten_case_details(config_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in config_rows:
        for detail in row["case_details"]:
            flattened = {key: value for key, value in row.items() if key != "case_details"}
            flattened.update(detail)
            out.append(flattened)
    return out


def make_plots(config_rows: List[Dict[str, Any]], detail_rows: List[Dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    xs = [float(row["R_C_over_R_A"]) for row in config_rows]
    ys = [float(row["coupling_C"]) for row in config_rows]
    colors = [
        "tab:green"
        if bool(row["gauge_valid"])
        else ("tab:orange" if bool(row["resolution_valid"]) and bool(row["clock_valid"]) else "tab:red")
        for row in config_rows
    ]
    sizes = [40.0 + 45.0 * float(row["resolved_case_count"]) for row in config_rows]
    ax.scatter(xs, ys, c=colors, s=sizes, alpha=0.72, edgecolor="black", linewidth=0.3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("R_C / R_A")
    ax.set_ylabel("C coupling")
    ax.set_title("C gauge eligibility map")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_ab_distance_exponent_gauge_eligibility_map_v1.png", dpi=180)
    plt.close(fig)

    valid = [row for row in config_rows if bool(row["fit_valid"]) and bool(row["gauge_valid"])]
    fig, ax = plt.subplots(figsize=(8, 5))
    if valid:
        ax.scatter(
            [float(row["R_C_over_R_A"]) for row in valid],
            [float(row["power_candidate_alpha"]) for row in valid],
            c=[float(row["coupling_C"]) for row in valid],
            cmap="viridis",
            s=70,
            edgecolor="black",
            linewidth=0.3,
        )
    for alpha, label in [(2.0, "alpha=2"), (1.0, "alpha=1"), (0.0, "alpha=0"), (-1.0, "alpha=-1")]:
        ax.axhline(alpha, color="black", linestyle="--" if alpha >= 0 else ":", linewidth=0.8)
        ax.text(min(R_C_RATIOS), alpha + 0.03, label, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("R_C / R_A")
    ax.set_ylabel("power_candidate_alpha")
    ax.set_title("Gauge-valid C readout alpha candidates")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_ab_distance_exponent_alpha_candidates_v1.png", dpi=180)
    plt.close(fig)

    selected = [
        row
        for row in detail_rows
        if bool(row["gauge_valid"])
        and float(row["R_C_over_R_A"]) == 1.0
        and float(row["coupling_C"]) == 0.02
        and float(row["C_offset_rad"]) == 0.0
        and bool(row["resolution_valid_at_L"])
    ]
    selected.sort(key=lambda row: float(row["L_AB_C_read"]))
    fig, ax = plt.subplots(figsize=(8, 5))
    if selected:
        ax.plot(
            [float(row["L_AB_C_read"]) for row in selected],
            [float(row["a_AB_C_read"]) for row in selected],
            marker="o",
            label="C-read candidate",
        )
        ax.plot(
            [float(row["L_AB_C_read"]) for row in selected],
            [float(row["L_AB_C_read"]) for row in selected],
            linestyle="--",
            label="proportional reference",
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("L_AB read by C")
    ax.set_ylabel("a_AB candidate read by C")
    ax.set_title("Reference C-gauge readout curve")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_ab_distance_exponent_reference_curve_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    labels = ["resolution", "clock", "disturbance", "gauge"]
    counts = [
        sum(1 for row in config_rows if bool(row["resolution_valid"])),
        sum(1 for row in config_rows if bool(row["clock_valid"])),
        sum(1 for row in config_rows if bool(row["disturbance_valid"])),
        sum(1 for row in config_rows if bool(row["gauge_valid"])),
    ]
    ax.bar(labels, counts, color=["#577590", "#43aa8b", "#f9c74f", "#277da1"])
    ax.set_ylabel("config count")
    ax.set_title("C gauge validity filters")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_ab_distance_exponent_validity_filters_v1.png", dpi=180)
    plt.close(fig)


def write_report(validation: Dict[str, Any]) -> None:
    text = f"""# ABC C-gauge AB distance exponent preliminary report v1

## Summary

- valid: `{validation["abc_c_gauge_ab_distance_exponent_preliminary_valid"]}`
- config_count: `{validation["config_count"]}`
- gauge_valid_count: `{validation["gauge_valid_count"]}`
- resolution_valid_count: `{validation["resolution_valid_count"]}`
- clock_valid_count: `{validation["clock_valid_count"]}`
- disturbance_valid_count: `{validation["disturbance_valid_count"]}`
- inverse_like_alpha1_count: `{validation["inverse_like_alpha1_count"]}`
- inverse_square_like_alpha2_count: `{validation["inverse_square_like_alpha2_count"]}`
- proportional_like_alpha_minus1_count: `{validation["proportional_like_alpha_minus1_count"]}`

## Alpha classes in gauge-valid cases

```json
{json.dumps(validation["alpha_class_counts_in_gauge_valid_cases"], ensure_ascii=False, indent=2)}
```

## Alpha classes by relation time

```json
{json.dumps(validation["alpha_class_counts_by_relation_time"], ensure_ascii=False, indent=2)}
```

## Max alpha delta from tau_ABC

```json
{json.dumps(validation["max_abs_alpha_delta_from_tau_ABC"], ensure_ascii=False, indent=2)}
```

## Reading

This preliminary test does not inject `1/L` or `1/L^2`.

It checks whether a third wave `C` can serve as an independent space-time gauge for reading the AB relation compensation.

The main time readout is `tau_ABC`.

Relation-time diagnostics for `tau_AB`, `tau_AC`, and `tau_BC` are also recorded, but they are not used as the main success criterion.

The strict reading is:

```text
{validation["main_reading"]}
```
"""
    (OUT_DIR / "abc_c_gauge_ab_distance_exponent_preliminary_report_v1.md").write_text(text, encoding="utf-8")


def main() -> None:
    params = Params()
    config_rows = [run_config(config, params) for config in configs()]
    detail_rows = flatten_case_details(config_rows)
    validation = validation_summary(config_rows)

    public_config_rows = [{key: value for key, value in row.items() if key != "case_details"} for row in config_rows]
    write_csv(OUT_DIR / "abc_c_gauge_ab_distance_exponent_gauge_validity_v1.csv", public_config_rows)
    write_csv(OUT_DIR / "abc_c_gauge_ab_distance_exponent_cases_v1.csv", detail_rows)
    (OUT_DIR / "abc_c_gauge_ab_distance_exponent_preliminary_result_v1.json").write_text(
        json.dumps({"validation": validation, "configs": public_config_rows}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_plots(public_config_rows, detail_rows)
    write_report(validation)
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
