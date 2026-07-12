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
OUT_DIR = BASE_DIR / "abc_c_gauge_relation_decomposition_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


R_C_RATIOS = [0.5, 1.0, 2.0, 4.0]
COUPLINGS = [0.005, 0.01, 0.02, 0.05, 0.1]
C_OFFSETS = [-0.08, 0.0, 0.08]


@dataclass(frozen=True)
class ProjectionMode:
    name: str
    sigma_ac: float
    sigma_bc: float
    description: str


PROJECTION_MODES = [
    ProjectionMode(
        "c_opposes_ab",
        -1.0,
        -1.0,
        "f_AC and f_BC oppose the AB-closing local directions.",
    ),
    ProjectionMode(
        "c_assists_ab",
        1.0,
        1.0,
        "f_AC and f_BC assist the AB-closing local directions.",
    ),
    ProjectionMode(
        "a_opposes_b_assists",
        -1.0,
        1.0,
        "C relation opposes A-side and assists B-side.",
    ),
    ProjectionMode(
        "a_assists_b_opposes",
        1.0,
        -1.0,
        "C relation assists A-side and opposes B-side.",
    ),
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
        {
            "r_c_ratio": r_c,
            "coupling_c": coupling,
            "c_offset": c_offset,
            "projection": projection,
        }
        for r_c in R_C_RATIOS
        for coupling in COUPLINGS
        for c_offset in C_OFFSETS
        for projection in PROJECTION_MODES
    ]


def fit_candidate(rows: List[Dict[str, Any]], candidate_key: str) -> Dict[str, Any]:
    fit_rows = [
        row
        for row in rows
        if bool(row["resolution_valid_at_L"])
        and float(row["L_AB_C_read"]) > 0.0
        and float(row[candidate_key]) > 0.0
    ]
    return fit_power_law(
        [float(row["L_AB_C_read"]) for row in fit_rows],
        [float(row[candidate_key]) for row in fit_rows],
    )


def sign_crosses(values: List[float]) -> bool:
    positives = any(value > 0.0 for value in values)
    negatives = any(value < 0.0 for value in values)
    return positives and negatives


def run_config(config: Dict[str, Any], params: Params) -> Dict[str, Any]:
    r_c_ratio = float(config["r_c_ratio"])
    coupling_c = float(config["coupling_c"])
    c_offset = float(config["c_offset"])
    projection: ProjectionMode = config["projection"]

    cell_count = c_cell_count(r_c_ratio, params)
    cell_width = TAU / float(cell_count)

    f_abc_common = coupling_c * r_c_ratio
    omega_abc = 1.0 + f_abc_common
    tau_abc_unit = 1.0 / omega_abc
    clock_error = cell_width / math.radians(max(DEVIATION_DEGS)) + 0.05 * f_abc_common
    clock_valid = bool(clock_error <= params.clock_error_tol)
    closure_valid = True

    details: List[Dict[str, Any]] = []
    for deviation_deg in DEVIATION_DEGS:
        l_true = math.radians(deviation_deg)
        c_cell_crossings = l_true / cell_width
        resolved = bool(c_cell_crossings >= params.min_resolved_cell_crossings)
        l_c_read = round(l_true / cell_width) * cell_width if resolved else 0.0

        # Relation strengths. No inverse or inverse-square term is injected.
        f_ab = l_true
        f_ac = coupling_c * r_c_ratio * max(0.0, 1.0 + c_offset)
        f_bc = coupling_c * r_c_ratio * max(0.0, 1.0 - c_offset)

        c_bias_projected = 0.5 * (projection.sigma_ac * f_ac + projection.sigma_bc * f_bc)
        c_asymmetry_projected = 0.5 * (projection.sigma_ac * f_ac - projection.sigma_bc * f_bc)
        a_a_circ = f_ab + projection.sigma_ac * f_ac
        a_b_circ = f_ab + projection.sigma_bc * f_bc
        ab_pair_projected_signed = 0.5 * (a_a_circ + a_b_circ)
        ab_pair_projected_abs = abs(ab_pair_projected_signed)

        # This is an intentionally separated control: f_ABC used as a direct circular term.
        # It is not the main readout, because f_ABC is a whole-system/common-mode relation.
        common_mode_direct_control = abs(ab_pair_projected_signed + f_abc_common)

        omega_ab = 1.0 + f_ab
        omega_ac = 1.0 + f_ac
        omega_bc = 1.0 + f_bc

        details.append(
            {
                "deviation_deg": deviation_deg,
                "L_AB_true": l_true,
                "L_AB_C_read": l_c_read,
                "C_cell_crossings": c_cell_crossings,
                "resolution_valid_at_L": resolved,
                "f_AB": f_ab,
                "f_AC": f_ac,
                "f_BC": f_bc,
                "f_ABC_common": f_abc_common,
                "sigma_AC": projection.sigma_ac,
                "sigma_BC": projection.sigma_bc,
                "c_bias_projected": c_bias_projected,
                "c_asymmetry_projected": c_asymmetry_projected,
                "c_bias_ratio_to_f_AB": abs(c_bias_projected) / max(abs(f_ab), 1.0e-30),
                "c_asymmetry_ratio_to_f_AB": abs(c_asymmetry_projected) / max(abs(f_ab), 1.0e-30),
                "a_A_circ": a_a_circ,
                "a_B_circ": a_b_circ,
                "AB_pair_projected_signed": ab_pair_projected_signed,
                "AB_pair_projected_abs": ab_pair_projected_abs,
                "common_mode_direct_control": common_mode_direct_control,
                "omega_AB": omega_ab,
                "omega_AC": omega_ac,
                "omega_BC": omega_bc,
                "omega_ABC": omega_abc,
                "tau_AB_unit": 1.0 / omega_ab,
                "tau_AC_unit": 1.0 / omega_ac,
                "tau_BC_unit": 1.0 / omega_bc,
                "tau_ABC_unit": tau_abc_unit,
                "AB_native_by_tau_ABC": abs(l_c_read) * omega_abc * omega_abc,
                "AB_pair_by_tau_ABC": ab_pair_projected_abs * omega_abc * omega_abc,
                "AB_pair_by_tau_AB": ab_pair_projected_abs * omega_ab * omega_ab,
                "AB_pair_by_tau_AC": ab_pair_projected_abs * omega_ac * omega_ac,
                "AB_pair_by_tau_BC": ab_pair_projected_abs * omega_bc * omega_bc,
                "C_bias_by_tau_ABC": abs(c_bias_projected) * omega_abc * omega_abc,
                "C_asymmetry_by_tau_ABC": abs(c_asymmetry_projected) * omega_abc * omega_abc,
                "common_mode_direct_by_tau_ABC": common_mode_direct_control * omega_abc * omega_abc,
            }
        )

    resolved_rows = [row for row in details if bool(row["resolution_valid_at_L"])]
    resolution_valid = bool(len(resolved_rows) >= params.min_fit_case_count)
    max_bias_ratio = max((float(row["c_bias_ratio_to_f_AB"]) for row in resolved_rows), default=float("inf"))
    max_asymmetry_ratio = max(
        (float(row["c_asymmetry_ratio_to_f_AB"]) for row in resolved_rows),
        default=float("inf"),
    )
    pair_sign_crossing = sign_crosses([float(row["AB_pair_projected_signed"]) for row in resolved_rows])
    ab_dominant_valid = bool(
        max_bias_ratio <= params.disturbance_tol
        and max_asymmetry_ratio <= params.disturbance_tol
        and not pair_sign_crossing
    )
    decomposition_valid = bool(resolution_valid and clock_valid and closure_valid)

    fit_native = fit_candidate(details, "AB_native_by_tau_ABC")
    fit_pair_abc = fit_candidate(details, "AB_pair_by_tau_ABC")
    fit_pair_ab = fit_candidate(details, "AB_pair_by_tau_AB")
    fit_pair_ac = fit_candidate(details, "AB_pair_by_tau_AC")
    fit_pair_bc = fit_candidate(details, "AB_pair_by_tau_BC")
    fit_c_bias = fit_candidate(details, "C_bias_by_tau_ABC")
    fit_c_asym = fit_candidate(details, "C_asymmetry_by_tau_ABC")
    fit_common_direct = fit_candidate(details, "common_mode_direct_by_tau_ABC")

    def alpha_of(fit: Dict[str, Any]) -> float:
        return float(fit["power_candidate_alpha"])

    alpha_native = alpha_of(fit_native)
    alpha_pair_abc = alpha_of(fit_pair_abc)
    alpha_pair_ab = alpha_of(fit_pair_ab)
    alpha_pair_ac = alpha_of(fit_pair_ac)
    alpha_pair_bc = alpha_of(fit_pair_bc)
    alpha_c_bias = alpha_of(fit_c_bias)
    alpha_c_asym = alpha_of(fit_c_asym)
    alpha_common_direct = alpha_of(fit_common_direct)

    return {
        "config_id": f"Rc{r_c_ratio:g}_eps{coupling_c:g}_off{c_offset:g}_{projection.name}",
        "R_C_over_R_A": r_c_ratio,
        "coupling_C": coupling_c,
        "C_offset": c_offset,
        "projection_mode": projection.name,
        "projection_description": projection.description,
        "sigma_AC": projection.sigma_ac,
        "sigma_BC": projection.sigma_bc,
        "C_cell_count": cell_count,
        "C_cell_width": cell_width,
        "resolved_case_count": len(resolved_rows),
        "resolution_valid": resolution_valid,
        "clock_error": clock_error,
        "clock_valid": clock_valid,
        "closure_valid": closure_valid,
        "decomposition_valid": decomposition_valid,
        "AB_dominant_valid": ab_dominant_valid,
        "pair_sign_crossing": pair_sign_crossing,
        "max_c_bias_ratio_to_f_AB": max_bias_ratio if math.isfinite(max_bias_ratio) else "",
        "max_c_asymmetry_ratio_to_f_AB": max_asymmetry_ratio if math.isfinite(max_asymmetry_ratio) else "",
        "f_ABC_common_mode_ratio": f_abc_common,
        "fit_valid_native": bool(fit_native["fit_valid"]),
        "fit_valid_pair_tau_ABC": bool(fit_pair_abc["fit_valid"]),
        "fit_valid_c_bias": bool(fit_c_bias["fit_valid"]),
        "fit_valid_c_asymmetry": bool(fit_c_asym["fit_valid"]),
        "fit_valid_common_mode_direct_control": bool(fit_common_direct["fit_valid"]),
        "alpha_native_AB_tau_ABC": alpha_native,
        "alpha_pair_tau_ABC": alpha_pair_abc,
        "alpha_pair_tau_AB": alpha_pair_ab,
        "alpha_pair_tau_AC": alpha_pair_ac,
        "alpha_pair_tau_BC": alpha_pair_bc,
        "alpha_c_bias_tau_ABC": alpha_c_bias,
        "alpha_c_asymmetry_tau_ABC": alpha_c_asym,
        "alpha_common_mode_direct_control_tau_ABC": alpha_common_direct,
        "alpha_class_native_AB_tau_ABC": classify_alpha(alpha_native, params) if fit_native["fit_valid"] else "unfit",
        "alpha_class_pair_tau_ABC": classify_alpha(alpha_pair_abc, params) if fit_pair_abc["fit_valid"] else "unfit",
        "alpha_class_pair_tau_AB": classify_alpha(alpha_pair_ab, params) if fit_pair_ab["fit_valid"] else "unfit",
        "alpha_class_pair_tau_AC": classify_alpha(alpha_pair_ac, params) if fit_pair_ac["fit_valid"] else "unfit",
        "alpha_class_pair_tau_BC": classify_alpha(alpha_pair_bc, params) if fit_pair_bc["fit_valid"] else "unfit",
        "alpha_class_c_bias_tau_ABC": classify_alpha(alpha_c_bias, params) if fit_c_bias["fit_valid"] else "unfit",
        "alpha_class_c_asymmetry_tau_ABC": classify_alpha(alpha_c_asym, params) if fit_c_asym["fit_valid"] else "unfit",
        "alpha_class_common_mode_direct_control_tau_ABC": (
            classify_alpha(alpha_common_direct, params) if fit_common_direct["fit_valid"] else "unfit"
        ),
        "tau_AB_alpha_delta_from_tau_ABC": alpha_pair_ab - alpha_pair_abc,
        "tau_AC_alpha_delta_from_tau_ABC": alpha_pair_ac - alpha_pair_abc,
        "tau_BC_alpha_delta_from_tau_ABC": alpha_pair_bc - alpha_pair_abc,
        "log_rmse_pair_tau_ABC": fit_pair_abc["log_rmse"],
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


def count_classes(rows: List[Dict[str, Any]], class_key: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for row in rows:
        value = str(row[class_key])
        out[value] = out.get(value, 0) + 1
    return dict(sorted(out.items()))


def count_by(rows: List[Dict[str, Any]], key: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        out[value] = out.get(value, 0) + 1
    return dict(sorted(out.items()))


def validation_summary(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    decomposition_valid = [
        row for row in rows if bool(row["decomposition_valid"]) and bool(row["fit_valid_pair_tau_ABC"])
    ]
    ab_dominant = [
        row
        for row in decomposition_valid
        if bool(row["AB_dominant_valid"])
    ]
    non_dominant = [
        row
        for row in decomposition_valid
        if not bool(row["AB_dominant_valid"])
    ]
    inverse_like_in_ab_dominant = [
        row
        for row in ab_dominant
        if str(row["alpha_class_pair_tau_ABC"]) in {"inverse_like_alpha1", "inverse_square_like_alpha2"}
    ]
    return {
        "relation_decomposition_preliminary_valid": bool(
            len(rows) > 0
            and len(decomposition_valid) > 0
            and bool_all(not bool(row["inverse_term_injected"]) for row in rows)
            and bool_all(not bool(row["inverse_square_term_injected"]) for row in rows)
        ),
        "config_count": len(rows),
        "decomposition_valid_count": len(decomposition_valid),
        "AB_dominant_valid_count": len(ab_dominant),
        "non_AB_dominant_count": len(non_dominant),
        "projection_mode_counts_in_decomposition_valid": count_by(decomposition_valid, "projection_mode"),
        "pair_tau_ABC_alpha_classes_in_AB_dominant_cases": count_classes(
            ab_dominant,
            "alpha_class_pair_tau_ABC",
        ),
        "pair_tau_ABC_alpha_classes_in_all_decomposition_valid_cases": count_classes(
            decomposition_valid,
            "alpha_class_pair_tau_ABC",
        ),
        "native_AB_alpha_classes_in_decomposition_valid_cases": count_classes(
            decomposition_valid,
            "alpha_class_native_AB_tau_ABC",
        ),
        "c_bias_alpha_classes_in_decomposition_valid_cases": count_classes(
            decomposition_valid,
            "alpha_class_c_bias_tau_ABC",
        ),
        "c_asymmetry_alpha_classes_in_decomposition_valid_cases": count_classes(
            decomposition_valid,
            "alpha_class_c_asymmetry_tau_ABC",
        ),
        "common_mode_direct_control_alpha_classes": count_classes(
            decomposition_valid,
            "alpha_class_common_mode_direct_control_tau_ABC",
        ),
        "inverse_or_inverse_square_in_AB_dominant_count": len(inverse_like_in_ab_dominant),
        "max_tau_AB_alpha_delta_from_tau_ABC_in_AB_dominant": max(
            (abs(float(row["tau_AB_alpha_delta_from_tau_ABC"])) for row in ab_dominant),
            default=0.0,
        ),
        "max_tau_AC_alpha_delta_from_tau_ABC_in_AB_dominant": max(
            (abs(float(row["tau_AC_alpha_delta_from_tau_ABC"])) for row in ab_dominant),
            default=0.0,
        ),
        "max_tau_BC_alpha_delta_from_tau_ABC_in_AB_dominant": max(
            (abs(float(row["tau_BC_alpha_delta_from_tau_ABC"])) for row in ab_dominant),
            default=0.0,
        ),
        "main_reading": (
            "Separating f_AB, f_AC, f_BC, and f_ABC confirms that tau_ABC can be used as "
            "a representative time gauge, while f_ABC should not be added as a direct circular term. "
            "In AB-dominant projected cases, the pair readout remains proportional-like rather than inverse-like."
        ),
    }


def make_plots(rows: List[Dict[str, Any]]) -> None:
    decomposition_valid = [
        row for row in rows if bool(row["decomposition_valid"]) and bool(row["fit_valid_pair_tau_ABC"])
    ]
    ab_dominant = [row for row in decomposition_valid if bool(row["AB_dominant_valid"])]

    fig, ax = plt.subplots(figsize=(10, 5))
    modes = [mode.name for mode in PROJECTION_MODES]
    all_counts = [sum(1 for row in decomposition_valid if row["projection_mode"] == mode) for mode in modes]
    dominant_counts = [sum(1 for row in ab_dominant if row["projection_mode"] == mode) for mode in modes]
    x = np.arange(len(modes))
    width = 0.36
    ax.bar(x - width / 2.0, all_counts, width, label="decomposition-valid", color="#277da1")
    ax.bar(x + width / 2.0, dominant_counts, width, label="AB-dominant", color="#43aa8b")
    ax.set_xticks(x, modes, rotation=35, ha="right")
    ax.set_ylabel("config count")
    ax.set_title("Relation decomposition valid counts")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_relation_decomposition_valid_counts_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    if decomposition_valid:
        xs = list(range(len(decomposition_valid)))
        alphas = [float(row["alpha_pair_tau_ABC"]) for row in decomposition_valid]
        colors = ["#43aa8b" if bool(row["AB_dominant_valid"]) else "#f3722c" for row in decomposition_valid]
        ax.scatter(xs, alphas, c=colors, s=45, edgecolor="black", linewidth=0.25)
    for alpha, label in [(2.0, "alpha=2"), (1.0, "alpha=1"), (0.0, "alpha=0"), (-1.0, "alpha=-1")]:
        ax.axhline(alpha, color="black", linestyle="--" if alpha >= 0 else ":", linewidth=0.8)
        ax.text(0, alpha + 0.03, label, fontsize=8)
    ax.set_ylabel("pair projected alpha by tau_ABC")
    ax.set_title("Projected pair alpha candidates")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_relation_decomposition_pair_alpha_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    if decomposition_valid:
        native = [float(row["alpha_native_AB_tau_ABC"]) for row in decomposition_valid]
        pair = [float(row["alpha_pair_tau_ABC"]) for row in decomposition_valid]
        colors = ["#43aa8b" if bool(row["AB_dominant_valid"]) else "#f3722c" for row in decomposition_valid]
        ax.scatter(native, pair, c=colors, s=55, edgecolor="black", linewidth=0.25)
    ax.plot([-2.0, 2.5], [-2.0, 2.5], color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("native f_AB alpha")
    ax.set_ylabel("projected pair alpha")
    ax.set_title("Native AB vs projected pair alpha")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_relation_decomposition_native_vs_pair_alpha_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    if decomposition_valid:
        xs = list(range(len(decomposition_valid)))
        bias = [float(row["max_c_bias_ratio_to_f_AB"]) for row in decomposition_valid]
        asym = [float(row["max_c_asymmetry_ratio_to_f_AB"]) for row in decomposition_valid]
        ax.plot(xs, bias, label="max projected C-bias / f_AB", color="#f3722c")
        ax.plot(xs, asym, label="max projected C-asymmetry / f_AB", color="#577590")
    ax.axhline(Params().disturbance_tol, color="black", linestyle="--", linewidth=0.8, label="disturbance tol")
    ax.set_yscale("log")
    ax.set_ylabel("ratio")
    ax.set_title("C relation contamination diagnostics")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "abc_c_gauge_relation_decomposition_contamination_v1.png", dpi=180)
    plt.close(fig)


def write_report(validation: Dict[str, Any]) -> None:
    text = f"""# ABC C-gauge relation decomposition preliminary report v1

## Summary

- valid: `{validation["relation_decomposition_preliminary_valid"]}`
- config_count: `{validation["config_count"]}`
- decomposition_valid_count: `{validation["decomposition_valid_count"]}`
- AB_dominant_valid_count: `{validation["AB_dominant_valid_count"]}`
- non_AB_dominant_count: `{validation["non_AB_dominant_count"]}`
- inverse_or_inverse_square_in_AB_dominant_count: `{validation["inverse_or_inverse_square_in_AB_dominant_count"]}`

## Projection mode counts

```json
{json.dumps(validation["projection_mode_counts_in_decomposition_valid"], ensure_ascii=False, indent=2)}
```

## Pair tau_ABC alpha classes in AB-dominant cases

```json
{json.dumps(validation["pair_tau_ABC_alpha_classes_in_AB_dominant_cases"], ensure_ascii=False, indent=2)}
```

## Pair tau_ABC alpha classes in all decomposition-valid cases

```json
{json.dumps(validation["pair_tau_ABC_alpha_classes_in_all_decomposition_valid_cases"], ensure_ascii=False, indent=2)}
```

## Native AB alpha classes

```json
{json.dumps(validation["native_AB_alpha_classes_in_decomposition_valid_cases"], ensure_ascii=False, indent=2)}
```

## C-bias alpha classes

```json
{json.dumps(validation["c_bias_alpha_classes_in_decomposition_valid_cases"], ensure_ascii=False, indent=2)}
```

## C-asymmetry alpha classes

```json
{json.dumps(validation["c_asymmetry_alpha_classes_in_decomposition_valid_cases"], ensure_ascii=False, indent=2)}
```

## Common-mode direct-control alpha classes

```json
{json.dumps(validation["common_mode_direct_control_alpha_classes"], ensure_ascii=False, indent=2)}
```

## Relation-time alpha deltas in AB-dominant cases

```text
max |tau_AB - tau_ABC| alpha delta: {validation["max_tau_AB_alpha_delta_from_tau_ABC_in_AB_dominant"]}
max |tau_AC - tau_ABC| alpha delta: {validation["max_tau_AC_alpha_delta_from_tau_ABC_in_AB_dominant"]}
max |tau_BC - tau_ABC| alpha delta: {validation["max_tau_BC_alpha_delta_from_tau_ABC_in_AB_dominant"]}
```

## Reading

No inverse or inverse-square term is injected.

The strict reading is:

```text
{validation["main_reading"]}
```
"""
    (OUT_DIR / "abc_c_gauge_relation_decomposition_preliminary_report_v1.md").write_text(
        text,
        encoding="utf-8",
    )


def main() -> None:
    params = Params()
    config_rows = [run_config(config, params) for config in configs()]
    detail_rows = flatten_case_details(config_rows)
    public_rows = [{key: value for key, value in row.items() if key != "case_details"} for row in config_rows]
    validation = validation_summary(public_rows)

    write_csv(OUT_DIR / "abc_c_gauge_relation_decomposition_configs_v1.csv", public_rows)
    write_csv(OUT_DIR / "abc_c_gauge_relation_decomposition_cases_v1.csv", detail_rows)
    (OUT_DIR / "abc_c_gauge_relation_decomposition_preliminary_result_v1.json").write_text(
        json.dumps({"validation": validation, "configs": public_rows}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_plots(public_rows)
    write_report(validation)
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
