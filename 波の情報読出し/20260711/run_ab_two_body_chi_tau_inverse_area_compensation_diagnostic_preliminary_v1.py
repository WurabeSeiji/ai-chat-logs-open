from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR / "ab_two_body_c1_internal_calibration_chi_tau_area_sweep_preliminary_result_v1"
OUT_DIR = BASE_DIR / "ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SERIES_PATH = INPUT_DIR / "ab_two_body_c1_internal_calibration_chi_tau_area_sweep_series_v1.csv"
SUMMARY_PATH = INPUT_DIR / "ab_two_body_c1_internal_calibration_chi_tau_area_sweep_case_summary_v1.csv"

AREA_TOL = 1.0e-12
ALPHA_TOL = 0.12


def bool_all(values: Iterable[bool]) -> bool:
    return all(bool(value) for value in values)


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in keys})


def as_float(value: str, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    return float(value)


def group_key(row: Dict[str, str]) -> Tuple[str, str, str, str]:
    return (row["case_id"], row["protocol"], row["tau_mode"], row["readout_mode"])


def fit_power(rows: List[Dict[str, Any]], candidate: str) -> Dict[str, Any]:
    selected = [row for row in rows if float(row.get(candidate, 0.0)) > 1.0e-30]
    selected.sort(key=lambda row: float(row["initial_deviation_rad"]))
    if len(selected) < 3:
        return {
            "candidate": candidate,
            "fit_valid": False,
            "loglog_slope": 0.0,
            "power_candidate_alpha": 0.0,
            "case_count": len(selected),
        }
    x = np.array([float(row["initial_deviation_rad"]) for row in selected], dtype=float)
    y = np.array([float(row[candidate]) for row in selected], dtype=float)
    slope, intercept = np.polyfit(np.log(x), np.log(y), 1)
    predicted = slope * np.log(x) + intercept
    residual = np.log(y) - predicted
    rmse = float(np.sqrt(np.mean(residual * residual)))
    return {
        "candidate": candidate,
        "fit_valid": True,
        "loglog_slope": float(slope),
        "power_candidate_alpha": -float(slope),
        "loglog_intercept": float(intercept),
        "log_rmse": rmse,
        "case_count": len(selected),
        "min_value": float(np.min(y)),
        "max_value": float(np.max(y)),
    }


def build_case_diagnostics(series_rows: List[Dict[str, str]], summary_rows: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    grouped: Dict[Tuple[str, str, str, str], List[Dict[str, str]]] = {}
    for row in series_rows:
        grouped.setdefault(group_key(row), []).append(row)

    diagnostics: List[Dict[str, Any]] = []
    for summary in summary_rows:
        key = group_key(summary)
        rows = grouped[key]
        max_area = as_float(summary["max_abs_A_chi_tau"])
        area_valid = bool(max_area > AREA_TOL and int(summary["rank_chi_tau"]) >= 2)
        max_f = max(abs(as_float(row["f_AB"])) for row in rows)
        max_q_raw = max(abs(as_float(row["Q_raw"])) for row in rows)
        max_closure_relaxation = max(abs(as_float(row["closure_relaxation"])) for row in rows)
        max_envelope = max(abs(as_float(row["envelope_V_AB"])) for row in rows)
        max_v = max(abs(as_float(row["V_AB"])) for row in rows)
        max_epsilon_c = as_float(summary["max_epsilon_c_abs"])

        if area_valid:
            inv_area = 1.0 / max_area
            f_over_area = max_f / max_area
            q_over_area = max_q_raw / max_area
            closure_over_area = max_closure_relaxation / max_area
            envelope_over_area = max_envelope / max_area
        else:
            inv_area = 0.0
            f_over_area = 0.0
            q_over_area = 0.0
            closure_over_area = 0.0
            envelope_over_area = 0.0

        diagnostics.append(
            {
                "case_id": summary["case_id"],
                "protocol": summary["protocol"],
                "tau_mode": summary["tau_mode"],
                "readout_mode": summary["readout_mode"],
                "initial_deviation_deg": as_float(summary["initial_deviation_deg"]),
                "initial_deviation_rad": math.radians(as_float(summary["initial_deviation_deg"])),
                "rank_chi_tau": int(summary["rank_chi_tau"]),
                "area_valid": area_valid,
                "max_abs_A_chi_tau": max_area,
                "inv_A_chi_tau_constructed_control": inv_area,
                "native_max_f_AB": max_f,
                "native_max_Q_raw": max_q_raw,
                "native_max_closure_relaxation": max_closure_relaxation,
                "native_max_envelope_V_AB": max_envelope,
                "native_max_V_AB": max_v,
                "native_max_epsilon_c": max_epsilon_c,
                "derived_f_AB_over_A": f_over_area,
                "derived_Q_raw_over_A": q_over_area,
                "derived_closure_relaxation_over_A": closure_over_area,
                "derived_envelope_over_A": envelope_over_area,
                "tau_is_step_used": summary["tau_is_step_used"],
                "external_c_used": summary["external_c_used"],
                "f_A_or_f_B_used": summary["f_A_or_f_B_used"],
            }
        )
    return diagnostics


def candidate_fits(diagnostics: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    candidates = [
        ("native_max_f_AB", "native"),
        ("native_max_Q_raw", "native"),
        ("native_max_closure_relaxation", "native"),
        ("native_max_envelope_V_AB", "native"),
        ("native_max_V_AB", "native"),
        ("native_max_epsilon_c", "native"),
        ("max_abs_A_chi_tau", "area"),
        ("inv_A_chi_tau_constructed_control", "constructed_reciprocal_control"),
        ("derived_f_AB_over_A", "derived_ratio"),
        ("derived_Q_raw_over_A", "derived_ratio"),
        ("derived_closure_relaxation_over_A", "derived_ratio"),
        ("derived_envelope_over_A", "derived_ratio"),
    ]
    fit_rows: List[Dict[str, Any]] = []
    for protocol in sorted({row["protocol"] for row in diagnostics}):
        for tau_mode in sorted({row["tau_mode"] for row in diagnostics}):
            for readout_mode in sorted({row["readout_mode"] for row in diagnostics}):
                selected = [
                    row
                    for row in diagnostics
                    if row["protocol"] == protocol
                    and row["tau_mode"] == tau_mode
                    and row["readout_mode"] == readout_mode
                    and bool(row["area_valid"])
                ]
                for candidate, candidate_kind in candidates:
                    fit = fit_power(selected, candidate)
                    fit.update(
                        {
                            "protocol": protocol,
                            "tau_mode": tau_mode,
                            "readout_mode": readout_mode,
                            "candidate_kind": candidate_kind,
                            "alpha_near_positive_2": bool(
                                fit["fit_valid"] and abs(float(fit["power_candidate_alpha"]) - 2.0) <= ALPHA_TOL
                            ),
                        }
                    )
                    fit_rows.append(fit)
    return fit_rows


def validation_summary(diagnostics: List[Dict[str, Any]], fits: List[Dict[str, Any]]) -> Dict[str, Any]:
    native_fits = [row for row in fits if row["candidate_kind"] == "native" and bool(row["fit_valid"])]
    constructed = [
        row
        for row in fits
        if row["candidate"] == "inv_A_chi_tau_constructed_control"
        and row["tau_mode"] == "tau_independent_c1"
        and row["readout_mode"] == "readout_off"
        and bool(row["fit_valid"])
    ]
    derived = [row for row in fits if row["candidate_kind"] == "derived_ratio" and bool(row["fit_valid"])]
    native_positive2 = [row for row in native_fits if bool(row["alpha_near_positive_2"])]
    constructed_positive2 = [row for row in constructed if bool(row["alpha_near_positive_2"])]
    best_native_alpha = max(native_fits, key=lambda row: float(row["power_candidate_alpha"])) if native_fits else {}
    c1_area_rows = [
        row
        for row in diagnostics
        if row["tau_mode"] == "tau_independent_c1" and row["readout_mode"] == "readout_off" and bool(row["area_valid"])
    ]
    return {
        "inverse_area_compensation_diagnostic_preliminary_valid": bool(
            len(constructed_positive2) >= 1
            and len(native_positive2) == 0
            and bool_all(str(row["tau_is_step_used"]).lower() == "false" for row in diagnostics)
            and bool_all(str(row["external_c_used"]).lower() == "false" for row in diagnostics)
            and bool_all(str(row["f_A_or_f_B_used"]).lower() == "false" for row in diagnostics)
        ),
        "diagnostic_case_count": len(diagnostics),
        "area_valid_case_count": sum(1 for row in diagnostics if bool(row["area_valid"])),
        "fit_count": len(fits),
        "native_fit_count": len(native_fits),
        "derived_fit_count": len(derived),
        "native_positive2_count": len(native_positive2),
        "constructed_reciprocal_positive2_count": len(constructed_positive2),
        "best_native_alpha_candidate": best_native_alpha,
        "c1_readout_off_area_min": min(float(row["max_abs_A_chi_tau"]) for row in c1_area_rows),
        "c1_readout_off_area_max": max(float(row["max_abs_A_chi_tau"]) for row in c1_area_rows),
        "main_reading": (
            "The reciprocal diagnostic 1/A gives alpha≈2 by construction, "
            "but no native readout in this preliminary dataset naturally shows alpha≈2."
        ),
    }


def make_plots(diagnostics: List[Dict[str, Any]], fits: List[Dict[str, Any]]) -> None:
    selected = [
        row
        for row in diagnostics
        if row["protocol"] == "Protocol_B"
        and row["tau_mode"] == "tau_independent_c1"
        and row["readout_mode"] == "readout_off"
        and bool(row["area_valid"])
    ]
    selected.sort(key=lambda row: float(row["initial_deviation_deg"]))
    x = [float(row["initial_deviation_rad"]) for row in selected]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, [float(row["max_abs_A_chi_tau"]) for row in selected], marker="o", label="A_chi_tau")
    ax.plot(
        x,
        [float(row["inv_A_chi_tau_constructed_control"]) for row in selected],
        marker="o",
        label="1/A_chi_tau constructed",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("initial deviation rad")
    ax.set_ylabel("diagnostic value")
    ax.set_title("area and constructed reciprocal")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_chi_tau_inverse_area_constructed_control_v1.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    native_candidates = [
        ("native_max_f_AB", "f_AB"),
        ("native_max_Q_raw", "Q_raw"),
        ("native_max_envelope_V_AB", "envelope"),
        ("native_max_V_AB", "V_AB"),
    ]
    for key, label in native_candidates:
        values = [float(row[key]) for row in selected]
        if max(values) > 1.0e-30:
            ax.plot(x, values, marker="o", label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("initial deviation rad")
    ax.set_ylabel("native readout value")
    ax.set_title("native readouts on c1 chi-tau surface")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_chi_tau_inverse_area_native_candidates_v1.png", dpi=180)
    plt.close(fig)

    fit_selected = [
        row
        for row in fits
        if row["protocol"] == "Protocol_B"
        and row["tau_mode"] == "tau_independent_c1"
        and row["readout_mode"] in {"readout_off", "readout_normal"}
        and row["candidate"]
        in {
            "native_max_f_AB",
            "native_max_Q_raw",
            "native_max_envelope_V_AB",
            "max_abs_A_chi_tau",
            "inv_A_chi_tau_constructed_control",
            "derived_f_AB_over_A",
            "derived_Q_raw_over_A",
        }
        and bool(row["fit_valid"])
    ]
    labels = [f'{row["candidate"]}\n{row["readout_mode"]}' for row in fit_selected]
    alphas = [float(row["power_candidate_alpha"]) for row in fit_selected]
    colors = [
        "tab:red" if row["candidate_kind"] == "constructed_reciprocal_control" else "tab:blue"
        for row in fit_selected
    ]
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(range(len(alphas)), alphas, color=colors)
    ax.axhline(2.0, color="black", linewidth=1.0, linestyle="--", label="alpha=2")
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xticks(range(len(labels)), labels, rotation=75, ha="right", fontsize=8)
    ax.set_ylabel("power_candidate_alpha")
    ax.set_title("native vs derived inverse-area diagnostic exponents")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "ab_two_body_chi_tau_inverse_area_alpha_comparison_v1.png", dpi=180)
    plt.close(fig)


def write_report(validation: Dict[str, Any]) -> None:
    best = validation["best_native_alpha_candidate"]
    text = f"""# AB二体 chi-tau 面積逆数補償診断予備実験レポート v1

## Summary

- valid: `{validation["inverse_area_compensation_diagnostic_preliminary_valid"]}`
- diagnostic_case_count: `{validation["diagnostic_case_count"]}`
- area_valid_case_count: `{validation["area_valid_case_count"]}`
- native_positive2_count: `{validation["native_positive2_count"]}`
- constructed_reciprocal_positive2_count: `{validation["constructed_reciprocal_positive2_count"]}`
- c1_readout_off_area_min: `{validation["c1_readout_off_area_min"]:.16e}`
- c1_readout_off_area_max: `{validation["c1_readout_off_area_max"]:.16e}`

## Best native alpha candidate

```json
{json.dumps(best, ensure_ascii=False, indent=2)}
```

## Reading

This diagnostic deliberately separates native readouts from constructed reciprocal controls.

`1/A_chi_tau` gives an alpha near `+2` by construction. This is not counted as a native discovery.

In the present dataset, no native readout candidate naturally gives alpha near `+2`.

Therefore the strict reading is:

```text
chi-tau area exists.
1 / chi-tau area has inverse-square scaling by construction.
native inverse-area compensation has not yet been detected.
```
"""
    (OUT_DIR / "ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_report_v1.md").write_text(
        text,
        encoding="utf-8",
    )


def main() -> None:
    series_rows = read_csv(SERIES_PATH)
    summary_rows = read_csv(SUMMARY_PATH)
    diagnostics = build_case_diagnostics(series_rows, summary_rows)
    fits = candidate_fits(diagnostics)
    validation = validation_summary(diagnostics, fits)

    write_csv(OUT_DIR / "ab_two_body_chi_tau_inverse_area_compensation_diagnostic_cases_v1.csv", diagnostics)
    write_csv(OUT_DIR / "ab_two_body_chi_tau_inverse_area_compensation_diagnostic_fits_v1.csv", fits)
    result = {
        "input_series": str(SERIES_PATH.relative_to(BASE_DIR)),
        "input_summary": str(SUMMARY_PATH.relative_to(BASE_DIR)),
        "validation": validation,
        "diagnostics": diagnostics,
        "fits": fits,
        "notes": [
            "Constructed reciprocal controls are not native readouts.",
            "A_chi_tau=0 or rank<2 cases are excluded from reciprocal diagnostics.",
            "Positive alpha near 2 is accepted only as native evidence if it appears in native candidates.",
        ],
    }
    (OUT_DIR / "ab_two_body_chi_tau_inverse_area_compensation_diagnostic_preliminary_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    make_plots(diagnostics, fits)
    write_report(validation)
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
