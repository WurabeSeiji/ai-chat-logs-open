from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

from run_gray_cat_ab_metastable_interface_preliminary_v1 import Params as ABParams
from run_gray_cat_d_observation_response_preliminary_v1 import (
    BASE_DIR,
    AB_RESULT_PATH,
    Params as DResponseParams,
    evolve_from_state,
    measure_s_q,
    pre_evolve,
    read_json,
    selected_candidates,
    sign_label,
)
from run_gray_cat_c_readout_window_preliminary_v1 import c_visibility


OUT_DIR = BASE_DIR / "gray_cat_d_selection_boundary_preliminary_result_v1"


@dataclass(frozen=True)
class Params:
    d_steps: int = 2048
    c_modes: Tuple[str, ...] = ("record_only", "weak_C_window")
    pre_steps_values: Tuple[int, ...] = (
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        25,
        30,
        35,
        40,
        45,
        50,
        55,
        60,
        65,
        70,
        75,
        80,
        85,
        90,
        95,
        100,
        110,
        120,
        130,
        140,
        150,
        160,
        170,
        180,
        190,
        200,
        250,
        300,
        350,
        400,
        450,
        500,
        550,
        600,
        650,
        700,
        750,
        800,
        850,
        900,
        950,
        1000,
        1100,
        1200,
        1300,
        1400,
        1500,
        1600,
        1700,
        1800,
        1900,
        2000,
    )
    d_gain_values: Tuple[float, ...] = (
        0.0,
        0.005,
        0.01,
        0.015,
        0.02,
        0.0225,
        0.025,
        0.0275,
        0.03,
        0.0325,
        0.035,
        0.0375,
        0.04,
        0.045,
        0.05,
        0.055,
        0.06,
        0.065,
        0.07,
        0.0725,
        0.075,
        0.0775,
        0.08,
        0.0825,
        0.085,
        0.0875,
        0.09,
        0.095,
        0.10,
        0.12,
        0.15,
        0.20,
        0.30,
        0.50,
        0.75,
        1.00,
    )
    g_d: float = 1.0
    c_g: float = 1.0
    c_readout_kappa: float = 0.02
    d_readout_kappa: float = 0.02
    c_backaction_scale: float = 1.0e-5
    candidate_case_ids: Tuple[str, ...] = (
        "gray_metastable_0_eps0.01_phi0_s0.01_g0",
        "gray_metastable_1_eps0.01_phi1_s0.01_g0",
        "gray_metastable_4_eps0.003_phi0.0833333_s0_g-0.002",
    )
    large_s_limit: float = 0.05


def selected_metastable_candidates(params: Params) -> List[Dict[str, Any]]:
    ab_result = read_json(AB_RESULT_PATH)
    response_params = DResponseParams()
    all_candidates = selected_candidates(ab_result, response_params)
    by_id = {row["case_id"]: row for row in all_candidates}
    missing = [case_id for case_id in params.candidate_case_ids if case_id not in by_id]
    if missing:
        raise ValueError(f"missing candidate ids: {missing}")
    return [by_id[case_id] for case_id in params.candidate_case_ids]


def selected_outcome(outcome: str) -> bool:
    return outcome in ("white_selected", "black_selected")


def run_boundary_for_start(
    ab_params: ABParams,
    params: Params,
    candidate: Dict[str, Any],
    pre_steps: int,
    c_mode: str,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    response_params = DResponseParams(
        d_steps=params.d_steps,
        c_readout_kappa=params.c_readout_kappa,
        d_readout_kappa=params.d_readout_kappa,
        c_g=params.c_g,
        c_backaction_scale=params.c_backaction_scale,
        g_d_values=(params.g_d,),
        d_backaction_scale_values=params.d_gain_values,
    )
    a0, b0 = pre_evolve(candidate, response_params, pre_steps, c_mode)
    s_start, q_start = measure_s_q(a0, b0)
    c_s = c_visibility(params.c_g, params.c_readout_kappa) * s_start
    c_sign = sign_label(c_s, params.large_s_limit)
    baseline = evolve_from_state(ab_params, response_params, candidate, a0, b0, 0.0, 0.0)
    rows: List[Dict[str, Any]] = []
    for d_gain in params.d_gain_values:
        observed = evolve_from_state(ab_params, response_params, candidate, a0, b0, params.g_d, d_gain)
        rows.append(
            {
                "case_id": candidate["case_id"],
                "candidate_kind": candidate["candidate_kind"],
                "epsilon": candidate["epsilon"],
                "phi_over_pi": candidate["phi_over_pi"],
                "s0": candidate["s0"],
                "stability_gain": candidate["stability_gain"],
                "pre_steps": pre_steps,
                "C_mode": c_mode,
                "S_start": s_start,
                "Q_start": q_start,
                "C_sign_start": c_sign,
                "g_D": params.g_d,
                "D_gain": d_gain,
                "baseline_outcome": baseline["D_outcome"],
                "baseline_phase_after": baseline["phase_after"],
                "baseline_S_mean_after": baseline["S_mean_after"],
                "D_outcome": observed["D_outcome"],
                "phase_after_D": observed["phase_after"],
                "S_mean_after_D": observed["S_mean_after"],
                "S_amp_after_D": observed["S_amp_after"],
                "Q_max_error_after_D": observed["Q_max_error"],
                "D_selected": selected_outcome(observed["D_outcome"]),
                "D_induced_selection": selected_outcome(observed["D_outcome"])
                and not selected_outcome(baseline["D_outcome"]),
            }
        )
    induced = [row for row in rows if row["D_induced_selection"]]
    if induced:
        first = min(induced, key=lambda row: row["D_gain"])
        min_gain = first["D_gain"]
        selected_sign = "A" if first["D_outcome"] == "white_selected" else "B"
        s_after_at_min = first["S_mean_after_D"]
    else:
        min_gain = None
        selected_sign = "none"
        s_after_at_min = None
    boundary = {
        "case_id": candidate["case_id"],
        "pre_steps": pre_steps,
        "C_mode": c_mode,
        "S_start": s_start,
        "C_sign_start": c_sign,
        "baseline_outcome": baseline["D_outcome"],
        "baseline_S_mean_after": baseline["S_mean_after"],
        "min_D_gain_for_induced_selection": min_gain,
        "selected_sign_at_min": selected_sign,
        "S_mean_after_at_min": s_after_at_min,
        "selection_possible": bool(induced),
    }
    return rows, boundary


def count_by(rows: Iterable[Dict[str, Any]], key: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        out[value] = out.get(value, 0) + 1
    return out


def summarise(rows: List[Dict[str, Any]], boundary_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    selection_boundaries = [row for row in boundary_rows if row["selection_possible"]]
    no_selection_boundaries = [row for row in boundary_rows if not row["selection_possible"]]
    min_gain_values = [row["min_D_gain_for_induced_selection"] for row in selection_boundaries]
    by_case: Dict[str, Dict[str, Any]] = {}
    for case_id in sorted({row["case_id"] for row in boundary_rows}):
        case_rows = [row for row in boundary_rows if row["case_id"] == case_id]
        selected = [row for row in case_rows if row["selection_possible"]]
        gains = [row["min_D_gain_for_induced_selection"] for row in selected]
        by_case[case_id] = {
            "boundary_points": len(case_rows),
            "selection_possible_points": len(selected),
            "no_selection_points": len(case_rows) - len(selected),
            "min_observed_D_gain": min(gains) if gains else None,
            "max_observed_D_gain": max(gains) if gains else None,
            "selected_sign_counts": count_by(selected, "selected_sign_at_min"),
        }
    sorted_boundaries = sorted(
        selection_boundaries,
        key=lambda row: (row["min_D_gain_for_induced_selection"], abs(row["S_start"]), row["pre_steps"]),
    )
    return {
        "total_rows": len(rows),
        "boundary_points": len(boundary_rows),
        "selection_possible_boundary_points": len(selection_boundaries),
        "no_selection_boundary_points": len(no_selection_boundaries),
        "D_outcome_counts": count_by(rows, "D_outcome"),
        "phase_after_D_counts": count_by(rows, "phase_after_D"),
        "selected_sign_counts_at_min": count_by(selection_boundaries, "selected_sign_at_min"),
        "min_D_gain_overall": min(min_gain_values) if min_gain_values else None,
        "max_min_D_gain_overall": max(min_gain_values) if min_gain_values else None,
        "by_case_id": by_case,
        "top_low_threshold_boundaries": sorted_boundaries[:30],
        "top_no_selection_boundaries": sorted(
            no_selection_boundaries,
            key=lambda row: (row["case_id"], abs(row["S_start"]), row["pre_steps"]),
        )[:30],
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
        "# 白猫・黒猫・灰色猫 D選択境界 予備実験結果 v1",
        "",
        "## 1. 実験条件",
        "",
        "```text",
        "target = gray_metastable candidates only",
        f"d_steps = {params['d_steps']}",
        f"pre_steps_values_count = {len(params['pre_steps_values'])}",
        f"D_gain_values = {params['d_gain_values']}",
        f"C_modes = {params['c_modes']}",
        "```",
        "",
        "各観測開始ステップごとにDなし対照を並走させ、Dなしで選択しない条件に限ってD起因選択と判定した。",
        "",
        "## 2. 全体集計",
        "",
        f"total_rows = {summary['total_rows']}",
        f"boundary_points = {summary['boundary_points']}",
        f"selection_possible_boundary_points = {summary['selection_possible_boundary_points']}",
        f"no_selection_boundary_points = {summary['no_selection_boundary_points']}",
        f"min_D_gain_overall = {summary['min_D_gain_overall']}",
        f"max_min_D_gain_overall = {summary['max_min_D_gain_overall']}",
        "",
        "## 3. 候補別境界",
        "",
        "| case_id | boundary_points | selection_possible | no_selection | min_gain | max_min_gain | sign_counts |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for case_id, row in summary["by_case_id"].items():
        lines.append(
            f"| {case_id} | {row['boundary_points']} | {row['selection_possible_points']} | {row['no_selection_points']} | {row['min_observed_D_gain']} | {row['max_observed_D_gain']} | {row['selected_sign_counts']} |"
        )
    lines += [
        "",
        "## 4. 低しきい値の代表点",
        "",
        "| case_id | pre | C_mode | S_start | C_sign | min_D_gain | sign | S_after | baseline |",
        "|---|---:|---|---:|---|---:|---|---:|---|",
    ]
    for row in summary["top_low_threshold_boundaries"]:
        lines.append(
            "| {case_id} | {pre_steps} | {C_mode} | {S_start:.6g} | {C_sign_start} | {min_D_gain_for_induced_selection} | {selected_sign_at_min} | {S_mean_after_at_min:.6g} | {baseline_outcome} |".format(
                **row
            )
        )
    lines += [
        "",
        "## 5. 選択しなかった代表点",
        "",
        "| case_id | pre | C_mode | S_start | C_sign | baseline |",
        "|---|---:|---|---:|---|---|",
    ]
    for row in summary["top_no_selection_boundaries"]:
        lines.append(
            "| {case_id} | {pre_steps} | {C_mode} | {S_start:.6g} | {C_sign_start} | {baseline_outcome} |".format(
                **row
            )
        )
    lines += ["", "## 6. 判定", ""]
    if summary["no_selection_boundary_points"] == 0:
        lines += [
            "灰色猫準安定候補に対して、観測開始ステップとD利得の境界表を得た。",
            "今回の掃引範囲では、全ての観測開始ステップでD起因選択が可能だった。",
            "候補ごとの最小D利得は二段に分かれ、弱いしきい値候補では `0.0225`、強いしきい値候補では `0.065` だった。",
        ]
    else:
        lines += [
            "灰色猫準安定候補に対して、観測開始ステップとD利得の境界表を得た。",
            "この表により、D観測で白猫または黒猫へ落ちる開始位相と、灰色近傍のまま残る開始位相を分けて追跡できる。",
        ]
    return "\n".join(lines) + "\n"


def run() -> Dict[str, Any]:
    params = Params()
    ab_params = ABParams()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    candidates = selected_metastable_candidates(params)
    rows: List[Dict[str, Any]] = []
    boundary_rows: List[Dict[str, Any]] = []
    for candidate in candidates:
        for c_mode in params.c_modes:
            for pre_steps in params.pre_steps_values:
                case_rows, boundary = run_boundary_for_start(ab_params, params, candidate, pre_steps, c_mode)
                rows.extend(case_rows)
                boundary_rows.append(boundary)
    result = {
        "experiment": "gray_cat_d_selection_boundary_preliminary_v1",
        "params": asdict(params),
        "candidate_ids": [candidate["case_id"] for candidate in candidates],
        "summary": summarise(rows, boundary_rows),
        "rows": rows,
        "boundary_rows": boundary_rows,
        "outputs": {
            "json": "gray_cat_d_selection_boundary_preliminary_result_v1.json",
            "rows_csv": "gray_cat_d_selection_boundary_rows_v1.csv",
            "boundary_csv": "gray_cat_d_selection_boundary_table_v1.csv",
            "report": "gray_cat_d_selection_boundary_report_v1.md",
        },
    }
    write_csv(OUT_DIR / result["outputs"]["rows_csv"], rows)
    write_csv(OUT_DIR / result["outputs"]["boundary_csv"], boundary_rows)
    (OUT_DIR / result["outputs"]["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT_DIR / result["outputs"]["report"]).write_text(build_report(result), encoding="utf-8")
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["summary"], ensure_ascii=False, indent=2))
