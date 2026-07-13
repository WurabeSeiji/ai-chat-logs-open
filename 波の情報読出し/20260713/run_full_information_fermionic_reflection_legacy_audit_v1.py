from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

from run_full_information_fermionic_localization_transfer_preliminary_v1 import (
    Params as FullParams,
    copy_distance,
    eta_distribution,
    inner,
    make_state,
    pure_harmonic_distribution,
    rho_components,
    rho_expectation,
    op_inner_p_chi,
)


BASE_DIR = Path(__file__).resolve().parent
ROOT_DIR = BASE_DIR.parent
OUT_DIR = BASE_DIR / "full_information_fermionic_reflection_legacy_audit_result_v1"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


P_TOL = 1.0e-2
MODE_TOL = 1.0e-2
COPY_DISTANCE_TOL = 1.0e-2


LEGACY_RESULTS: List[Dict[str, str]] = [
    {
        "series": "20260710",
        "kind": "copy_reflection",
        "name": "elastic_collision_simulation",
        "path": "20260710/elastic_collision_simulation_result_v1/elastic_collision_result_v1.json",
        "valid_key": "judgement.elastic_collision_map_valid",
    },
    {
        "series": "20260710",
        "kind": "fermionic_ab_c_replacement",
        "name": "fermionic_interference_reflection",
        "path": "20260710/fermionic_interference_reflection_result_v1/fermionic_interference_reflection_result_v1.json",
        "valid_key": "verdict.mechanism_valid_minimal",
    },
    {
        "series": "20260710",
        "kind": "copy_reflection",
        "name": "elastic_collision_multi_collision",
        "path": "20260710/elastic_collision_multi_collision_result_v1/multi_collision_result_v1.json",
        "valid_key": "summary.multi_collision_valid",
    },
    {
        "series": "20260710",
        "kind": "control_maps",
        "name": "elastic_collision_control_maps",
        "path": "20260710/elastic_collision_control_maps_result_v1/control_maps_result_v1.json",
        "valid_key": "summary.reflection_valid_maps",
    },
    {
        "series": "20260710",
        "kind": "copy_reflection_cases",
        "name": "elastic_collision_cell_resolution_sweep",
        "path": "20260710/elastic_collision_cell_resolution_sweep_result_v1/cell_resolution_sweep_result_v1.json",
        "valid_key": "case_valid",
    },
    {
        "series": "20260710",
        "kind": "copy_reflection_cases",
        "name": "elastic_collision_label_robustness",
        "path": "20260710/elastic_collision_label_robustness_result_v1/label_robustness_result_v1.json",
        "valid_key": "case_valid",
    },
    {
        "series": "20260710",
        "kind": "copy_reflection_cases",
        "name": "elastic_collision_eta_resolution_sweep",
        "path": "20260710/elastic_collision_eta_resolution_sweep_result_v1/eta_resolution_sweep_result_v1.json",
        "valid_key": "case_valid",
    },
    {
        "series": "20260710",
        "kind": "copy_reflection_cases",
        "name": "elastic_collision_observer_sweep",
        "path": "20260710/elastic_collision_observer_sweep_result_v1/observer_sweep_result_v1.json",
        "valid_key": "model_valid",
    },
    {
        "series": "20260710",
        "kind": "copy_reflection_cases",
        "name": "elastic_collision_observation_perturbation",
        "path": "20260710/elastic_collision_observation_perturbation_result_v1/observation_perturbation_result_v1.json",
        "valid_key": "model_valid",
    },
    {
        "series": "20260710",
        "kind": "copy_reflection_cases",
        "name": "elastic_collision_asymmetry_sweep",
        "path": "20260710/elastic_collision_asymmetry_sweep_result_v1/asymmetry_sweep_result_v1.json",
        "valid_key": "case_valid",
    },
    {
        "series": "20260711",
        "kind": "abc_simple_reflection",
        "name": "abc_multigauge_interference_readout",
        "path": "20260711/abc_multigauge_interference_readout_result_v1/abc_multigauge_interference_readout_result_v1.json",
        "valid_key": "verdicts.multigauge_measurement_valid",
    },
    {
        "series": "20260711",
        "kind": "abc_simple_reflection",
        "name": "abc_multigauge_interference_readout_multi_collision",
        "path": "20260711/abc_multigauge_interference_readout_multi_collision_result_v1/abc_multigauge_interference_readout_multi_collision_result_v1.json",
        "valid_key": "verdicts.multi_collision_multigauge_valid",
    },
    {
        "series": "20260711",
        "kind": "abc_simple_reflection_cases",
        "name": "abc_multigauge_interference_readout_asymmetric_amplitude_sweep",
        "path": "20260711/abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1/abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v1.json",
        "valid_key": "individual_multigauge_valid",
    },
    {
        "series": "20260711",
        "kind": "abc_simple_reflection_cases",
        "name": "abc_multigauge_interference_readout_robustness_sweep",
        "path": "20260711/abc_multigauge_interference_readout_robustness_sweep_result_v1/abc_multigauge_interference_readout_robustness_sweep_result_v1.json",
        "valid_key": "case_valid",
    },
    {
        "series": "20260711",
        "kind": "abc_generalized_elastic_cases",
        "name": "abc_multigauge_generalized_elastic_collision_readout",
        "path": "20260711/abc_multigauge_generalized_elastic_collision_readout_result_v1/abc_multigauge_generalized_elastic_collision_readout_result_v1.json",
        "valid_key": "individual_readout_valid",
    },
    {
        "series": "20260711",
        "kind": "abc_generalized_elastic_cases",
        "name": "abc_multigauge_generalized_elastic_collision_velocity_sweep",
        "path": "20260711/abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1/abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v1.json",
        "valid_key": "individual_readout_valid",
    },
    {
        "series": "20260711",
        "kind": "abc_generalized_elastic_cases",
        "name": "abc_multigauge_generalized_elastic_collision_multi_collision",
        "path": "20260711/abc_multigauge_generalized_elastic_collision_multi_collision_result_v1/abc_multigauge_generalized_elastic_collision_multi_collision_result_v1.json",
        "valid_key": "case_valid",
    },
    {
        "series": "20260711",
        "kind": "abc_generalized_elastic_cases",
        "name": "abc_multigauge_generalized_elastic_collision_noise_robustness",
        "path": "20260711/abc_multigauge_generalized_elastic_collision_noise_robustness_result_v1/abc_multigauge_generalized_elastic_collision_noise_robustness_result_v1.json",
        "valid_key": "noise_row_valid",
    },
    {
        "series": "20260711",
        "kind": "abc_generalized_elastic_cases",
        "name": "abc_multigauge_generalized_elastic_collision_extreme_R_sweep",
        "path": "20260711/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v1/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v1.json",
        "valid_key": "boundary_case_valid",
    },
]


def read_json(path: str) -> Dict[str, Any]:
    return json.loads((ROOT_DIR / path).read_text(encoding="utf-8"))


def nested_get(obj: Dict[str, Any], dotted: str, default: Any = None) -> Any:
    current: Any = obj
    for key in dotted.split("."):
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def bool_from_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, list):
        return len(value) > 0
    if isinstance(value, (int, float)):
        return bool(value)
    return False


def pick_params(data: Dict[str, Any], case: Dict[str, Any]) -> Dict[str, Any]:
    params = dict(data.get("parameters") or data.get("params") or {})
    case_params = case.get("parameters") if isinstance(case.get("parameters"), dict) else {}
    params.update(case_params)
    for key in [
        "A_A",
        "A_B",
        "q_A0",
        "q_B0",
        "m_A",
        "m_B",
        "Nh_chi_A",
        "Nh_chi_B",
        "Nh_tau_A",
        "Nh_tau_B",
    ]:
        if key in case:
            params[key] = case[key]
    return params


def generalized_velocity(r_a: float, r_b: float, u_a: float, u_b: float) -> Tuple[float, float]:
    denom = r_a + r_b
    return (
        ((r_a - r_b) / denom) * u_a + (2.0 * r_b / denom) * u_b,
        (2.0 * r_a / denom) * u_a + ((r_b - r_a) / denom) * u_b,
    )


def case_rows(data: Dict[str, Any], meta: Dict[str, str]) -> List[Dict[str, Any]]:
    if meta["name"] == "elastic_collision_control_maps":
        return [
            row
            for row in data.get("cases", [])
            if str(row.get("map")) in {"reflection", "label_exchange_reflection"}
        ]
    if meta["name"] == "abc_multigauge_generalized_elastic_collision_noise_robustness":
        return list(data.get("summary_rows", []))
    if isinstance(data.get("cases"), list):
        return list(data["cases"])
    if isinstance(data.get("case_summaries"), list):
        return list(data["case_summaries"])
    if isinstance(data.get("case_results"), list):
        out: List[Dict[str, Any]] = []
        for entry in data["case_results"]:
            if isinstance(entry, dict) and isinstance(entry.get("case_summary"), dict):
                merged = dict(entry["case_summary"])
                if "parameters" in entry and isinstance(entry["parameters"], dict):
                    merged["parameters"] = entry["parameters"]
                out.append(merged)
        if out:
            return out
    if meta["name"] == "fermionic_interference_reflection":
        summary = data["ab_c_replacement_summary"]
        return [
            {
                "case": "ab_c_replacement_pi",
                "q_A_after": summary["q_readouts"]["A_after"],
                "q_B_after": summary["q_readouts"]["B_after"],
                "case_valid": summary["verdict"]["integrated_ab_c_replacement_valid"],
            }
        ]
    return [{"case": meta["name"]}]


def legacy_valid(data: Dict[str, Any], case: Dict[str, Any], meta: Dict[str, str]) -> bool:
    key = meta["valid_key"]
    if key == "noise_row_valid":
        if str(case.get("noise_mode")) == "zero_mean_gauge_noise":
            return bool(case.get("zero_mean_multigauge_valid", False))
        if str(case.get("noise_mode")) == "common_bias_control":
            noise_level = float(case.get("noise_level", 0.0))
            return bool(noise_level == 0.0 or case.get("common_bias_detected", False))
        return False
    if "." in key:
        return bool_from_value(nested_get(data, key, False))
    if key in case:
        return bool_from_value(case[key])
    return bool_from_value(nested_get(data, key, False))


def target_q(meta: Dict[str, str], params: Dict[str, Any], case: Dict[str, Any]) -> Tuple[float, float]:
    q_a0 = float(params.get("q_A0", 1.0))
    q_b0 = float(params.get("q_B0", -1.0))
    if "q_A_after" in case and "q_B_after" in case:
        return float(case["q_A_after"]), float(case["q_B_after"])
    if meta["kind"].startswith("abc_generalized"):
        r_a = float(params.get("A_A", 1.0)) ** 2
        r_b = float(params.get("A_B", 1.0)) ** 2
        return generalized_velocity(r_a, r_b, q_a0, q_b0)
    return -q_a0, -q_b0


def finite_int(value: Any, default: int) -> int:
    try:
        number = int(round(float(value)))
    except (TypeError, ValueError):
        return default
    return number


def audit_full_information(
    params_raw: Dict[str, Any],
    q_target_a: float,
    q_target_b: float,
) -> Dict[str, Any]:
    n_a = finite_int(params_raw.get("Nh_chi_A"), 99)
    n_b = finite_int(params_raw.get("Nh_chi_B"), 99)
    q_a0 = float(params_raw.get("q_A0", 1.0))
    q_b0 = float(params_raw.get("q_B0", -1.0))
    m_a = finite_int(params_raw.get("m_A"), 1)
    m_b = finite_int(params_raw.get("m_B"), 2)
    a_a = float(params_raw.get("A_A", 1.0))
    a_b = float(params_raw.get("A_B", 1.0))

    params = replace(
        FullParams(),
        high_n=max(abs(n_a), abs(n_b), 3),
        q_A=q_a0,
        q_B=q_b0,
        m_A=m_a,
        m_B=m_b,
        A_A=a_a,
        A_B=a_b,
        chi_grid_n=max(256, 2 * (max(abs(n_a), abs(n_b), 3) + 3)),
    )
    alpha = make_state(params, n_a, q_a0, m_a, True, a_a)
    beta = make_state(params, n_b, q_b0, m_b, True, a_b)
    comp = rho_components(alpha, beta, math.pi)

    p_read = float(rho_expectation(alpha, beta, comp, lambda v, w: op_inner_p_chi(params, v, w)).real)
    eta = eta_distribution(params, alpha, beta, comp, [m_a, m_b])
    target_a = make_state(params, n_a, q_target_a, m_a, True, a_a)
    target_b = make_state(params, n_b, q_target_b, m_b, True, a_b)
    copy_distance_a = copy_distance(alpha, beta, comp, target_a)
    copy_distance_b = copy_distance(alpha, beta, comp, target_b)
    p_error_a = abs(p_read - q_target_a)
    p_error_b = abs(p_read - q_target_b)
    mode_prob_a = float(eta.get(m_a, 0.0))
    mode_prob_b = float(eta.get(m_b, 0.0))
    reproduced = bool(
        p_error_a <= P_TOL
        and p_error_b <= P_TOL
        and mode_prob_a >= 1.0 - MODE_TOL
        and mode_prob_b >= 1.0 - MODE_TOL
        and copy_distance_a <= COPY_DISTANCE_TOL
        and copy_distance_b <= COPY_DISTANCE_TOL
    )
    return {
        "N_A": n_a,
        "N_B": n_b,
        "A_A": a_a,
        "A_B": a_b,
        "q_A0": q_a0,
        "q_B0": q_b0,
        "m_A": m_a,
        "m_B": m_b,
        "q_target_A": q_target_a,
        "q_target_B": q_target_b,
        "full_p_read": p_read,
        "full_p_error_A": p_error_a,
        "full_p_error_B": p_error_b,
        "full_P_m_A": mode_prob_a,
        "full_P_m_B": mode_prob_b,
        "full_copy_distance_A": copy_distance_a,
        "full_copy_distance_B": copy_distance_b,
        "full_reproduces_legacy": reproduced,
        "significant_difference": not reproduced,
    }


def audit_rows() -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    rows: List[Dict[str, Any]] = []
    experiment_rows: List[Dict[str, Any]] = []
    for meta in LEGACY_RESULTS:
        data = read_json(meta["path"])
        cases = case_rows(data, meta)
        exp_start = len(rows)
        for index, case in enumerate(cases):
            params_raw = pick_params(data, case)
            q_a, q_b = target_q(meta, params_raw, case)
            full = audit_full_information(params_raw, q_a, q_b)
            old_valid = legacy_valid(data, case, meta)
            row = {
                "series": meta["series"],
                "experiment": meta["name"],
                "kind": meta["kind"],
                "case_index": index,
                "case": str(case.get("case") or case.get("map") or case.get("name") or index),
                "legacy_valid": old_valid,
                **full,
            }
            rows.append(row)
        exp_rows = rows[exp_start:]
        legacy_valid_count = sum(1 for row in exp_rows if row["legacy_valid"])
        full_reproduced_count = sum(1 for row in exp_rows if row["full_reproduces_legacy"])
        significant_count = sum(1 for row in exp_rows if row["legacy_valid"] and row["significant_difference"])
        experiment_rows.append(
            {
                "series": meta["series"],
                "experiment": meta["name"],
                "kind": meta["kind"],
                "case_count": len(exp_rows),
                "legacy_valid_count": legacy_valid_count,
                "full_reproduced_count": full_reproduced_count,
                "significant_difference_among_legacy_valid": significant_count,
                "all_legacy_valid_cases_reproduced": bool(
                    legacy_valid_count > 0 and significant_count == 0 and full_reproduced_count >= legacy_valid_count
                ),
            }
        )
    return rows, experiment_rows


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def make_plots(rows: List[Dict[str, Any]], experiment_rows: List[Dict[str, Any]]) -> Dict[str, str]:
    labels = [row["experiment"].replace("abc_multigauge_", "abc_").replace("elastic_collision_", "ec_") for row in experiment_rows]
    legacy_valid = [int(row["legacy_valid_count"]) for row in experiment_rows]
    full_reproduced = [int(row["full_reproduced_count"]) for row in experiment_rows]
    significant = [int(row["significant_difference_among_legacy_valid"]) for row in experiment_rows]

    fig, ax = plt.subplots(figsize=(13, 6), constrained_layout=True)
    x = np.arange(len(labels))
    ax.bar(x - 0.25, legacy_valid, width=0.25, label="legacy valid cases")
    ax.bar(x, full_reproduced, width=0.25, label="full-info reproduced cases")
    ax.bar(x + 0.25, significant, width=0.25, label="significant differences")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=8)
    ax.set_ylabel("case count")
    ax.set_title("Legacy fermionic reflection series vs full-information exchange audit")
    ax.legend()
    path1 = OUT_DIR / "full_information_fermionic_reflection_legacy_audit_counts_v1.png"
    fig.savefig(path1, dpi=160)
    plt.close(fig)

    valid_rows = [row for row in rows if row["legacy_valid"]]
    p_errors = [max(float(row["full_p_error_A"]), float(row["full_p_error_B"])) for row in valid_rows]
    mode_mins = [min(float(row["full_P_m_A"]), float(row["full_P_m_B"])) for row in valid_rows]
    copy_max = [max(float(row["full_copy_distance_A"]), float(row["full_copy_distance_B"])) for row in valid_rows]
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), constrained_layout=True)
    axes[0].hist(p_errors, bins=30)
    axes[0].axvline(P_TOL, color="black", linestyle="--")
    axes[0].set_xlabel("max p error")
    axes[1].hist(mode_mins, bins=30)
    axes[1].axvline(1.0 - MODE_TOL, color="black", linestyle="--")
    axes[1].set_xlabel("min target mode probability")
    axes[2].hist(copy_max, bins=30)
    axes[2].axvline(COPY_DISTANCE_TOL, color="black", linestyle="--")
    axes[2].set_xlabel("max copy distance")
    path2 = OUT_DIR / "full_information_fermionic_reflection_legacy_audit_metric_histograms_v1.png"
    fig.savefig(path2, dpi=160)
    plt.close(fig)
    return {"count_plot": path1.name, "metric_plot": path2.name}


def build_report(result: Dict[str, Any]) -> str:
    aggregate = result["aggregate"]
    experiment_lines = []
    for row in result["experiment_rows"]:
        experiment_lines.append(
            f"| {row['series']} | {row['experiment']} | {row['case_count']} | "
            f"{row['legacy_valid_count']} | {row['full_reproduced_count']} | "
            f"{row['significant_difference_among_legacy_valid']} | "
            f"`{str(row['all_legacy_valid_cases_reproduced']).lower()}` |"
        )

    sorted_worst_rows = sorted(
        [row for row in result["case_rows"] if row["legacy_valid"]],
        key=lambda row: max(float(row["full_p_error_A"]), float(row["full_p_error_B"]), float(row["full_copy_distance_A"]), float(row["full_copy_distance_B"])),
        reverse=True,
    )
    worst_rows = []
    seen_worst_keys = set()
    for row in sorted_worst_rows:
        key = (row["experiment"], row["case"])
        if key in seen_worst_keys:
            continue
        seen_worst_keys.add(key)
        worst_rows.append(row)
        if len(worst_rows) >= 12:
            break
    worst_lines = []
    for row in worst_rows:
        worst_lines.append(
            f"| {row['experiment']} | {row['case']} | {row['q_target_A']:.6g} | {row['q_target_B']:.6g} | "
            f"{row['full_p_read']:.6g} | {max(row['full_p_error_A'], row['full_p_error_B']):.6g} | "
            f"{min(row['full_P_m_A'], row['full_P_m_B']):.6g} | "
            f"{max(row['full_copy_distance_A'], row['full_copy_distance_B']):.6g} |"
        )

    return f"""# 全情報交換干渉によるフェルミオン反射系列差分監査 v1

## 目的

20260710 および 20260711 のフェルミオン反射または弾性衝突読出し系列について、旧実験の同一条件を読み直し、保存コピー近似ではなく全情報交換干渉写像を適用した場合に、旧結果を再現できるかを検査した。

ここでは推測で判定しない。旧 JSON の有効ケースを読み込み、同じ `A_A`, `A_B`, `q_A0`, `q_B0`, `m_A`, `m_B`, `Nh_chi_A`, `Nh_chi_B` を用いて、全情報交換縮約から `p`、識別モード純度、保存コピー距離を再計算した。

## 総合判定

| 指標 | 値 |
|---|---:|
| 監査対象実験数 | `{aggregate['experiment_count']}` |
| 監査対象ケース数 | `{aggregate['case_count']}` |
| 旧有効ケース数 | `{aggregate['legacy_valid_case_count']}` |
| 全情報交換で旧条件を再現したケース数 | `{aggregate['full_reproduced_case_count']}` |
| 旧有効ケース中の有意差分ケース数 | `{aggregate['significant_difference_case_count']}` |
| 全旧有効ケースを再現 | `{str(aggregate['all_legacy_valid_cases_reproduced']).lower()}` |

![count plot]({result['outputs']['count_plot']})

![metric plot]({result['outputs']['metric_plot']})

## 実験別監査

| series | experiment | cases | legacy valid | full-info reproduced | significant diff among legacy valid | all reproduced |
|---|---|---:|---:|---:|---:|---|
{chr(10).join(experiment_lines)}

## 最大差分例

| experiment | case | q target A | q target B | full p read | max p error | min mode prob | max copy distance |
|---|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(worst_lines)}

## 解釈

旧系列は、衝突点で `q` または一般化速度を更新し、識別モード、振幅、倍音構造を保存する実効的な保存コピー近似で動いていた。

全情報交換干渉写像では、A と B の一体情報をまとめて交換合成した後に縮約するため、旧実験が前提にしていた個別スロットの保存コピーが成立しない。等振幅の基本条件でも、識別モードはおおむね A/B 半々に混合し、同じ縮約状態から A 用と B 用の二つの反射後状態を同時に復元できない。

したがって、過去のフェルミオン反射系列、ABC 多ゲージ読出し系列、一般化弾性衝突系列は、旧目的に対しては有効な近似実験だったが、局在性移乗や倍音移乗を検査するための全情報干渉実験としては、そのままでは再利用できない。

## 結論

同一条件で監査した結果、保存コピー近似と全情報交換干渉写像の差は有意である。

よって、エネルギー系・加速度系など、旧フェルミオン反射写像に依存する後続実験へ全情報交換干渉を持ち込む場合は、旧実験を単純に流用せず、全情報交換版として再実装し直した上で差分を評価する必要がある。

## 出力

| 種別 | ファイル |
|---|---|
| JSON | `{result['outputs']['json']}` |
| case CSV | `{result['outputs']['case_csv']}` |
| experiment CSV | `{result['outputs']['experiment_csv']}` |
| count plot | `{result['outputs']['count_plot']}` |
| metric plot | `{result['outputs']['metric_plot']}` |
| report | `{result['outputs']['report']}` |
"""


def run() -> Dict[str, Any]:
    rows, experiment_rows = audit_rows()
    outputs = {
        "json": "full_information_fermionic_reflection_legacy_audit_result_v1.json",
        "case_csv": "full_information_fermionic_reflection_legacy_audit_cases_v1.csv",
        "experiment_csv": "full_information_fermionic_reflection_legacy_audit_experiments_v1.csv",
        "report": "full_information_fermionic_reflection_legacy_audit_report_v1.md",
    }
    outputs.update(make_plots(rows, experiment_rows))
    aggregate = {
        "experiment_count": len(experiment_rows),
        "case_count": len(rows),
        "legacy_valid_case_count": sum(1 for row in rows if row["legacy_valid"]),
        "full_reproduced_case_count": sum(1 for row in rows if row["full_reproduces_legacy"]),
        "significant_difference_case_count": sum(1 for row in rows if row["legacy_valid"] and row["significant_difference"]),
    }
    aggregate["all_legacy_valid_cases_reproduced"] = bool(
        aggregate["legacy_valid_case_count"] > 0
        and aggregate["significant_difference_case_count"] == 0
        and aggregate["full_reproduced_case_count"] >= aggregate["legacy_valid_case_count"]
    )
    result = {
        "experiment": "full_information_fermionic_reflection_legacy_audit_v1",
        "purpose": "Audit every 20260710/20260711 fermionic-reflection or elastic-collision readout series against the full-information exchange map under the same recorded conditions.",
        "thresholds": {
            "p_tol": P_TOL,
            "mode_tol": MODE_TOL,
            "copy_distance_tol": COPY_DISTANCE_TOL,
        },
        "legacy_results": LEGACY_RESULTS,
        "aggregate": aggregate,
        "experiment_rows": experiment_rows,
        "case_rows": rows,
        "outputs": outputs,
    }
    write_csv(OUT_DIR / outputs["case_csv"], rows)
    write_csv(OUT_DIR / outputs["experiment_csv"], experiment_rows)
    (OUT_DIR / outputs["json"]).write_text(json.dumps(json_safe(result), ensure_ascii=False, indent=2), encoding="utf-8")
    report = build_report(json_safe(result))
    (OUT_DIR / outputs["report"]).write_text(report, encoding="utf-8")
    (BASE_DIR / "全情報交換干渉によるフェルミオン反射系列差分監査予備実験検証メモ_v1.md").write_text(
        report,
        encoding="utf-8",
    )
    return result


if __name__ == "__main__":
    print(json.dumps(json_safe(run()["aggregate"]), ensure_ascii=False, indent=2))
