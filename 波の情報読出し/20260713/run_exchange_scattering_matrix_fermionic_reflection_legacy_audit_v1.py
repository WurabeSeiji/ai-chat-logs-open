from __future__ import annotations

import csv
import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

from run_full_information_fermionic_reflection_legacy_audit_v1 import (
    LEGACY_RESULTS,
    bool_from_value,
    case_rows,
    finite_int,
    json_safe,
    legacy_valid,
    pick_params,
    read_json,
    target_q,
)
from run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1 import (
    Params,
    make_state,
    normalized_distance,
    norm2,
    pure_expect_p,
    scattering_coefficients,
    scattering_outputs,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "exchange_scattering_matrix_fermionic_reflection_legacy_audit_result_v1"
OUT_DIR.mkdir(exist_ok=True)

P_TOL = 1.0e-2
COPY_DISTANCE_TOL = 1.0e-2


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def as_params(params_raw: Dict[str, Any]) -> Tuple[Params, Dict[str, Any]]:
    n_a = finite_int(params_raw.get("Nh_chi_A"), 99)
    n_b = finite_int(params_raw.get("Nh_chi_B"), 99)
    q_a0 = float(params_raw.get("q_A0", 1.0))
    q_b0 = float(params_raw.get("q_B0", -1.0))
    m_a = finite_int(params_raw.get("m_A"), 1)
    m_b = finite_int(params_raw.get("m_B"), 2)
    a_a = float(params_raw.get("A_A", 1.0))
    a_b = float(params_raw.get("A_B", 1.0))
    max_n = max(abs(n_a), abs(n_b), 3)
    params = replace(
        Params(),
        high_n=max_n,
        q_A=q_a0,
        q_B=q_b0,
        m_A=m_a,
        m_B=m_b,
        A_A=a_a,
        A_B=a_b,
        chi_grid_n=max(256, 2 * (max_n + 3)),
    )
    scalars = {
        "N_A": n_a,
        "N_B": n_b,
        "q_A0": q_a0,
        "q_B0": q_b0,
        "m_A": m_a,
        "m_B": m_b,
        "A_A": a_a,
        "A_B": a_b,
    }
    return params, scalars


def nearly_same(x: float, y: float, tol: float = P_TOL) -> bool:
    return abs(float(x) - float(y)) <= tol


def weighted_p(
    r_prob: float,
    q_a: float,
    q_b: float,
    a_a: float,
    a_b: float,
    channel: str,
) -> float:
    t_prob = 1.0 - r_prob
    if channel == "minus_out":
        n = r_prob * (a_a**2) * (-q_a) + t_prob * (a_b**2) * q_b
        d = r_prob * (a_a**2) + t_prob * (a_b**2)
    elif channel == "plus_out":
        n = t_prob * (a_a**2) * q_a + r_prob * (a_b**2) * (-q_b)
        d = t_prob * (a_a**2) + r_prob * (a_b**2)
    else:
        raise ValueError(channel)
    if d <= 0.0:
        return float("nan")
    return n / d


def best_probability_fit(
    q_a: float,
    q_b: float,
    a_a: float,
    a_b: float,
    q_target_a: float,
    q_target_b: float,
) -> Dict[str, float]:
    best = {
        "best_R": float("nan"),
        "best_T": float("nan"),
        "best_p_minus": float("nan"),
        "best_p_plus": float("nan"),
        "best_max_p_error": float("inf"),
        "best_l2_p_error": float("inf"),
    }
    for r_prob in np.linspace(0.0, 1.0, 2001):
        p_minus = weighted_p(float(r_prob), q_a, q_b, a_a, a_b, "minus_out")
        p_plus = weighted_p(float(r_prob), q_a, q_b, a_a, a_b, "plus_out")
        e_a = abs(p_minus - q_target_a)
        e_b = abs(p_plus - q_target_b)
        max_error = max(e_a, e_b)
        l2_error = math.sqrt(e_a * e_a + e_b * e_b)
        if (max_error, l2_error) < (best["best_max_p_error"], best["best_l2_p_error"]):
            best.update(
                {
                    "best_R": float(r_prob),
                    "best_T": float(1.0 - r_prob),
                    "best_p_minus": float(p_minus),
                    "best_p_plus": float(p_plus),
                    "best_max_p_error": float(max_error),
                    "best_l2_p_error": float(l2_error),
                }
            )
    return best


def audit_scattering_matrix(
    params_raw: Dict[str, Any],
    q_target_a: float,
    q_target_b: float,
) -> Dict[str, Any]:
    params, scalars = as_params(params_raw)
    n_a = int(scalars["N_A"])
    n_b = int(scalars["N_B"])
    q_a0 = float(scalars["q_A0"])
    q_b0 = float(scalars["q_B0"])
    a_a = float(scalars["A_A"])
    a_b = float(scalars["A_B"])

    complete_reflection_target = nearly_same(q_target_a, -q_a0) and nearly_same(q_target_b, -q_b0)
    out = scattering_outputs(params, n_a, n_b, True, math.pi)
    target_a = make_state(params, n_a, q_target_a, int(scalars["m_A"]), True, a_a)
    target_b = make_state(params, n_b, q_target_b, int(scalars["m_B"]), True, a_b)

    p_minus = pure_expect_p(params, out["out_minus"])
    p_plus = pure_expect_p(params, out["out_plus"])
    p_error_a = abs(p_minus - q_target_a)
    p_error_b = abs(p_plus - q_target_b)
    copy_distance_a = normalized_distance(out["out_minus"], target_a)
    copy_distance_b = normalized_distance(out["out_plus"], target_b)
    norm_minus = norm2(out["out_minus"])
    norm_plus = norm2(out["out_plus"])
    p_reproduced = p_error_a <= P_TOL and p_error_b <= P_TOL
    state_reproduced = copy_distance_a <= COPY_DISTANCE_TOL and copy_distance_b <= COPY_DISTANCE_TOL

    best = best_probability_fit(q_a0, q_b0, a_a, a_b, q_target_a, q_target_b)
    t, r, T, R = scattering_coefficients(math.pi)
    return {
        **scalars,
        "q_target_A": q_target_a,
        "q_target_B": q_target_b,
        "target_class": "complete_reflection" if complete_reflection_target else "non_complete_or_generalized",
        "complete_R": R,
        "complete_T": T,
        "complete_t_real": float(t.real),
        "complete_t_imag": float(t.imag),
        "complete_r_real": float(r.real),
        "complete_r_imag": float(r.imag),
        "p_minus": p_minus,
        "p_plus": p_plus,
        "p_error_A": p_error_a,
        "p_error_B": p_error_b,
        "copy_distance_A": copy_distance_a,
        "copy_distance_B": copy_distance_b,
        "norm_minus": norm_minus,
        "norm_plus": norm_plus,
        "p_reproduced": bool(p_reproduced),
        "state_reproduced": bool(state_reproduced),
        "scattering_reproduces_legacy": bool(p_reproduced and state_reproduced),
        "needs_nontrivial_scattering_or_carrier_update": bool(not (p_reproduced and state_reproduced)),
        **best,
        "best_probability_fit_reproduces_p": bool(best["best_max_p_error"] <= P_TOL),
    }


def audit_rows() -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    rows: List[Dict[str, Any]] = []
    experiment_rows: List[Dict[str, Any]] = []
    for meta in LEGACY_RESULTS:
        data = read_json(meta["path"])
        cases = case_rows(data, meta)
        start = len(rows)
        for index, case in enumerate(cases):
            params_raw = pick_params(data, case)
            q_a, q_b = target_q(meta, params_raw, case)
            scattering = audit_scattering_matrix(params_raw, q_a, q_b)
            old_valid = legacy_valid(data, case, meta)
            row = {
                "series": meta["series"],
                "experiment": meta["name"],
                "kind": meta["kind"],
                "case_index": index,
                "case": str(case.get("case") or case.get("map") or case.get("name") or index),
                "legacy_valid": old_valid,
                **scattering,
            }
            rows.append(row)
        exp_rows = rows[start:]
        legacy_valid_rows = [row for row in exp_rows if row["legacy_valid"]]
        complete_rows = [row for row in legacy_valid_rows if row["target_class"] == "complete_reflection"]
        non_complete_rows = [row for row in legacy_valid_rows if row["target_class"] != "complete_reflection"]
        reproduced_rows = [row for row in legacy_valid_rows if row["scattering_reproduces_legacy"]]
        p_best_rows = [row for row in legacy_valid_rows if row["best_probability_fit_reproduces_p"]]
        experiment_rows.append(
            {
                "series": meta["series"],
                "experiment": meta["name"],
                "kind": meta["kind"],
                "case_count": len(exp_rows),
                "legacy_valid_count": len(legacy_valid_rows),
                "complete_reflection_target_count": len(complete_rows),
                "non_complete_target_count": len(non_complete_rows),
                "scattering_reproduced_count": len(reproduced_rows),
                "best_probability_fit_p_count": len(p_best_rows),
                "needs_v2_recalculation_count": len(legacy_valid_rows) - len(reproduced_rows),
                "all_legacy_valid_cases_reproduced": bool(
                    len(legacy_valid_rows) > 0 and len(reproduced_rows) == len(legacy_valid_rows)
                ),
            }
        )
    return rows, experiment_rows


def build_report(result: Dict[str, Any]) -> str:
    aggregate = result["aggregate"]
    experiment_lines = []
    for row in result["experiment_rows"]:
        experiment_lines.append(
            f"| {row['series']} | {row['experiment']} | {row['legacy_valid_count']} | "
            f"{row['complete_reflection_target_count']} | {row['non_complete_target_count']} | "
            f"{row['scattering_reproduced_count']} | {row['best_probability_fit_p_count']} | "
            f"{row['needs_v2_recalculation_count']} | `{str(row['all_legacy_valid_cases_reproduced']).lower()}` |"
        )

    valid_rows = [row for row in result["case_rows"] if row["legacy_valid"]]
    worst_rows = sorted(
        valid_rows,
        key=lambda row: max(float(row["p_error_A"]), float(row["p_error_B"]), float(row["copy_distance_A"]), float(row["copy_distance_B"])),
        reverse=True,
    )[:16]
    worst_lines = []
    for row in worst_rows:
        worst_lines.append(
            f"| {row['experiment']} | {row['case']} | {row['target_class']} | "
            f"{row['q_target_A']:.6g} | {row['q_target_B']:.6g} | "
            f"{row['p_minus']:.6g} | {row['p_plus']:.6g} | "
            f"{max(row['p_error_A'], row['p_error_B']):.6g} | "
            f"{max(row['copy_distance_A'], row['copy_distance_B']):.6g} | "
            f"{row['best_R']:.4g} | {row['best_max_p_error']:.6g} |"
        )

    return f"""# 交換干渉散乱行列によるフェルミオン反射系列差分監査 v1

## 目的

20260710 および 20260711 のフェルミオン反射・弾性衝突読出し系列について、旧 JSON の同一条件を読み直し、交換干渉散乱行列版で旧主読出しが再現されるかを監査した。

ここでの主判定は、名前の毛ではない。旧有効ケースに対して、二つの出力チャネルが旧ターゲットの運動量位相 `q` と状態形状を再現するかで判定する。

## 総合判定

| 指標 | 値 |
|---|---:|
| 監査対象実験数 | `{aggregate['experiment_count']}` |
| 監査対象ケース数 | `{aggregate['case_count']}` |
| 旧有効ケース数 | `{aggregate['legacy_valid_case_count']}` |
| 完全反射ターゲットの旧有効ケース数 | `{aggregate['complete_reflection_target_count']}` |
| 非完全反射または一般化ターゲットの旧有効ケース数 | `{aggregate['non_complete_target_count']}` |
| 散乱行列完全反射で旧主読出しを再現したケース数 | `{aggregate['scattering_reproduced_case_count']}` |
| 確率比のみなら最適近似できるケース数 | `{aggregate['best_probability_fit_p_count']}` |
| V2再計算検討が必要な旧有効ケース数 | `{aggregate['needs_v2_recalculation_case_count']}` |
| 全旧有効ケースを再現 | `{str(aggregate['all_legacy_valid_cases_reproduced']).lower()}` |

## 実験別監査

| series | experiment | legacy valid | complete target | non-complete target | scattering reproduced | best probability p-fit | needs V2 | all reproduced |
|---|---|---:|---:|---:|---:|---:|---:|---|
{chr(10).join(experiment_lines)}

## 最大差分例

| experiment | case | target class | q target A | q target B | p minus | p plus | max p error | max copy distance | best R | best max p error |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(worst_lines)}

## 判定

完全反射ターゲット、すなわち旧ターゲットが `q_A -> -q_A`, `q_B -> -q_B` であるケースでは、交換干渉散乱行列の `R=1` 極限が旧主読出しを再現する。

一方、一般化弾性衝突や非完全反射ターゲットでは、単純な `R=1` 反射だけでは旧ターゲットを再現しない。さらに、`R:T` の確率比だけを調整しても再現できないケースは、散乱確率ではなく、出力キャリア位相そのものを更新する一般化散乱写像が必要である。

したがって、20260710/20260711 側を V2 化する必要性は一律ではない。完全反射の主結論だけに依存する論文は、散乱行列版でも主読出しが保持される。一般化速度、非等振幅、ABC 多ゲージ一般化弾性衝突に依存する論文は、V2 再計算の対象である。

## 出力

| 種別 | ファイル |
|---|---|
| JSON | `{result['outputs']['json']}` |
| case CSV | `{result['outputs']['case_csv']}` |
| experiment CSV | `{result['outputs']['experiment_csv']}` |
| report | `{result['outputs']['report']}` |
"""


def run() -> Dict[str, Any]:
    rows, experiment_rows = audit_rows()
    aggregate = {
        "experiment_count": len(experiment_rows),
        "case_count": len(rows),
        "legacy_valid_case_count": sum(1 for row in rows if row["legacy_valid"]),
        "complete_reflection_target_count": sum(
            1 for row in rows if row["legacy_valid"] and row["target_class"] == "complete_reflection"
        ),
        "non_complete_target_count": sum(
            1 for row in rows if row["legacy_valid"] and row["target_class"] != "complete_reflection"
        ),
        "scattering_reproduced_case_count": sum(
            1 for row in rows if row["legacy_valid"] and row["scattering_reproduces_legacy"]
        ),
        "best_probability_fit_p_count": sum(
            1 for row in rows if row["legacy_valid"] and row["best_probability_fit_reproduces_p"]
        ),
    }
    aggregate["needs_v2_recalculation_case_count"] = (
        aggregate["legacy_valid_case_count"] - aggregate["scattering_reproduced_case_count"]
    )
    aggregate["all_legacy_valid_cases_reproduced"] = bool(
        aggregate["legacy_valid_case_count"] > 0
        and aggregate["scattering_reproduced_case_count"] == aggregate["legacy_valid_case_count"]
    )
    outputs = {
        "json": "exchange_scattering_matrix_fermionic_reflection_legacy_audit_result_v1.json",
        "case_csv": "exchange_scattering_matrix_fermionic_reflection_legacy_audit_cases_v1.csv",
        "experiment_csv": "exchange_scattering_matrix_fermionic_reflection_legacy_audit_experiments_v1.csv",
        "report": "exchange_scattering_matrix_fermionic_reflection_legacy_audit_report_v1.md",
    }
    result = {
        "experiment": "exchange_scattering_matrix_fermionic_reflection_legacy_audit_v1",
        "purpose": "Audit 20260710/20260711 legacy fermionic reflection series against the exchange-interference scattering-matrix map under recorded conditions.",
        "thresholds": {
            "p_tol": P_TOL,
            "copy_distance_tol": COPY_DISTANCE_TOL,
        },
        "aggregate": aggregate,
        "experiment_rows": experiment_rows,
        "case_rows": rows,
        "outputs": outputs,
    }
    safe = json_safe(result)
    write_csv(OUT_DIR / outputs["case_csv"], safe["case_rows"])
    write_csv(OUT_DIR / outputs["experiment_csv"], safe["experiment_rows"])
    (OUT_DIR / outputs["json"]).write_text(json.dumps(safe, ensure_ascii=False, indent=2), encoding="utf-8")
    report = build_report(safe)
    (OUT_DIR / outputs["report"]).write_text(report, encoding="utf-8")
    (BASE_DIR / "交換干渉散乱行列によるフェルミオン反射系列差分監査予備実験検証メモ_v1.md").write_text(
        report,
        encoding="utf-8",
    )
    return safe


if __name__ == "__main__":
    print(json.dumps(run()["aggregate"], ensure_ascii=False, indent=2))
