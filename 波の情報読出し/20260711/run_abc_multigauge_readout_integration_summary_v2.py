from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, List


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "abc_multigauge_readout_integration_summary_result_v2"
OUT_DIR.mkdir(exist_ok=True)


EXPERIMENTS: List[Dict[str, str]] = [
    {
        "name": "single_collision_multigauge_readout",
        "purpose": "単回ABC衝突で p/E/R を多ゲージ干渉読出しする",
        "path": "abc_multigauge_interference_readout_result_v2/abc_multigauge_interference_readout_result_v2.json",
        "section": "verdicts",
        "valid_key": "multigauge_measurement_valid",
        "case_key": "",
    },
    {
        "name": "multi_collision_multigauge_readout",
        "purpose": "対称ABC衝突の反復で p/E/R 読出しを維持する",
        "path": "abc_multigauge_interference_readout_multi_collision_result_v2/abc_multigauge_interference_readout_multi_collision_result_v2.json",
        "section": "verdicts",
        "valid_key": "multi_collision_multigauge_valid",
        "case_key": "ab_collision_count",
    },
    {
        "name": "readout_robustness_sweep",
        "purpose": "複数の読出し器構成で p/E/R 再構成が安定する",
        "path": "abc_multigauge_interference_readout_robustness_sweep_result_v2/abc_multigauge_interference_readout_robustness_sweep_result_v2.json",
        "section": "aggregate_verdict",
        "valid_key": "robustness_sweep_valid",
        "case_key": "case_count",
    },
    {
        "name": "asymmetric_amplitude_diagnostic",
        "purpose": "非対称Rで単純反転が保存を破ることを検出する",
        "path": "abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v2/abc_multigauge_interference_readout_asymmetric_amplitude_sweep_result_v2.json",
        "section": "aggregate_verdict",
        "valid_key": "asymmetric_amplitude_diagnostic_valid",
        "case_key": "case_count",
    },
    {
        "name": "generalized_elastic_collision_readout",
        "purpose": "非対称Rで R*p と R*p^2 を保存する一般化写像を読む",
        "path": "abc_multigauge_generalized_elastic_collision_readout_result_v2/abc_multigauge_generalized_elastic_collision_readout_result_v2.json",
        "section": "aggregate_verdict",
        "valid_key": "generalized_elastic_collision_readout_valid",
        "case_key": "case_count",
    },
    {
        "name": "generalized_velocity_sweep",
        "purpose": "非単位・非対称位相勾配でも一般化写像が成立する",
        "path": "abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v2/abc_multigauge_generalized_elastic_collision_velocity_sweep_result_v2.json",
        "section": "aggregate_verdict",
        "valid_key": "velocity_sweep_generalized_collision_valid",
        "case_key": "case_count",
    },
    {
        "name": "generalized_multi_collision",
        "purpose": "一般化写像を複数回AB衝突へ反復適用する",
        "path": "abc_multigauge_generalized_elastic_collision_multi_collision_result_v2/abc_multigauge_generalized_elastic_collision_multi_collision_result_v2.json",
        "section": "aggregate_verdict",
        "valid_key": "generalized_multi_collision_valid",
        "case_key": "case_count",
    },
    {
        "name": "generalized_noise_robustness",
        "purpose": "ゼロ平均読出しノイズの相殺と共通バイアス検出を確認する",
        "path": "abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2/abc_multigauge_generalized_elastic_collision_noise_robustness_result_v2.json",
        "section": "aggregate_verdict",
        "valid_key": "noise_robustness_valid",
        "case_key": "case_count",
    },
    {
        "name": "generalized_extreme_R_sweep",
        "purpose": "極端なR比でも一般化写像と読出しが維持されるか調べる",
        "path": "abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2/abc_multigauge_generalized_elastic_collision_extreme_R_sweep_result_v2.json",
        "section": "aggregate_verdict",
        "valid_key": "extreme_R_sweep_valid",
        "case_key": "case_count",
    },
]


def load_experiment(row: Dict[str, str]) -> Dict[str, Any]:
    path = BASE_DIR / row["path"]
    data = json.loads(path.read_text(encoding="utf-8"))
    section = data[row["section"]]
    case_count = section.get(row["case_key"], "") if row["case_key"] else ""
    return {
        "name": row["name"],
        "purpose": row["purpose"],
        "path": row["path"],
        "section": row["section"],
        "valid_key": row["valid_key"],
        "valid": bool(section[row["valid_key"]]),
        "case_count_or_collision_count": case_count,
        "single_gauge_only_used": bool(section.get("single_gauge_only_used", False)),
        "max_p_error": section.get("max_p_abs_error", section.get("max_p_abs_error_all_cases", "")),
        "max_E_error": section.get("max_E_abs_error", section.get("max_E_abs_error_all_cases", "")),
        "max_R_error": section.get("max_R_abs_error", section.get("max_R_abs_error_all_cases", "")),
        "raw_verdict": section,
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    keys = [
        "name",
        "purpose",
        "path",
        "valid_key",
        "valid",
        "case_count_or_collision_count",
        "single_gauge_only_used",
        "max_p_error",
        "max_E_error",
        "max_R_error",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in keys})


def write_report(result: Dict[str, Any]) -> None:
    lines = [
        "# ABC Multigauge Readout Integration Summary v1",
        "",
        "## Purpose",
        "",
        "This summary collects the executed ABC multigauge readout experiments for the 2026-07-11 series.",
        "It is an index of executed numerical evidence, not an additional physical assumption.",
        "",
        "## Aggregate Verdict",
        "",
    ]
    for key, value in result["aggregate_verdict"].items():
        lines.append(f"- {key}: `{value}`")
    lines.extend(
        [
            "",
            "## Experiment Table",
            "",
            "| experiment | purpose | count | valid | single gauge only |",
            "|---|---|---:|---|---|",
        ]
    )
    for row in result["summary_rows"]:
        lines.append(
            f"| {row['name']} | {row['purpose']} | {row['case_count_or_collision_count']} | "
            f"`{row['valid']}` | `{row['single_gauge_only_used']}` |"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "| kind | file |",
            "|---|---|",
            "| JSON | `abc_multigauge_readout_integration_summary_result_v2.json` |",
            "| CSV | `abc_multigauge_readout_integration_summary_v2.csv` |",
        ]
    )
    (OUT_DIR / "abc_multigauge_readout_integration_summary_report_v2.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def run() -> Dict[str, Any]:
    rows = [load_experiment(row) for row in EXPERIMENTS]
    return {
        "experiment": "abc_multigauge_readout_integration_summary_v2",
        "summary_rows": rows,
        "aggregate_verdict": {
            "experiment_count": len(rows),
            "all_experiments_valid": all(bool(row["valid"]) for row in rows),
            "single_gauge_only_used_any": any(bool(row["single_gauge_only_used"]) for row in rows),
            "integration_summary_valid": all(bool(row["valid"]) for row in rows)
            and not any(bool(row["single_gauge_only_used"]) for row in rows),
        },
    }


def write_outputs(result: Dict[str, Any]) -> None:
    (OUT_DIR / "abc_multigauge_readout_integration_summary_result_v2.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(OUT_DIR / "abc_multigauge_readout_integration_summary_v2.csv", result["summary_rows"])
    write_report(result)


def main() -> None:
    result = run()
    write_outputs(result)
    print(json.dumps(result["aggregate_verdict"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
