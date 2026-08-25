#!/usr/bin/env python3
"""Stage A2a固定報告書を生成し、後続Stageへ進まず停止する。"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
RAW = HERE / "raw"
PROCESSED = HERE / "processed"
FIGURES = HERE / "figures"
REPORTS = HERE / "reports"
LOGS = HERE / "logs"
CONFIG = json.loads((HERE / "config_locked.json").read_text(encoding="utf-8"))
EXPECTED = json.loads((HERE / "expected_hashes.json").read_text(encoding="utf-8"))
RUN1, RUN2 = CONFIG["run_ids"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def md_table(columns: list[str], data: list[dict], labels: dict[str, str] | None = None) -> str:
    labels = labels or {}
    out = [
        "| " + " | ".join(labels.get(c, c) for c in columns) + " |",
        "|" + "|".join(["---"] * len(columns)) + "|",
    ]
    for row in data:
        out.append("| " + " | ".join(str(row.get(c, "")) for c in columns) + " |")
    return "\n".join(out)


def top_adjacent_changes(by_level: list[dict], field: str, count: int = 3) -> list[tuple[str, str, float]]:
    valid = [r for r in by_level if r.get(f"seedless_{field}_exact", "") != ""]
    changes = []
    for a, b in zip(valid, valid[1:]):
        va = float(a[f"seedless_{field}_exact"])
        vb = float(b[f"seedless_{field}_exact"])
        changes.append((a["level_label"], b["level_label"], abs(vb - va)))
    return sorted(changes, key=lambda x: x[2], reverse=True)[:count]


def main() -> None:
    required_gates = [
        LOGS / "source_verification.json",
        LOGS / "execution_manifest.json",
        PROCESSED / "exec_comparison_summary.json",
        PROCESSED / "seeded_comparison_summary.json",
        FIGURES / "figure_manifest.json",
    ]
    if not all(p.is_file() for p in required_gates):
        raise SystemExit("EXECUTION_FAILED: 報告書作成に必要な前工程記録が不足")
    source_gate = json.loads(required_gates[0].read_text(encoding="utf-8"))
    execution = json.loads(required_gates[1].read_text(encoding="utf-8"))
    exec_compare = json.loads(required_gates[2].read_text(encoding="utf-8"))
    seeded_compare = json.loads(required_gates[3].read_text(encoding="utf-8"))
    figure_manifest = json.loads(required_gates[4].read_text(encoding="utf-8"))

    complete = (
        source_gate.get("status") == "VERIFIED"
        and execution.get("status") == "COMPLETED"
        and exec_compare.get("status") == "COMPARED"
        and exec_compare.get("execs_bitwise_identical") is True
        and exec_compare.get("numerical_health_passed") is True
        and seeded_compare.get("status") == "COMPARISON_TABLES_COMPLETE"
        and figure_manifest.get("status") == "FIGURES_COMPLETE"
        and figure_manifest.get("figure_count") == 14
    )
    status = "SEEDLESS_COMPARISON_COMPLETE" if complete else "SEEDLESS_COMPARISON_INCOMPLETE"

    summary1 = json.loads((RAW / RUN1 / "run_summary.json").read_text(encoding="utf-8"))
    summary2 = json.loads((RAW / RUN2 / "run_summary.json").read_text(encoding="utf-8"))
    first_passage = rows(PROCESSED / "seedless_first_passage_levels.csv")
    growth = rows(PROCESSED / "seeded_vs_seedless_growth_rate.csv")
    by_level = rows(PROCESSED / "seeded_vs_seedless_by_f_level.csv")
    health = rows(PROCESSED / "numerical_health.csv")
    source_rows = []
    for group in ("sources", "dependencies"):
        for name, item in EXPECTED[group].items():
            path = REPO / item["path"]
            source_rows.append({
                "file": name,
                "role": item["role"],
                "sha256": sha256(path),
                "absolute_path": str(path),
            })

    found_steps = [int(r["first_passage_step"]) for r in first_passage if r["status"] == "found"]
    rate_values = np.asarray([
        float(r["seedless_mean_exponential_rate_per_step"])
        for r in growth if r["seedless_mean_exponential_rate_per_step"]
    ])
    top_change_text = []
    for field, label in (
        ("q3_over_q1", "q3/q1"),
        ("q4_over_q1", "q4/q1"),
        ("direction_3_occupation", "direction 3 occupation"),
        ("direction_4_occupation", "direction 4 occupation"),
        ("kernel_occupation", "kernel occupation"),
    ):
        parts = [
            f"{a}→{b}: |Δ|={delta:.6e}"
            for a, b, delta in top_adjacent_changes(by_level, field)
        ]
        top_change_text.append(f"- {label}: " + ("; ".join(parts) if parts else "比較可能な連続水準なし"))

    source_table = md_table(["file", "role", "sha256", "absolute_path"], source_rows)
    first_table = md_table(
        ["level_label", "status", "first_passage_step", "f_at_first_passage",
         "q3_over_q1", "q4_over_q1", "direction_3_occupation",
         "direction_4_occupation", "kernel_occupation"],
        first_passage,
    )
    growth_table = md_table(
        ["lower_level", "upper_level", "seeded_step_difference", "seedless_step_difference",
         "seeded_mean_exponential_rate_per_step", "seedless_mean_exponential_rate_per_step"],
        growth,
    )
    commands = "\n".join(f"{i}. `python3 {name}`" for i, name in enumerate(CONFIG["execution_order"], 1))
    elapsed_lines = "\n".join(
        f"- `{r['run_id']}`: {r['elapsed_seconds']:.6f} s" for r in execution["runs"]
    )
    figure_lines = "\n".join(
        f"- `{name}`" for name in figure_manifest["files"] if name.endswith(".png")
    )
    report = f"""# Stage A2a seedless N=5 報告書

## 総合状態

**{status}**

対象は `N=5`、`float64`、`step=0..5000` のみである。高精度親、Delta掃引、`N=40`、`N=300`、Stage B/Cは実行していない。H1/H2/H0、下位階層の存在、増幅停止、方向成立、上下位の同型性は判定していない。

## 固定原本とSHA-256

{source_table}

全原本、Stage A0入力、Stage A1b入力は `verify_sources.py` で固定SHA-256と照合し、`{source_gate['status']}` となった。原本は編集していない。

## 実行仕様と初期状態

- 親: `make_parent(LowRankSystem(5), np.random.default_rng(40265722), iters=1200, tol=1e-12)`
- exec 1 親残差: `{summary1['parent_residual']:.17e}`
- exec 2 親残差: `{summary2['parent_residual']:.17e}`
- 初期状態: `Z0 = v.copy()`。親とのbitwise一致: `{summary1['initial_state_bitwise_equal_parent']}`
- 明示状態seed追加: `{summary1['explicit_state_seed_added']}`
- `zero_closure_kernel_seed` 呼出し: `{summary1['zero_closure_kernel_seed_called']}`
- delta、kick、noise、量子化、retraction、高精度化: なし
- warm-start: 親生成後の同一PRNGから `rng.normal` を生成し、既存 `sigma_max_power` のみに使用した。状態には加えていない。
- f・基本誤差: 毎step、q: 5 stepごと、占有: 25 stepごと。指定f水準初回通過stepではq/占有を追加実測した。
- 測定値・支配平面は状態更新およびwarm-startへフィードバックしていない。
- `norm_error = |Z†Z-1|`、`closure_error = |Z^T Z|`、`conservation_error = |Z†Z-(Z0†Z0)|`。

## Python環境と実行

- Python: `{sys.version.replace(chr(10), ' ')}`
- NumPy: `{np.__version__}`
- OS: `{platform.platform()}`

{commands}

実行時間:

{elapsed_lines}

## exec 1/2再現確認

- bitwise一致: `{exec_compare['execs_bitwise_identical']}`
- 数値健全性: `{exec_compare['numerical_health_passed']}`
- exec 1 f行数/q行数/占有行数: `{summary1['f_row_count']}` / `{summary1['q_row_count']}` / `{summary1['occupation_row_count']}`
- exec 2 f行数/q行数/占有行数: `{summary2['f_row_count']}` / `{summary2['q_row_count']}` / `{summary2['occupation_row_count']}`
- exec 1最大 norm/closure/conservation error: `{summary1['maximum_norm_error']:.6e}` / `{summary1['maximum_closure_error']:.6e}` / `{summary1['maximum_conservation_error']:.6e}`
- 最大projection closure error: `{summary1['maximum_projection_closure_error']:.6e}`

詳細: [exec1_vs_exec2.md](../processed/exec1_vs_exec2.md)、[numerical_health.md](../processed/numerical_health.md)

## seedなし初期fと全f初回到達

- seedなし初期f: `{summary1['initial_f']:.17e}`
- seedあり基準初期f: `{seeded_compare['seeded_initial_f']:.17e}`
- seedなし到達水準数: `{seeded_compare['seedless_levels_found']}/31`
- 見つかった初回到達step範囲: `{min(found_steps) if found_steps else 'not_found'}..{max(found_steps) if found_steps else 'not_found'}`

{first_table}

CSV/Markdown正本: [seedless_first_passage_levels.csv](../processed/seedless_first_passage_levels.csv)、[seedless_first_passage_levels.md](../processed/seedless_first_passage_levels.md)

## decade間step差と平均指数率

平均指数率は、指定された隣接f座標間の `ln(level_upper/level_lower) / Δstep` であり、イベント定義ではない。

{growth_table}

seedなしで計算可能な全隣接区間の率分布: minimum=`{float(np.min(rate_values)) if rate_values.size else float('nan'):.6e}`、median=`{float(np.median(rate_values)) if rate_values.size else float('nan'):.6e}`、maximum=`{float(np.max(rate_values)) if rate_values.size else float('nan'):.6e}`。

詳細: [seedless_decade_growth_rates.md](../processed/seedless_decade_growth_rates.md)、[seeded_vs_seedless_growth_rate.md](../processed/seeded_vs_seedless_growth_rate.md)

## seedあり基準との時間ずれ

表示規則 `f>=1e-12` の初回到達を0とした。

- seedあり: absolute step `{seeded_compare['seeded_alignment_step']}`
- seedなし: absolute step `{seeded_compare['seedless_alignment_step']}`
- seedなし − seedあり: `{seeded_compare['seedless_minus_seeded_alignment_step']}` step

これは表示用平行移動であり、物理イベント時刻の採用ではない。

## 時間平行移動後のf曲線差

- 共通relative-step範囲: `{seeded_compare['aligned_relative_step_range']}`
- 最大 `|Δf|`: `{seeded_compare['aligned_max_absolute_f_difference']:.6e}`
- 最大 `|Δlog10(f)|`: `{seeded_compare['aligned_max_absolute_log10_f_difference']:.6e}`
- `Δlog10(f)` RMSE: `{seeded_compare['aligned_log10_f_rmse']:.6e}`

正本: [time_aligned_f_comparison.csv](../processed/time_aligned_f_comparison.csv)

## f水準ごとのq3/q4・direction 3/4・kernel比較

seedなし側は各f水準初回到達stepでの追加実測値である。seedあり側はStage A1bの直前/直後の実保存レコードを両方保持した。補間は使用していないため、差も「seedなし実測 − seedあり直前実測」と「seedなし実測 − seedあり直後実測」の2本を保存した。

正本: [seeded_vs_seedless_by_f_level.csv](../processed/seeded_vs_seedless_by_f_level.csv)、[seeded_vs_seedless_by_f_level.md](../processed/seeded_vs_seedless_by_f_level.md)

## 単一の指数系列として説明可能な範囲

このStageでは「率が同一とみなせる許容幅」が固定されていないため、単一指数系列の範囲を自動採用しない。上の全隣接f座標の率、図5、およびCSVが、人間による範囲固定の資料である。したがって、データから直接提示できるのは各座標対の率とその変動までであり、単一指数系列という物理解釈の採用ではない。

## 増幅率または方向占有の構造変化が大きい固定f座標間

以下は供給された隣接f水準間の絶対差を機械的に降順表示した上位3組であり、閾値・イベント・方向成立の判定ではない。

{chr(10).join(top_change_text)}

## データから直接言えること

- 2実行は保存CSVと支配平面NPYを含めbitwise一致した。
- `Z0=v` のfloat64自然床から、指定31水準のうち `{seeded_compare['seedless_levels_found']}` 水準について5000 step内の初回到達と同stepのq/占有を実測できた。
- seedあり・seedなしの絶対step差、`f>=1e-12` 表示平行移動後の曲線差、各f水準でのq3/q4・direction 3/4・kernelの差が保存された。
- `rank_q`、q3/q4の有限値、direction 3/4占有、f増幅は別々の列として保存され、同一イベントとは扱っていない。

## データだけでは言えないこと

- 下位階層の存在、増幅停止時刻、方向成立時刻
- H1/H2/H0のいずれが正しいか
- 上位遷移と下位遷移の同型性
- seedあり参照の未保存stepにおけるq/占有の値（補間も再実行もしていない）
- 率変動のどこを単一指数系列の境界とするか

## 必須表

- [seedless_first_passage_levels](../processed/seedless_first_passage_levels.md)
- [seedless_decade_growth_rates](../processed/seedless_decade_growth_rates.md)
- [seeded_vs_seedless_by_f_level](../processed/seeded_vs_seedless_by_f_level.md)
- [seeded_vs_seedless_growth_rate](../processed/seeded_vs_seedless_growth_rate.md)
- [exec1_vs_exec2](../processed/exec1_vs_exec2.md)
- [numerical_health](../processed/numerical_health.md)

## 必須図

{figure_lines}

## 最終停止

Stage A2aの報告書作成をもって停止する。高精度親、Delta掃引、`N=40`、`N=300`、Stage B/Cへ進まない。
"""
    REPORTS.mkdir(parents=True, exist_ok=True)
    out = REPORTS / "stage_A2a_seedless_N5_report.md"
    if out.exists():
        raise SystemExit("EXECUTION_FAILED: 既存報告書を上書きしない")
    out.write_text(report, encoding="utf-8")
    print(status)
    print(out)


if __name__ == "__main__":
    main()
