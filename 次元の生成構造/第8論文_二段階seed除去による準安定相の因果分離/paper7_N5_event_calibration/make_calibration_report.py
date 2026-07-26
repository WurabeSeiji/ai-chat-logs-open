#!/usr/bin/env python3
"""Create the Stage A1 calibration report without adopting a single event."""

from __future__ import annotations

import csv
import json
import math
import platform
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parent
CONFIG_PATH = PACKAGE_ROOT / "config_candidates.json"
VERIFY_PATH = PACKAGE_ROOT / "logs" / "input_verification.json"
GROWTH_MANIFEST_PATH = PACKAGE_ROOT / "logs" / "growth_analysis_manifest.json"
RANK_MANIFEST_PATH = PACKAGE_ROOT / "logs" / "rank_analysis_manifest.json"
FIGURE_MANIFEST_PATH = PACKAGE_ROOT / "logs" / "figures_manifest.json"
PROCESSED_DIR = PACKAGE_ROOT / "processed"
REPORT_DIR = PACKAGE_ROOT / "reports"
REPORT_PATH = REPORT_DIR / "paper7_N5_event_calibration_report.md"
TEXT_LOG_PATH = PACKAGE_ROOT / "logs" / "make_calibration_report.log"


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def distribution(values: list[int | float]) -> dict:
    if not values:
        return {
            "count": 0,
            "minimum": None,
            "q25": None,
            "median": None,
            "q75": None,
            "maximum": None,
            "unique": 0,
        }
    array = np.asarray(values, dtype=float)
    return {
        "count": len(values),
        "minimum": float(np.min(array)),
        "q25": float(np.quantile(array, 0.25)),
        "median": float(np.median(array)),
        "q75": float(np.quantile(array, 0.75)),
        "maximum": float(np.max(array)),
        "unique": len(set(float(value) for value in values)),
    }


def fmt(value) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        if value.is_integer():
            return str(int(value))
        return f"{value:.10g}"
    return str(value)


def distribution_row(name: str, stats: dict) -> str:
    return (
        f"| {name} | {stats['count']} | {fmt(stats['minimum'])} | "
        f"{fmt(stats['q25'])} | {fmt(stats['median'])} | {fmt(stats['q75'])} | "
        f"{fmt(stats['maximum'])} | {stats['unique']} |"
    )


def file_size(path: Path) -> str:
    size = path.stat().st_size
    units = ["B", "KiB", "MiB", "GiB"]
    value = float(size)
    unit = units[0]
    for unit in units:
        if value < 1024 or unit == units[-1]:
            break
        value /= 1024
    return f"{value:.2f} {unit}"


def main() -> int:
    started = time.perf_counter()
    if REPORT_PATH.exists() or TEXT_LOG_PATH.exists():
        raise RuntimeError("校正報告書の上書きを避けて停止")
    required_manifests = [
        VERIFY_PATH,
        GROWTH_MANIFEST_PATH,
        RANK_MANIFEST_PATH,
        FIGURE_MANIFEST_PATH,
    ]
    if any(not path.is_file() for path in required_manifests):
        raise RuntimeError("Stage A1先行工程の記録が欠落")

    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    verification = json.loads(VERIFY_PATH.read_text(encoding="utf-8"))
    growth_manifest = json.loads(GROWTH_MANIFEST_PATH.read_text(encoding="utf-8"))
    rank_manifest = json.loads(RANK_MANIFEST_PATH.read_text(encoding="utf-8"))
    figure_manifest = json.loads(FIGURE_MANIFEST_PATH.read_text(encoding="utf-8"))

    table_stems = [
        "growth_intervals_all_candidates",
        "growth_end_all_candidates",
        "rank4_onset_all_candidates",
        "growth_end_vs_rank4_onset_all_pairs",
        "candidate_summary",
    ]
    required_tables = [
        PROCESSED_DIR / f"{stem}.{extension}"
        for stem in table_stems
        for extension in ("csv", "md")
    ]
    all_inputs_ok = verification.get("success") is True
    all_analysis_ok = (
        growth_manifest.get("success") is True
        and rank_manifest.get("success") is True
    )
    all_tables_present = all(path.is_file() for path in required_tables)
    all_figures_present = (
        figure_manifest.get("success") is True
        and figure_manifest.get("figure_count_by_format") == {"png": 10, "svg": 10}
        and all(Path(path).is_file() for path in figure_manifest["generated_files"])
    )
    if not all_inputs_ok:
        overall = "INPUT_MISMATCH"
    elif all_analysis_ok and all_tables_present and all_figures_present:
        overall = "CALIBRATION_DATA_COMPLETE"
    else:
        overall = "CALIBRATION_DATA_INCOMPLETE"

    intervals = load_csv(PROCESSED_DIR / "growth_intervals_all_candidates.csv")
    growth_ends = load_csv(PROCESSED_DIR / "growth_end_all_candidates.csv")
    rank_onsets = load_csv(PROCESSED_DIR / "rank4_onset_all_candidates.csv")
    found_ends = [row for row in growth_ends if row["status"] == "found"]
    found_onsets = [row for row in rank_onsets if row["status"] == "found"]

    interval_starts = [int(row["interval_start"]) for row in intervals]
    interval_ends = [int(row["interval_end"]) for row in intervals]
    growth_end_steps = [int(row["growth_end_candidate"]) for row in found_ends]
    rank_onset_steps = [int(row["rank4_onset_candidate"]) for row in found_onsets]
    pair_differences = [
        onset - ending for ending in growth_end_steps for onset in rank_onset_steps
    ]

    start_dist = distribution(interval_starts)
    interval_end_dist = distribution(interval_ends)
    growth_end_dist = distribution(growth_end_steps)
    rank_onset_dist = distribution(rank_onset_steps)
    pair_dist = distribution(pair_differences)

    input_checks = verification["input_checks"]
    time_axis = verification["time_axis"]
    lines = [
        "# 第7論文 N=5 イベント定義校正報告（Stage A1）",
        "",
        "## 1. 範囲と総合状態",
        "",
        "- 実施内容: Stage A0で完全再現されたN=5 CSV 3件の後処理だけ",
        "- 原本力学コードの実行: なし",
        "- 新しい軌道生成: なし",
        "- 第8論文本実験、高精度親、Series 1〜3、N=40、N=300: 未実行",
        "- 単一イベント時刻の採用: なし",
        "- 仮説H1/H2/H0の判定: なし",
        f"- **総合状態: `{overall}`**",
        "",
        "本報告は候補集合の校正資料であり、指数成長開始、指数成長終了、rank_Q=4持続開始のいずれについても最終時刻を決定しない。",
        "",
        "## 2. 入力SHA-256",
        "",
        "| 入力 | 行数 | 絶対パス | SHA-256 | Stage A0照合 |",
        "|---|---:|---|---|---|",
    ]
    for name in ("fcurve", "q_svd", "paper7_long"):
        check = input_checks[name]
        lines.append(
            f"| `{name}` | {check['row_count']} | `{check['path']}` | "
            f"`{check['actual_sha256']}` | "
            f"{'一致' if check['success'] else '不一致'} |"
        )

    lines.extend(
        [
            "",
            "照合はStage A0報告書の `REPRODUCED_EXACTLY` とCSV別 `exact` 行、Stage A0 CSV比較記録の `reproduced_sha256`、実ファイルSHA-256の四者で行った。",
            "",
            "## 3. 使用列と時間軸",
            "",
            "- `fcurve`: `tau`, `f`",
            "- `q_svd`: `step`, `time`, `relative_time`, `q1`, `q2`, `q3`, `q4`, `rank_q`",
            "- `paper7_long_timeseries`: `step`, `time`, `crossing_flag`, `splitting_fraction`",
            "",
            f"- 共通軸: `{time_axis['shared_axis']}`",
            f"- f範囲: `{time_axis['f_step_range']}`、1 step連続: `{time_axis['f_contiguous_unit_steps']}`",
            f"- q範囲: `{time_axis['q_step_range']}`、`step-relative_time`: `{time_axis['q_step_minus_relative_time_unique']}`",
            f"- 長時間CSV範囲: `{time_axis['paper7_step_range']}`",
            f"- 既存crossing: `{config['existing_crossing']}`",
            "- qの未保存stepは補間していない。rank持続長は連続する保存レコード数であり、実step spanも候補表に併記した。",
            "",
            "## 4. 全候補パラメータ",
            "",
            f"- 回帰窓: `{config['regression_windows']}`",
            f"- R²閾値: `{config['r2_thresholds']}`",
            f"- 成長区間最小持続長: `{config['growth_minimum_durations']}`",
            f"- 成長終了条件: `{config['growth_end_conditions']}`",
            f"- 成長終了持続長: `{config['growth_end_persistence']}`",
            f"- rank相対閾値: `{config['rank_relative_thresholds']}`",
            f"- rank=4持続長（保存レコード数）: `{config['rank_persistence_records']}`",
            "",
            "回帰対象は自然対数 `log_f` である。中心窓の通常最小二乗を用い、残差標準偏差は `sqrt(SSE/(window-2))` とした。`f<=0` は正数へ置換せず、対数不能点はNaNとした。",
            "",
            "## 5. 各窓の成長区間候補数",
            "",
            "| window | 候補行数 |",
            "|---:|---:|",
        ]
    )
    for window in config["regression_windows"]:
        lines.append(
            f"| {window} | {growth_manifest['growth_interval_counts_by_window'][str(window)]} |"
        )
    lines.extend(
        [
            f"| **全窓** | **{growth_manifest['growth_interval_candidate_count']}** |",
            "",
            "同じ最大連続区間が複数のR²閾値・最小持続長を満たす場合も、各パラメータ候補行を削除していない。",
            "",
            "## 6. 候補分布",
            "",
            "| 候補集合 | 件数 | 最小step | Q25 | 中央値 | Q75 | 最大step | 異なるstep数 |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
            distribution_row("成長区間開始", start_dist),
            distribution_row("成長区間終了端", interval_end_dist),
            distribution_row("成長終了候補（found）", growth_end_dist),
            distribution_row("rank=4持続開始候補（found）", rank_onset_dist),
            distribution_row("rank開始 - 成長終了（全組合せ）", pair_dist),
            "",
            f"- 成長終了候補表全行: `{growth_manifest['growth_end_row_count']}`",
            f"- 成長終了 `found`: `{growth_manifest['growth_end_found_count']}`",
            f"- 成長終了 `not_found`: `{growth_manifest['growth_end_not_found_count']}`",
            f"- rank候補表全行: `{rank_manifest['rank4_onset_row_count']}`",
            f"- rank候補 `found`: `{rank_manifest['rank4_onset_found_count']}`",
            f"- 成長終了×rank開始の全ペア: `{rank_manifest['all_pair_count']}`",
            "",
            "## 7. パラメータに対して反復する候補群",
            "",
            "以下は、異なる候補パラメータで同じstepが何回現れたかの頻度である。頻度順は採用順位ではない。",
            "",
        ]
    )
    recurring_groups = [
        ("成長区間開始", interval_starts),
        ("成長終了", growth_end_steps),
        ("rank=4持続開始", rank_onset_steps),
    ]
    for label, values in recurring_groups:
        lines.extend(
            [
                f"### {label}",
                "",
                "| 頻度順位 | step | 出現候補行数 |",
                "|---:|---:|---:|",
            ]
        )
        for rank, (candidate_step, count) in enumerate(
            sorted(Counter(values).items(), key=lambda item: (-item[1], item[0]))[:10],
            1,
        ):
            lines.append(f"| {rank} | {candidate_step} | {count} |")
        lines.append("")

    lines.extend(
        [
            "これらは「パラメータに対して安定して反復する候補群」を抽出するための記述統計に限られ、どのstepも採用していない。",
            "",
            "## 8. パラメータに強く依存する候補群",
            "",
            "| パラメータ群 | 件数 | 最小 | 中央値 | 最大 | 異なるstep数 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for window in config["regression_windows"]:
        values = [
            int(row["interval_start"])
            for row in intervals
            if int(row["window"]) == int(window)
        ]
        stats = distribution(values)
        lines.append(
            f"| 成長開始 window={window} | {stats['count']} | {fmt(stats['minimum'])} | "
            f"{fmt(stats['median'])} | {fmt(stats['maximum'])} | {stats['unique']} |"
        )
    for condition_name in config["growth_end_conditions"]:
        values = [
            int(row["growth_end_candidate"])
            for row in found_ends
            if row["end_condition"] == condition_name
        ]
        stats = distribution(values)
        lines.append(
            f"| 成長終了 condition={condition_name} | {stats['count']} | "
            f"{fmt(stats['minimum'])} | {fmt(stats['median'])} | "
            f"{fmt(stats['maximum'])} | {stats['unique']} |"
        )
    for threshold in config["rank_relative_thresholds"]:
        values = [
            int(row["rank4_onset_candidate"])
            for row in found_onsets
            if float(row["relative_threshold"]) == float(threshold)
        ]
        stats = distribution(values)
        lines.append(
            f"| rank開始 threshold={threshold:.0e} | {stats['count']} | "
            f"{fmt(stats['minimum'])} | {fmt(stats['median'])} | "
            f"{fmt(stats['maximum'])} | {stats['unique']} |"
        )
    lines.extend(
        [
            "",
            "候補範囲・異なるstep数が広い群は、窓、R²、区間持続長、終了条件、終了持続長、rank閾値、rank保存レコード持続長への依存を人間が確認すべき群である。ここでも許容範囲や代表値は設定していない。",
            "",
            "## 9. 解析不能または曖昧な範囲",
            "",
            "- 各回帰窓の半窓より端側では中心回帰を定義できず、該当列はNaNである。",
            f"- `fcurve` はstep {time_axis['f_step_range'][1]}で終了するため、それ以後の成長slope・成長終了候補は解析できない。",
            "- q系列は可変間隔で保存されている。連続保存レコード間の未観測stepでrank=4が維持されたかは、この入力だけでは判定できない。",
            "- rank持続長を「保存レコード数」ではなく「物理step数」とする最終仕様は未決定である。",
            "- `paper7_long_timeseries` は25 step間隔であり、crossing=1167そのものは保存点ではない。保存点のflag整合だけを検査した。",
            "- 成長区間候補が複数ある系列について、どの区間を主区間と呼ぶかは未決定である。",
            "- 成長終了候補が見つからないパラメータ組合せも削除せず `not_found` として保存した。",
            "- 成長終了候補とrank開始候補の対応付け規則がないため、全直積だけを保存した。",
            "",
            "## 10. q比と既存rank定義の検証",
            "",
            "- 既存定義: `rank_Q = count(q_j > 1e-8 q1), j=1,...,4`",
            f"- 既存`rank_q`との不一致行数: `{rank_manifest['existing_rank_mismatch_count']}`",
            "- 比較用rank閾値は頑健性確認だけであり、既存定義を変更していない。",
            "- `q3/q1`, `q4/q1`, `min(q3,q4)/q1`, `q3-q4`, `q3/q4` を保存した。",
            "- ゼロ除算はNaNであり、比に固定イベント閾値を置いていない。",
            "",
            "## 11. 出力表",
            "",
            "| 表 | CSV | Markdown |",
            "|---|---|---|",
        ]
    )
    for stem in table_stems:
        csv_path = PROCESSED_DIR / f"{stem}.csv"
        md_path = PROCESSED_DIR / f"{stem}.md"
        lines.append(
            f"| `{stem}` | `{csv_path}` ({file_size(csv_path)}) | "
            f"`{md_path}` ({file_size(md_path)}) |"
        )

    lines.extend(
        [
            "",
            "候補数が多い表も削除・間引きしていない。全ペア表は候補IDを使って正規化し、個別パラメータは成長終了表とrank開始表から参照できる。",
            "",
            "## 12. 出力図",
            "",
        ]
    )
    for path in figure_manifest["generated_files"]:
        lines.append(f"- `{path}`")

    lines.extend(
        [
            "",
            "全図は候補を同等に表示し、採用候補の色・線幅・注釈による強調を行っていない。",
            "",
            "## 13. 人間が決定すべき事項",
            "",
            "1. 指数成長回帰に採用する窓幅、R²閾値、最小持続長。",
            "2. 複数の連続成長区間がある場合の対象区間。",
            "3. 成長終了条件A/B/Cと終了持続長。",
            "4. `rank_Q=4` 頑健性確認で重視する相対閾値。",
            "5. rank持続を保存レコード数で定義するか、別の毎stepデータを要求するか。",
            "6. 成長終了候補とrank開始候補を対応付ける規則。",
            "7. 既存crossing=1167と新たに校正するイベント群の関係。",
            "",
            "上記事項が承認されるまで、単一のイベント時刻を固定してはならない。",
            "",
            "## 14. 実行環境と停止",
            "",
            f"- Python: `{sys.version}`",
            f"- NumPy: `{np.__version__}`",
            f"- OS: `{platform.platform()}`",
            f"- 報告生成日時（UTC）: `{datetime.now(timezone.utc).isoformat()}`",
            f"- 報告生成時間（秒）: `{time.perf_counter() - started:.6f}`",
            f"- **最終状態: `{overall}`**",
            "- Stage A1はここで停止する。Stage A2へ進まない。",
        ]
    )

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    TEXT_LOG_PATH.write_text(
        f"report={REPORT_PATH}\noverall={overall}\nseconds={time.perf_counter() - started:.6f}\n",
        encoding="utf-8",
    )
    print(f"[REPORT] {REPORT_PATH}")
    print(f"[OVERALL] {overall}")
    return 0 if overall == "CALIBRATION_DATA_COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
