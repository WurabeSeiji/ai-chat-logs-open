#!/usr/bin/env python3
"""Create the Stage A1b transition-anatomy report without selecting events."""

from __future__ import annotations

import bisect
import csv
import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


PACKAGE_ROOT = Path(__file__).resolve().parent
PROCESSED_DIR = PACKAGE_ROOT / "processed"
LOG_DIR = PACKAGE_ROOT / "logs"
REPORT_DIR = PACKAGE_ROOT / "reports"
REPORT_PATH = REPORT_DIR / "paper7_N5_transition_anatomy_report.md"
TEXT_LOG_PATH = LOG_DIR / "make_transition_report.log"
VERIFY_PATH = LOG_DIR / "input_verification.json"
ANALYSIS_MANIFEST_PATH = LOG_DIR / "analysis_manifest.json"
FIGURE_MANIFEST_PATH = LOG_DIR / "figures_manifest.json"
CROSSING = 1167
FIXED_WINDOWS = [
    ("0-500", 0, 500),
    ("500-1000", 500, 1000),
    ("800-1400", 800, 1400),
    ("1000-1800", 1000, 1800),
    ("1400-2500", 1400, 2500),
]
REQUIRED_TABLE_STEMS = [
    "f_first_passage_levels",
    "f_decade_growth_rates",
    "first_passage_nearest_occupation_records",
    "first_passage_nearest_q_records",
    "transition_window_descriptive_statistics",
]


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str | float | int | None) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def fmt(value, digits: int = 8) -> str:
    if value is None:
        return "—"
    if isinstance(value, str):
        return value
    number = float(value)
    if math.isnan(number):
        return "NaN"
    if math.isinf(number):
        return "Infinity" if number > 0 else "-Infinity"
    if number == 0:
        return "0"
    if abs(number) >= 1e4 or abs(number) < 1e-4:
        return f"{number:.{digits}e}"
    if number.is_integer():
        return str(int(number))
    return f"{number:.{digits}g}"


def bracket(rows: list[dict], target: int) -> tuple[dict | None, dict | None]:
    steps = [int(row["step"]) for row in rows]
    before_index = bisect.bisect_right(steps, target) - 1
    after_index = bisect.bisect_left(steps, target)
    before = rows[before_index] if before_index >= 0 else None
    after = rows[after_index] if after_index < len(rows) else None
    return before, after


def stats_lookup(
    rows: list[dict[str, str]],
    window_name: str,
    variable: str,
) -> dict[str, str]:
    for row in rows:
        if row["window_name"] == window_name and row["variable"] == variable:
            return row
    raise KeyError((window_name, variable))


def file_size(path: Path) -> str:
    size = path.stat().st_size
    value = float(size)
    unit = "B"
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024 or unit == "GiB":
            break
        value /= 1024
    return f"{value:.2f} {unit}"


def main() -> int:
    started = time.perf_counter()
    if REPORT_PATH.exists() or TEXT_LOG_PATH.exists():
        raise RuntimeError("Stage A1b報告書の上書きを避けて停止")
    for path in (VERIFY_PATH, ANALYSIS_MANIFEST_PATH, FIGURE_MANIFEST_PATH):
        if not path.is_file():
            raise RuntimeError(f"先行工程記録が欠落: {path}")
    verification = json.loads(VERIFY_PATH.read_text(encoding="utf-8"))
    analysis = json.loads(ANALYSIS_MANIFEST_PATH.read_text(encoding="utf-8"))
    figures = json.loads(FIGURE_MANIFEST_PATH.read_text(encoding="utf-8"))
    direct = json.loads(
        (PROCESSED_DIR / "transition_direct_observations.json").read_text(
            encoding="utf-8"
        )
    )

    required_tables = [
        PROCESSED_DIR / f"{stem}.{extension}"
        for stem in REQUIRED_TABLE_STEMS
        for extension in ("csv", "md")
    ]
    input_ok = verification.get("success") is True
    outputs_ok = (
        analysis.get("success") is True
        and all(path.is_file() for path in required_tables)
        and analysis.get("full_time_candidate_cartesian_product_generated") is False
    )
    figures_ok = (
        figures.get("success") is True
        and figures.get("figure_count_by_format") == {"png": 13, "svg": 13}
        and all(Path(path).is_file() for path in figures["generated_files"])
    )
    if not input_ok:
        overall = "INPUT_MISMATCH"
    elif outputs_ok and figures_ok:
        overall = "TRANSITION_ANATOMY_COMPLETE"
    else:
        overall = "TRANSITION_ANATOMY_INCOMPLETE"

    first_passage = load_csv(PROCESSED_DIR / "f_first_passage_levels.csv")
    rates = load_csv(PROCESSED_DIR / "f_decade_growth_rates.csv")
    stats = load_csv(
        PROCESSED_DIR / "transition_window_descriptive_statistics.csv"
    )
    f_rows = load_csv(PROCESSED_DIR / "transition_f_metrics_0_3000.csv")
    occ_rows = load_csv(PROCESSED_DIR / "occupation_actual_records_0_3000.csv")
    q_rows = load_csv(PROCESSED_DIR / "q_actual_records_0_3000.csv")
    f_by_step = {int(row["step"]): row for row in f_rows}

    valid_decade_rates = [
        as_float(row["mean_exponential_rate_per_step"])
        for row in rates
        if row["pair_type"] == "decade_to_decade" and row["status"] == "found"
    ]
    rate_array = np.asarray(valid_decade_rates, dtype=float)
    rate_summary = {
        "count": len(rate_array),
        "min": float(np.min(rate_array)) if len(rate_array) else None,
        "q25": float(np.quantile(rate_array, 0.25)) if len(rate_array) else None,
        "median": float(np.median(rate_array)) if len(rate_array) else None,
        "q75": float(np.quantile(rate_array, 0.75)) if len(rate_array) else None,
        "max": float(np.max(rate_array)) if len(rate_array) else None,
        "mean": float(np.mean(rate_array)) if len(rate_array) else None,
        "std": float(np.std(rate_array)) if len(rate_array) else None,
    }

    first_rank4 = direct["first_rank_q_4_saved_record"]
    first_rank4_step = int(first_rank4["step"]) if first_rank4 else None
    first_rank4_f = (
        as_float(f_by_step[first_rank4_step]["f"])
        if first_rank4_step in f_by_step
        else None
    )
    first_rank4_occ_before, first_rank4_occ_after = bracket(
        occ_rows,
        first_rank4_step,
    )

    input_checks = verification["input_checks"]
    time_axis = verification["time_axis"]
    lines = [
        "# 第7論文 N=5 最初の主成長エピソード構造観察報告（Stage A1b）",
        "",
        "## 1. 範囲と総合状態",
        "",
        "- 入力: Stage A0で完全再現されたN=5 CSV 3件だけ",
        "- 基本観察範囲: absolute step 0〜3000",
        f"- 固定拡大範囲: `{FIXED_WINDOWS}`",
        "- 原本力学コード、新しい軌道、第8論文本実験、高精度計算、N=40、N=300: 未実行",
        "- 全時間域候補の直積: 未生成",
        "- 単一の成長開始・終了・方向成立時刻: 未採用",
        "- H1/H2/H0: 未判定",
        f"- **総合状態: `{overall}`**",
        "",
        "## 2. 入力SHA-256",
        "",
        "| 入力 | 行数 | 絶対パス | SHA-256 | Stage A0照合 |",
        "|---|---:|---|---|---|",
    ]
    for name in ("fcurve", "q_svd", "paper7_long"):
        item = input_checks[name]
        lines.append(
            f"| `{name}` | {item['row_count']} | `{item['path']}` | "
            f"`{item['actual_sha256']}` | "
            f"{'一致' if item['success'] else '不一致'} |"
        )

    lines.extend(
        [
            "",
            "## 3. 時間軸と保存間隔",
            "",
            f"- 共通時間軸: `{time_axis['axis']}`",
            f"- f: 全範囲 `{time_axis['f_full_range']}`、観察窓内1 step間隔",
            f"- q: 全範囲 `{time_axis['q_full_range']}`、観察窓内保存間隔 `{time_axis['q_intervals_in_observation_range']}` step",
            f"- 状態占有: 全範囲 `{time_axis['paper7_full_range']}`、観察窓内保存間隔 `{time_axis['paper7_intervals_in_observation_range']}` step",
            f"- 全q保存行で `step-relative_time={time_axis['q_step_minus_relative_time']}`",
            "- q未保存stepは補間していない。",
            "- 状態占有の解析表・記述統計・初回到達対応は実保存値だけを使用した。",
            "- 状態占有の線形補間は表示専用ファイルへ分離し、図の線にだけ使用した。",
            "",
            "## 4. fの各水準初回到達step",
            "",
            "| level | source | status | first step | f at step | step-crossing |",
            "|---:|---|---|---:|---:|---:|",
        ]
    )
    for row in first_passage:
        lines.append(
            f"| `{row['level_label']}` | `{row['level_source']}` | "
            f"`{row['status']}` | {row['first_passage_step'] or '—'} | "
            f"{fmt(as_float(row['f_at_first_passage']))} | "
            f"{row['relative_to_crossing'] or '0' if row['relative_to_crossing'] == '0' else row['relative_to_crossing'] or '—'} |"
        )

    lines.extend(
        [
            "",
            "これらは増幅過程の座標であり、採用イベント閾値ではない。正の最小fを含むdecadeの下端から0.1までを列挙したため、最初の複数水準が同じ保存stepに到達する場合も削除していない。",
            "",
            "## 5. decade間の増幅率の安定性",
            "",
            "| decade pair | lower step | upper step | Δstep | ln(level ratio)/Δstep |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in rates:
        if row["pair_type"] != "decade_to_decade":
            continue
        lines.append(
            f"| `{fmt(as_float(row['lower_level']))} → {fmt(as_float(row['upper_level']))}` | "
            f"{row['lower_first_passage_step'] or '—'} | "
            f"{row['upper_first_passage_step'] or '—'} | "
            f"{row['step_difference'] or '0' if row['step_difference'] == '0' else row['step_difference'] or '—'} | "
            f"{fmt(as_float(row['mean_exponential_rate_per_step']))} |"
        )
    lines.extend(
        [
            "",
            f"- 有効なdecade-to-decade平均率数: `{rate_summary['count']}`",
            f"- 最小 / Q25 / 中央値 / Q75 / 最大: `{fmt(rate_summary['min'])} / {fmt(rate_summary['q25'])} / {fmt(rate_summary['median'])} / {fmt(rate_summary['q75'])} / {fmt(rate_summary['max'])}`",
            f"- 平均 / 標準偏差: `{fmt(rate_summary['mean'])} / {fmt(rate_summary['std'])}`",
            "",
            "平均率は水準間で完全一定ではない。分布と各行を提示するだけで、安定区間や主指数を自動選択していない。",
            "",
            "## 6. direction 3/4占有が増加する固定観察step帯",
            "",
            "次表は人間が指定した固定窓内の最初と最後の実保存値、および差である。増加帯の境界を新たに推定していない。",
            "",
            "| 固定窓 | d3 first→last | d3差 | d4 first→last | d4差 |",
            "|---|---|---:|---|---:|",
        ]
    )
    for label, _, _ in FIXED_WINDOWS:
        window_name = "zoom_" + label.replace("-", "_")
        d3 = stats_lookup(stats, window_name, "direction_3_occupation")
        d4 = stats_lookup(stats, window_name, "direction_4_occupation")
        lines.append(
            f"| `{label}` | {fmt(as_float(d3['first_value']))} → "
            f"{fmt(as_float(d3['last_value']))} | {fmt(as_float(d3['net_change']))} | "
            f"{fmt(as_float(d4['first_value']))} → {fmt(as_float(d4['last_value']))} | "
            f"{fmt(as_float(d4['net_change']))} |"
        )

    lines.extend(
        [
            "",
            "## 7. q3/q4が増加する固定観察step帯",
            "",
            "| 固定窓 | q3/q1 first→last | 差 | q4/q1 first→last | 差 |",
            "|---|---|---:|---|---:|",
        ]
    )
    for label, _, _ in FIXED_WINDOWS:
        window_name = "zoom_" + label.replace("-", "_")
        q3 = stats_lookup(stats, window_name, "q3_over_q1")
        q4 = stats_lookup(stats, window_name, "q4_over_q1")
        lines.append(
            f"| `{label}` | {fmt(as_float(q3['first_value']))} → "
            f"{fmt(as_float(q3['last_value']))} | {fmt(as_float(q3['net_change']))} | "
            f"{fmt(as_float(q4['first_value']))} → {fmt(as_float(q4['last_value']))} | "
            f"{fmt(as_float(q4['net_change']))} |"
        )

    lines.extend(
        [
            "",
            "## 8. A〜Eの重要な区別",
            "",
            "| 区別 | データ上の記述 | 同一視しない理由 |",
            "|---|---|---|",
            f"| A. rank_q=4 | 最初の実保存rank_q=4はstep `{first_rank4_step}`。その行のq3/q1=`{fmt(as_float(first_rank4['q3_over_q1']))}`、q4/q1=`{fmt(as_float(first_rank4['q4_over_q1']))}` | 既存相対閾値への応答であり、状態占有増加やf増幅の成立を意味しない |",
            f"| B. q3/q4の有限・非ゼロ値 | q3は初期保存行から有限値。q4の最初の正の実保存値はstep `{direct['first_q4_positive_saved_record']['step']}` | 有限値・非ゼロ値と大きさの増加は別である |",
            "| C. direction 3/4占有増加 | 25 step間隔の実保存占有を固定窓ごとに上表へ記載 | qのrankや有限値とは異なる状態占有である |",
            "| D. fの指数増幅 | decade初回到達のstep間隔と平均指数率として記載 | 単一rank行や単一占有行とは異なる区間的挙動である |",
            f"| E. f>0.05 | 既存crossing=`{CROSSING}` | 第7論文既存閾値であり、A〜Dを定義しない |",
            "",
            "### rank_qの早期数値床応答",
            "",
            f"- 最初のrank_q=4保存行: step `{first_rank4_step}`、f=`{fmt(first_rank4_f)}`。",
            f"- その時点を挟む状態占有の実保存step: `{first_rank4_occ_before['step']}` と `{first_rank4_occ_after['step']}`。",
            f"- beforeのdirection 3/4: `{fmt(as_float(first_rank4_occ_before['direction_3_occupation']))}` / `{fmt(as_float(first_rank4_occ_before['direction_4_occupation']))}`。",
            f"- afterのdirection 3/4: `{fmt(as_float(first_rank4_occ_after['direction_3_occupation']))}` / `{fmt(as_float(first_rank4_occ_after['direction_4_occupation']))}`。",
            f"- crossing前のrank_q=4保存行数: `{direct['rank_q_4_saved_records_before_crossing']}`。",
            "",
            "rank_q=4がfの大域増幅およびdirection 3/4の有限占有より大幅に早く現れるため、既存rank閾値が早期の数値床に反応していることと整合する。ただし、この後処理だけで数値誤差の原因や物理的無効性までは証明しない。",
            "",
            "## 9. crossing=1167前後の実保存値",
            "",
            "fは毎step、qは5 step、占有は25 step保存なので、同一stepへ補間せず各系列の実レコードで挟む。",
            "",
            "| 系列 | before step | before値 | after/crossing step | after/crossing値 |",
            "|---|---:|---|---:|---|",
        ]
    )
    crossing_f = direct["crossing_f_records"]
    crossing_occ = direct["crossing_occupation_actual_bracket"]
    crossing_q = direct["crossing_q_actual_bracket"]
    lines.extend(
        [
            f"| f | {crossing_f['previous']['step']} | f={fmt(crossing_f['previous']['f'])} | "
            f"{crossing_f['crossing']['step']} | f={fmt(crossing_f['crossing']['f'])} |",
            f"| occupation | {crossing_occ['before_or_at']['step']} | "
            f"d3={fmt(crossing_occ['before_or_at']['direction_3_occupation'])}, "
            f"d4={fmt(crossing_occ['before_or_at']['direction_4_occupation'])}, "
            f"kernel={fmt(crossing_occ['before_or_at']['kernel_occupation'])} | "
            f"{crossing_occ['after_or_at']['step']} | "
            f"d3={fmt(crossing_occ['after_or_at']['direction_3_occupation'])}, "
            f"d4={fmt(crossing_occ['after_or_at']['direction_4_occupation'])}, "
            f"kernel={fmt(crossing_occ['after_or_at']['kernel_occupation'])} |",
            f"| q | {crossing_q['before_or_at']['step']} | "
            f"q3/q1={fmt(crossing_q['before_or_at']['q3_over_q1'])}, "
            f"q4/q1={fmt(crossing_q['before_or_at']['q4_over_q1'])}, "
            f"rank={crossing_q['before_or_at']['rank_q']} | "
            f"{crossing_q['after_or_at']['step']} | "
            f"q3/q1={fmt(crossing_q['after_or_at']['q3_over_q1'])}, "
            f"q4/q1={fmt(crossing_q['after_or_at']['q4_over_q1'])}, "
            f"rank={crossing_q['after_or_at']['rank_q']} |",
            "",
            "crossingの1 step前からcrossingまでにfが既存0.05水準を超える。一方、qと状態占有は異なる保存間隔で既に連続的に変化しており、同一stepで新たに全量が同時成立したとは読まない。",
            "",
            "## 10. crossing後に遅れて変化する量",
            "",
            "| 固定step | f | d3実保存 | d4実保存 | kernel実保存 | q3/q1実保存 | q4/q1実保存 | rank_q |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for step in (1200, 1400, 1800, 2500):
        occ_before, occ_after = bracket(occ_rows, step)
        q_before, q_after = bracket(q_rows, step)
        occ = occ_before if int(occ_before["step"]) == step else occ_after
        q = q_before if int(q_before["step"]) == step else q_after
        lines.append(
            f"| {step} | {fmt(as_float(f_by_step[step]['f']))} | "
            f"{fmt(as_float(occ['direction_3_occupation']))} | "
            f"{fmt(as_float(occ['direction_4_occupation']))} | "
            f"{fmt(as_float(occ['kernel_occupation']))} | "
            f"{fmt(as_float(q['q3_over_q1']))} | "
            f"{fmt(as_float(q['q4_over_q1']))} | {q['rank_q']} |"
        )
    lines.extend(
        [
            "",
            "この固定step表と図9・10は、crossing後もdirection 3/4占有、kernel、q比が同じ速度では変化しないことを示す。遅れの開始・終了時刻は採用していない。",
            "",
            "## 11. 増幅から飽和・振動へ移る固定観察帯",
            "",
        ]
    )
    f_800_1400 = stats_lookup(stats, "zoom_800_1400", "f")
    f_1400_2500 = stats_lookup(stats, "zoom_1400_2500", "f")
    run_1400_2500 = stats_lookup(stats, "zoom_1400_2500", "running_max_f")
    lines.extend(
        [
            f"- 800〜1400のf: min `{fmt(as_float(f_800_1400['minimum']))}`、max `{fmt(as_float(f_800_1400['maximum']))}`、正差分 `{f_800_1400['positive_differences']}`、負差分 `{f_800_1400['negative_differences']}`。",
            f"- 1400〜2500のf: min `{fmt(as_float(f_1400_2500['minimum']))}`、max `{fmt(as_float(f_1400_2500['maximum']))}`、正差分 `{f_1400_2500['positive_differences']}`、負差分 `{f_1400_2500['negative_differences']}`。",
            f"- 1400〜2500のrunning maximum: `{fmt(as_float(run_1400_2500['first_value']))}` → `{fmt(as_float(run_1400_2500['last_value']))}`。",
            "",
            "固定窓1400〜2500ではfの正負両方向差分が多数あり、単調なdecade通過だけでは記述できない飽和・振動的挙動が見える。これは移行帯の観察であり、成長終了時刻の決定ではない。",
            "",
            "## 12. データから直接言えること",
            "",
            "1. fはstep 0〜3000の範囲で多数の10進水準を順次通過し、既存crossing=1167で初めて0.05を超える。",
            "2. decade間のstep差と平均指数率は保存データから直接計算できるが、全decadeで完全一定ではない。",
            "3. 既存rank_q=4はcrossingや大域的状態占有増加より非常に早い保存行で現れる。",
            "4. direction 3/4占有、q3/q4比、kernel、fは、共通absolute step上で異なる推移を示す。",
            "5. crossing前後の各系列は保存間隔が異なるため、同時性は実保存レコードの範囲でしか言えない。",
            "6. 1400〜2500の固定窓ではfに増減が共存し、初期の単調な増幅だけとは異なる。",
            "",
            "## 13. データだけでは言えないこと",
            "",
            "1. 単一の指数成長開始時刻、終了時刻、方向成立時刻。",
            "2. rank_q=4の早期応答が物理方向成立を意味するかどうか。",
            "3. q未保存stepまたは占有未保存stepでの厳密な同時性。",
            "4. direction 3とdirection 4のどちらを先行方向と解釈すべきか。",
            "5. 増幅から飽和・振動へ移る境界の一意性。",
            "6. H1/H2/H0のいずれが正しいか。",
            "",
            "## 14. 人間が次に固定すべき最小定義",
            "",
            "1. 「主増幅」をrunning maximumのdecade通過で定義するか、log(f)回帰で定義するか。",
            "2. 増幅の終了をrunning maximumの停滞、局所slope、振動幅のどれで記述するか。",
            "3. q3/q4について数値床を除外する最小振幅・持続・保存レコード規則。",
            "4. direction 3/4状態占有について、増加を認定する振幅と持続の規則。",
            "5. 異なる保存間隔の系列間で、同時・先行・遅延を記述するbracketing規則。",
            "6. 既存crossing=1167を参照座標のまま使うか、他イベントの定義へ組み込むか。",
            "",
            "これらを人間が固定するまで、単一イベント時刻を採用しない。",
            "",
            "## 15. 出力表",
            "",
            "| 表 | CSV | Markdown |",
            "|---|---|---|",
        ]
    )
    for stem in REQUIRED_TABLE_STEMS:
        csv_path = PROCESSED_DIR / f"{stem}.csv"
        md_path = PROCESSED_DIR / f"{stem}.md"
        lines.append(
            f"| `{stem}` | `{csv_path}` ({file_size(csv_path)}) | "
            f"`{md_path}` ({file_size(md_path)}) |"
        )

    lines.extend(
        [
            "",
            "## 16. 出力図",
            "",
        ]
    )
    for path in figures["generated_files"]:
        lines.append(f"- `{path}`")
    lines.extend(
        [
            "",
            "図の単一イベント参照線は既存crossing=1167だけである。状態占有の線形補間は表示専用、qは実保存点だけを使用した。",
            "",
            "## 17. 実行環境と停止",
            "",
            f"- Python: `{sys.version}`",
            f"- NumPy: `{np.__version__}`",
            f"- OS: `{platform.platform()}`",
            f"- 報告生成日時（UTC）: `{datetime.now(timezone.utc).isoformat()}`",
            f"- 報告生成時間（秒）: `{time.perf_counter() - started:.6f}`",
            f"- **最終状態: `{overall}`**",
            "- Stage A1bはここで停止する。Stage A2へ進まない。",
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
    return 0 if overall == "TRANSITION_ANATOMY_COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
