#!/usr/bin/env python3
"""Create the final Stage A0 Paper 7 N=5 reproduction report."""

from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
CONFIG_PATH = PACKAGE_ROOT / "config_locked.json"
EXPECTED_PATH = PACKAGE_ROOT / "expected_hashes.json"
VERIFY_PATH = PACKAGE_ROOT / "comparison" / "source_verification.json"
RUN_MANIFEST_PATH = PACKAGE_ROOT / "logs" / "run_manifest.json"
CSV_RESULT_PATH = PACKAGE_ROOT / "comparison" / "csv_comparison.json"
JSON_RESULT_PATH = PACKAGE_ROOT / "comparison" / "json_comparison.json"
FIGURE_RESULT_PATH = PACKAGE_ROOT / "comparison" / "figure_comparison.json"
SUMMARY_PATH = PACKAGE_ROOT / "comparison" / "comparison_summary.json"
REPORT_DIR = PACKAGE_ROOT / "reports"
REPORT_PATH = REPORT_DIR / "paper7_N5_reproduction_report.md"
LOG_PATH = PACKAGE_ROOT / "logs" / "make_reproduction_report.log"


def fmt_number(value) -> str:
    if value is None:
        return "—"
    if isinstance(value, str):
        return value
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        return f"{value:.17g}"
    return str(value)


def md_path(path: str | Path) -> str:
    return f"`{path}`"


def main() -> int:
    started = time.perf_counter()
    if REPORT_PATH.exists():
        raise RuntimeError(f"報告書の上書きを避けて停止: {REPORT_PATH}")

    required = [
        CONFIG_PATH,
        EXPECTED_PATH,
        VERIFY_PATH,
        RUN_MANIFEST_PATH,
        CSV_RESULT_PATH,
        JSON_RESULT_PATH,
        FIGURE_RESULT_PATH,
        SUMMARY_PATH,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError("報告書入力欠落: " + ", ".join(missing))

    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    expected = json.loads(EXPECTED_PATH.read_text(encoding="utf-8"))
    verification = json.loads(VERIFY_PATH.read_text(encoding="utf-8"))
    run_manifest = json.loads(RUN_MANIFEST_PATH.read_text(encoding="utf-8"))
    csv_results = json.loads(CSV_RESULT_PATH.read_text(encoding="utf-8"))
    json_results = json.loads(JSON_RESULT_PATH.read_text(encoding="utf-8"))
    figure_results = json.loads(FIGURE_RESULT_PATH.read_text(encoding="utf-8"))
    comparison_summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))

    if verification.get("success") is not True:
        raise RuntimeError("原本SHA-256検証が成功していない")
    if run_manifest.get("success") is not True or run_manifest.get("N_values_executed") != [5]:
        raise RuntimeError("N=5再生成が成功していない")
    if comparison_summary.get("success") is not True:
        raise RuntimeError("比較工程が完了していない")

    numeric = comparison_summary["numeric_comparison"]
    figures_ok = (
        comparison_summary["all_required_figures_exist"]
        and comparison_summary["all_figure_dimensions_match"]
    )
    if numeric == "exact" and figures_ok:
        overall = "REPRODUCED_EXACTLY"
    elif numeric == "numerical_match" and figures_ok:
        overall = "REPRODUCED_NUMERICALLY"
    else:
        overall = "NOT_REPRODUCED"

    env = run_manifest["environment"]
    lines = [
        "# 第7論文 N=5 再現報告（Stage A0）",
        "",
        "## 1. 実行範囲と総合判定",
        "",
        f"- 対象: 第7論文 N=5 の既存結果だけ",
        f"- 実行したN: `{run_manifest['N_values_executed']}`",
        "- 第8論文本実験、高精度親、Series 1〜3、N=40、N=300: 未実行",
        f"- 数値比較: `{numeric}`",
        f"- 必須PNG生成: `{'一致' if comparison_summary['all_required_figures_exist'] else '欠落あり'}`",
        f"- PNG画像サイズ: `{'一致' if comparison_summary['all_figure_dimensions_match'] else '不一致あり'}`",
        f"- **総合判定: `{overall}`**",
        "",
        "画像SHA-256の不一致だけでは総合判定を失敗にしていない。数値CSV・JSON、PNG生成有無、画像サイズを優先した。",
        "",
        "## 2. 使用した既存コードとSHA-256",
        "",
        "| ファイル | 役割 | 絶対パス | SHA-256 | 検証 |",
        "|---|---|---|---|---|",
    ]
    verification_by_name = {
        (item["group"], item["name"]): item for item in verification["checks"]
    }
    for name, item in expected["sources"].items():
        check = verification_by_name[("sources", name)]
        lines.append(
            f"| `{name}` | {item['role']} | {md_path(check['path'])} | `{check['actual_sha256']}` | "
            f"{'一致' if check['sha256_match'] else '不一致'} |"
        )

    lines.extend(
        [
            "",
            "### import・固定hash依存",
            "",
            "| ファイル | 役割 | 絶対パス | SHA-256 | 検証 |",
            "|---|---|---|---|---|",
        ]
    )
    for name, item in expected["dependencies"].items():
        check = verification_by_name[("dependencies", name)]
        lines.append(
            f"| `{name}` | {item['role']} | {md_path(check['path'])} | `{check['actual_sha256']}` | "
            f"{'一致' if check['sha256_match'] else '不一致'} |"
        )

    lines.extend(
        [
            "",
            "原本コードは編集していない。出力先モジュール変数だけをラッパー内で `reproduced/` 配下へ差し替えた。"
            "`make_paper7_figures.py` と `make_saturation_comparison.py` の `NS` は実行時に `[5]` へ固定した。",
            "",
            "## 3. 実行環境",
            "",
            f"- Python executable: {md_path(env['python_executable'])}",
            f"- Python: `{env['python_version']}`",
            f"- Python implementation: `{env['python_implementation']}`",
            f"- NumPy: `{env['numpy_version']}`",
            f"- Matplotlib: `{env['matplotlib_version']}`",
            f"- platform: `{env['platform']}`",
            f"- OS: `{env['os']['system']} {env['os']['release']} {env['os']['machine']}`",
            "",
            "## 4. 実行コマンドと実行時間",
            "",
        ]
    )
    commands = [
        PACKAGE_ROOT / "verify_sources.py",
        PACKAGE_ROOT / "run_reproduction.py",
        PACKAGE_ROOT / "compare_outputs.py",
        PACKAGE_ROOT / "make_reproduction_report.py",
    ]
    for index, command in enumerate(commands, 1):
        lines.append(f"{index}. `PYTHONDONTWRITEBYTECODE=1 python3 {command}`")
    lines.extend(
        [
            "",
            "| 工程 | 実行時間（秒） |",
            "|---|---:|",
            f"| `verify_sources.py` | {verification['duration_seconds']:.6f} |",
            f"| `run_reproduction.py` 合計 | {run_manifest['duration_seconds']:.6f} |",
        ]
    )
    for item in run_manifest["timings"]:
        lines.append(f"| └ {item['step']} | {item['seconds']:.6f} |")
    lines.extend(
        [
            f"| `compare_outputs.py` | {comparison_summary['duration_seconds']:.6f} |",
            f"| `make_reproduction_report.py`（報告生成時点） | {time.perf_counter() - started:.6f} |",
            "",
            "## 5. 再生成ファイル一覧",
            "",
        ]
    )
    for path in run_manifest["reproduced_files"]:
        lines.append(f"- {md_path(PACKAGE_ROOT / path)}")

    lines.extend(
        [
            "",
            "## 6. 比較対象一覧",
            "",
            "### CSV",
            "",
        ]
    )
    for item in expected["csv_baselines"].values():
        lines.append(f"- {md_path(REPO_ROOT / item['expected'])}")
    lines.extend(["", "### JSON", ""])
    for item in expected["json_baselines"].values():
        lines.append(f"- {md_path(REPO_ROOT / item['expected'])}")
    lines.extend(["", "### PNG", ""])
    for item in expected["png_baselines"].values():
        lines.append(f"- {md_path(REPO_ROOT / item['expected'])}")

    lines.extend(
        [
            "",
            "## 7. CSV比較",
            "",
            "| CSV | 基準行数 | 再現行数 | 列数 | 列名 | 判定 |",
            "|---|---:|---:|---:|---|---|",
        ]
    )
    for name, result in csv_results.items():
        column_count = len(result.get("expected_columns", []))
        lines.append(
            f"| `{name}` | {result.get('expected_row_count', '—')} | "
            f"{result.get('reproduced_row_count', '—')} | {column_count} | "
            f"{'一致' if result.get('column_names_match') else '不一致'} | `{result['classification']}` |"
        )

    for name, result in csv_results.items():
        lines.extend(
            [
                "",
                f"### {name}",
                "",
                "| 列 | 種別 | 最大絶対誤差 | 最大相対誤差 | NaN位置 | bitwise/完全一致 | 許容差内 |",
                "|---|---|---:|---:|---|---|---|",
            ]
        )
        for column, detail in result.get("column_results", {}).items():
            if detail["kind"] == "float":
                nan_status = (
                    f"一致 ({len(detail['expected_nan_rows'])}箇所)"
                    if detail["nan_positions_match"]
                    else "不一致"
                )
                lines.append(
                    f"| `{column}` | float | {fmt_number(detail['max_absolute_error'])} | "
                    f"{fmt_number(detail['max_relative_error'])} | {nan_status} | "
                    f"{'一致' if detail['bitwise_match'] else '不一致'} | "
                    f"{'はい' if detail['within_tolerance'] else 'いいえ'} |"
                )
            else:
                lines.append(
                    f"| `{column}` | {detail['kind']} | — | — | — | "
                    f"{'一致' if detail['exact_match'] else '不一致'} | "
                    f"{'はい' if detail['exact_match'] else 'いいえ'} |"
                )
        nan_columns = [
            (column, detail["expected_nan_rows"], detail["reproduced_nan_rows"])
            for column, detail in result.get("column_results", {}).items()
            if detail["kind"] == "float" and detail["expected_nan_rows"]
        ]
        if nan_columns:
            lines.extend(["", "NaN位置（0始まりのデータ行番号）:"])
            for column, expected_nan, reproduced_nan in nan_columns:
                lines.append(
                    f"- `{column}`: 基準 `{expected_nan}` / 再現 `{reproduced_nan}`"
                )

    lines.extend(
        [
            "",
            "## 8. JSON比較",
            "",
            "| JSON | キー | 型 | 非float値 | float | 判定 |",
            "|---|---|---|---|---|---|",
        ]
    )
    for name, result in json_results.items():
        lines.append(
            f"| `{name}` | {'一致' if result.get('keys_match') else '不一致'} | "
            f"{'一致' if result.get('types_match') else '不一致'} | "
            f"{'一致' if result.get('exact_nonfloat_values_match') else '不一致'} | "
            f"{'bitwise一致' if result.get('all_float_bitwise') else '許容差内' if result.get('all_float_numerical') else '不一致'} | "
            f"`{result['classification']}` |"
        )
        differences = (
            result.get("missing_keys", [])
            + result.get("extra_keys", [])
            + [item["path"] for item in result.get("length_mismatches", [])]
            + [item["path"] for item in result.get("type_mismatches", [])]
            + [item["path"] for item in result.get("value_differences", [])]
        )
        lines.append(
            f"- `{name}` JSON差分: "
            + ("なし" if not differences else ", ".join(f"`{item}`" for item in differences))
        )

    lines.extend(
        [
            "",
            "## 9. 図差分",
            "",
            "| PNG | 生成 | 基準サイズ | 再現サイズ | サイズ | SHA-256 | 判定 |",
            "|---|---|---|---|---|---|---|",
        ]
    )
    for name, result in figure_results.items():
        lines.append(
            f"| `{name}` | {'あり' if result['reproduced_exists'] else 'なし'} | "
            f"`{result.get('expected_dimensions', '—')}` | "
            f"`{result.get('reproduced_dimensions', '—')}` | "
            f"{'一致' if result.get('dimensions_match') else '不一致'} | "
            f"{'一致' if result.get('sha256_match') else '不一致'} | "
            f"`{result['classification']}` |"
        )
        if result.get("expected_sha256"):
            lines.append(f"- `{name}` 基準SHA-256: `{result['expected_sha256']}`")
            lines.append(f"- `{name}` 再現SHA-256: `{result['reproduced_sha256']}`")

    p7_summary = json.loads(
        (
            PACKAGE_ROOT
            / "reproduced/exact_lowN_eigenspectrum_v2/paper7_longtime/summary/N00005_5color_meta.json"
        ).read_text(encoding="utf-8")
    )
    transverse_summary = json.loads(
        (
            PACKAGE_ROOT
            / "reproduced/exact_lowN_eigenspectrum_v2/paper7_longtime/summary/N00005_transverse_meta.json"
        ).read_text(encoding="utf-8")
    )
    saturation_summary = json.loads(
        (
            PACKAGE_ROOT
            / "reproduced/exact_lowN_eigenspectrum_v2/diagnostics/N00005_saturation.json"
        ).read_text(encoding="utf-8")
    )
    metastable_summary = json.loads(
        (
            PACKAGE_ROOT
            / "reproduced/metastable_series_result_v1/metastable_N00005_delta1e-15_seed0.json"
        ).read_text(encoding="utf-8")
    )
    lines.extend(
        [
            "",
            "## 10. 主要な既存出力の再現値",
            "",
            f"- f(t) crossing: `{metastable_summary['base_summary']['crossing_tau']}`",
            f"- 五成分時系列 crossing: `{p7_summary['crossing']}`",
            f"- q/rank診断 crossing: `{saturation_summary['crossing']}`",
            f"- q3（準安定）: `{saturation_summary['q3_meta']}`",
            f"- q4（準安定）: `{saturation_summary['q4_meta']}`",
            f"- q3（最終）: `{saturation_summary['q3_final']}`",
            f"- q4（最終）: `{saturation_summary['q4_final']}`",
            f"- rank_Q（最終）: `{saturation_summary['rank_q_final']}`",
            f"- 横摂動 t0: `{transverse_summary['t0']}`",
            f"- lambda_transverse max: `{transverse_summary['lambda_max_for_N']}`",
            f"- lambda_transverse分類: `{transverse_summary['classification']}`",
            "",
            "これらは既存コードが返す既存量であり、新しい物理量またはイベント定義は追加していない。",
            "",
            "## 11. 不一致の原因候補",
            "",
        ]
    )
    causes = []
    for name, result in csv_results.items():
        if result["classification"] == "mismatch":
            causes.append(f"CSV `{name}` の数値・構造不一致")
    for name, result in json_results.items():
        if result["classification"] == "mismatch":
            causes.append(f"JSON `{name}` のキー・値不一致")
    for name, result in figure_results.items():
        if not result.get("reproduced_exists"):
            causes.append(f"PNG `{name}` の生成欠落")
        elif not result.get("dimensions_match"):
            causes.append(f"PNG `{name}` の画像サイズ不一致")
        elif not result.get("sha256_match"):
            causes.append(f"PNG `{name}` は描画環境差候補（画像サイズ一致、SHA-256のみ不一致）")
    if causes:
        for cause in causes:
            lines.append(f"- {cause}")
    else:
        lines.append("- 不一致なし")

    lines.extend(
        [
            "",
            "数値不一致が存在する場合も、本パッケージは原因候補を報告するだけで、自動修正・代替実装・再実行を行わない。",
            "",
            "## 12. 最終判定",
            "",
            f"**{overall}**",
            "",
            f"- 報告生成日時（UTC）: `{datetime.now(timezone.utc).isoformat()}`",
            "- Stage A0はここで停止する。",
        ]
    )

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    LOG_PATH.write_text(
        f"report={REPORT_PATH}\noverall={overall}\nseconds={time.perf_counter() - started:.6f}\n",
        encoding="utf-8",
    )
    print(f"[REPORT] {REPORT_PATH}")
    print(f"[OVERALL] {overall}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
