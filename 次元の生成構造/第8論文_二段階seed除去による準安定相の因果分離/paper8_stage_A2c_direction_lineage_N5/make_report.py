#!/usr/bin/env python3
"""Stage A2c報告書を生成し、後続Stageへ進まず停止する。"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RAW = HERE / "raw"
PROCESSED = HERE / "processed"
FIGURES = HERE / "figures"
REPORTS = HERE / "reports"
LOGS = HERE / "logs"
CFG = json.loads((HERE / "config_locked.json").read_text(encoding="utf-8"))
EXPECTED = json.loads((HERE / "expected_hashes.json").read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def md_table(fields: list[str], data: list[dict]) -> str:
    out = [
        "| " + " | ".join(fields) + " |",
        "|" + "|".join(["---"] * len(fields)) + "|",
    ]
    for row in data:
        out.append("| " + " | ".join(str(row.get(f, "")).replace("|", "\\|") for f in fields) + " |")
    return "\n".join(out)


def main() -> None:
    gates = {
        "source": LOGS / "source_verification.json",
        "replay": LOGS / "replay_summary.json",
        "transverse": LOGS / "transverse_reconstruction_summary.json",
        "analysis": PROCESSED / "lineage_analysis_summary.json",
        "figures": FIGURES / "figure_manifest.json",
    }
    if not all(p.is_file() for p in gates.values()):
        raise SystemExit("EXECUTION_FAILED: report gates missing")
    source = json.loads(gates["source"].read_text(encoding="utf-8"))
    replay = json.loads(gates["replay"].read_text(encoding="utf-8"))
    transverse = json.loads(gates["transverse"].read_text(encoding="utf-8"))
    analysis = json.loads(gates["analysis"].read_text(encoding="utf-8"))
    figures = json.loads(gates["figures"].read_text(encoding="utf-8"))
    if not (
        source.get("status") == "VERIFIED"
        and replay.get("status") == "REPLAY_COMPLETE"
        and replay.get("trajectory_f_bitwise_matches_stage_a0")
        and replay.get("q_diagnostics_bitwise_match_stage_a0")
        and transverse.get("status") == "TRANSVERSE_RECONSTRUCTED"
        and transverse.get("stage_a0_first_record_strings_all_match")
        and analysis.get("status") == "ANALYSIS_COMPLETE"
        and analysis.get("numerical_health_passed")
        and figures.get("status") == "FIGURES_COMPLETE"
        and figures.get("figure_count") == 15
    ):
        raise SystemExit("EXECUTION_FAILED: one or more report gates failed")
    if REPORTS.exists() and any(REPORTS.iterdir()):
        raise SystemExit("EXECUTION_FAILED: reports/が空ではないため上書きを拒否")

    classification = analysis["classification"]
    if classification not in CFG["classification"]["labels"]:
        raise SystemExit("EXECUTION_FAILED: classification outside locked labels")
    flevels = rows(PROCESSED / "lineage_by_f_level.csv")
    qbands = rows(PROCESSED / "lineage_by_q_resolution_band.csv")
    quality = rows(PROCESSED / "direction_projector_quality.csv")
    late_trans = rows(PROCESSED / "late_direction_vs_transverse_overlap.csv")
    continuity = rows(PROCESSED / "consecutive_subspace_continuity.csv")
    late = rows(PROCESSED / "early_vs_late_direction_overlap.csv")
    trans_rows = rows(PROCESSED / "transverse_direction_reconstruction.csv")

    source_rows = []
    for group in ("sources", "dependencies"):
        for name, item in EXPECTED[group].items():
            path = REPO / item["path"]
            source_rows.append({
                "file": name,
                "sha256": sha256(path),
                "absolute_path": str(path),
            })

    pre_quality = [r for r in quality if int(r["step"]) < replay["crossing"]]
    pre_counts = {
        band["id"]: sum(r["q_resolution_band"] == band["id"] for r in pre_quality)
        for band in CFG["q_resolution_bands"]
    }
    late_window_rows = [
        r for r in late
        if CFG["late_reference"]["start_step"] <= int(r["step"]) <= CFG["late_reference"]["end_step"]
    ]
    late_window_overlap_median = float(np.median([float(r["overlap"]) for r in late_window_rows]))
    late_window_angle_median = float(np.median([float(r["maximum_principal_angle_rad"]) for r in late_window_rows]))
    late_window_overlap_min = min(float(r["overlap"]) for r in late_window_rows)
    late_window_angle_max = max(float(r["maximum_principal_angle_rad"]) for r in late_window_rows)
    crossing_cont = [
        r for r in continuity
        if CFG["display_windows"]["crossing_zoom"][0] < int(r["step_to"]) <= CFG["display_windows"]["crossing_zoom"][1]
    ]
    crossing_rotation_sum = sum(float(r["maximum_principal_angle_rad"]) for r in crossing_cont)

    if classification == "FIRST_DIRECTIONS_CONTINUOUS":
        answer1 = (
            "固定された記述的規則では「連続」と分類された。"
            "ただし、q解像度床未満の最早期点まで物理的に同一方向だったとは断定しない。"
        )
        answer2 = "いいえ。固定規則上、比較可能範囲のD34は既存Tperp群と低overlapであり、late D34側へ対応した。"
    elif classification == "MATCHES_LATE_TRANSVERSE_GERM":
        answer1 = "いいえ。固定規則上、比較可能なD34はD34_lateには対応しなかった。"
        answer2 = "はい。固定規則上、既存Tperpの少なくとも一つに対応した。"
    elif classification == "ROTATING_OR_MIXED_LINEAGE":
        answer1 = "単純な連続とは分類されなかった。比較可能範囲は回転または混合した系譜として記述された。"
        answer2 = "既存Tperpへの単純な一致とも分類されなかった。"
    else:
        answer1 = "解像度または既存横方向再構成の条件が不足し、判定不能と分類された。"
        answer2 = "判定不能である。"
    answer3 = (
        f"固定step 900〜1400で、連続step間最大主角の最大値は "
        f"`{analysis['crossing_window_max_consecutive_angle_rad']:.6e}` rad、"
        f"projector distance最大値は `{analysis['crossing_window_max_projector_distance']:.6e}`。"
        f"同区間の列swapは `{analysis['crossing_window_basis_swap_count']}` 回、"
        f"sign flip総数は `{analysis['crossing_window_sign_flip_total']}`。"
        "列交換・符号反転と2次元部分空間回転は別表で分離した。"
    )

    f_fields = [
        "level_label", "status", "step", "q_resolution_band",
        "D34_vs_late_overlap", "D34_vs_late_max_angle_rad",
        "D34_vs_Tperp_seed0_overlap", "D34_vs_Tperp_seed1_overlap",
        "D34_vs_Tperp_seed2_overlap",
    ]
    q_fields = [
        "q_resolution_band", "point_count", "first_step", "last_step",
        "median_D34_vs_late_overlap", "median_D34_vs_late_max_angle_rad",
        "median_max_D34_vs_Tperp_overlap", "maximum_consecutive_max_angle_rad",
        "basis_column_swap_count", "sign_flip_total",
    ]
    late_trans_fields = [
        "seed", "epsilons_sharing_this_direction", "overlap",
        "theta_1_rad", "theta_2_rad", "projector_distance",
        "same_t0_S4_vs_Tperp_orthogonality_error",
    ]
    source_table = md_table(["file", "sha256", "absolute_path"], source_rows)
    f_table = md_table(f_fields, flevels)
    q_table = md_table(q_fields, qbands)
    lt_table = md_table(late_trans_fields, late_trans)
    pre_band_text = "\n".join(f"- `{k}`: {v} step" for k, v in pre_counts.items())
    trans_medians = analysis["classification_D34_Tperp_overlap_median_by_seed"]
    trans_median_text = ", ".join(f"seed {i}: `{v:.6e}`" for i, v in enumerate(trans_medians))
    figure_list = "\n".join(f"- `{name}`" for name in figures["files"] if name.endswith(".png"))
    table_list = "\n".join(
        f"- [{stem}](../processed/{stem}.md)"
        for stem in (
            "direction_basis_snapshots",
            "direction_projector_quality",
            "consecutive_subspace_continuity",
            "early_vs_late_direction_overlap",
            "early_vs_transverse_overlap",
            "late_direction_vs_transverse_overlap",
            "lineage_by_f_level",
            "lineage_by_q_resolution_band",
            "transverse_direction_reconstruction",
            "numerical_health",
        )
    )

    report = f"""# Stage A2c N=5 方向基底の系譜追跡 報告書

## 実行状態

**STAGE_A2C_COMPLETE**

記述的分類: **{classification}**

この分類は方向系譜だけを対象とする。H1/H2/H0、三方向成立step、物理的追加次元の存在は判定していない。

## 使用原本とSHA-256

{source_table}

固定原本、依存原本、Stage A0成果物は実行前にSHA-256照合した。原本コードと既存成果物は編集・上書きしていない。

## Python環境

- Python: `{sys.version.replace(chr(10), " ")}`
- NumPy: `{np.__version__}`
- OS: `{platform.platform()}`

## 軌道再現

- 対象: `N=5`, `float64`, Stage A0の `delta=1e-15`, `seed_index=0`
- 再実行範囲: step `0..5000`
- fのStage A0とのbitwise一致: `{replay['trajectory_f_bitwise_matches_stage_a0']}`（`{replay['trajectory_f_comparison_count']}`点）
- f最大絶対誤差: `{replay['trajectory_f_max_absolute_error']:.17e}`
- q1〜q4、rank_q、gram_rank、dominant eigenvalueのStage A0とのbitwise一致: `{replay['q_diagnostics_bitwise_match_stage_a0']}`（`{replay['q_comparison_count']}`保存点）
- q最大絶対誤差: `{replay['q_max_absolute_error']:.17e}`
- existing crossing: `{replay['crossing']}`
- 親残差: `{replay['parent_residual']:.17e}`

したがって、以下は新しい軌道探索ではなく、Stage A0軌道内部で既存関数が生成する基底の抽出である。

## D34(t)の既存定義

各stepで既存 `gram_reduce` と `dominant_plane` により `Bdom(t)` を得て、既存

`D34(t) = s4_new_dirs(B0, Bdom(t))`

をそのまま使用した。比較の正本は列ではなく

`P34(t) = D34(t) D34(t)^T`

である。個別direction 3/4表示は、既存5色コードと同様に固定other空間へ射影・QR後、既存 `align_2d` で連続化した。

`S4(t)` は既存 `run_paper7_transverse.py:s4_basis` で再構成し、`PS4(t)=S4(t)S4(t)^T` を保存した。

## D34_lateの比較用構成

step `{CFG['late_reference']['start_step']}..{CFG['late_reference']['end_step']}` の全 `{analysis['D34_late_member_count']}` 個の `P34(t)` を平均し、平均射影行列の上位2固有ベクトルを `D34_late` とした。これは比較用代表部分空間であり、新しい物理方向定義ではない。

- 平均射影行列固有値: `{analysis['D34_late_mean_projector_eigenvalues']}`
- late window内のD34対D34_late overlap中央値/最小値: `{late_window_overlap_median:.6e}` / `{late_window_overlap_min:.6e}`
- late window内の最大主角中央値/最大値: `{late_window_angle_median:.6e}` / `{late_window_angle_max:.6e}` rad

## 横摂動方向の既存コードからの再構成

既存 `run_paper7_transverse.py` の固定値をそのまま用いた。

- `t0 = crossing + 3000 = {transverse['t0']}`
- 方向PRNG seed: `{transverse['direction_prng_seed']}`
- 既存方向seed数: `{transverse['unique_seed_directions']}`
- epsilon: `{transverse['epsilons']}`
- 規則: `eta.real` と `eta.imag` をそれぞれ既存 `S4(t0)^perp` へ射影し、複素ベクトル全体を正規化
- 保存正本: 原本が使った複素 `eta`
- 2次元比較表現: `span(eta.real, eta.imag)` のQR基底
- 任意方向の追加: `{transverse['new_arbitrary_direction_added']}`
- Stage A0の最初の保存step 4200における全12レコード文字列一致: `{transverse['stage_a0_first_record_strings_all_match']}`
- 基準軌道状態のbitwise一致: `{transverse['baseline_states_bitwise_match_replay']}`

同じt0のS4とTperpの最大直交誤差は `{analysis['t0_S4_vs_Tperp_max_orthogonality_error']:.6e}` であり、構成上の直交を数値誤差内で確認した。

## 急拡大前のD34の数値解像度

crossing前step 0〜1166のq解像度帯別点数:

{pre_band_text}

- crossing前で `min(q3,q4)/q1 >= 1e-6` の比較点数: `{analysis['pre_crossing_resolved_point_count']}`
- その比較点でのD34対D34_late overlap中央値: `{analysis['pre_crossing_resolved_late_overlap_median']}`
- その比較点でのD34対Tperp最大overlap中央値: `{analysis['pre_crossing_resolved_max_transverse_overlap_median']}`

`min(q3,q4)/q1 < 1e-8` の早期基底も削除していないが、数値解像度不足帯として分離した。この帯の列方向を物理的に確定した方向とは扱わない。

## 急拡大中のD34部分空間回転

固定step 900〜1400について:

- 連続step間最大主角の最大値: `{analysis['crossing_window_max_consecutive_angle_rad']:.6e}` rad
- 連続step間projector distance最大値: `{analysis['crossing_window_max_projector_distance']:.6e}`
- 最大主角の区間内総和（経路長の記述量）: `{crossing_rotation_sum:.6e}` rad
- basis column swap: `{analysis['crossing_window_basis_swap_count']}` 回
- sign flip: `{analysis['crossing_window_sign_flip_total']}` 列

総和はイベント閾値ではなく、固定観察窓内の連続step回転量の単純和である。

## 急拡大後D34との主角・overlap

分類評価は `min(q3,q4)/q1 >= 1e-6` を初めて満たすstep `{analysis['classification_evaluation_start_step']}` からstep 5000までとした。

- D34(t)対D34_late overlap中央値: `{analysis['classification_late_overlap_median']:.6e}`
- 最大主角中央値: `{analysis['classification_late_max_angle_median_rad']:.6e}` rad

## 横摂動方向との主角・overlap

評価範囲でのD34(t)対Tperp overlap中央値:

{trans_median_text}

固定集約値（3 seed中央値の最大）: `{analysis['classification_transverse_aggregate_max_seed_median']:.6e}`

D34_late対各Tperp:

{lt_table}

epsilonは同一seedの方向を変えないため、方向部分空間比較では3個のunique seed方向を用い、12個のepsilon条件は再構成検証表に保持した。

## 個別direction列の交換と部分空間系譜

`consecutive_subspace_continuity` は、各連続stepについて以下を別々に保存した。

- raw列割当からのcolumn swap
- 対応列のsign flip
- signed permutationから残るD34内部回転
- `align_2d` 後の個別列内積
- P34主角・overlap・projector distanceによるambient部分空間回転

したがって、direction 3/4の色交換や符号反転だけを、物理的2次元部分空間の選び直しとは解釈していない。

## f水準別の方向系譜

{f_table}

## q解像度帯別の方向系譜

{q_table}

これらのq帯は数値解像度別表示であり、物理的方向成立閾値ではない。

## 記述的分類

**{classification}**

固定判定量:

- evaluation start: `{analysis['classification_evaluation_start_step']}`
- D34対late overlap中央値: `{analysis['classification_late_overlap_median']}`
- D34対late最大主角中央値: `{analysis['classification_late_max_angle_median_rad']}` rad
- D34対Tperp seed別overlap中央値: `{analysis['classification_D34_Tperp_overlap_median_by_seed']}`
- transverse集約: `{analysis['classification_transverse_aggregate_max_seed_median']}`

## 三つの問いへの回答

1. 急拡大前に微小に見えた方向部分空間は、急拡大後に最初に成立する三方向の部分空間へ連続しているか。

   {answer1}

2. その微小方向部分空間は、準安定後の追加方向萌芽の横摂動部分空間に近いか。

   {answer2}

3. 急拡大途中に、方向部分空間の交換・回転・選び直しが起きているか。

   {answer3}

## データから直接言えること

- Stage A0軌道のfとq診断をbitwise一致で再実行できた。
- P34によるD34部分空間系譜、D34_late、既存3個のTperp方向との主角・overlap・射影距離を全stepで比較できた。
- t0のS4とTperpは既存生成規則どおり数値誤差内で直交していた。
- 個別列の交換・符号反転と、2次元部分空間自体のambient回転を分離できた。

## データだけでは言えないこと

- q解像度床未満の最早期基底が物理的に確定した方向であるか
- 単一の三方向成立step
- H1/H2/H0の判定
- 追加方向萌芽が自然軌道上に実在する時刻
- 高精度親や別Delta/Nにおける同じ系譜

## 必須表

{table_list}

## 必須図

{figure_list}

## 最終停止

Stage A2cの報告書完成をもって停止する。高精度、Delta掃引、`N=40`、`N=300`、Stage B/Cへ進まない。
"""
    REPORTS.mkdir(parents=True, exist_ok=True)
    output = REPORTS / "stage_A2c_direction_lineage_N5_report.md"
    output.write_text(report, encoding="utf-8")
    print("STAGE_A2C_COMPLETE")
    print(output)


if __name__ == "__main__":
    main()
