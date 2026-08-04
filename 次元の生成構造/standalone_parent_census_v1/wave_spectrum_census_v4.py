#!/usr/bin/env python3
"""Hを置かない関係波・倍音スペクトル読出し器 v4。

入力は closed_wave_trajectory_v1。連続する mathcal N 標本を
X[e,t] (shape M x mathcal N) として読み、倍音は事後DFTだけで得る。
周波数セル、関係成分、絶対位相、相対位相をまとめて消さない。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA = "phase_resolved_wave_spectrum_v4"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Hを置かず、M×mathcal N標本から全周波数セルの振幅・位相を読む。"
    )
    parser.add_argument("--input", type=Path, required=True, help="closed_wave_trajectory_v1")
    parser.add_argument("--output", type=Path, required=True, help="新規出力ディレクトリ")
    parser.add_argument(
        "--window-start",
        type=int,
        default=0,
        help="読み始める保存標本番号。既定は0。",
    )
    return parser.parse_args()


def load_input(input_dir: Path, window_start: int) -> dict[str, Any]:
    manifest_path = input_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "closed_wave_trajectory_v1":
        raise SystemExit(f"未対応schema: {manifest.get('schema')!r}")
    if manifest.get("source") != "white_null":
        raise SystemExit("本v4は white_null 入力だけを受理する")

    theory = manifest["theory_inputs"]
    resolution = int(theory["resolution_mathcal_N"])
    n_body = int(theory["n_body"])
    expected_relations = n_body * (n_body - 1) // 2
    if int(theory["m_relations"]) != expected_relations:
        raise SystemExit("M != N(N-1)/2")

    trajectory = np.load(input_dir / "trajectory.npy", mmap_mode="r")
    edges = np.load(input_dir / "edges.npy")
    steps = np.load(input_dir / "steps.npy")
    if trajectory.ndim != 2 or trajectory.shape[1] != expected_relations:
        raise SystemExit(f"trajectory shape不正: {trajectory.shape}")
    if edges.shape != (expected_relations, 2):
        raise SystemExit(f"edges shape不正: {edges.shape}")
    if window_start < 0 or window_start + resolution > trajectory.shape[0]:
        raise SystemExit(
            f"mathcal N={resolution}連続標本が不足: start={window_start}, "
            f"stored={trajectory.shape[0]}"
        )
    selected_steps = np.asarray(steps[window_start : window_start + resolution], dtype=np.int64)
    if len(selected_steps) > 1 and not np.all(np.diff(selected_steps) == 1):
        raise SystemExit("倍音読出しには record_stride=1 の連続標本が必要")

    # 正本の向き: 行=関係、列=分解能標本。
    waves = np.asarray(
        trajectory[window_start : window_start + resolution, :], dtype=np.complex128
    ).T.copy()
    return {
        "manifest": manifest,
        "manifest_path": manifest_path,
        "n_body": n_body,
        "resolution": resolution,
        "m_relations": expected_relations,
        "edges": np.asarray(edges, dtype=np.int64),
        "selected_steps": selected_steps,
        "waves": waves,
    }


def signed_frequency(index: int, resolution: int) -> int:
    return index if index <= resolution // 2 else index - resolution


def frequency_role(k_signed: int) -> str:
    if k_signed == 0:
        return "直流成分"
    if abs(k_signed) == 1:
        return "基底波"
    return "倍音"


def phase_cell(phase: np.ndarray, resolution: int) -> tuple[np.ndarray, np.ndarray]:
    cell_width = 2.0 * math.pi / resolution
    wrapped = np.mod(phase, 2.0 * math.pi)
    cells = np.mod(np.rint(wrapped / cell_width).astype(np.int64), resolution)
    centres = cells * cell_width
    error = np.angle(np.exp(1j * (wrapped - centres)))
    return cells, error


def fmt(value: float, digits: int = 8) -> str:
    if not math.isfinite(value):
        return "—"
    return f"{value:.{digits}g}"


def edge_label(edge: np.ndarray) -> str:
    return f"({int(edge[0]) + 1},{int(edge[1]) + 1})"


def analyse(data: dict[str, Any], reader_sha256: str) -> dict[str, Any]:
    waves = data["waves"]
    edges = data["edges"]
    resolution = data["resolution"]
    m_relations = data["m_relations"]
    spectrum = np.fft.fft(waves, axis=1, norm="ortho")

    eps = np.finfo(np.float64).eps
    closure_tolerance = float(100.0 * eps * max(resolution, m_relations))
    phase_amplitude_floor = float(
        100.0 * eps * max(float(np.max(np.abs(spectrum))), np.finfo(float).tiny)
    )

    time_power = np.sum(np.abs(waves) ** 2, axis=0)
    time_closure = np.sum(waves * waves, axis=0)
    time_closure_relative = np.abs(time_closure) / np.maximum(time_power, np.finfo(float).tiny)
    initial = waves[:, 0]
    initial_real = initial.real
    initial_imag = initial.imag

    total_spectral_power = float(np.sum(np.abs(spectrum) ** 2))
    singular_values = np.linalg.svd(spectrum.T, compute_uv=False)
    rank_tolerance = float(max(spectrum.T.shape) * eps * singular_values[0])
    independent_direction_count = int(np.count_nonzero(singular_values > rank_tolerance))

    frequency_rows: list[dict[str, Any]] = []
    component_rows: list[dict[str, Any]] = []
    phase_signatures: set[tuple[int, tuple[int, ...]]] = set()

    for index in range(resolution):
        k_signed = signed_frequency(index, resolution)
        vector = spectrum[:, index]
        amplitudes = np.abs(vector)
        phases = np.angle(vector)
        cells, cell_errors = phase_cell(phases, resolution)
        reference = int(np.argmax(amplitudes))
        relative_phases = np.angle(vector * np.exp(-1j * phases[reference]))
        relative_cells, relative_cell_errors = phase_cell(relative_phases, resolution)
        power = float(np.sum(amplitudes**2))
        closure = complex(np.sum(vector * vector))
        closure_relative = float(abs(closure) / max(power, np.finfo(float).tiny))
        period = None if k_signed == 0 else resolution // math.gcd(abs(k_signed), resolution)
        nonzero = amplitudes > phase_amplitude_floor
        absolute_signature = tuple(int(x) if nz else -1 for x, nz in zip(cells, nonzero))
        relative_signature = tuple(
            int(x) if nz else -1 for x, nz in zip(relative_cells, nonzero)
        )
        phase_signatures.add((k_signed, absolute_signature))

        frequency_rows.append(
            {
                "frequency_cell": index,
                "signed_multiple_k": k_signed,
                "role": frequency_role(k_signed),
                "orientation": (
                    "静止" if k_signed == 0 else ("正方向" if k_signed > 0 else "逆方向")
                ),
                "period_in_samples": period,
                "power": power,
                "power_fraction": power / total_spectral_power,
                "closure_abs": float(abs(closure)),
                "closure_relative": closure_relative,
                "individually_zero_closed": closure_relative <= closure_tolerance,
                "reference_relation": reference + 1,
                "reference_edge": edge_label(edges[reference]),
                "absolute_phase_cells": list(absolute_signature),
                "relative_phase_cells": list(relative_signature),
                "max_phase_cell_error_deg": float(
                    np.max(np.abs(np.degrees(cell_errors[nonzero]))) if np.any(nonzero) else 0.0
                ),
                "max_relative_phase_cell_error_deg": float(
                    np.max(np.abs(np.degrees(relative_cell_errors[nonzero])))
                    if np.any(nonzero)
                    else 0.0
                ),
            }
        )

        for relation in range(m_relations):
            component_rows.append(
                {
                    "relation_id": relation + 1,
                    "edge": edge_label(edges[relation]),
                    "frequency_cell": index,
                    "signed_multiple_k": k_signed,
                    "role": frequency_role(k_signed),
                    "amplitude": float(amplitudes[relation]),
                    "power": float(amplitudes[relation] ** 2),
                    "power_fraction_within_relation": float(
                        amplitudes[relation] ** 2
                        / max(float(np.sum(np.abs(spectrum[relation]) ** 2)), np.finfo(float).tiny)
                    ),
                    "absolute_phase_deg": float(np.degrees(phases[relation])),
                    "absolute_phase_cell": int(cells[relation]) if nonzero[relation] else None,
                    "absolute_phase_cell_error_deg": float(np.degrees(cell_errors[relation]))
                    if nonzero[relation]
                    else None,
                    "relative_phase_to_bin_reference_deg": float(
                        np.degrees(relative_phases[relation])
                    )
                    if nonzero[relation]
                    else None,
                    "relative_phase_cell": int(relative_cells[relation])
                    if nonzero[relation]
                    else None,
                }
            )

    relation_rows: list[dict[str, Any]] = []
    for relation in range(m_relations):
        samples = waves[relation]
        coefficients = spectrum[relation]
        powers = np.abs(coefficients) ** 2
        dominant_index = int(np.argmax(powers))
        dominant_signed = signed_frequency(dominant_index, resolution)
        sample_closure = complex(np.sum(samples * samples))
        sample_power = float(np.sum(np.abs(samples) ** 2))
        probability = powers / max(float(np.sum(powers)), np.finfo(float).tiny)
        nonzero_probability = probability[probability > 0.0]
        entropy = float(-np.sum(nonzero_probability * np.log(nonzero_probability)))
        effective_cells = float(np.exp(entropy))
        sample_cells, _ = phase_cell(np.angle(samples), resolution)
        relation_rows.append(
            {
                "relation_id": relation + 1,
                "edge": edge_label(edges[relation]),
                "sample_count": resolution,
                "rms_amplitude": float(np.sqrt(np.mean(np.abs(samples) ** 2))),
                "min_amplitude": float(np.min(np.abs(samples))),
                "max_amplitude": float(np.max(np.abs(samples))),
                "time_phase_cell_count": int(len(np.unique(sample_cells))),
                "resolved_frequency_cell_count": int(np.count_nonzero(np.abs(coefficients) > phase_amplitude_floor)),
                "dominant_frequency_cell": dominant_index,
                "dominant_signed_multiple_k": dominant_signed,
                "dominant_role": frequency_role(dominant_signed),
                "dominant_power_fraction": float(powers[dominant_index] / np.sum(powers)),
                "effective_frequency_cell_count": effective_cells,
                "relation_series_closure_abs": float(abs(sample_closure)),
                "relation_series_closure_relative": float(
                    abs(sample_closure) / max(sample_power, np.finfo(float).tiny)
                ),
            }
        )

    result = {
        "schema": SCHEMA,
        "reader": Path(__file__).name,
        "reader_sha256": reader_sha256,
        "source_manifest_sha256": sha256_file(data["manifest_path"]),
        "input_contract": {
            "n_body": data["n_body"],
            "m_relations": m_relations,
            "resolution_mathcal_N": resolution,
            "lambda0": 2.0 * math.pi / resolution,
            "initial_source": data["manifest"]["source"],
            "initial_construction": data["manifest"]["source_info"]["construction"],
            "window_steps": data["selected_steps"].tolist(),
            "state_shape_relation_by_sample": [m_relations, resolution],
            "frequency_axis_is_posthoc_dft": True,
            "preassigned_frequency_layers": False,
            "preassigned_frequency_weights": False,
        },
        "numerical_rules": {
            "dft": "numpy.fft.fft(axis=sample, norm=ortho)",
            "phase_cell_width_deg": 360.0 / resolution,
            "phase_cell_assignment": "nearest cell modulo resolution",
            "phase_amplitude_floor": phase_amplitude_floor,
            "closure_relative_tolerance": closure_tolerance,
            "rank_tolerance": rank_tolerance,
        },
        "audits": {
            "initial_zero_closure_abs": float(abs(np.sum(initial * initial))),
            "initial_norm": float(np.sum(np.abs(initial) ** 2)),
            "initial_real_norm_squared": float(initial_real @ initial_real),
            "initial_imag_norm_squared": float(initial_imag @ initial_imag),
            "initial_real_imag_dot": float(initial_real @ initial_imag),
            "max_time_zero_closure_abs": float(np.max(np.abs(time_closure))),
            "max_time_zero_closure_relative": float(np.max(time_closure_relative)),
            "time_sample_zero_closed_count": int(
                np.count_nonzero(time_closure_relative <= closure_tolerance)
            ),
            "time_sample_count": resolution,
            "parseval_time_power": float(np.sum(np.abs(waves) ** 2)),
            "parseval_spectrum_power": total_spectral_power,
            "parseval_relative_error": float(
                abs(np.sum(np.abs(waves) ** 2) - total_spectral_power)
                / max(total_spectral_power, np.finfo(float).tiny)
            ),
            "phase_resolved_frequency_component_count": len(phase_signatures),
            "independent_spectrum_direction_count": independent_direction_count,
            "independent_direction_upper_bound_M": m_relations,
            "individually_zero_closed_frequency_count": int(
                sum(row["individually_zero_closed"] for row in frequency_rows)
            ),
        },
        "relation_rows": relation_rows,
        "frequency_rows": frequency_rows,
        "component_rows": component_rows,
    }
    result["_arrays"] = {"waves": waves, "spectrum": spectrum}
    return result


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[label for _, label in fields])
        writer.writeheader()
        for row in rows:
            writer.writerow({label: row[key] for key, label in fields})


def write_outputs(output: Path, result: dict[str, Any]) -> None:
    output.mkdir(parents=True, exist_ok=False)
    arrays = result.pop("_arrays")
    np.save(output / "relation_waves.npy", arrays["waves"])
    np.save(output / "posthoc_spectrum.npy", arrays["spectrum"])
    (output / "census.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    relation_fields = [
        ("relation_id", "関係ID"),
        ("edge", "関係辺"),
        ("sample_count", "標本数"),
        ("rms_amplitude", "実効振幅"),
        ("min_amplitude", "最小振幅"),
        ("max_amplitude", "最大振幅"),
        ("time_phase_cell_count", "時間位相セル数"),
        ("resolved_frequency_cell_count", "検出周波数セル数"),
        ("dominant_frequency_cell", "支配周波数セル"),
        ("dominant_signed_multiple_k", "支配符号付き倍数k"),
        ("dominant_role", "支配成分種別"),
        ("dominant_power_fraction", "支配強度比"),
        ("effective_frequency_cell_count", "実効周波数セル数"),
        ("relation_series_closure_abs", "関係標本列の閉塞絶対値"),
        ("relation_series_closure_relative", "関係標本列の閉塞相対値"),
    ]
    component_fields = [
        ("relation_id", "関係ID"),
        ("edge", "関係辺"),
        ("frequency_cell", "周波数セル"),
        ("signed_multiple_k", "符号付き倍数k"),
        ("role", "成分種別"),
        ("amplitude", "振幅"),
        ("power", "強度"),
        ("power_fraction_within_relation", "関係内強度比"),
        ("absolute_phase_deg", "絶対位相度"),
        ("absolute_phase_cell", "絶対位相セル"),
        ("absolute_phase_cell_error_deg", "絶対位相セル誤差度"),
        ("relative_phase_to_bin_reference_deg", "同周波数基準関係への相対位相度"),
        ("relative_phase_cell", "相対位相セル"),
    ]
    frequency_fields = [
        ("frequency_cell", "周波数セル"),
        ("signed_multiple_k", "符号付き倍数k"),
        ("role", "成分種別"),
        ("orientation", "回転方向"),
        ("period_in_samples", "離散周期"),
        ("power", "全関係強度"),
        ("power_fraction", "全体強度比"),
        ("closure_abs", "周波数成分の閉塞絶対値"),
        ("closure_relative", "周波数成分の閉塞相対値"),
        ("individually_zero_closed", "周波数成分が単独零閉塞"),
        ("reference_relation", "位相基準関係ID"),
        ("reference_edge", "位相基準辺"),
        ("absolute_phase_cells", "全関係の絶対位相セル署名"),
        ("relative_phase_cells", "全関係の相対位相セル署名"),
        ("max_phase_cell_error_deg", "最大絶対位相セル誤差度"),
        ("max_relative_phase_cell_error_deg", "最大相対位相セル誤差度"),
    ]
    write_csv(output / "relation_waves_ja.csv", result["relation_rows"], relation_fields)
    write_csv(output / "spectral_components_ja.csv", result["component_rows"], component_fields)
    write_csv(output / "frequency_cells_ja.csv", result["frequency_rows"], frequency_fields)

    contract = result["input_contract"]
    audit = result["audits"]
    lines = [
        f"# N={contract['n_body']} 倍音段数を事前指定しない位相分解波一覧 v4",
        "",
        "## 1. 実装契約",
        "",
        f"- 状態配列は $X\\in\\mathbb{{C}}^{{{contract['m_relations']}\\times{contract['resolution_mathcal_N']}}}$。行は関係、列は分解能標本。",
        "- 周波数成分は生成前に置かず、この標本列へ事後DFTを施して全セルを読む。",
        "- 各周波数セルについて全関係成分の振幅、絶対位相、相対位相を保存し、位相の異なる成分を統合しない。",
        "- `relation_waves.npy` が時間領域の正本、`posthoc_spectrum.npy` が事後読出し結果。",
        "",
        "## 2. 検証結果",
        "",
        "| 項目 | 結果 |",
        "|---|---:|",
        f"| N体数 $N$ | {contract['n_body']} |",
        f"| 関係数 $M=N(N-1)/2$ | {contract['m_relations']} |",
        f"| 分解能 $\\mathcal N$ | {contract['resolution_mathcal_N']} |",
        f"| 状態配列形状 | {contract['m_relations']} × {contract['resolution_mathcal_N']} |",
        f"| 初期入力 | {contract['initial_source']} |",
        f"| 初期入力の構成 | `{contract['initial_construction']}` |",
        f"| 初期実部ノルム二乗 | {audit['initial_real_norm_squared']:.16f} |",
        f"| 初期虚部ノルム二乗 | {audit['initial_imag_norm_squared']:.16f} |",
        f"| 初期実部・虚部内積 | {audit['initial_real_imag_dot']:.3e} |",
        f"| 初期零閉塞絶対残差 | {audit['initial_zero_closure_abs']:.3e} |",
        f"| 零閉塞した時点 | {audit['time_sample_zero_closed_count']} / {audit['time_sample_count']} |",
        f"| 最大時点零閉塞相対残差 | {audit['max_time_zero_closure_relative']:.3e} |",
        f"| 位相を保存した周波数成分 | {audit['phase_resolved_frequency_component_count']} |",
        f"| 独立スペクトル方向数 | {audit['independent_spectrum_direction_count']}（上限M={audit['independent_direction_upper_bound_M']}） |",
        f"| 単独零閉塞した周波数成分 | {audit['individually_zero_closed_frequency_count']} |",
        f"| Parseval相対誤差 | {audit['parseval_relative_error']:.3e} |",
        "",
        f"## 3. 関係波{contract['m_relations']}本の概要",
        "",
        "| 関係ID | 辺 | 実効振幅 | 時間位相セル数 | 検出周波数セル数 | 支配k | 種別 | 支配強度比 | 実効周波数セル数 | 関係標本列の閉塞相対値 |",
        "|---:|:---:|---:|---:|---:|---:|:---|---:|---:|---:|",
    ]
    for row in result["relation_rows"]:
        lines.append(
            f"| {row['relation_id']} | {row['edge']} | {fmt(row['rms_amplitude'])} | "
            f"{row['time_phase_cell_count']} | {row['resolved_frequency_cell_count']} | "
            f"{row['dominant_signed_multiple_k']:+d} | {row['dominant_role']} | "
            f"{100.0 * row['dominant_power_fraction']:.4f}% | "
            f"{row['effective_frequency_cell_count']:.4f} | "
            f"{row['relation_series_closure_relative']:.3e} |"
        )

    lines.extend(
        [
            "",
            "## 4. 全周波数セル",
            "",
            "以下は1セル1行であり、基底波と倍音を一行へ合算していない。各セルのM成分位相は `frequency_cells_ja.csv`、各関係まで展開した全成分は `spectral_components_ja.csv` に保存した。",
            "",
            "| セル | k | 種別 | 方向 | 周期 | 全体強度比 | 単独零閉塞相対値 | 単独零閉塞 | 位相基準辺 | 相対位相セル署名 |",
            "|---:|---:|:---|:---|---:|---:|---:|:---:|:---:|:---|",
        ]
    )
    for row in result["frequency_rows"]:
        period = "—" if row["period_in_samples"] is None else str(row["period_in_samples"])
        closed = "はい" if row["individually_zero_closed"] else "いいえ"
        signature = ",".join(str(x) for x in row["relative_phase_cells"])
        lines.append(
            f"| {row['frequency_cell']} | {row['signed_multiple_k']:+d} | {row['role']} | "
            f"{row['orientation']} | {period} | {100.0 * row['power_fraction']:.6f}% | "
            f"{row['closure_relative']:.3e} | {closed} | {row['reference_edge']} | `{signature}` |"
        )

    lines.extend(
        [
            "",
            "## 5. 解釈上の境界",
            "",
            "- 時点ごとの零閉塞と、DFTで分離した各周波数成分の単独零閉塞は別の検査である。後者が成立しなければ、そのセルだけを一つの閉じた波とは呼ばない。",
            "- 周波数セル数は分解能による読出し住所の数であり、物理的な波数を事前指定した値ではない。",
            "- 位相署名は絶対位相セルと、同一周波数内の最大振幅関係を基準にした相対位相セルの両方を保存した。",
            "- 本表は波と倍音の実測表であり、粒子名・電荷・質量・スピン・寿命を未導出のまま付与しない。",
            "",
        ]
    )
    (output / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    input_dir = args.input.resolve()
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"出力先が既に存在するため上書きしない: {output}")
    data = load_input(input_dir, args.window_start)
    result = analyse(data, sha256_file(Path(__file__).resolve()))
    write_outputs(output, result)
    audit = result["audits"]
    print(
        json.dumps(
            {
                "output": str(output),
                "state_shape": result["input_contract"]["state_shape_relation_by_sample"],
                "time_zero_closed": f"{audit['time_sample_zero_closed_count']}/{audit['time_sample_count']}",
                "phase_resolved_frequency_components": audit[
                    "phase_resolved_frequency_component_count"
                ],
                "independent_spectrum_directions": audit[
                    "independent_spectrum_direction_count"
                ],
                "individually_zero_closed_frequencies": audit[
                    "individually_zero_closed_frequency_count"
                ],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
