#!/usr/bin/env python3
"""倍音対応 make_parent の位相分解・関係波一覧表 v3。

状態 Z in C^(M x H) を、生成器自身の契約どおり

    行 e : 一つの関係波（完全グラフ K_N の辺）
    列 n : 同じ関係波が持つ n 倍音

として読む。一つの代表位相へ潰さず、各関係波について全倍音の振幅と位相を
保存する。有限分解能 mathcal N では、位相を幅 2*pi/mathcal N のセルへ丸め、
全倍音位相セル列が異なる関係波を別の位相分解波として数える。

本器は初期生成状態だけを読む。緩和、塊閾値、粒子名の手置きは行わない。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
GENERATOR = REPO / "次元の生成構造" / "make_parent_harmonic_unit_v1" / "make_parent_harmonic_v1.py"
SCHEMA = "phase_resolved_relation_wave_census_v3"
DEFAULT_CASES = [(5, 8, 40260801), (40, 4, 40260802)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Z[M,H] を一行一関係波として、全倍音の振幅・位相を一覧化する。"
    )
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        metavar="N:H:SEED",
        help="計算条件。省略時は 5:8:40260801 と 40:4:40260802。",
    )
    parser.add_argument(
        "--resolution",
        type=int,
        default=144,
        help="位相分解能 mathcal N。位相セル幅は 2*pi/mathcal N。",
    )
    parser.add_argument("--output", type=Path, required=True, help="新規出力ディレクトリ")
    parser.add_argument(
        "--summary-waves",
        type=int,
        default=20,
        help="summary.md に強度順で表示する関係波数。全件は各Nのwave_table.md/CSV。",
    )
    return parser.parse_args()


def parse_cases(raw_cases: list[str]) -> list[tuple[int, int, int]]:
    if not raw_cases:
        return DEFAULT_CASES.copy()
    parsed = []
    for raw in raw_cases:
        try:
            n_text, h_text, seed_text = raw.split(":")
            case = (int(n_text), int(h_text), int(seed_text))
        except ValueError as exc:
            raise SystemExit(f"--case は N:H:SEED 形式です: {raw!r}") from exc
        if case[0] < 2 or case[1] < 1:
            raise SystemExit(f"N>=2, H>=1 が必要です: {raw!r}")
        parsed.append(case)
    if len({case[0] for case in parsed}) != len(parsed):
        raise SystemExit("同じNを複数回指定できません")
    return parsed


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_generator():
    spec = importlib.util.spec_from_file_location("make_parent_harmonic_census_v3", GENERATOR)
    if spec is None or spec.loader is None:
        raise SystemExit(f"生成器を読み込めません: {GENERATOR}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def wrap_phase(value: np.ndarray | float) -> np.ndarray | float:
    return (value + np.pi) % (2.0 * np.pi) - np.pi


def phase_cells(phases: np.ndarray, resolution: int) -> tuple[np.ndarray, np.ndarray]:
    delta = 2.0 * np.pi / resolution
    cells = np.floor((np.mod(phases, 2.0 * np.pi) / delta) + 0.5).astype(int) % resolution
    centers = cells * delta
    errors = wrap_phase(phases - centers)
    return cells, np.asarray(errors, dtype=float)


def phase_network_readout(values: np.ndarray, cells: np.ndarray, resolution: int) -> dict[str, Any]:
    amplitudes = np.abs(values)
    phases = np.angle(values)
    amplitude_sum = float(np.sum(amplitudes))
    coherence = float(abs(np.sum(values)) / max(amplitude_sum, np.finfo(float).tiny))

    pair_weight = 0.0
    cancellation_weight = 0.0
    for left in range(len(values)):
        for right in range(left + 1, len(values)):
            weight = float(amplitudes[left] * amplitudes[right])
            pair_weight += weight
            cancellation_weight += weight * math.sin(
                float(wrap_phase(phases[left] - phases[right])) / 2.0
            ) ** 2
    cancellation = cancellation_weight / pair_weight if pair_weight else 0.0

    if len(values) == 1:
        network_type = "単一倍音（位相差なし）"
        network_label = "—"
    elif np.all(cells == cells[0]):
        network_type = "同相型（ボゾン的極限）"
        network_label = "B"
    elif len(values) == 2 and resolution % 2 == 0 and (
        int(cells[0] - cells[1]) % resolution == resolution // 2
    ):
        network_type = "逆相二成分型（フェルミオン的極限）"
        network_label = "F"
    else:
        network_type = "中間位相差型（エルミオン的）"
        network_label = "E"
    return {
        "phase_network_BFE": network_label,
        "phase_network_type": network_type,
        "coherence_rho": coherence,
        "weighted_cancellation_F": float(cancellation),
    }


def relation_wave_rows(
    z_matrix: np.ndarray,
    n_body: int,
    resolution: int,
) -> list[dict[str, Any]]:
    m_relations, harmonic_count = z_matrix.shape
    edge_a, edge_b = np.triu_indices(n_body, k=1)
    if len(edge_a) != m_relations:
        raise RuntimeError("Mと完全グラフの辺数が一致しません")

    amplitudes = np.abs(z_matrix)
    phases = np.angle(z_matrix)
    cells, errors = phase_cells(phases, resolution)
    signatures = [tuple(int(x) for x in row) for row in cells]
    unique_signatures = {signature: index + 1 for index, signature in enumerate(sorted(set(signatures)))}

    rows: list[dict[str, Any]] = []
    for edge_index in range(m_relations):
        values = np.asarray(z_matrix[edge_index], dtype=np.complex128)
        amps = amplitudes[edge_index]
        power = amps**2
        total_power = float(np.sum(power))
        fractions = power / max(total_power, np.finfo(float).tiny)
        dominant = int(np.argmax(power)) + 1
        closure = complex(np.sum(values * values))
        closure_relative = float(abs(closure) / max(total_power, np.finfo(float).tiny))
        network = phase_network_readout(values, cells[edge_index], resolution)
        response = complex(np.sum(values) / max(float(np.sum(amps)), np.finfo(float).tiny))
        odd_power = float(np.sum(power[np.arange(harmonic_count) % 2 == 0]))
        even_power = float(np.sum(power[np.arange(harmonic_count) % 2 == 1]))
        parity_total = odd_power + even_power
        odd_fraction = odd_power / parity_total if parity_total else 0.0
        even_fraction = even_power / parity_total if parity_total else 0.0
        signature = signatures[edge_index]
        internal_relative = [
            float(wrap_phase(phases[edge_index, harmonic - 1] - harmonic * phases[edge_index, 0]))
            for harmonic in range(1, harmonic_count + 1)
        ]
        rows.append(
            {
                "wave_id": f"W{edge_index + 1:04d}",
                "edge_index": edge_index,
                "vertex_i": int(edge_a[edge_index]) + 1,
                "vertex_j": int(edge_b[edge_index]) + 1,
                "phase_class_id": f"P{unique_signatures[signature]:04d}",
                "phase_signature": list(signature),
                "harmonic_count": harmonic_count,
                "all_harmonics_nonzero": bool(np.all(amps > 0.0)),
                "amplitudes": [float(x) for x in amps],
                "power_fractions": [float(x) for x in fractions],
                "phases_rad": [float(x) for x in phases[edge_index]],
                "phases_deg": [float(np.degrees(x)) for x in phases[edge_index]],
                "phase_cells": [int(x) for x in cells[edge_index]],
                "phase_cell_errors_rad": [float(x) for x in errors[edge_index]],
                "time_origin_invariant_phases_rad": internal_relative,
                "dominant_harmonic_n": dominant,
                "dominant_power_fraction": float(fractions[dominant - 1]),
                "instantaneous_amplitude_squared_R0": total_power,
                "mass_like_status": "単時点の振幅二乗量。多ゲージ安定性未測定",
                "row_zero_closure_abs": float(abs(closure)),
                "row_zero_closure_relative": closure_relative,
                "row_zero_closed": bool(closure_relative <= 1e-10),
                "odd_harmonic_power_fraction": odd_fraction,
                "even_harmonic_power_fraction": even_fraction,
                "harmonic_parity_balance": float(odd_fraction - even_fraction),
                "external_axis_response_S0": float(response.real),
                "external_axis_response_Spi_over_2": float(response.imag),
                "external_axis_response_strength": float(abs(response)),
                "spin_status": "外部軸相関応答。内在スピン固有値ではない",
                "charge_readout": None,
                "charge_status": "初期生成契約に相互作用・観測波がなく未読",
                "lifetime_steps": None,
                "lifetime_status": "倍音段間を結合する時間発展がなく未測定",
                "conjugate_phase_candidate_id": None,
                "conjugate_amplitude_residual": None,
                **network,
            }
        )

    signature_to_indices: dict[tuple[int, ...], list[int]] = {}
    for index, signature in enumerate(signatures):
        signature_to_indices.setdefault(signature, []).append(index)
    for index, row in enumerate(rows):
        target = tuple((-cell) % resolution for cell in signatures[index])
        candidates = [candidate for candidate in signature_to_indices.get(target, []) if candidate != index]
        if not candidates:
            continue
        source_amplitude = amplitudes[index] / max(float(np.linalg.norm(amplitudes[index])), np.finfo(float).tiny)
        best = min(
            candidates,
            key=lambda candidate: float(
                np.linalg.norm(
                    source_amplitude
                    - amplitudes[candidate]
                    / max(float(np.linalg.norm(amplitudes[candidate])), np.finfo(float).tiny)
                )
            ),
        )
        best_amplitude = amplitudes[best] / max(float(np.linalg.norm(amplitudes[best])), np.finfo(float).tiny)
        row["conjugate_phase_candidate_id"] = rows[best]["wave_id"]
        row["conjugate_amplitude_residual"] = float(np.linalg.norm(source_amplitude - best_amplitude))
    return rows


def harmonic_layer_rows(
    z_matrix: np.ndarray,
    info: dict[str, Any],
    resolution: int,
) -> list[dict[str, Any]]:
    cells, _ = phase_cells(np.angle(z_matrix), resolution)
    rows = []
    for harmonic_index in range(z_matrix.shape[1]):
        values = np.asarray(z_matrix[:, harmonic_index], dtype=np.complex128)
        level = info["levels"][harmonic_index]
        power = float(np.sum(np.abs(values) ** 2))
        closure = complex(np.sum(values * values))
        harmonic_n = harmonic_index + 1
        rows.append(
            {
                "harmonic_n": harmonic_n,
                "omega_over_omega0": harmonic_n,
                "angular_step_rad": float(2.0 * np.pi * harmonic_n / resolution),
                "discrete_period_steps": int(resolution // math.gcd(harmonic_n, resolution)),
                "generator_branch": level["branch"],
                "sigma1": float(level["sigma1"]),
                "family": (
                    "N-1"
                    if abs(float(level["sigma1"]) - (info["n_vertices"] - 1)) < 1e-9
                    else "broken"
                ),
                "layer_power": power,
                "layer_closure_abs": float(abs(closure)),
                "layer_closure_relative": float(abs(closure) / max(power, np.finfo(float).tiny)),
                "phase_resolved_component_count": int(len(set(int(x) for x in cells[:, harmonic_index]))),
                "occupied_phase_cells": sorted(set(int(x) for x in cells[:, harmonic_index])),
            }
        )
    return rows


def analyze_case(generator, n_body: int, harmonic_count: int, seed: int, resolution: int) -> dict[str, Any]:
    z_matrix, info = generator.make_parent_harmonic(
        n_body,
        harmonic_count,
        seed,
        iters=2000,
        restarts=10,
        tol=1e-12,
    )
    waves = relation_wave_rows(z_matrix, n_body, resolution)
    layers = harmonic_layer_rows(z_matrix, info, resolution)
    signatures = {tuple(row["phase_signature"]) for row in waves}
    phase_class_sizes: dict[str, int] = {}
    for row in waves:
        phase_class_sizes[row["phase_class_id"]] = phase_class_sizes.get(row["phase_class_id"], 0) + 1
    return {
        "schema": SCHEMA,
        "n_body": n_body,
        "m_relations": int(z_matrix.shape[0]),
        "harmonic_register_count_H": harmonic_count,
        "seed": seed,
        "resolution_mathcal_N": resolution,
        "phase_cell_width_rad": float(2.0 * np.pi / resolution),
        "phase_cell_width_deg": float(360.0 / resolution),
        "base_angular_step_omega0": float(2.0 * np.pi / resolution),
        "generator_info": info,
        "summary": {
            "relation_wave_slot_count_M": int(z_matrix.shape[0]),
            "phase_resolved_relation_wave_count_W": len(signatures),
            "phase_signature_collision_count": int(z_matrix.shape[0] - len(signatures)),
            "harmonic_components_total_M_times_H": int(z_matrix.size),
            "harmonics_per_relation_wave": harmonic_count,
            "harmonic_register_coverage_fraction": float(harmonic_count / resolution),
            "all_relation_waves_have_all_harmonics": bool(
                all(row["all_harmonics_nonzero"] for row in waves)
            ),
            "individually_row_zero_closed_wave_count": sum(int(row["row_zero_closed"]) for row in waves),
            "zero_closed_harmonic_layer_count": sum(
                int(row["layer_closure_relative"] <= 1e-10) for row in layers
            ),
            "conjugate_phase_candidate_wave_count": sum(
                int(row["conjugate_phase_candidate_id"] is not None) for row in waves
            ),
            "phase_class_sizes": phase_class_sizes,
        },
        "harmonic_layers": layers,
        "relation_waves": waves,
    }


def fmt(value: float | None, digits: int = 6) -> str:
    if value is None:
        return "—"
    return f"{float(value):.{digits}g}"


def spectrum_text(row: dict[str, Any]) -> str:
    return "; ".join(
        f"n={index + 1}: {row['amplitudes'][index]:.5f}∠{row['phases_deg'][index]:+.2f}°"
        f"[q={row['phase_cells'][index]}]"
        for index in range(row["harmonic_count"])
    )


def case_markdown(result: dict[str, Any]) -> str:
    n_body = result["n_body"]
    summary = result["summary"]
    lines = [
        f"# 位相分解・関係波一覧 — N={n_body}",
        "",
        "## 数え上げ",
        "",
        "| N体数 | 関係波枠M | 位相分解波W | 位相署名の重複 | 倍音段H | 倍音成分総数M×H | 行ごとに零閉鎖する波 | 零閉鎖する倍音段 |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        f"| {n_body} | {result['m_relations']} | {summary['phase_resolved_relation_wave_count_W']} | "
        f"{summary['phase_signature_collision_count']} | {result['harmonic_register_count_H']} | "
        f"{summary['harmonic_components_total_M_times_H']} | "
        f"{summary['individually_row_zero_closed_wave_count']} | "
        f"{summary['zero_closed_harmonic_layer_count']} |",
        "",
        f"位相セル幅は $2\\pi/{result['resolution_mathcal_N']}$ = "
        f"{result['phase_cell_width_deg']:.6g}°。全倍音の位相セル列が違えば別の波として数えた。",
        "`W` は関係波枠ごとの全倍音位相署名の種類数なので、必ず $W\\le M$ である。",
        f"この計算で実装された倍音段は H={result['harmonic_register_count_H']} で、"
        f"分解能 $\\mathcal N={result['resolution_mathcal_N']}$ の "
        f"{100.0 * summary['harmonic_register_coverage_fraction']:.3f}% にすぎない。",
        "",
        "## 倍音段の閉鎖表",
        "",
        "| 倍音n | 周波数nω0 | 離散周期 | 族 | 生成枝 | σ1 | 位相セル数 | 段の零閉鎖相対残差 |",
        "| ---: | ---: | ---: | --- | --- | ---: | ---: | ---: |",
    ]
    for row in result["harmonic_layers"]:
        lines.append(
            f"| {row['harmonic_n']} | {row['omega_over_omega0']} | "
            f"{row['discrete_period_steps']} | {row['family']} | {row['generator_branch']} | "
            f"{row['sigma1']:.9g} | {row['phase_resolved_component_count']} | "
            f"{row['layer_closure_relative']:.3e} |"
        )
    lines.extend(
        [
            "",
            "## 関係波の全件表",
            "",
            "一行が一つの関係波である。同じ行の `スペクトル` 欄に、全倍音の振幅・位相・位相セルを列挙する。",
            "代表位相一つへの平均は行っていない。",
            "",
            "| 波ID | 関係辺 | 位相類 | 倍音数 | スペクトル A∠φ[q] | 支配倍音 | 支配強度 | 振幅二乗R0 | 行零閉鎖残差 | B/F/E | 位相差型 | ρ | F | S(0) | S(π/2) | 共役位相候補 |",
            "| --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in result["relation_waves"]:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["wave_id"],
                    f"({row['vertex_i']},{row['vertex_j']})",
                    row["phase_class_id"],
                    str(row["harmonic_count"]),
                    spectrum_text(row),
                    str(row["dominant_harmonic_n"]),
                    f"{100.0 * row['dominant_power_fraction']:.3f}%",
                    fmt(row["instantaneous_amplitude_squared_R0"], 8),
                    f"{row['row_zero_closure_relative']:.3e}",
                    row["phase_network_BFE"],
                    row["phase_network_type"],
                    fmt(row["coherence_rho"], 6),
                    fmt(row["weighted_cancellation_F"], 6),
                    fmt(row["external_axis_response_S0"], 6),
                    fmt(row["external_axis_response_Spi_over_2"], 6),
                    row["conjugate_phase_candidate_id"] or "—",
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## 読出し境界",
            "",
            "- `振幅二乗R0` は単時点の $\\sum_n|Z_{e,n}|^2$。質量型候補に必要な安定性・多ゲージ再構成は未測定。",
            "- `S(0), S(π/2)` は外部軸との相関応答であり、内在スピン固有値ではない。",
            "- 電荷と寿命は、現生成器に倍音段間相互作用を含む時間発展がないため、この初期表からは読めない。",
            "- `行零閉鎖残差` は各関係波について $|\\sum_n Z_{e,n}^2|/\\sum_n|Z_{e,n}|^2$ を直接測った値。",
            "",
        ]
    )
    return "\n".join(lines)


def summary_markdown(results: list[dict[str, Any]], summary_waves: int) -> str:
    lines = [
        "# N=5・N=40 位相分解波一覧 v3",
        "",
        "## 表の単位",
        "",
        "倍音対応生成器の $Z\\in\\mathbb C^{M\\times H}$ を、行＝関係波、列＝同じ波の倍音として読んだ。",
        "位相が異なる波を一つへ平均せず、各関係波の全倍音位相を保存している。",
        "",
        "| N体数 | M | 位相分解波W | H | 一波あたり倍音数 | M×H | 行零閉鎖波 | 零閉鎖倍音段 | 共役位相候補波 |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for result in results:
        summary = result["summary"]
        lines.append(
            f"| {result['n_body']} | {result['m_relations']} | "
            f"{summary['phase_resolved_relation_wave_count_W']} | "
            f"{result['harmonic_register_count_H']} | {summary['harmonics_per_relation_wave']} | "
            f"{summary['harmonic_components_total_M_times_H']} | "
            f"{summary['individually_row_zero_closed_wave_count']} | "
            f"{summary['zero_closed_harmonic_layer_count']} | "
            f"{summary['conjugate_phase_candidate_wave_count']} |"
        )
    lines.extend(
        [
            "",
            "`W` は分解能 $\\mathcal N=144$ で量子化した全倍音位相署名の種類数で、$W\\le M$。",
            "`M×H` は波数ではなく、関係波が持つ倍音成分の総数である。",
            "",
            "## 生成器について直接分かったこと",
            "",
            "- 現生成器は白色雑音を入力していない。各倍音段を別々の自己無撞着円偏波閉包として作り、段振幅を $1/\\sqrt H$ に揃えている。",
            "- 分解能は $\\mathcal N=144$ だが、実際のレジスタはN=5でH=8、N=40でH=4の打切りであり、許容倍音全体ではない。",
            "- 各倍音列は零閉鎖している。一方、行を一つの多倍音関係波として $\\sum_n Z_{e,n}^2=0$ を要求すると、該当は両Nとも0本である。",
            "- したがって現生成器は『零閉鎖した倍音段の等振幅集合』は作るが、『各関係波がそれ自身で零閉鎖した白色多倍音波の集合』は作っていない。",
            "",
            f"## 各Nの強度上位{summary_waves}波",
            "",
            "| N | 波ID | 関係辺 | 位相類 | 全倍音位相セル | 支配倍音 | 振幅二乗R0 | 行零閉鎖残差 | B/F/E | 位相差型 |",
            "| --- | --- | --- | --- | --- | ---: | ---: | ---: | --- | --- |",
        ]
    )
    for result in results:
        order = sorted(
            result["relation_waves"],
            key=lambda row: row["instantaneous_amplitude_squared_R0"],
            reverse=True,
        )[:summary_waves]
        for row in order:
            lines.append(
                f"| {result['n_body']} | {row['wave_id']} | ({row['vertex_i']},{row['vertex_j']}) | "
                f"{row['phase_class_id']} | {','.join(str(x) for x in row['phase_signature'])} | "
                f"{row['dominant_harmonic_n']} | {row['instantaneous_amplitude_squared_R0']:.8g} | "
                f"{row['row_zero_closure_relative']:.3e} | {row['phase_network_BFE']} | "
                f"{row['phase_network_type']} |"
            )
    lines.extend(
        [
            "",
            "全関係波と全倍音位相は、各Nの `wave_table.md` と `relation_waves_ja.csv` に保存した。",
            "",
        ]
    )
    return "\n".join(lines)


def write_csv(result: dict[str, Any], path: Path) -> None:
    harmonic_count = result["harmonic_register_count_H"]
    fixed_fields = [
        "N体数",
        "波ID",
        "辺番号",
        "頂点i",
        "頂点j",
        "位相類ID",
        "倍音数H",
        "支配倍音n",
        "支配倍音強度比",
        "振幅二乗R0",
        "行零閉鎖絶対残差",
        "行零閉鎖相対残差",
        "行零閉鎖",
        "B/F/E",
        "位相差型",
        "位相整列度rho",
        "打消し度F",
        "奇数倍音強度比",
        "偶数倍音強度比",
        "外部軸応答S0",
        "外部軸応答S(pi/2)",
        "共役位相候補",
        "共役振幅残差",
        "質量型判定",
        "電荷判定",
        "スピン判定",
        "寿命判定",
    ]
    harmonic_fields = []
    for harmonic in range(1, harmonic_count + 1):
        harmonic_fields.extend(
            [
                f"倍音{harmonic}振幅",
                f"倍音{harmonic}強度比",
                f"倍音{harmonic}位相度",
                f"倍音{harmonic}位相セル",
                f"倍音{harmonic}セル誤差rad",
                f"倍音{harmonic}内部相対位相rad",
            ]
        )
    fields = fixed_fields + harmonic_fields
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in result["relation_waves"]:
            output: dict[str, Any] = {
                "N体数": result["n_body"],
                "波ID": row["wave_id"],
                "辺番号": row["edge_index"],
                "頂点i": row["vertex_i"],
                "頂点j": row["vertex_j"],
                "位相類ID": row["phase_class_id"],
                "倍音数H": row["harmonic_count"],
                "支配倍音n": row["dominant_harmonic_n"],
                "支配倍音強度比": row["dominant_power_fraction"],
                "振幅二乗R0": row["instantaneous_amplitude_squared_R0"],
                "行零閉鎖絶対残差": row["row_zero_closure_abs"],
                "行零閉鎖相対残差": row["row_zero_closure_relative"],
                "行零閉鎖": "はい" if row["row_zero_closed"] else "いいえ",
                "B/F/E": row["phase_network_BFE"],
                "位相差型": row["phase_network_type"],
                "位相整列度rho": row["coherence_rho"],
                "打消し度F": row["weighted_cancellation_F"],
                "奇数倍音強度比": row["odd_harmonic_power_fraction"],
                "偶数倍音強度比": row["even_harmonic_power_fraction"],
                "外部軸応答S0": row["external_axis_response_S0"],
                "外部軸応答S(pi/2)": row["external_axis_response_Spi_over_2"],
                "共役位相候補": row["conjugate_phase_candidate_id"] or "",
                "共役振幅残差": row["conjugate_amplitude_residual"],
                "質量型判定": row["mass_like_status"],
                "電荷判定": row["charge_status"],
                "スピン判定": row["spin_status"],
                "寿命判定": row["lifetime_status"],
            }
            for index in range(harmonic_count):
                harmonic = index + 1
                output[f"倍音{harmonic}振幅"] = row["amplitudes"][index]
                output[f"倍音{harmonic}強度比"] = row["power_fractions"][index]
                output[f"倍音{harmonic}位相度"] = row["phases_deg"][index]
                output[f"倍音{harmonic}位相セル"] = row["phase_cells"][index]
                output[f"倍音{harmonic}セル誤差rad"] = row["phase_cell_errors_rad"][index]
                output[f"倍音{harmonic}内部相対位相rad"] = row[
                    "time_origin_invariant_phases_rad"
                ][index]
            writer.writerow(output)


def main() -> None:
    args = parse_args()
    cases = parse_cases(args.case)
    if args.resolution < 2:
        raise SystemExit("--resolution は2以上が必要です")
    if args.summary_waves < 1:
        raise SystemExit("--summary-waves は1以上が必要です")
    if args.output.exists():
        raise SystemExit(f"出力先は既に存在します。上書きしません: {args.output}")

    generator = load_generator()
    results = [
        analyze_case(generator, n_body, harmonic_count, seed, args.resolution)
        for n_body, harmonic_count, seed in cases
    ]
    results.sort(key=lambda result: result["n_body"])
    args.output.mkdir(parents=True, exist_ok=False)

    for result in results:
        case_dir = args.output / f"N{result['n_body']}"
        case_dir.mkdir()
        (case_dir / "census.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        (case_dir / "wave_table.md").write_text(
            case_markdown(result) + "\n", encoding="utf-8"
        )
        write_csv(result, case_dir / "relation_waves_ja.csv")

    combined = {
        "schema": SCHEMA,
        "reader": Path(__file__).name,
        "reader_sha256": sha256(Path(__file__)),
        "generator": str(GENERATOR.resolve()),
        "generator_sha256": sha256(GENERATOR),
        "results": results,
    }
    (args.output / "summary.json").write_text(
        json.dumps(combined, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (args.output / "summary.md").write_text(
        summary_markdown(results, args.summary_waves) + "\n", encoding="utf-8"
    )
    print(
        " / ".join(
            f"N={result['n_body']}: M={result['m_relations']}, "
            f"W={result['summary']['phase_resolved_relation_wave_count_W']}, "
            f"H={result['harmonic_register_count_H']}"
            for result in results
        )
    )
    print(f"saved: {args.output.resolve()}")


if __name__ == "__main__":
    main()
