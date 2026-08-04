#!/usr/bin/env python3
"""n_only_multifrequency_parent_v1 の独立・位相分解一覧器。

生成器をimportせず、保存契約だけを読む。1行は一つの周波数波であり、
その波を構成する全M関係成分の振幅・絶対位相・相対位相を別CSVへ全件保存する。
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Nだけから生成された親波の位相分解一覧")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def phase_cells(phases: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    width = 2.0 * math.pi / n
    wrapped = np.mod(phases, 2.0 * math.pi)
    cells = np.mod(np.rint(wrapped / width).astype(np.int64), n)
    error = np.angle(np.exp(1j * (wrapped - cells * width)))
    return cells, error


def edge_label(edge: np.ndarray) -> str:
    return f"({int(edge[0]) + 1},{int(edge[1]) + 1})"


def read_attempt_modes(manifest: dict[str, Any]) -> dict[int, dict[str, Any]]:
    successful = [attempt for attempt in manifest["attempts"] if attempt["success"]]
    if len(successful) != 1:
        raise SystemExit(f"成功試行が一意でない: {len(successful)}")
    return {int(mode["frequency_multiple"]): mode for mode in successful[0]["modes"]}


def analyse(input_dir: Path) -> dict[str, Any]:
    manifest = json.loads((input_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != "n_only_multifrequency_parent_v1":
        raise SystemExit(f"未対応schema: {manifest.get('schema')!r}")
    if manifest.get("status") != "success":
        raise SystemExit("成功したmake_parent出力ではない")
    n = int(manifest["function_contract"]["N"])
    m = int(manifest["derived"]["M"])
    if m != n * (n - 1) // 2:
        raise SystemExit("M != N(N-1)/2")
    parent = np.load(input_dir / "parent_modes.npy")
    white = np.load(input_dir / "white_input.npy")
    weights = np.load(input_dir / "mode_weights.npy")
    edges = np.load(input_dir / "edges.npy")
    if parent.shape != (m, n) or white.shape != (m, n):
        raise SystemExit(f"波行列shape不正: parent={parent.shape}, white={white.shape}")
    if weights.shape != (n,) or edges.shape != (m, 2):
        raise SystemExit("weight/edge shape不正")

    mode_audit = read_attempt_modes(manifest)
    total_power = float(np.sum(np.abs(parent) ** 2))
    tolerance = float(manifest["fixed_solver_contract"]["tolerance"])
    waves: list[dict[str, Any]] = []
    components: list[dict[str, Any]] = []

    for column in range(n):
        multiple = column + 1
        vector = parent[:, column]
        amplitudes = np.abs(vector)
        phases = np.angle(vector)
        cells, cell_errors = phase_cells(phases, n)
        reference = int(np.argmax(amplitudes))
        relative = np.angle(vector * np.exp(-1j * phases[reference]))
        relative_cells, relative_errors = phase_cells(relative, n)
        power = float(np.sum(amplitudes**2))
        closure = float(abs(complex(np.sum(vector * vector))))
        closure_relative = closure / max(power, np.finfo(float).tiny)
        phase_increment = 2.0 * math.pi * multiple / n
        order = n // math.gcd(n, multiple)
        audit = mode_audit[multiple]
        wave_id = f"W{multiple:04d}"
        waves.append(
            {
                "wave_id": wave_id,
                "frequency_multiple": multiple,
                "kind": "基底波" if multiple == 1 else "倍音",
                "phase_increment_rad": phase_increment,
                "phase_increment_deg": math.degrees(phase_increment),
                "finite_order": order,
                "amplitude": math.sqrt(power),
                "power": power,
                "power_fraction": power / total_power,
                "white_input_weight": float(weights[column]),
                "solver_iterations": int(audit["iterations"]),
                "solver_residual": float(audit["residual"]),
                "sigma_max": float(audit["sigma_max"]),
                "closure_abs": closure,
                "closure_relative": closure_relative,
                "zero_closed": bool(closure_relative < tolerance),
                "relation_component_count": m,
                "nonzero_relation_component_count": int(np.count_nonzero(amplitudes > 0.0)),
                "reference_relation_id": reference + 1,
                "reference_edge": edge_label(edges[reference]),
                "absolute_phase_cell_signature": [int(x) for x in cells],
                "relative_phase_cell_signature": [int(x) for x in relative_cells],
            }
        )
        for relation in range(m):
            components.append(
                {
                    "wave_id": wave_id,
                    "frequency_multiple": multiple,
                    "relation_id": relation + 1,
                    "edge": edge_label(edges[relation]),
                    "amplitude": float(amplitudes[relation]),
                    "power": float(amplitudes[relation] ** 2),
                    "power_fraction_within_wave": float(amplitudes[relation] ** 2 / power),
                    "absolute_phase_deg": float(math.degrees(phases[relation])),
                    "absolute_phase_cell": int(cells[relation]),
                    "absolute_phase_cell_error_deg": float(math.degrees(cell_errors[relation])),
                    "relative_phase_deg": float(math.degrees(relative[relation])),
                    "relative_phase_cell": int(relative_cells[relation]),
                    "relative_phase_cell_error_deg": float(math.degrees(relative_errors[relation])),
                }
            )

    phase_signatures = {
        (row["frequency_multiple"], tuple(row["absolute_phase_cell_signature"])) for row in waves
    }
    return {
        "schema": "n_only_parent_wave_table_v1",
        "source": {
            "input_directory": str(input_dir),
            "accepted_seed": manifest["accepted_seed"],
            "seed_was_explicit": manifest["seed_was_explicit"],
            "prng": manifest["prng"],
            "numpy_version": manifest["numpy_version"],
            "attempts": manifest["attempts"],
        },
        "summary": {
            "N": n,
            "M": m,
            "lambda0": float(manifest["derived"]["lambda0"]),
            "stable_wave_count": len(waves),
            "phase_distinguished_wave_count": len(phase_signatures),
            "wave_count_upper_bound_M": m,
            "all_waves_zero_closed": all(row["zero_closed"] for row in waves),
            "max_wave_closure_relative": max(row["closure_relative"] for row in waves),
            "max_solver_residual": max(row["solver_residual"] for row in waves),
            "total_closure_abs": float(abs(complex(np.sum(parent * parent)))),
            "frobenius_norm": float(np.linalg.norm(parent)),
        },
        "waves": waves,
        "components": components,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[label for _, label in fields])
        writer.writeheader()
        for row in rows:
            writer.writerow({label: row[key] for key, label in fields})


def write_outputs(output: Path, result: dict[str, Any]) -> None:
    output.mkdir(parents=True, exist_ok=False)
    (output / "census.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    wave_fields = [
        ("wave_id", "波ID"),
        ("frequency_multiple", "振動数倍数n"),
        ("kind", "成分種別"),
        ("phase_increment_deg", "1刻み位相増分度"),
        ("finite_order", "有限回帰次数"),
        ("amplitude", "波全体振幅"),
        ("power", "波強度"),
        ("power_fraction", "全体強度比"),
        ("white_input_weight", "白色入力由来強度重み"),
        ("solver_iterations", "収束反復数"),
        ("solver_residual", "自己無撞着残差"),
        ("sigma_max", "生成子最大特異値"),
        ("closure_abs", "零閉塞絶対残差"),
        ("closure_relative", "零閉塞相対残差"),
        ("zero_closed", "単独零閉塞"),
        ("relation_component_count", "関係成分数"),
        ("nonzero_relation_component_count", "非零関係成分数"),
        ("reference_relation_id", "位相基準関係ID"),
        ("reference_edge", "位相基準辺"),
        ("absolute_phase_cell_signature", "絶対位相セル署名"),
        ("relative_phase_cell_signature", "相対位相セル署名"),
    ]
    component_fields = [
        ("wave_id", "波ID"),
        ("frequency_multiple", "振動数倍数n"),
        ("relation_id", "関係ID"),
        ("edge", "関係辺"),
        ("amplitude", "振幅"),
        ("power", "強度"),
        ("power_fraction_within_wave", "波内強度比"),
        ("absolute_phase_deg", "絶対位相度"),
        ("absolute_phase_cell", "絶対位相セル"),
        ("absolute_phase_cell_error_deg", "絶対位相セル誤差度"),
        ("relative_phase_deg", "基準関係への相対位相度"),
        ("relative_phase_cell", "相対位相セル"),
        ("relative_phase_cell_error_deg", "相対位相セル誤差度"),
    ]
    write_csv(output / "waves_ja.csv", result["waves"], wave_fields)
    write_csv(output / "wave_components_ja.csv", result["components"], component_fields)

    summary = result["summary"]
    source = result["source"]
    lines = [
        f"# N={summary['N']} 白色雑音起源・位相分解親波一覧",
        "",
        "## 1. 条件と生成結果",
        "",
        "| 項目 | 値 |",
        "|---|---:|",
        f"| 唯一の理論入力N | {summary['N']} |",
        f"| 関係数M=N(N−1)/2 | {summary['M']} |",
        f"| 基底位相刻みλ0=2π/N | {summary['lambda0']:.12g} |",
        f"| 採用seed | `{source['accepted_seed']}` |",
        f"| seedを外部指定したか | {'はい' if source['seed_was_explicit'] else 'いいえ（OS生成）'} |",
        f"| 成功までの試行数 | {len(source['attempts'])} |",
        f"| 安定波数W | {summary['stable_wave_count']} |",
        f"| 位相を区別した波数 | {summary['phase_distinguished_wave_count']} |",
        f"| 波数上限M | {summary['wave_count_upper_bound_M']} |",
        f"| 全波が単独零閉塞 | {'はい' if summary['all_waves_zero_closed'] else 'いいえ'} |",
        f"| 最大自己無撞着残差 | {summary['max_solver_residual']:.3e} |",
        f"| 最大零閉塞相対残差 | {summary['max_wave_closure_relative']:.3e} |",
        f"| 全波合成の零閉塞絶対残差 | {summary['total_closure_abs']:.3e} |",
        "",
        "## 2. 波一覧",
        "",
        "一行が一つの周波数・位相波である。関係成分を合計して消さず、全成分は `wave_components_ja.csv` に保存した。",
        "",
        "| 波ID | n | 種別 | 位相増分 | 回帰次数 | 振幅 | 強度比 | 反復 | 自己無撞着残差 | 零閉塞相対残差 | 基準辺 | 相対位相セル署名 |",
        "|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|:---:|:---|",
    ]
    for wave in result["waves"]:
        signature = ",".join(str(x) for x in wave["relative_phase_cell_signature"])
        lines.append(
            f"| {wave['wave_id']} | {wave['frequency_multiple']} | {wave['kind']} | "
            f"{wave['phase_increment_deg']:.6f}° | {wave['finite_order']} | "
            f"{wave['amplitude']:.9f} | {100.0 * wave['power_fraction']:.6f}% | "
            f"{wave['solver_iterations']} | {wave['solver_residual']:.3e} | "
            f"{wave['closure_relative']:.3e} | {wave['reference_edge']} | `{signature}` |"
        )
    lines.extend(
        [
            "",
            "## 3. 読み方",
            "",
            "- 各波は全M関係成分にまたがり、それ自身で零閉塞する。",
            "- 異なるnは別行であり、各行の絶対位相・相対位相も成分ごとに保存する。",
            "- 強度重みは白色入力の列ノルムから引き継ぎ、等振幅を手で置いていない。",
            "- 収束反復数はソルバー監査値であり、物理的な寿命や生成確率ではない。",
            "",
        ]
    )
    (output / "wave_table.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    input_dir = args.input.resolve()
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"output already exists; refusing to overwrite: {output}")
    result = analyse(input_dir)
    write_outputs(output, result)
    print(
        json.dumps(
            {
                "output": str(output),
                "N": result["summary"]["N"],
                "M": result["summary"]["M"],
                "waves": result["summary"]["stable_wave_count"],
                "components": len(result["components"]),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
