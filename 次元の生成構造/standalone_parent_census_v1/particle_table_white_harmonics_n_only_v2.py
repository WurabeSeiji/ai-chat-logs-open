#!/usr/bin/env python3
"""make_parent保存直後の波を、既存の安定波長分類表で読む。

このプログラムは make_parent を import ・変更せず、保存契約
``n_only_white_closed_harmonic_parent_v2`` だけを読む。分類の正本は
``analytic_particle_tables_*/N*/粒子属性族一覧.csv`` である。

分類原則
--------
* 主キーは波長と倍音波長の組であり、位相ではない。
* 位相は実数値を保存し、0/180度へ丸めない。
* N点DFTの循環次数 |k| に対し、|k| が N を割る場合だけ
  波長 q*lambda0, q=N/|k| を安定候補とする。k=0 はN点節上で
  k=Nと区別できない lambda0 表現として別記する。
* 一つの保存波形に分類表外の波長が一本でも含まれれば、
  安定族へ強制投影せず「分類外」とする。
* 分類表内の波長集合でも、最長の基底 q に対して他の d が
  d|q を満たさなければ「単一安定波でない」とする。
* 零閉塞は make_parent の監査量であり、安定波分類の代用にしない。
* M=N(N-1)/2 は入力に保存された関係波数として監査するが、
  Mを安定粒子数とは解釈しない。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np


INPUT_SCHEMA = "n_only_white_closed_harmonic_parent_v2"
OUTPUT_SCHEMA = "stable_wavelength_classification_after_parent_v1"
HERE = Path(__file__).resolve().parent


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(HERE))
    except ValueError:
        return resolved.name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="make_parent直後の波を安定波長分類表で分類する"
    )
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def signed_order(bin_index: int, n: int) -> int:
    return bin_index if bin_index <= n // 2 else bin_index - n


def edge_label(edge: np.ndarray) -> str:
    return f"({int(edge[0]) + 1},{int(edge[1]) + 1})"


def parse_base_q(text: str) -> int:
    suffix = "lambda0"
    if not text.endswith(suffix):
        raise ValueError(f"基底波長表記が不正です: {text!r}")
    return int(text[: -len(suffix)])


def find_catalog(n: int) -> Path:
    candidates = sorted(
        HERE.glob(f"analytic_particle_tables_*/N{n}/粒子属性族一覧.csv")
    )
    if len(candidates) != 1:
        raise SystemExit(
            f"N={n} の安定波分類表を1個に特定できません: "
            f"{[str(path) for path in candidates]}"
        )
    return candidates[0]


def load_catalog(n: int) -> tuple[Path, dict[tuple[int, int, int], dict[str, str]]]:
    path = find_catalog(n)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    lookup: dict[tuple[int, int, int], dict[str, str]] = {}
    for row in rows:
        if int(row["分解能N"]) != n:
            raise SystemExit(f"分類表のNが一致しません: {path}")
        key = (
            parse_base_q(row["基底波長"]),
            int(row["選択する奇数倍音数"]),
            int(row["選択する偶数倍音数"]),
        )
        if key in lookup:
            raise SystemExit(f"分類表のキーが重複しています: {key}")
        lookup[key] = row
    return path, lookup


def order_wavelength(order: int, n: int) -> tuple[str, int | None, bool, bool]:
    """DFT次数を波長倍率 q へ変換する。

    戻り値は (表示, q, 安定候補, lambda0エイリアス) 。
    """

    absolute = abs(order)
    if absolute == 0:
        return "1lambda0（N点節上のk=Nエイリアス）", 1, True, True
    quotient, remainder = divmod(n, absolute)
    if remainder == 0:
        return f"{quotient}lambda0", quotient, True, False
    divisor = math.gcd(n, absolute)
    return f"{n // divisor}/{absolute // divisor}lambda0", None, False, False


def family_attributes(family: dict[str, str] | None) -> dict[str, Any]:
    if family is None:
        return {
            "family_id": None,
            "reflection_equal_amplitude": None,
            "closure_address": None,
            "conjugate_address": None,
            "charge_from_address": None,
            "particle_antiparticle": None,
            "spin_from_winding_cover": None,
        }
    charge_text = family["住所電荷量sin2(pi*m/n)"]
    return {
        "family_id": family["族ID"],
        "reflection_equal_amplitude": family["反射率R（等振幅）"],
        "closure_address": family["閉鎖住所"],
        "conjugate_address": family["共役住所"],
        "charge_from_address": None if charge_text == "—" else float(charge_text),
        "particle_antiparticle": family["粒子・反粒子"],
        "spin_from_winding_cover": family["巻数・被覆から読むスピン"],
    }


def classify_spectrum(
    spectrum: np.ndarray,
    n: int,
    m: int,
    catalog: dict[tuple[int, int, int], dict[str, str]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    power = np.abs(spectrum) ** 2
    total_power = float(np.sum(power))
    eps = float(np.finfo(spectrum.real.dtype).eps)
    numerical_floor = (
        max(n, m) * eps * max(float(np.max(power)), np.finfo(float).tiny)
    )
    present = power > numerical_floor

    components: list[dict[str, Any]] = []
    wavelength_power: defaultdict[int, float] = defaultdict(float)
    allowed_wavelengths: set[int] = set()
    outside_orders: list[int] = []
    outside_power = 0.0

    for bin_index in range(n):
        order = signed_order(bin_index, n)
        wavelength, q, allowed, alias = order_wavelength(order, n)
        is_present = bool(present[bin_index])
        if is_present and allowed:
            assert q is not None
            allowed_wavelengths.add(q)
            wavelength_power[q] += float(power[bin_index])
        elif is_present:
            outside_orders.append(order)
            outside_power += float(power[bin_index])
        components.append(
            {
                "DFT_bin": bin_index,
                "cyclic_order": order,
                "wavelength": wavelength,
                "wavelength_multiple": q,
                "stationary_wavelength_allowed": allowed,
                "lambda0_node_alias": alias,
                "present_above_numerical_floor": is_present,
                "amplitude": float(abs(spectrum[bin_index])),
                "power": float(power[bin_index]),
                "power_fraction": float(power[bin_index] / total_power),
                "phase_deg": float(math.degrees(np.angle(spectrum[bin_index]))),
            }
        )

    base_q: int | None = None
    harmonic_orders: list[int] = []
    odd_count = 0
    even_count = 0
    family: dict[str, str] | None = None
    classification_status = "分類外"
    if outside_orders:
        classification_reason = "安定波分類表外の波長を含む"
    elif not allowed_wavelengths:
        classification_reason = "数値誤差床より上の波長成分がない"
    else:
        base_q = max(allowed_wavelengths)
        incompatible = sorted(q for q in allowed_wavelengths if base_q % q != 0)
        if incompatible:
            classification_reason = "単一基底波の倍音束ではない"
        else:
            harmonic_orders = sorted(
                base_q // q for q in allowed_wavelengths if q != base_q
            )
            odd_count = sum(order % 2 == 1 for order in harmonic_orders)
            even_count = len(harmonic_orders) - odd_count
            family = catalog.get((base_q, odd_count, even_count))
            if family is None:
                classification_reason = "波長構成に対応する安定族が分類表にない"
            else:
                classification_status = "安定波分類表に一致"
                classification_reason = "基底波長と倍音波長の組が一致"

    actual_reflection: float | None = None
    if family is not None and base_q is not None:
        odd_power = 0.0
        even_power = 0.0
        for q, component_power in wavelength_power.items():
            if q == base_q:
                continue
            harmonic_order = base_q // q
            if harmonic_order % 2:
                odd_power += component_power
            else:
                even_power += component_power
        harmonic_power = odd_power + even_power
        actual_reflection = 0.0 if harmonic_power == 0.0 else 0.5 * odd_power / harmonic_power

    result = {
        "classification_status": classification_status,
        "classification_reason": classification_reason,
        "base_wavelength_multiple": base_q,
        "base_wavelength": None if base_q is None else f"{base_q}lambda0",
        "present_allowed_wavelengths": sorted(allowed_wavelengths),
        "present_harmonic_orders": harmonic_orders,
        "selected_odd_harmonic_count": odd_count,
        "selected_even_harmonic_count": even_count,
        "outside_stationary_order_count": len(outside_orders),
        "outside_stationary_orders": outside_orders,
        "outside_stationary_power_fraction": outside_power / total_power,
        "actual_reflection_from_power": actual_reflection,
        "phase_is_classification_key": False,
        "phase_rounding_to_0_or_180": False,
        "numerical_floor": numerical_floor,
        **family_attributes(family),
    }
    return result, components


def analyse(input_dir: Path) -> dict[str, Any]:
    manifest_path = input_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != INPUT_SCHEMA:
        raise SystemExit(f"未対応schema: {manifest.get('schema')!r}")
    if manifest.get("status") != "success":
        raise SystemExit(f"成功した親ではありません: {input_dir}")

    n = int(manifest["function_contract"]["N"])
    m = int(manifest["derived_only_from_N"]["M"])
    if m != n * (n - 1) // 2:
        raise SystemExit("M != N(N-1)/2")

    waves = np.load(
        input_dir / manifest["arrays"]["relation_waves"]["file"], mmap_mode="r"
    )
    parent = np.load(
        input_dir / manifest["arrays"]["parent_vector"]["file"], mmap_mode="r"
    )
    edges = np.load(input_dir / manifest["arrays"]["edges"]["file"], mmap_mode="r")
    if waves.shape != (m, n) or parent.shape != (m,) or edges.shape != (m, 2):
        raise SystemExit(
            f"配列shape不正: waves={waves.shape}, parent={parent.shape}, edges={edges.shape}"
        )

    catalog_path, catalog = load_catalog(n)
    rows: list[dict[str, Any]] = []
    component_aggregate: dict[int, dict[str, Any]] = {}
    closure_abs_values: list[float] = []
    closure_relative_values: list[float] = []
    nested_closure = 0.0 + 0.0j
    reconstruction_max = 0.0

    for relation in range(m):
        samples = np.asarray(waves[relation])
        spectrum = np.fft.fft(samples) / math.sqrt(n)
        reconstructed = np.fft.ifft(spectrum) * math.sqrt(n)
        reconstruction_max = max(
            reconstruction_max, float(np.max(np.abs(reconstructed - samples)))
        )
        power_sum = float(np.sum(np.abs(samples) ** 2))
        closure_value = complex(np.sum(samples**2))
        nested_closure += closure_value
        closure_abs = float(abs(closure_value))
        closure_relative = closure_abs / max(power_sum, np.finfo(float).tiny)
        closure_abs_values.append(closure_abs)
        closure_relative_values.append(closure_relative)

        classification, components = classify_spectrum(spectrum, n, m, catalog)
        outside_orders = classification.pop("outside_stationary_orders")
        outside_order_bytes = json.dumps(
            outside_orders, separators=(",", ":")
        ).encode("ascii")
        if len(outside_orders) <= 32:
            outside_orders_preview: list[int | str] = outside_orders
        else:
            outside_orders_preview = [
                *outside_orders[:16], "...", *outside_orders[-16:]
            ]
        classified_phases = [
            {
                "cyclic_order": item["cyclic_order"],
                "wavelength": item["wavelength"],
                "amplitude": item["amplitude"],
                "phase_deg": item["phase_deg"],
            }
            for item in components
            if item["present_above_numerical_floor"]
            and item["stationary_wavelength_allowed"]
            and (n <= 40 or classification["family_id"] is not None)
        ]
        rows.append(
            {
                "wave_id": f"W{relation + 1:04d}",
                "relation_id": relation + 1,
                "edge": edge_label(edges[relation]),
                "wave_amplitude": math.sqrt(power_sum),
                "parent_phase_deg": float(math.degrees(np.angle(parent[relation]))),
                "classified_component_phases": classified_phases,
                "outside_stationary_orders_preview": outside_orders_preview,
                "outside_stationary_orders_sha256": hashlib.sha256(
                    outside_order_bytes
                ).hexdigest(),
                "zero_closure_abs": closure_abs,
                "zero_closure_relative": closure_relative,
                "zero_closed": closure_relative < 1e-12,
                **classification,
            }
        )

        for item in components:
            order = int(item["cyclic_order"])
            aggregate = component_aggregate.setdefault(
                order,
                {
                    "cyclic_order": order,
                    "wavelength": item["wavelength"],
                    "wavelength_multiple": item["wavelength_multiple"],
                    "stationary_wavelength_allowed": item[
                        "stationary_wavelength_allowed"
                    ],
                    "lambda0_node_alias": item["lambda0_node_alias"],
                    "present_wave_count": 0,
                    "total_power": 0.0,
                    "complex_amplitude_sum": 0.0 + 0.0j,
                    "absolute_amplitude_sum": 0.0,
                },
            )
            if item["present_above_numerical_floor"]:
                aggregate["present_wave_count"] += 1
            aggregate["total_power"] += item["power"]
            value = spectrum[item["DFT_bin"]]
            aggregate["complex_amplitude_sum"] += value
            aggregate["absolute_amplitude_sum"] += abs(value)

    spectrum_rows: list[dict[str, Any]] = []
    total_spectral_power = sum(item["total_power"] for item in component_aggregate.values())
    for order in sorted(component_aggregate):
        item = component_aggregate[order]
        amplitude_sum = float(item.pop("absolute_amplitude_sum"))
        complex_sum = complex(item.pop("complex_amplitude_sum"))
        spectrum_rows.append(
            {
                **item,
                "power_fraction": item["total_power"] / total_spectral_power,
                "circular_mean_phase_deg": float(math.degrees(np.angle(complex_sum))),
                "phase_coherence": abs(complex_sum) / max(amplitude_sum, np.finfo(float).tiny),
            }
        )

    stable_rows = [row for row in rows if row["family_id"] is not None]
    reason_counts = Counter(row["classification_reason"] for row in rows)
    family_counts = Counter(row["family_id"] for row in stable_rows)
    family_power = defaultdict(float)
    for row in stable_rows:
        family_power[row["family_id"]] += row["wave_amplitude"] ** 2
    stable_families = []
    for family_id in sorted(family_counts):
        example = next(row for row in stable_rows if row["family_id"] == family_id)
        stable_families.append(
            {
                "family_id": family_id,
                "base_wavelength": example["base_wavelength"],
                "odd_harmonic_count": example["selected_odd_harmonic_count"],
                "even_harmonic_count": example["selected_even_harmonic_count"],
                "occurrence_count": family_counts[family_id],
                "total_power": family_power[family_id],
                "reflection_equal_amplitude": example["reflection_equal_amplitude"],
                "closure_address": example["closure_address"],
                "charge_from_address": example["charge_from_address"],
                "spin_from_winding_cover": example["spin_from_winding_cover"],
            }
        )

    allowed_power = sum(
        item["total_power"]
        for item in spectrum_rows
        if item["stationary_wavelength_allowed"]
    )
    outside_power = total_spectral_power - allowed_power
    return {
        "schema": OUTPUT_SCHEMA,
        "source": {
            "input_directory": portable_path(input_dir),
            "manifest_sha256": sha256_file(manifest_path),
            "accepted_seed": manifest["accepted_seed"],
            "generator": manifest["generator"],
            "generator_sha256": manifest["generator_sha256"],
            "make_parent_modified_by_reader": False,
            "universal_interaction_used_or_modified": False,
        },
        "classification_catalog": {
            "path": portable_path(catalog_path),
            "sha256": sha256_file(catalog_path),
            "family_count": len(catalog),
            "phase_used_as_key": False,
        },
        "summary": {
            "N": n,
            "M": m,
            "lambda0": float(manifest["derived_only_from_N"]["lambda0"]),
            "generated_relation_wave_count": m,
            "stable_classified_wave_count": len(stable_rows),
            "unclassified_wave_count": m - len(stable_rows),
            "classification_reason_counts": dict(reason_counts),
            "observed_stable_family_count": len(stable_families),
            "allowed_stationary_wavelength_power": allowed_power,
            "outside_stationary_wavelength_power": outside_power,
            "allowed_stationary_wavelength_power_fraction": (
                allowed_power / total_spectral_power
            ),
            "outside_stationary_wavelength_power_fraction": (
                outside_power / total_spectral_power
            ),
            "all_relation_waves_zero_closed": all(row["zero_closed"] for row in rows),
            "max_relation_wave_closure_abs": max(closure_abs_values),
            "max_relation_wave_closure_relative": max(closure_relative_values),
            "nested_total_closure_abs": float(abs(nested_closure)),
            "parent_vector_closure_abs": float(abs(complex(np.asarray(parent) @ np.asarray(parent)))),
            "DFT_reconstruction_max_error": reconstruction_max,
            "classification_is_readout_only": True,
            "unclassified_is_preserved": True,
        },
        "stable_families_observed": stable_families,
        "waves": rows,
        "spectrum_by_cyclic_order": spectrum_rows,
    }


def serialise_cell(value: Any) -> Any:
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "はい" if value else "いいえ"
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return value


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[label for _, label in fields])
        writer.writeheader()
        for row in rows:
            writer.writerow({label: serialise_cell(row[key]) for key, label in fields})


WAVE_FIELDS = [
    ("wave_id", "生成波ID"),
    ("relation_id", "関係ID"),
    ("edge", "関係辺"),
    ("classification_status", "安定波分類"),
    ("classification_reason", "分類理由"),
    ("family_id", "安定族ID"),
    ("base_wavelength", "基底波長"),
    ("present_allowed_wavelengths", "存在する許容波長倍率"),
    ("present_harmonic_orders", "倍音次数"),
    ("outside_stationary_order_count", "分類表外波長数"),
    ("outside_stationary_orders_preview", "分類表外循環次数_最大32件"),
    ("outside_stationary_orders_sha256", "分類表外循環次数_SHA256"),
    ("outside_stationary_power_fraction", "分類表外強度比"),
    ("selected_odd_harmonic_count", "奇数倍音数"),
    ("selected_even_harmonic_count", "偶数倍音数"),
    ("reflection_equal_amplitude", "反射率R_等振幅参考"),
    ("actual_reflection_from_power", "反射率R_実強度"),
    ("closure_address", "閉鎖住所"),
    ("conjugate_address", "共役住所"),
    ("charge_from_address", "住所電荷量"),
    ("particle_antiparticle", "粒子反粒子読出し"),
    ("spin_from_winding_cover", "巻数被覆スピン読出し"),
    ("wave_amplitude", "生成波振幅"),
    ("parent_phase_deg", "親位相度_連続値"),
    ("classified_component_phases", "許容波長成分の振幅位相_連続値"),
    ("zero_closure_relative", "零閉塞相対残差_監査のみ"),
]

SPECTRUM_FIELDS = [
    ("cyclic_order", "循環次数k"),
    ("wavelength", "波長"),
    ("stationary_wavelength_allowed", "安定波長候補"),
    ("lambda0_node_alias", "lambda0節上エイリアス"),
    ("present_wave_count", "存在する関係波数"),
    ("total_power", "合計強度"),
    ("power_fraction", "全強度比"),
    ("circular_mean_phase_deg", "円平均位相度"),
    ("phase_coherence", "位相コヒーレンス"),
]


def write_run(output: Path, result: dict[str, Any]) -> None:
    output.mkdir(parents=True, exist_ok=False)
    census = {key: value for key, value in result.items() if key != "waves"}
    census["wave_table"] = {
        "file": "make_parent直後の波分類.csv",
        "row_count": len(result["waves"]),
        "all_rows_written": True,
    }
    (output / "census.json").write_text(
        json.dumps(census, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(output / "make_parent直後の波分類.csv", result["waves"], WAVE_FIELDS)
    write_csv(
        output / "循環次数別波長集計.csv",
        result["spectrum_by_cyclic_order"],
        SPECTRUM_FIELDS,
    )

    summary = result["summary"]
    source = result["source"]
    catalog = result["classification_catalog"]
    lines = [
        f"# N={summary['N']} make_parent直後の安定波分類",
        "",
        "make_parentと万能相互作用演算は変更・実行していない。",
        "保存されたmake_parent直後の波を、既存の安定波分類表で読んだ結果である。",
        "位相は分類キーにせず、0/180度に丸めていない。",
        "",
        "## 直後分類の集計",
        "",
        "| 項目 | 値 |",
        "|---|---:|",
        f"| N | {summary['N']} |",
        f"| M=N(N−1)/2（生成関係波数。粒子数ではない） | {summary['M']} |",
        f"| 安定波分類表に一致 | {summary['stable_classified_wave_count']} |",
        f"| 分類外 | {summary['unclassified_wave_count']} |",
        f"| 観測された安定族数 | {summary['observed_stable_family_count']} |",
        f"| 安定波長候補の強度比 | {summary['allowed_stationary_wavelength_power_fraction']:.12g} |",
        f"| 分類表外波長の強度比 | {summary['outside_stationary_wavelength_power_fraction']:.12g} |",
        f"| 零閉塞成立（分類とは別の監査） | {'はい' if summary['all_relation_waves_zero_closed'] else 'いいえ'} |",
        f"| 最大零閉塞相対残差 | {summary['max_relation_wave_closure_relative']:.3e} |",
        f"| 採用seed | `{source['accepted_seed']}` |",
        f"| 使用した安定波分類族数 | {catalog['family_count']} |",
        "",
        "## 分類外の理由",
        "",
        "| 理由 | 波数 |",
        "|---|---:|",
    ]
    for reason, count in sorted(summary["classification_reason_counts"].items()):
        lines.append(f"| {reason} | {count} |")

    lines.extend(
        [
            "",
            "## 実際に観測された安定族",
            "",
            "| 族ID | 基底波長 | 奇数倍音数 | 偶数倍音数 | 出現数 | 合計強度 | 反射率R（等振幅参考） | 閉鎖住所 | 電荷量 | スピン読出し |",
            "|---|---:|---:|---:|---:|---:|---:|---|---:|---|",
        ]
    )
    if not result["stable_families_observed"]:
        lines.append("| — | — | — | — | 0 | 0 | — | — | — | 現時点で安定族なし |")
    else:
        for row in result["stable_families_observed"]:
            charge = "—" if row["charge_from_address"] is None else f"{row['charge_from_address']:.12g}"
            lines.append(
                f"| {row['family_id']} | {row['base_wavelength']} | "
                f"{row['odd_harmonic_count']} | {row['even_harmonic_count']} | "
                f"{row['occurrence_count']} | {row['total_power']:.12g} | "
                f"{row['reflection_equal_amplitude']} | {row['closure_address']} | "
                f"{charge} | {row['spin_from_winding_cover']} |"
            )
    lines.extend(
        [
            "",
            "全生成波の判定、分類表外循環次数、連続位相（N<=40）は",
            "`make_parent直後の波分類.csv` に保存した。",
            "循環次数ごとの波長、許容判定、強度は `循環次数別波長集計.csv` にある。",
            "",
            "## 分類の境界",
            "",
            "- 分類表外の波長を含む波を、最近の安定族へ投影していない。",
            "- 零閉塞していても、波長構成が分類表外なら「分類外」である。",
            "- k=0成分はN点節上で最短波長lambda0と区別できないため、エイリアスと明記した。",
            "",
        ]
    )
    (output / "make_parent直後の分類表.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def combined_markdown(results: list[dict[str, Any]]) -> str:
    lines = [
        "# make_parent直後の安定波分類まとめ",
        "",
        "make_parent本体と万能相互作用演算は変更していない。",
        "既存の安定波分類表を、make_parent保存直後の波に対して使用した。",
        "",
        "| N | M（粒子数ではない） | 安定分類に一致 | 分類外 | 観測安定族数 | 安定波長候補強度比 | 分類表外強度比 | 零閉塞 | seed |",
        "|---:|---:|---:|---:|---:|---:|---:|---|---:|",
    ]
    for result in results:
        summary = result["summary"]
        lines.append(
            f"| {summary['N']} | {summary['M']} | "
            f"{summary['stable_classified_wave_count']} | {summary['unclassified_wave_count']} | "
            f"{summary['observed_stable_family_count']} | "
            f"{summary['allowed_stationary_wavelength_power_fraction']:.12g} | "
            f"{summary['outside_stationary_wavelength_power_fraction']:.12g} | "
            f"{'はい' if summary['all_relation_waves_zero_closed'] else 'いいえ'} | "
            f"{result['source']['accepted_seed']} |"
        )
    lines.extend(
        [
            "",
            "「分類外」は失敗や削除対象ではなく、安定波分類表にまだ入らない",
            "make_parent直後の波長構成をそのまま数えたものである。",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    if output.exists():
        raise SystemExit(f"output already exists; refusing to overwrite: {output}")
    output.mkdir(parents=True, exist_ok=False)

    results: list[dict[str, Any]] = []
    seen: set[int] = set()
    for input_path in args.input:
        result = analyse(input_path.resolve())
        n = result["summary"]["N"]
        if n in seen:
            raise SystemExit(f"同じNが重複しています: {n}")
        seen.add(n)
        results.append(result)
        write_run(output / f"N{n}", result)
        print(
            f"N={n} M={result['summary']['M']} "
            f"stable={result['summary']['stable_classified_wave_count']} "
            f"unclassified={result['summary']['unclassified_wave_count']}"
        )

    results.sort(key=lambda item: item["summary"]["N"])
    (output / "summary.md").write_text(combined_markdown(results), encoding="utf-8")
    (output / "summary.json").write_text(
        json.dumps(
            {
                "schema": OUTPUT_SCHEMA,
                "reader": Path(__file__).name,
                "reader_sha256": sha256_file(Path(__file__).resolve()),
                "results": [
                    {
                        "source": result["source"],
                        "classification_catalog": result["classification_catalog"],
                        "summary": result["summary"],
                        "stable_families_observed": result[
                            "stable_families_observed"
                        ],
                        "wave_table": {
                            "file": f"N{result['summary']['N']}/make_parent直後の波分類.csv",
                            "row_count": len(result["waves"]),
                            "all_rows_written": True,
                        },
                    }
                    for result in results
                ],
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
