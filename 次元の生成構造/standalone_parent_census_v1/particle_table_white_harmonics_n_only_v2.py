#!/usr/bin/env python3
"""白色零閉塞親から、M本の関係波と全倍音を独立に一覧化する。

生成器をimportしない。``n_only_white_closed_harmonic_parent_v2`` の保存契約
だけを読む。主表の一行は一つの関係波であり、N=5なら10行、N=40なら
780行になる。倍音は保存波形のN点DFTからだけ読み、全M*N成分の次数・
振幅・位相を別CSVへ保存する。
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


INPUT_SCHEMA = "n_only_white_closed_harmonic_parent_v2"
OUTPUT_SCHEMA = "n_only_white_closed_particle_table_v2"
CONSTANT_WATCH = [0.3, 0.7, 0.302822, 0.697178]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Nだけの白色零閉塞親から粒子表を作る")
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def signed_order(bin_index: int, n: int) -> int:
    return bin_index if bin_index <= n // 2 else bin_index - n


def edge_label(edge: np.ndarray) -> str:
    return f"({int(edge[0]) + 1},{int(edge[1]) + 1})"


def address_for(k: int, n: int) -> tuple[int, int, str]:
    if k == 0:
        return 0, 1, "0/1"
    divisor = math.gcd(abs(k), n)
    numerator = k // divisor
    denominator = n // divisor
    return numerator, denominator, f"{numerator}/{denominator}"


def harmonic_kind(k: int) -> str:
    if k == 0:
        return "直流"
    if k == 1:
        return "基本波"
    if k == -1:
        return "逆回転基本波"
    if k > 1:
        return f"第{k}倍音"
    return f"逆回転第{abs(k)}倍音"


def parity_type(odd_power: float, even_power: float, floor: float) -> str:
    if odd_power <= floor and even_power <= floor:
        return "非直流倍音なし"
    if even_power <= floor:
        return "純奇数倍音（F型）"
    if odd_power <= floor:
        return "純偶数倍音（B型）"
    if odd_power >= even_power:
        return "奇数・偶数倍音混合（奇数優勢）"
    return "奇数・偶数倍音混合（偶数優勢）"


def phase_network_readout(
    values: np.ndarray, tolerance: float
) -> tuple[str, str, float]:
    """DFT成分の相対位相から B/F/E 位相網を読む。

    位相セル数や外部分解能は置かない。数値許容差は浮動小数誤差の
    判定にだけ用い、B=全成分が同相、F=二成分が逆相、E=それ以外とする。
    """

    if len(values) <= 1:
        return "—", "単一成分（位相網なし）", 0.0
    units = values / np.abs(values)
    relative = units * np.conjugate(units[0])
    same_phase_residual = float(np.max(np.abs(relative - 1.0)))
    if same_phase_residual <= tolerance:
        return "B", "同相網（ボゾン型極限）", same_phase_residual
    if len(values) == 2:
        opposite_phase_residual = float(abs(relative[1] + 1.0))
        if opposite_phase_residual <= tolerance:
            return "F", "逆相二成分網（フェルミオン型極限）", opposite_phase_residual
    pairwise = units[:, None] * np.conjugate(units[None, :])
    off_diagonal = ~np.eye(len(values), dtype=bool)
    intermediate_residual = float(
        np.min(
            np.minimum(np.abs(pairwise - 1.0), np.abs(pairwise + 1.0))[
                off_diagonal
            ]
        )
    )
    return "E", "中間位相網（エルミオン型）", intermediate_residual


def cyclic_return_orders(
    samples: np.ndarray, tolerance: float = 1e-10
) -> tuple[int, int]:
    """状態 w と二次量 w^2 の最小循環回帰次数を別々に読む。"""

    state_norm = float(np.linalg.norm(samples))
    quadratic = samples * samples
    quadratic_norm = float(np.linalg.norm(quadratic))
    state_order = len(samples)
    quadratic_order = len(samples)
    for lag in range(1, len(samples) + 1):
        state_residual = float(
            np.linalg.norm(np.roll(samples, -lag) - samples)
            / max(state_norm, np.finfo(float).tiny)
        )
        if state_residual <= tolerance:
            state_order = lag
            break
    for lag in range(1, len(samples) + 1):
        quadratic_residual = float(
            np.linalg.norm(np.roll(quadratic, -lag) - quadratic)
            / max(quadratic_norm, np.finfo(float).tiny)
        )
        if quadratic_residual <= tolerance:
            quadratic_order = lag
            break
    return state_order, quadratic_order


def half_turn_cover(
    samples: np.ndarray, tolerance: float = 1e-10
) -> tuple[complex | None, float | None, float | None, float | None, float | None]:
    """半周期で状態が +w / -w のどちらへ戻るかと二次量の回帰を測る。"""

    n = len(samples)
    if n % 2:
        return None, None, None, None, None
    shifted = np.roll(samples, -n // 2)
    state_norm = max(float(np.linalg.norm(samples)), np.finfo(float).tiny)
    quadratic = samples * samples
    quadratic_norm = max(float(np.linalg.norm(quadratic)), np.finfo(float).tiny)
    overlap = complex(np.vdot(samples, shifted) / (state_norm * state_norm))
    plus_residual = float(np.linalg.norm(shifted - samples) / state_norm)
    minus_residual = float(np.linalg.norm(shifted + samples) / state_norm)
    quadratic_residual = float(
        np.linalg.norm(shifted * shifted - quadratic) / quadratic_norm
    )
    if plus_residual <= tolerance and quadratic_residual <= tolerance:
        cover_ratio: float | None = 1.0
    elif minus_residual <= tolerance and quadratic_residual <= tolerance:
        cover_ratio = 2.0
    else:
        cover_ratio = None
    return overlap, plus_residual, minus_residual, quadratic_residual, cover_ratio


def spin_type(
    parity: str,
    half_turn_ratio: float | None,
) -> str:
    """倍音パリティと半周期の符号/二次量回帰からスピン型候補を読む。"""

    if parity.startswith("純奇数"):
        if half_turn_ratio == 2.0:
            return "2:1被覆確認（半整数スピン型）"
        return "奇数倍音・2:1被覆候補（半整数スピン型）"
    if parity.startswith("純偶数"):
        if half_turn_ratio == 1.0:
            return "1:1被覆確認（整数スピン型）"
        return "偶数倍音・整数スピン型候補"
    if parity.startswith("奇数・偶数"):
        return "奇数・偶数倍音の混合スピン型"
    return "非直流倍音なし"


def cyclic_correlation_lifetime(samples: np.ndarray) -> tuple[int | None, list[float]]:
    power = float(np.vdot(samples, samples).real)
    if power == 0.0:
        return None, []
    correlations = []
    first_crossing: int | None = None
    for lag in range(1, len(samples)):
        value = float(abs(np.vdot(samples, np.roll(samples, -lag))) / power)
        correlations.append(value)
        if first_crossing is None and value < math.exp(-1.0):
            first_crossing = lag
    return first_crossing, correlations


def branch_gram(spectrum: np.ndarray, n: int) -> tuple[float, float, float, float, complex]:
    max_pair = (n - 1) // 2
    if max_pair == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0 + 0.0j
    positive = spectrum[np.arange(1, max_pair + 1)]
    negative = spectrum[(-np.arange(1, max_pair + 1)) % n]
    norm_positive = float(np.vdot(positive, positive).real)
    norm_negative = float(np.vdot(negative, negative).real)
    overlap = complex(np.vdot(positive, negative))
    determinant = float(max(0.0, norm_positive * norm_negative - abs(overlap) ** 2))
    return determinant, math.sqrt(determinant), norm_positive, norm_negative, overlap


def phase_distinguished_count(waves: np.ndarray) -> int:
    # 位相を商で消さない。同じ複素波形だけを同一とする。
    return len({np.ascontiguousarray(row).tobytes() for row in waves})


def analyse(input_dir: Path) -> dict[str, Any]:
    manifest = json.loads((input_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != INPUT_SCHEMA:
        raise SystemExit(f"未対応schema: {manifest.get('schema')!r}")
    if manifest.get("status") != "success":
        raise SystemExit(f"成功した親ではありません: {input_dir}")

    n = int(manifest["function_contract"]["N"])
    m = int(manifest["derived_only_from_N"]["M"])
    if m != n * (n - 1) // 2:
        raise SystemExit("M != N(N-1)/2")
    waves = np.load(input_dir / manifest["arrays"]["relation_waves"]["file"])
    parent = np.load(input_dir / manifest["arrays"]["parent_vector"]["file"])
    edges = np.load(input_dir / manifest["arrays"]["edges"]["file"])
    if waves.shape != (m, n) or parent.shape != (m,) or edges.shape != (m, 2):
        raise SystemExit(
            f"配列shape不正: waves={waves.shape}, parent={parent.shape}, edges={edges.shape}"
        )

    spectra = np.fft.fft(waves, axis=1) / math.sqrt(n)
    reconstructed = np.fft.ifft(spectra, axis=1) * math.sqrt(n)
    reconstruction_error = float(np.max(np.abs(reconstructed - waves)))
    row_power = np.sum(np.abs(waves) ** 2, axis=1)
    row_closure = np.abs(np.sum(waves**2, axis=1))
    row_closure_relative = row_closure / np.maximum(row_power, np.finfo(float).tiny)
    orders = np.array([signed_order(k, n) for k in range(n)], dtype=int)
    eps = float(np.finfo(waves.real.dtype).eps)

    rows: list[dict[str, Any]] = []
    harmonics: list[dict[str, Any]] = []
    all_charge_records: list[dict[str, Any]] = []

    for relation in range(m):
        wave_id = f"W{relation + 1:04d}"
        samples = waves[relation]
        spectrum = spectra[relation]
        power = np.abs(spectrum) ** 2
        total_power = float(np.sum(power))
        numerical_floor = max(n, m) * eps * max(float(np.max(power)), np.finfo(float).tiny)
        present = power > numerical_floor
        non_dc = orders != 0
        present_non_dc = present & non_dc
        base_mask = non_dc & (np.abs(orders) == 1)
        higher_harmonic_mask = non_dc & (np.abs(orders) >= 2)
        present_base = present & base_mask
        present_harmonics = present & higher_harmonic_mask
        if np.any(present_non_dc):
            candidate_bins = np.flatnonzero(present_non_dc)
            dominant_bin = int(candidate_bins[np.argmax(power[candidate_bins])])
        else:
            dominant_bin = 0
        dominant_k = int(orders[dominant_bin])
        addr_m, addr_n, address = address_for(dominant_k, n)
        charge_magnitude = float(math.sin(math.pi * addr_m / addr_n) ** 2)

        non_dc_power = float(np.sum(power[non_dc]))
        odd_mask = non_dc & (np.abs(orders) % 2 == 1)
        even_mask = non_dc & (np.abs(orders) % 2 == 0)
        odd_power = float(np.sum(power[odd_mask]))
        even_power = float(np.sum(power[even_mask]))
        parity_denom = odd_power + even_power
        odd_fraction = odd_power / parity_denom if parity_denom else 0.0
        even_fraction = even_power / parity_denom if parity_denom else 0.0
        parity = parity_type(odd_power, even_power, numerical_floor * n)

        phase_tolerance = 100.0 * max(n, m) * eps
        bfe_label, bfe_type, bfe_residual = phase_network_readout(
            spectrum[present_non_dc], phase_tolerance
        )

        effective_harmonics = 0.0
        if non_dc_power:
            fractions = power[non_dc] / non_dc_power
            effective_harmonics = float(1.0 / np.sum(fractions * fractions))

        gram_det, gram_sqrt, positive_power, negative_power, overlap = branch_gram(
            spectrum, n
        )
        branch_denominator = positive_power + negative_power
        rotation_polarization = (
            (positive_power - negative_power) / branch_denominator
            if branch_denominator
            else 0.0
        )
        if rotation_polarization > 100.0 * eps:
            particle_side = "粒子側（正周回）"
        elif rotation_polarization < -100.0 * eps:
            particle_side = "反粒子側（負周回）"
        else:
            particle_side = "自己共役側（正負周回同強度）"
        charge_signed = (
            float(math.copysign(charge_magnitude, rotation_polarization))
            if abs(rotation_polarization) > 100.0 * eps
            else 0.0
        )

        (
            half_turn,
            half_turn_plus_residual,
            half_turn_minus_residual,
            half_turn_quadratic_residual,
            cover_ratio,
        ) = half_turn_cover(samples)
        lifetime, correlations = cyclic_correlation_lifetime(samples)
        state_order, quadratic_order = cyclic_return_orders(samples)
        present_non_dc_orders = [int(k) for k in orders[present_non_dc]]
        present_harmonic_orders = [int(k) for k in orders[present_harmonics]]
        base_count = int(np.count_nonzero(present_base))
        harmonic_count = int(np.count_nonzero(present_harmonics))
        if base_count and harmonic_count:
            wave_composition = "基本波＋倍音"
        elif base_count:
            wave_composition = "基本波のみ"
        elif harmonic_count:
            wave_composition = "倍音のみ（基本波なし）"
        else:
            wave_composition = "直流のみ"
        fundamental_bin = 1 if n > 1 else 0

        rows.append(
            {
                "wave_id": wave_id,
                "relation_id": relation + 1,
                "edge": edge_label(edges[relation]),
                "wave_amplitude": math.sqrt(total_power),
                "parent_phase_deg": float(math.degrees(np.angle(parent[relation]))),
                "sample_count_N": n,
                "wave_composition": wave_composition,
                "detected_non_dc_component_count": int(
                    np.count_nonzero(present_non_dc)
                ),
                "detected_non_dc_orders": present_non_dc_orders,
                "base_wave_component_count": base_count,
                "detected_harmonic_count": harmonic_count,
                "detected_harmonic_orders": present_harmonic_orders,
                "fundamental_amplitude": float(abs(spectrum[fundamental_bin])),
                "fundamental_phase_deg": float(
                    math.degrees(np.angle(spectrum[fundamental_bin]))
                ),
                "dominant_harmonic_order": dominant_k,
                "dominant_harmonic_amplitude": float(abs(spectrum[dominant_bin])),
                "dominant_harmonic_phase_deg": float(
                    math.degrees(np.angle(spectrum[dominant_bin]))
                ),
                "effective_harmonic_count": effective_harmonics,
                "odd_harmonic_power_fraction": odd_fraction,
                "even_harmonic_power_fraction": even_fraction,
                "BF_readout": parity,
                "BFE_phase_network_readout": bfe_label,
                "BFE_phase_network_type": bfe_type,
                "BFE_phase_network_residual": bfe_residual,
                "BFE_numerical_tolerance": phase_tolerance,
                "positive_rotation_power": positive_power,
                "negative_rotation_power": negative_power,
                "rotation_polarization": rotation_polarization,
                "particle_antiparticle_readout": particle_side,
                "dominant_address": address,
                "finite_return_order": addr_n,
                "charge_magnitude_from_address": charge_magnitude,
                "charge_signed_from_rotation": charge_signed,
                "mass_squared_type_gram_det": gram_det,
                "mass_type_sqrt_gram_det": gram_sqrt,
                "branch_overlap_real": float(overlap.real),
                "branch_overlap_imag": float(overlap.imag),
                "half_turn_overlap_real": None if half_turn is None else float(half_turn.real),
                "half_turn_overlap_imag": None if half_turn is None else float(half_turn.imag),
                "half_turn_plus_residual": half_turn_plus_residual,
                "half_turn_minus_residual": half_turn_minus_residual,
                "half_turn_quadratic_residual": half_turn_quadratic_residual,
                "state_return_order": state_order,
                "quadratic_return_order": quadratic_order,
                "cover_ratio": cover_ratio,
                "spin_type_readout": spin_type(parity, cover_ratio),
                "correlation_lifetime_samples": lifetime,
                "correlation_lifetime_definition": "最初に循環自己相関絶対値が1/e未満となる標本遅れ",
                "cyclic_correlations": correlations,
                "zero_closure_abs": float(row_closure[relation]),
                "zero_closure_relative": float(row_closure_relative[relation]),
                "zero_closed": bool(row_closure_relative[relation] < 1e-12),
                "DFT_reconstruction_error": float(
                    np.max(np.abs(reconstructed[relation] - samples))
                ),
            }
        )

        for bin_index in range(n):
            k = int(orders[bin_index])
            component_m, component_n, component_address = address_for(k, n)
            component_charge = float(
                math.sin(math.pi * component_m / component_n) ** 2
            )
            item = {
                "wave_id": wave_id,
                "relation_id": relation + 1,
                "edge": edge_label(edges[relation]),
                "DFT_bin": bin_index,
                "harmonic_order": k,
                "harmonic_kind": harmonic_kind(k),
                "amplitude": float(abs(spectrum[bin_index])),
                "power": float(power[bin_index]),
                "power_fraction_within_wave": float(power[bin_index] / total_power),
                "phase_deg": float(math.degrees(np.angle(spectrum[bin_index]))),
                "present_above_numerical_floor": bool(present[bin_index]),
                "finite_address": component_address,
                "finite_return_order": component_n,
                "BF_of_harmonic": (
                    "F型（奇数倍音）"
                    if k and abs(k) % 2 == 1
                    else "B型（偶数倍音）" if k else "直流"
                ),
                "rotation_branch": "正周回" if k > 0 else "負周回" if k < 0 else "直流",
                "charge_magnitude_from_address": component_charge,
                "zero_closure_is_wave_bundle_property": True,
            }
            harmonics.append(item)
            if k != 0:
                all_charge_records.append(item)

    watch = []
    for target in CONSTANT_WATCH:
        closest = min(
            all_charge_records,
            key=lambda item: abs(item["charge_magnitude_from_address"] - target),
        )
        watch.append(
            {
                "target": target,
                "closest_value": closest["charge_magnitude_from_address"],
                "absolute_difference": abs(
                    closest["charge_magnitude_from_address"] - target
                ),
                "wave_id": closest["wave_id"],
                "harmonic_order": closest["harmonic_order"],
                "address": closest["finite_address"],
            }
        )

    return {
        "schema": OUTPUT_SCHEMA,
        "source": {
            "input_directory": str(input_dir),
            "accepted_seed": manifest["accepted_seed"],
            "seed_selection_protocol": manifest.get("seed_selection_protocol"),
            "prng": manifest["prng"],
            "numpy_version": manifest["numpy_version"],
            "generator": manifest["generator"],
            "generator_sha256": manifest["generator_sha256"],
        },
        "summary": {
            "N": n,
            "M": m,
            "lambda0": float(manifest["derived_only_from_N"]["lambda0"]),
            "relation_wave_state_count": m,
            "phase_distinguished_wave_count": phase_distinguished_count(waves),
            "total_DFT_components": m * n,
            "detected_nonzero_harmonic_components": int(
                sum(row["detected_harmonic_count"] for row in rows)
            ),
            "detected_nonzero_base_wave_components": int(
                sum(row["base_wave_component_count"] for row in rows)
            ),
            "all_relation_waves_zero_closed": all(row["zero_closed"] for row in rows),
            "max_relation_wave_closure_abs": float(np.max(row_closure)),
            "max_relation_wave_closure_relative": float(np.max(row_closure_relative)),
            "nested_total_closure_abs": float(abs(complex(np.sum(waves**2)))),
            "parent_vector_closure_abs": float(abs(complex(parent @ parent))),
            "DFT_reconstruction_max_error": reconstruction_error,
            "base_wave_is_abs_k_1": True,
            "harmonics_were_read_after_generation": True,
            "harmonics_were_preassigned": False,
            "BFE_uses_external_phase_resolution": False,
        },
        "dimensionless_constant_watch": watch,
        "waves": rows,
        "harmonics": harmonics,
    }


def fmt(value: Any, digits: int = 6) -> str:
    if value is None:
        return "—"
    if isinstance(value, (list, tuple)):
        return ",".join(str(x) for x in value)
    number = float(value)
    if not np.isfinite(number):
        return "—"
    return f"{number:.{digits}g}"


def write_csv(
    path: Path, rows: list[dict[str, Any]], fields: list[tuple[str, str]]
) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[label for _, label in fields])
        writer.writeheader()
        for row in rows:
            writer.writerow({label: row[key] for key, label in fields})


def write_run(output: Path, result: dict[str, Any]) -> None:
    output.mkdir(parents=True, exist_ok=False)
    (output / "census.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    wave_fields = [
        ("wave_id", "波ID"),
        ("relation_id", "関係ID"),
        ("edge", "関係辺"),
        ("wave_amplitude", "波振幅"),
        ("parent_phase_deg", "親位相度"),
        ("sample_count_N", "標本数N"),
        ("wave_composition", "波の構成"),
        ("detected_non_dc_component_count", "非直流成分数"),
        ("detected_non_dc_orders", "非直流次数一覧"),
        ("base_wave_component_count", "基本波成分数"),
        ("detected_harmonic_count", "倍音数"),
        ("detected_harmonic_orders", "倍音次数一覧"),
        ("fundamental_amplitude", "基本波振幅"),
        ("fundamental_phase_deg", "基本波位相度"),
        ("dominant_harmonic_order", "最大強度倍音次数"),
        ("dominant_harmonic_amplitude", "最大強度倍音振幅"),
        ("dominant_harmonic_phase_deg", "最大強度倍音位相度"),
        ("effective_harmonic_count", "実効非直流成分数"),
        ("odd_harmonic_power_fraction", "奇数倍音強度比"),
        ("even_harmonic_power_fraction", "偶数倍音強度比"),
        ("BFE_phase_network_readout", "B/F/E位相網"),
        ("BFE_phase_network_type", "B/F/E位相網型"),
        ("BFE_phase_network_residual", "B/F/E判定残差"),
        ("BF_readout", "B/F読出し"),
        ("particle_antiparticle_readout", "粒子反粒子読出し"),
        ("rotation_polarization", "正負周回偏極"),
        ("dominant_address", "最大強度倍音住所"),
        ("finite_return_order", "有限回帰次数"),
        ("charge_magnitude_from_address", "住所電荷量"),
        ("charge_signed_from_rotation", "符号付き住所電荷量"),
        ("mass_squared_type_gram_det", "質量二乗型Gram量"),
        ("mass_type_sqrt_gram_det", "質量型sqrtGram量"),
        ("half_turn_overlap_real", "半周重なり実部"),
        ("half_turn_plus_residual", "半周同符号残差"),
        ("half_turn_minus_residual", "半周反符号残差"),
        ("half_turn_quadratic_residual", "半周二次量残差"),
        ("state_return_order", "状態回帰次数"),
        ("quadratic_return_order", "二次量回帰次数"),
        ("cover_ratio", "半周期被覆比"),
        ("spin_type_readout", "スピン型読出し"),
        ("correlation_lifetime_samples", "相関寿命型量標本"),
        ("zero_closure_abs", "零閉塞絶対残差"),
        ("zero_closure_relative", "零閉塞相対残差"),
        ("zero_closed", "零閉塞成立"),
    ]
    harmonic_fields = [
        ("wave_id", "波ID"),
        ("relation_id", "関係ID"),
        ("edge", "関係辺"),
        ("DFT_bin", "DFTビン"),
        ("harmonic_order", "倍音次数"),
        ("harmonic_kind", "倍音種別"),
        ("amplitude", "倍音振幅"),
        ("power", "倍音強度"),
        ("power_fraction_within_wave", "波内強度比"),
        ("phase_deg", "倍音位相度"),
        ("present_above_numerical_floor", "数値誤差床より上"),
        ("finite_address", "有限住所"),
        ("finite_return_order", "有限回帰次数"),
        ("BF_of_harmonic", "倍音B/F"),
        ("rotation_branch", "回転枝"),
        ("charge_magnitude_from_address", "住所電荷量"),
    ]
    write_csv(output / "particle_waves_ja.csv", result["waves"], wave_fields)
    write_csv(output / "harmonics_ja.csv", result["harmonics"], harmonic_fields)

    summary = result["summary"]
    source = result["source"]
    lines = [
        f"# N={summary['N']} 白色零閉塞・粒子表",
        "",
        "## 1. 状態数と閉塞",
        "",
        "| 項目 | 値 |",
        "|---|---:|",
        f"| N | {summary['N']} |",
        f"| M=N(N−1)/2 | {summary['M']} |",
        f"| 関係波状態数 | {summary['relation_wave_state_count']} |",
        f"| 位相を区別した波数 | {summary['phase_distinguished_wave_count']} |",
        f"| 採用seed | `{source['accepted_seed']}` |",
        f"| 基本波長 λ0=2π/N | {summary['lambda0']:.12g} |",
        f"| 全DFT成分数 | {summary['total_DFT_components']} |",
        f"| 数値誤差床より上の基本波成分 | {summary['detected_nonzero_base_wave_components']} |",
        f"| 数値誤差床より上の倍音成分（|k|≥2） | {summary['detected_nonzero_harmonic_components']} |",
        f"| 全関係波が零閉塞 | {'はい' if summary['all_relation_waves_zero_closed'] else 'いいえ'} |",
        f"| 最大零閉塞相対残差 | {summary['max_relation_wave_closure_relative']:.3e} |",
        f"| 全体系の零閉塞絶対残差 | {summary['nested_total_closure_abs']:.3e} |",
        f"| DFT逆変換最大誤差 | {summary['DFT_reconstruction_max_error']:.3e} |",
        "",
        "倍音は生成時に置いていない。各関係波の生成後にN点DFTを行い、基本波 $|k|=1$ と",
        "倍音 $|k|\\ge2$ を分離した。正負次数の振幅・位相を含む全成分は `harmonics_ja.csv` にある。",
        "",
        "## 2. 粒子主表（一波一行）",
        "",
        "この表が主成果である。B/F/E、倍音偶奇、電荷型量、質量型量、スピン型量を同じ波の一行で読む。",
        "",
        "| 波ID | 関係 | 波の構成 | 倍音数 | 支配次数 | 支配位相 | B/F/E | 倍音偶奇 | 粒子／反粒子 | 電荷型 q | 質量型 μ | スピン型 | 寿命型 | 零閉塞残差 |",
        "|---|---|---|---:|---:|---:|---|---|---|---:|---:|---|---:|---:|",
    ]
    for row in result["waves"]:
        lines.append(
            f"| {row['wave_id']} | {row['edge']} | {row['wave_composition']} | "
            f"{row['detected_harmonic_count']} | {row['dominant_harmonic_order']:+d} | "
            f"{row['dominant_harmonic_phase_deg']:+.3f}° | "
            f"{row['BFE_phase_network_readout']} | {row['BF_readout']} | "
            f"{row['particle_antiparticle_readout']} | "
            f"{row['charge_signed_from_rotation']:+.9g} | "
            f"{row['mass_type_sqrt_gram_det']:.7g} | {row['spin_type_readout']} | "
            f"{fmt(row['correlation_lifetime_samples'], 7)} | "
            f"{row['zero_closure_relative']:.3e} |"
        )
    lines.extend(
        [
            "",
            "B/F/Eは位相網の読出しであり、B=全成分同相、F=逆相二成分、E=中間位相網である。",
            "外部の位相分解能や144セルは使わず、N点DFTの相対位相そのものを数値誤差内で判定した。",
            "",
            "## 3. 読出し量の式",
            "",
            "- 質量二乗型量: $\\det\\Gamma=N_+N_- - |\\langle c_+|c_-\\rangle|^2$。",
            "- 住所電荷型量: 既約住所 $m/n$ に対して $q_{\\mathrm{addr}}=\\sin^2(\\pi m/n)$。符号は正負周回偏極から読む。",
            "- スピン型量: 半周移動で $T^{N/2}w=+w$ なら1:1、$T^{N/2}w=-w$ かつ二次量が戻れば2:1被覆として読む。状態と二次量の最小回帰次数も併記する。",
            "- 寿命型量: 循環自己相関 $|\\langle w,T^\\ell w\\rangle|/\\lVert w\\rVert^2$ が初めて $1/e$ 未満になる標本遅れ。",
            "",
            "## 4. 波・倍音・閉塞",
            "",
            "| 波ID | 関係 | 波振幅 | 親位相 | 波の構成 | 基本波成分数 | 倍音数 | 倍音次数 | 基本波振幅 | 基本波位相 | 最大強度次数 | 同振幅 | 同位相 | 実効成分数 | 零閉塞相対残差 |",
            "|---|---|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in result["waves"]:
        lines.append(
            f"| {row['wave_id']} | {row['edge']} | {row['wave_amplitude']:.7g} | "
            f"{row['parent_phase_deg']:+.3f}° | {row['wave_composition']} | "
            f"{row['base_wave_component_count']} | {row['detected_harmonic_count']} | "
            f"{fmt(row['detected_harmonic_orders'])} | {row['fundamental_amplitude']:.7g} | "
            f"{row['fundamental_phase_deg']:+.3f}° | {row['dominant_harmonic_order']:+d} | "
            f"{row['dominant_harmonic_amplitude']:.7g} | {row['dominant_harmonic_phase_deg']:+.3f}° | "
            f"{row['effective_harmonic_count']:.5g} | {row['zero_closure_relative']:.3e} |"
        )
    lines.extend(
        [
            "",
            "## 5. 質量型・電荷型・粒子／反粒子",
            "",
            "| 波ID | 質量二乗型 detΓ | 質量型 √detΓ | 最大強度住所 | 有限回帰次数 | 住所電荷量 | 符号付き住所電荷量 | 正負周回偏極 | 粒子／反粒子 |",
            "|---|---:|---:|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in result["waves"]:
        lines.append(
            f"| {row['wave_id']} | {row['mass_squared_type_gram_det']:.7g} | "
            f"{row['mass_type_sqrt_gram_det']:.7g} | {row['dominant_address']} | "
            f"{row['finite_return_order']} | {row['charge_magnitude_from_address']:.9g} | "
            f"{row['charge_signed_from_rotation']:+.9g} | {row['rotation_polarization']:+.7f} | "
            f"{row['particle_antiparticle_readout']} |"
        )
    lines.extend(
        [
            "",
            "## 6. B/F/E・倍音偶奇・スピン型・相関寿命型量",
            "",
            "| 波ID | B/F/E | 位相網型 | 奇数倍音強度比 | 偶数倍音強度比 | 倍音偶奇読出し | 半周重なりRe | 状態回帰 | 二次量回帰 | 被覆比 | スピン型読出し | 相関寿命型量 |",
            "|---|---|---|---:|---:|---|---:|---:|---:|---:|---|---:|",
        ]
    )
    for row in result["waves"]:
        lines.append(
            f"| {row['wave_id']} | {row['BFE_phase_network_readout']} | "
            f"{row['BFE_phase_network_type']} | "
            f"{row['odd_harmonic_power_fraction']:.7f} | "
            f"{row['even_harmonic_power_fraction']:.7f} | {row['BF_readout']} | "
            f"{fmt(row['half_turn_overlap_real'], 7)} | {row['state_return_order']} | "
            f"{row['quadratic_return_order']} | {fmt(row['cover_ratio'], 7)} | "
            f"{row['spin_type_readout']} | "
            f"{fmt(row['correlation_lifetime_samples'], 7)} |"
        )
    lines.extend(
        [
            "",
            "相関寿命型量は、N点波形の循環自己相関の絶対値が初めて $1/e$ 未満になる標本遅れである。",
            "物理時間への換算は行っていない。半周重なりはNが偶数の場合だけ定義する。",
            "",
            "## 7. 無次元数監視",
            "",
            "| 監視値 | 最近接住所電荷量 | 差 | 波ID | 倍音次数 | 住所 |",
            "|---:|---:|---:|---|---:|---|",
        ]
    )
    for item in result["dimensionless_constant_watch"]:
        lines.append(
            f"| {item['target']:.6f} | {item['closest_value']:.12g} | "
            f"{item['absolute_difference']:.3e} | {item['wave_id']} | "
            f"{item['harmonic_order']:+d} | {item['address']} |"
        )
    (output / "particle_table.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def combined_markdown(results: list[dict[str, Any]]) -> str:
    lines = [
        "# N=5・N=40 白色零閉塞・粒子表まとめ",
        "",
        "| N | M | 関係波状態数 | 位相区別波数 | 全DFT成分 | 基本波成分 | 倍音成分 | 全波零閉塞 | 最大相対残差 | 全体系残差 | seed |",
        "|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|",
    ]
    for result in results:
        summary = result["summary"]
        lines.append(
            f"| {summary['N']} | {summary['M']} | {summary['relation_wave_state_count']} | "
            f"{summary['phase_distinguished_wave_count']} | {summary['total_DFT_components']} | "
            f"{summary['detected_nonzero_base_wave_components']} | "
            f"{summary['detected_nonzero_harmonic_components']} | "
            f"{'はい' if summary['all_relation_waves_zero_closed'] else 'いいえ'} | "
            f"{summary['max_relation_wave_closure_relative']:.3e} | "
            f"{summary['nested_total_closure_abs']:.3e} | {result['source']['accepted_seed']} |"
        )
    lines.extend(
        [
            "",
            "## 粒子主表",
            "",
            "| N | 波ID | 関係 | 波の構成 | 倍音数 | 支配次数 | 支配位相 | B/F/E | 倍音偶奇 | 粒子／反粒子 | 電荷型 q | 質量型 μ | スピン型 | 寿命型 | 零閉塞残差 |",
            "|---:|---|---|---|---:|---:|---:|---|---|---|---:|---:|---|---:|---:|",
        ]
    )
    for result in results:
        n = result["summary"]["N"]
        for row in result["waves"]:
            lines.append(
                f"| {n} | {row['wave_id']} | {row['edge']} | {row['wave_composition']} | "
                f"{row['detected_harmonic_count']} | {row['dominant_harmonic_order']:+d} | "
                f"{row['dominant_harmonic_phase_deg']:+.3f}° | "
                f"{row['BFE_phase_network_readout']} | {row['BF_readout']} | "
                f"{row['particle_antiparticle_readout']} | "
                f"{row['charge_signed_from_rotation']:+.9g} | "
                f"{row['mass_type_sqrt_gram_det']:.7g} | {row['spin_type_readout']} | "
                f"{fmt(row['correlation_lifetime_samples'], 7)} | "
                f"{row['zero_closure_relative']:.3e} |"
            )
    lines.extend(
        [
            "",
            "各Nの `particle_table.md` に同じM本全行と計算内訳を、`particle_waves_ja.csv` に同じ全行を、",
            "`harmonics_ja.csv` に全M×N倍音成分の次数・振幅・位相を保存した。",
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

    results = []
    seen = set()
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
            f"waves={result['summary']['relation_wave_state_count']} "
            f"harmonic_components={result['summary']['total_DFT_components']}"
        )

    results.sort(key=lambda item: item["summary"]["N"])
    (output / "summary.md").write_text(
        combined_markdown(results), encoding="utf-8"
    )
    (output / "summary.json").write_text(
        json.dumps(
            {
                "schema": OUTPUT_SCHEMA,
                "reader": Path(__file__).name,
                "reader_sha256": sha256_file(Path(__file__).resolve()),
                "results": results,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
