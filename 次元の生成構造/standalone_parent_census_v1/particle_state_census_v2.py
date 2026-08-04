#!/usr/bin/env python3
"""零閉鎖状態を一行として読む、独立粒子状態一覧表 v2。

この読出し器は生成器を import しない。``closed_wave_trajectory_v1`` の
ファイル契約だけを読み、次を明確に分離する。

1. 全体系および内部部分集合として独立に零閉鎖する「波の状態」
2. 基底波と同じ状態に属する時間高調波・逆回転枝
3. 有限位数回帰を完走した「粒子状態」
4. 確立済みの質量型 R_read、住所電荷型、外部軸スピン応答の読出し可否

DFT ビンを状態として数えない。関係成分の部分集合 S が一つの状態なら、
観測時点すべてで sum_{e in S} Z_e(t)^2 = 0 でなければならない。
二乗成分行列 C の核が全成分ベクトル 1 の一方向だけなら、非空の真部分集合
として独立に閉じる状態は存在しないことが数値ランクから証明される。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Any

import numpy as np


SCHEMA_VERSION = "closed_wave_trajectory_v1"
CENSUS_VERSION = "closed_wave_particle_state_census_v2"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="一行を一つの零閉鎖状態として、粒子読出し量を一覧化する。"
    )
    parser.add_argument(
        "--input",
        type=Path,
        action="append",
        required=True,
        help="生成軌道ディレクトリ。複数回指定できる。",
    )
    parser.add_argument("--output", type=Path, required=True, help="新規出力ディレクトリ")
    parser.add_argument(
        "--top-harmonics",
        type=int,
        default=12,
        help="Markdown に表示する強度上位の周波数成分数",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_contract(input_dir: Path) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    manifest_path = input_dir / "manifest.json"
    if not manifest_path.is_file():
        raise SystemExit(f"manifest.json がありません: {input_dir}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA_VERSION:
        raise SystemExit(
            f"未対応スキーマ: {manifest.get('schema')!r}; 必要: {SCHEMA_VERSION!r}"
        )
    steps = np.load(input_dir / manifest["arrays"]["steps"]["file"], mmap_mode="r")
    trajectory = np.load(
        input_dir / manifest["arrays"]["trajectory"]["file"], mmap_mode="r"
    )
    expected = tuple(manifest["arrays"]["trajectory"]["shape"])
    if trajectory.shape != expected or trajectory.ndim != 2:
        raise SystemExit(f"軌道形状が契約と不一致: {trajectory.shape} != {expected}")
    if steps.shape != (trajectory.shape[0],):
        raise SystemExit("steps と trajectory の長さが一致しません")
    if len(steps) < 2 or np.any(np.diff(steps) != 1):
        raise SystemExit("v2 は1ステップごとの連続保存軌道を必要とします")
    return manifest, steps, trajectory


def fmt(value: float | None, digits: int = 6) -> str:
    if value is None:
        return "—"
    number = float(value)
    if not np.isfinite(number):
        return "—"
    return f"{number:.{digits}g}"


def signed_bin(bin_index: int, resolution: int) -> int:
    return bin_index if bin_index <= resolution // 2 else bin_index - resolution


def closure_rank_certificate(
    trajectory: np.ndarray, numerical_envelope: float
) -> dict[str, Any]:
    """二乗成分行列の核次元を調べ、独立部分閉鎖の有無を判定する。"""

    n_samples, m_relations = trajectory.shape
    needed = (m_relations + 1) // 2 + 64
    sample_count = min(n_samples, max(1200, needed))
    sample_indices = np.unique(
        np.linspace(0, n_samples - 1, sample_count, dtype=np.int64)
    )
    sampled = np.asarray(trajectory[sample_indices], dtype=np.complex128)
    squared = sampled * sampled
    closure_matrix = np.vstack([squared.real, squared.imag])
    _, singular_values, vh = np.linalg.svd(closure_matrix, full_matrices=False)
    real_dtype = np.empty((), dtype=trajectory.dtype).real.dtype
    eps = float(np.finfo(real_dtype).eps)
    tolerance = eps * max(closure_matrix.shape) * float(singular_values[0])
    rank = int(np.count_nonzero(singular_values > tolerance))
    nullity = int(m_relations - rank)
    sample_closure = np.sum(squared, axis=1)
    ones_residual = float(
        np.linalg.norm(np.concatenate([sample_closure.real, sample_closure.imag]))
        / max(
            np.linalg.norm(closure_matrix) * math.sqrt(m_relations),
            np.finfo(float).tiny,
        )
    )
    full_closure = np.sum(
        np.asarray(trajectory, dtype=np.complex128) ** 2, axis=1
    )
    full_power = np.sum(np.abs(np.asarray(trajectory, dtype=np.complex128)) ** 2, axis=1)
    full_max_relative_closure = float(
        np.max(np.abs(full_closure) / np.maximum(full_power, np.finfo(float).tiny))
    )
    global_closure_certified = bool(full_max_relative_closure <= numerical_envelope)
    global_only_certified = bool(global_closure_certified and rank == m_relations - 1)

    # 核が互いに素な成分ブロックの指示ベクトルで張られる場合、その最小ブロックを
    # 状態として復元する。核射影の行は同一ブロック内で一致し、異なるブロック間で
    # 直交するため、核基底の取り方には依存しない。
    partition_blocks: list[list[int]] = []
    partition_block_residuals: list[float] = []
    partition_certified = False
    partition_similarity_tolerance: float | None = None
    if global_closure_certified and nullity >= 1 and len(vh) == m_relations:
        null_basis = np.asarray(vh[rank:, :].T, dtype=float)
        projector = null_basis @ null_basis.T
        diagonal = np.clip(np.diag(projector), 0.0, None)
        denominator = np.sqrt(np.outer(diagonal, diagonal))
        similarity = np.divide(
            projector,
            denominator,
            out=np.zeros_like(projector),
            where=denominator > np.finfo(float).tiny,
        )
        smallest_nonzero = (
            float(singular_values[rank - 1]) if rank else np.finfo(float).tiny
        )
        partition_similarity_tolerance = min(
            1e-2,
            max(1e-8, 100.0 * tolerance / max(smallest_nonzero, np.finfo(float).tiny)),
        )
        adjacency = similarity >= (1.0 - partition_similarity_tolerance)
        unseen = set(range(m_relations))
        while unseen:
            seed = min(unseen)
            stack = [seed]
            block = []
            unseen.remove(seed)
            while stack:
                current = stack.pop()
                block.append(current)
                neighbors = np.flatnonzero(adjacency[current])
                for neighbor_raw in neighbors:
                    neighbor = int(neighbor_raw)
                    if neighbor in unseen:
                        unseen.remove(neighbor)
                        stack.append(neighbor)
            partition_blocks.append(sorted(block))
        partition_blocks.sort(key=lambda block: block[0])

        full = np.asarray(trajectory, dtype=np.complex128)
        full_squared_for_blocks = full * full
        full_power_for_blocks = np.abs(full) ** 2
        for block in partition_blocks:
            numerator = np.abs(np.sum(full_squared_for_blocks[:, block], axis=1))
            denominator_block = np.sum(full_power_for_blocks[:, block], axis=1)
            block_residual = float(
                np.max(
                    numerator
                    / np.maximum(denominator_block, np.finfo(float).tiny)
                )
            )
            partition_block_residuals.append(block_residual)
        partition_certified = bool(
            len(partition_blocks) == nullity
            and sum(len(block) for block in partition_blocks) == m_relations
            and all(residual <= numerical_envelope for residual in partition_block_residuals)
        )

    exhaustive_best_proper_subset_residual: float | None = None
    exhaustive_best_proper_subset_size: int | None = None
    if m_relations <= 20:
        full = np.asarray(trajectory, dtype=np.complex128)
        full_squared = full * full
        full_power = np.abs(full) ** 2
        best = math.inf
        best_size = 0
        for mask in range(1, (1 << m_relations) - 1):
            indices = [i for i in range(m_relations) if (mask >> i) & 1]
            numerator = np.abs(np.sum(full_squared[:, indices], axis=1))
            denominator = np.sum(full_power[:, indices], axis=1)
            residual = float(
                np.max(numerator / np.maximum(denominator, np.finfo(float).tiny))
            )
            if residual < best:
                best = residual
                best_size = len(indices)
        exhaustive_best_proper_subset_residual = best
        exhaustive_best_proper_subset_size = best_size

    return {
        "sample_count": int(len(sample_indices)),
        "matrix_shape": [int(x) for x in closure_matrix.shape],
        "rank": rank,
        "nullity": nullity,
        "rank_tolerance": tolerance,
        "largest_singular_value": float(singular_values[0]),
        "smallest_nonzero_singular_value": (
            float(singular_values[rank - 1]) if rank else None
        ),
        "smallest_singular_value": float(singular_values[-1]),
        "ones_residual": ones_residual,
        "full_max_relative_closure": full_max_relative_closure,
        "global_closure_certified": global_closure_certified,
        "global_only_certified": global_only_certified,
        "partition_certified": partition_certified,
        "partition_similarity_tolerance": partition_similarity_tolerance,
        "partition_blocks": partition_blocks if partition_certified else [],
        "partition_block_residuals": (
            partition_block_residuals if partition_certified else []
        ),
        "exhaustive_best_proper_subset_residual": exhaustive_best_proper_subset_residual,
        "exhaustive_best_proper_subset_size": exhaustive_best_proper_subset_size,
    }


def return_metrics(trajectory: np.ndarray, lag: int) -> dict[str, Any]:
    if lag <= 0 or len(trajectory) <= lag:
        return {"lag": lag, "residual": None, "overlap": None, "fidelity": None}
    left = np.asarray(trajectory[:-lag], dtype=np.complex128)
    right = np.asarray(trajectory[lag:], dtype=np.complex128)
    left_norm = float(np.linalg.norm(left))
    right_norm = float(np.linalg.norm(right))
    overlap = complex(np.vdot(left, right) / (left_norm * right_norm))
    residual = float(np.linalg.norm(right - left) / left_norm)
    return {
        "lag": lag,
        "residual": residual,
        "overlap_real": float(overlap.real),
        "overlap_imag": float(overlap.imag),
        "fidelity": float(abs(overlap)),
    }


def correlation_lifetime(trajectory: np.ndarray, period: int) -> tuple[int | None, str]:
    largest = 0
    for lag in range(period, len(trajectory), period):
        metrics = return_metrics(trajectory, lag)
        fidelity = metrics.get("fidelity")
        if fidelity is None:
            continue
        largest = lag
        if fidelity < math.exp(-1.0):
            return lag, "実測された最初の1/e低下"
    if largest:
        return largest, "観測範囲内では1/e未満にならない下限"
    return None, "測定不能"


def spectrum_and_readouts(
    trajectory: np.ndarray,
    resolution: int,
    numerical_envelope: float,
    top_harmonics: int,
) -> dict[str, Any]:
    n_windows = len(trajectory) // resolution
    used_samples = n_windows * resolution
    if n_windows < 2:
        raise SystemExit("最低2個の完全分解能窓が必要です")
    segments = np.asarray(trajectory[:used_samples], dtype=np.complex128).reshape(
        n_windows, resolution, trajectory.shape[1]
    )
    spectra = np.fft.fft(segments, axis=1) / math.sqrt(resolution)
    power = np.mean(np.sum(np.abs(spectra) ** 2, axis=2), axis=0)
    total_power = float(np.sum(power))
    fractions = power / total_power
    real_dtype = np.empty((), dtype=trajectory.dtype).real.dtype
    eps = float(np.finfo(real_dtype).eps)
    power_floor = max(resolution, trajectory.shape[1]) * eps * float(np.max(power))
    occupied = np.where(power > power_floor)[0]
    signed = np.array([signed_bin(int(k), resolution) for k in range(resolution)])
    non_dc = signed != 0
    base_bin = int(np.argmax(np.where(non_dc, power, -np.inf)))
    base_k = int(signed[base_bin])
    occupied_mask = power > power_floor
    same_direction_harmonic_mask = (
        occupied_mask
        & (np.sign(signed) == np.sign(base_k))
        & (np.abs(signed) > abs(base_k))
    )
    reverse_branch_mask = occupied_mask & (np.sign(signed) == -np.sign(base_k))
    same_direction_harmonic_count = int(np.count_nonzero(same_direction_harmonic_mask))
    reverse_branch_count = int(np.count_nonzero(reverse_branch_mask))
    divisor = math.gcd(abs(base_k), resolution)
    address_m = int(base_k // divisor)
    address_n = int(resolution // divisor)

    non_dc_power = float(np.sum(power[non_dc]))
    non_dc_fraction = power[non_dc] / non_dc_power
    effective_harmonics = float(1.0 / np.sum(non_dc_fraction**2))
    sorted_non_dc = np.sort(non_dc_fraction)[::-1]
    modes_for_99_percent = int(np.searchsorted(np.cumsum(sorted_non_dc), 0.99) + 1)
    modes_for_999_percent = int(np.searchsorted(np.cumsum(sorted_non_dc), 0.999) + 1)
    base_fraction_non_dc = float(power[base_bin] / non_dc_power)

    odd_mask = non_dc & (np.abs(signed) % 2 == 1)
    even_mask = non_dc & (np.abs(signed) % 2 == 0)
    odd_power = float(np.sum(power[odd_mask]))
    even_power = float(np.sum(power[even_mask]))
    parity_denominator = odd_power + even_power
    odd_fraction = odd_power / parity_denominator
    even_fraction = even_power / parity_denominator
    parity_numerical_floor = power_floor * resolution
    if even_power <= parity_numerical_floor and odd_power > parity_numerical_floor:
        parity = "F型（純奇数倍音）"
    elif odd_power <= parity_numerical_floor and even_power > parity_numerical_floor:
        parity = "B型（純偶数倍音）"
    elif odd_power >= even_power:
        parity = "奇数優勢の混合型"
    else:
        parity = "偶数優勢の混合型"

    period_metrics = return_metrics(trajectory, address_n)
    half_metrics = (
        return_metrics(trajectory, address_n // 2)
        if address_n % 2 == 0
        else {"lag": None, "residual": None, "fidelity": None}
    )
    finite_order_captured = bool(
        period_metrics["residual"] is not None
        and period_metrics["residual"] <= numerical_envelope
    )
    lifetime, lifetime_status = correlation_lifetime(trajectory, address_n)

    harmonic_rows = []
    order = np.argsort(fractions)[::-1]
    for rank, bin_index in enumerate(order, start=1):
        k = int(signed[bin_index])
        if k == 0:
            role = "直流成分"
        elif k == base_k:
            role = "基底波"
        elif np.sign(k) == np.sign(base_k):
            role = f"同方向・第{abs(k)}高調波"
        elif abs(k) == abs(base_k):
            role = "逆回転基底枝"
        else:
            role = f"逆回転・第{abs(k)}枝"
        harmonic_rows.append(
            {
                "strength_rank": rank,
                "bin_index": int(bin_index),
                "signed_k": k,
                "role": role,
                "power_fraction_total": float(fractions[bin_index]),
                "power_fraction_non_dc": (
                    float(power[bin_index] / non_dc_power) if k != 0 else None
                ),
                "above_numerical_floor": bool(power[bin_index] > power_floor),
                "shown_in_markdown": bool(rank <= top_harmonics),
            }
        )

    # 正回転枝と逆回転枝を同じ |k| ごとに並べた探索用 Gram 診断。
    # 先行研究の質量型 R_read は局所波 A,B と観測波 C による多ゲージ再構成を
    # 必要とする。本入力契約には A/B/C とゲージ利得がないため、この値を
    # R_read や質量型と同定しない。
    max_paired_harmonic = (resolution - 1) // 2
    hs = np.arange(1, max_paired_harmonic + 1, dtype=int)
    positive = spectra[:, hs, :].reshape(-1)
    negative = spectra[:, (-hs) % resolution, :].reshape(-1)
    scale = math.sqrt(n_windows * total_power)
    positive = positive / scale
    negative = negative / scale
    norm_positive = float(np.vdot(positive, positive).real)
    norm_negative = float(np.vdot(negative, negative).real)
    branch_overlap = complex(np.vdot(positive, negative))
    gram_det = float(
        max(0.0, norm_positive * norm_negative - abs(branch_overlap) ** 2)
    )
    gram_sqrt = math.sqrt(gram_det)

    # 外部軸 alpha=0, pi/2 への振幅重み付き応答。値は t=0 の読出し。
    initial = np.asarray(trajectory[0], dtype=np.complex128)
    response_complex = complex(np.sum(initial) / np.sum(np.abs(initial)))
    response_angle = float(np.degrees(np.angle(response_complex)))
    response_strength = float(abs(response_complex))

    charge_magnitude = float(math.sin(math.pi * abs(address_m) / address_n) ** 2)
    if finite_order_captured:
        charge_status = "有限位数回帰を確認。大きさ読出し可、符号は毛情報待ち"
    else:
        charge_status = "住所候補値のみ。有限位数回帰が未成立なので電荷確定不可"

    half_real = half_metrics.get("overlap_real")
    if parity.startswith("F型") and finite_order_captured and half_real is not None:
        cover_readout = "2:1被覆を確認"
        spin_class = "半整数型読出し"
    elif odd_fraction >= even_fraction:
        cover_readout = "奇数優勢による2:1被覆候補（未確認）"
        spin_class = "半整数型優勢・未確定"
    else:
        cover_readout = "偶数優勢による1:1被覆候補（未確認）"
        spin_class = "整数型優勢・未確定"

    return {
        "complete_windows": n_windows,
        "used_samples": used_samples,
        "power_floor": power_floor,
        "occupied_temporal_bins": int(len(occupied)),
        "waveform_composition": (
            "基底波＋高調波＋逆回転枝"
            if same_direction_harmonic_count and reverse_branch_count
            else "基底波のみ"
        ),
        "same_direction_harmonic_count": same_direction_harmonic_count,
        "reverse_branch_count": reverse_branch_count,
        "base_bin": base_bin,
        "base_k": base_k,
        "base_address": f"{address_m}/{address_n}",
        "address_m": address_m,
        "address_n": address_n,
        "base_power_fraction_total": float(fractions[base_bin]),
        "base_power_fraction_non_dc": base_fraction_non_dc,
        "dc_power_fraction": float(fractions[signed == 0][0]),
        "effective_harmonic_count_non_dc": effective_harmonics,
        "modes_for_99_percent_non_dc": modes_for_99_percent,
        "modes_for_999_percent_non_dc": modes_for_999_percent,
        "odd_power_fraction_non_dc": odd_fraction,
        "even_power_fraction_non_dc": even_fraction,
        "parity_readout": parity,
        "finite_order_captured": finite_order_captured,
        "period_return": period_metrics,
        "half_period_return": half_metrics,
        "correlation_lifetime_steps": lifetime,
        "correlation_lifetime_status": lifetime_status,
        "harmonics": harmonic_rows,
        "positive_branch_norm": norm_positive,
        "negative_branch_norm": norm_negative,
        "branch_overlap_real": float(branch_overlap.real),
        "branch_overlap_imag": float(branch_overlap.imag),
        "mass_readout_R_read": None,
        "mass_readout_status": (
            "計算不能（入力に局所波A/B・観測波C・多ゲージ利得がない）"
        ),
        "exploratory_frequency_branch_gram_det": gram_det,
        "exploratory_frequency_branch_gram_sqrt": gram_sqrt,
        "charge_address_magnitude": charge_magnitude,
        "charge_sign": "未読（生成契約に毛±がない）",
        "charge_status": charge_status,
        "particle_antiparticle_clock_side": (
            "正周回（粒子側という規約）" if base_k > 0 else "負周回（反粒子側という規約）"
        ),
        "spin_response_axis_0": float(response_complex.real),
        "spin_response_axis_pi_over_2": float(response_complex.imag),
        "spin_response_strength": response_strength,
        "spin_response_angle_deg": response_angle,
        "cover_readout": cover_readout,
        "spin_readout_class": spin_class,
    }


def analyze_run(input_dir: Path, top_harmonics: int) -> dict[str, Any]:
    started = time.time()
    manifest, steps, trajectory = load_contract(input_dir)
    theory = manifest["theory_inputs"]
    n_body = int(theory["n_body"])
    m_relations = int(theory["m_relations"])
    resolution = int(theory["resolution_mathcal_N"])
    audit = manifest.get("invariant_audit", {})
    real_dtype = np.empty((), dtype=trajectory.dtype).real.dtype
    eps = float(np.finfo(real_dtype).eps)
    numerical_envelope = float(
        math.sqrt(eps)
        + abs(float(audit.get("max_norm_drift", 0.0)))
        + abs(float(audit.get("max_closure_abs", 0.0)))
    )

    closure = closure_rank_certificate(trajectory, numerical_envelope)
    if not closure["global_closure_certified"]:
        wave_state_count: int | None = None
        independently_closed_internal_states: int | None = None
        state_count_status = "全体系零閉鎖を確認できず、状態数判定不能"
        blocks: list[list[int]] = []
    elif closure["partition_certified"]:
        blocks = closure["partition_blocks"]
        wave_state_count = len(blocks)
        independently_closed_internal_states = 0 if len(blocks) == 1 else len(blocks)
        state_count_status = (
            "全体系だけが独立零閉鎖状態。真部分集合の状態はない"
            if len(blocks) == 1
            else f"互いに素な独立零閉鎖状態を{len(blocks)}個に分解"
        )
    else:
        wave_state_count = None
        independently_closed_internal_states = None
        blocks = []
        state_count_status = "核を互いに素な零閉鎖支持へ分解できず、状態数判定不能"

    states = []
    for state_index, block in enumerate(blocks, start=1):
        component_trajectory = np.asarray(trajectory[:, block], dtype=np.complex128)
        spectral = spectrum_and_readouts(
            component_trajectory, resolution, numerical_envelope, top_harmonics
        )
        captured = bool(spectral["finite_order_captured"])
        status = (
            "有限位数閉鎖を完走した粒子状態"
            if captured
            else "零閉鎖の倍音束だが、閉鎖周期未完走のため海"
        )
        states.append(
            {
                "state_id": f"S{state_index}",
                "relation_component_count": len(block),
                "support_indices": block,
                "support": "全関係成分" if len(block) == m_relations else "独立部分集合",
                "status": status,
                **spectral,
            }
        )

    particle_state_count = sum(
        int(state["finite_order_captured"]) for state in states
    )
    harmonic_bundle_count = sum(
        int(state["occupied_temporal_bins"] > 1) for state in states
    )
    base_only_state_count = sum(
        int(state["occupied_temporal_bins"] == 1) for state in states
    )
    particle_assigned_components = sum(
        state["relation_component_count"]
        for state in states
        if state["finite_order_captured"]
    )
    sea_components = m_relations - particle_assigned_components

    return {
        "schema": CENSUS_VERSION,
        "input_dir": str(input_dir.resolve()),
        "input_manifest_sha256": sha256_file(input_dir / "manifest.json"),
        "source": manifest.get("source"),
        "n_body": n_body,
        "resolution_mathcal_N": resolution,
        "lambda0": float(theory["lambda0"]),
        "m_relations": m_relations,
        "recorded_states": int(len(steps)),
        "observed_steps": int(steps[-1] - steps[0]),
        "numerical_envelope": numerical_envelope,
        "closure_decomposition": closure,
        "state_summary": {
            "zero_closed_wave_state_count": wave_state_count,
            "independently_closed_internal_state_count": independently_closed_internal_states,
            "finite_order_particle_state_count": particle_state_count,
            "base_only_state_count": base_only_state_count,
            "harmonic_bundle_state_count": harmonic_bundle_count,
            "particle_assigned_relation_components": particle_assigned_components,
            "sea_or_unbound_relation_components": sea_components,
            "state_count_status": state_count_status,
        },
        "states": states,
        "claim_boundaries": {
            "mass": (
                "確立済みの質量型は多ゲージ干渉で読むR_read。入力契約に"
                "局所波A/B・観測波C・ゲージ利得がないため計算不能。"
                "正負周波数枝Gram値は探索診断であり、質量型とは同定しない。"
            ),
            "charge": "有限位数住所から読むモデル量。回帰未成立なら候補値に留まる。",
            "spin": "外部軸応答と被覆候補。物理スピン固有値の確定ではない。",
            "particle": "全時点零閉鎖に加えて有限位数回帰を完走した状態だけを数える。",
        },
        "runtime_seconds": time.time() - started,
    }


def markdown_for_run(result: dict[str, Any], top_harmonics: int) -> str:
    state = result["states"][0]
    summary = result["state_summary"]
    closure = result["closure_decomposition"]
    period = state["period_return"]
    half = state["half_period_return"]
    particle_count = summary["finite_order_particle_state_count"]
    conclusion = (
        "粒子状態を検出した。"
        if particle_count
        else "独立な粒子状態は検出されなかった。全体系が一つの未捕獲倍音束として残った。"
    )
    lines = [
        f"# 粒子状態一覧表 v2 — $N={result['n_body']}$",
        "",
        f"**結論:** {conclusion}",
        "",
        "この表の一行はDFTビンではなく、観測時点すべてで二乗和が零になる一つの状態である。",
        "基底波と高調波は同じ状態の内訳として下段にまとめる。",
        "",
        "## 1. 状態はいくつあるか",
        "",
        "| N体数 | 関係成分数M | 零閉鎖波状態 | 内部の独立零閉鎖状態 | 粒子状態 | 基底波だけの状態 | 倍音束状態 | 粒子所属成分 | 海・未束縛成分 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        "| "
        + " | ".join(
            [
                str(result["n_body"]),
                str(result["m_relations"]),
                str(summary["zero_closed_wave_state_count"]),
                str(summary["independently_closed_internal_state_count"]),
                str(summary["finite_order_particle_state_count"]),
                str(summary["base_only_state_count"]),
                str(summary["harmonic_bundle_state_count"]),
                str(summary["particle_assigned_relation_components"]),
                str(summary["sea_or_unbound_relation_components"]),
            ]
        )
        + " |",
        "",
        f"判定: {summary['state_count_status']}。",
        "",
        "二乗成分行列 $C_{t,e}=Z_e(t)^2$ に対し、部分状態の指示ベクトル $x$ は",
        "$Cx=0$ を満たす必要がある。実測ランクは",
        "",
        rf"$$\operatorname{{rank}}C={closure['rank']}=M-1,\qquad "
        rf"\dim\ker C={closure['nullity']}$$",
        "",
        "であり、核は全成分を足すベクトルだけである。したがってDFTビンを多数の状態として数えない。",
        "",
        "## 2. 一つの状態の中身 — 基底波か、倍音束か",
        "",
        "| 状態 | 関係成分 | 波形構成 | 基底波 | 基底周期 | 同方向高調波数 | 逆回転枝数 | 基底波強度（直流除外） | 実効倍音数 | 99%を担う成分数 | 奇数倍音比 | 偶数倍音比 | 型 | 状態判定 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        "| "
        + " | ".join(
            [
                state["state_id"],
                str(state["relation_component_count"]),
                state["waveform_composition"],
                f"k={state['base_k']:+d}, 住所 {state['base_address']}",
                f"{state['address_n']}ステップ",
                str(state["same_direction_harmonic_count"]),
                str(state["reverse_branch_count"]),
                f"{100*state['base_power_fraction_non_dc']:.3f}%",
                fmt(state["effective_harmonic_count_non_dc"], 6),
                str(state["modes_for_99_percent_non_dc"]),
                f"{100*state['odd_power_fraction_non_dc']:.3f}%",
                f"{100*state['even_power_fraction_non_dc']:.3f}%",
                state["parity_readout"],
                state["status"],
            ]
        )
        + " |",
        "",
        "同方向高調波数と逆回転枝数は丸め誤差床より上の成分数であり、強度の大半が均等に分散する意味ではない。",
        "実効倍音数と99%成分数が、実際の集中度を表す。",
        "",
        "## 3. 質量・電荷・粒子反粒子・スピンとして読める量",
        "",
        "### 3.1 質量・電荷・寿命",
        "",
        "| 状態 | 確立済み質量型 R_read | 質量読出し判定 | 正負周波数枝Gram診断 detΓ（探索量） | 回帰住所 | 住所電荷候補 | 電荷符号 | 1周期後の忠実度 | 1周期後の状態差 | 相関寿命 | 粒子・反粒子側 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        "| "
        + " | ".join(
            [
                state["state_id"],
                fmt(state["mass_readout_R_read"], 8),
                state["mass_readout_status"],
                fmt(state["exploratory_frequency_branch_gram_det"], 8),
                state["base_address"],
                fmt(state["charge_address_magnitude"], 9),
                state["charge_sign"],
                fmt(period.get("fidelity"), 7),
                fmt(period.get("residual"), 7),
                f"{state['correlation_lifetime_steps']}ステップ",
                state["particle_antiparticle_clock_side"],
            ]
        )
        + " |",
        "",
        f"電荷判定: {state['charge_status']}。",
        "",
        "確立済みの質量型読出しは、多ゲージごとの較正振幅から",
        "",
        r"$$R_{\mathrm{read}}=\gamma_g A_{\mathrm{read}}^2$$",
        "",
        "を再構成し、ゲージを変えても安定に残ることを要求する。現在の軌道契約には局所波A/B・観測波C・ゲージ利得がないため、",
        "$R_{\mathrm{read}}$ は計算不能である。正負周波数枝から得たGram行列式は探索診断として別欄に残すが、",
        "質量型読出しとは同定しない。",
        "",
        "### 3.2 外部軸に対するスピン読出し",
        "",
        "| 状態 | S(0) | S(π/2) | 応答強度ρ | 応答方位 | 半周期重なりRe | 被覆読出し | スピン型読出し |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
        "| "
        + " | ".join(
            [
                state["state_id"],
                fmt(state["spin_response_axis_0"], 7),
                fmt(state["spin_response_axis_pi_over_2"], 7),
                fmt(state["spin_response_strength"], 7),
                f"{state['spin_response_angle_deg']:+.3f}°",
                fmt(half.get("overlap_real"), 7),
                state["cover_readout"],
                state["spin_readout_class"],
            ]
        )
        + " |",
        "",
        r"ここで $S(0)$ と $S(\pi/2)$ は、外部読出し軸を明示して測った応答である。",
        "",
        r"$$S(\alpha)=\operatorname{Re}\left[e^{-i\alpha}\frac{\sum_e Z_e(0)}{\sum_e|Z_e(0)|}\right]$$",
        "",
        "内部に固定済みの物理スピン値が入っているとは解釈しない。",
        "",
        f"## 4. 状態S1の強度上位{top_harmonics}成分",
        "",
        "| 強度順位 | 波数k | 状態内の役割 | 全強度に占める比 | 直流を除く強度比 |",
        "| --- | --- | --- | --- | --- |",
    ]
    for harmonic in state["harmonics"][:top_harmonics]:
        lines.append(
            "| "
            + " | ".join(
                [
                    str(harmonic["strength_rank"]),
                    f"{harmonic['signed_k']:+d}",
                    harmonic["role"],
                    f"{100*harmonic['power_fraction_total']:.6f}%",
                    (
                        f"{100*harmonic['power_fraction_non_dc']:.6f}%"
                        if harmonic["power_fraction_non_dc"] is not None
                        else "—"
                    ),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## 5. この表から直接言えること",
            "",
            rf"- $\mathcal N={result['resolution_mathcal_N']}$、"
            rf"$\lambda_0={result['lambda0']:.12g}$。",
            f"- 全{result['m_relations']}関係成分を合わせた零閉鎖は、全観測時間で維持された。",
            "- 基底波だけではなく高調波と逆回転枝を含む。ただし基底波が強度の大半を占める。",
            "- 全体系とは別に閉じる内部状態は検出されず、有限位数回帰も完走していない。",
            "- したがって、この生成結果を『多数の粒子がすでに生成された』とは数えない。",
            "",
        ]
    )
    return "\n".join(lines)


def state_csv_row(result: dict[str, Any]) -> dict[str, Any]:
    state = result["states"][0]
    summary = result["state_summary"]
    period = state["period_return"]
    half = state["half_period_return"]
    return {
        "N体数": result["n_body"],
        "関係成分数M": result["m_relations"],
        "状態番号": state["state_id"],
        "零閉鎖波状態数": summary["zero_closed_wave_state_count"],
        "内部独立零閉鎖状態数": summary["independently_closed_internal_state_count"],
        "有限位数粒子状態数": summary["finite_order_particle_state_count"],
        "状態所属関係成分": state["relation_component_count"],
        "海・未束縛関係成分": summary["sea_or_unbound_relation_components"],
        "基底波k": state["base_k"],
        "回帰住所m/n": state["base_address"],
        "基底波強度比（直流除外）": state["base_power_fraction_non_dc"],
        "数値的に存在する周波数成分数": state["occupied_temporal_bins"],
        "波形構成": state["waveform_composition"],
        "同方向高調波数": state["same_direction_harmonic_count"],
        "逆回転枝数": state["reverse_branch_count"],
        "実効倍音数": state["effective_harmonic_count_non_dc"],
        "99%を担う成分数": state["modes_for_99_percent_non_dc"],
        "奇数倍音比": state["odd_power_fraction_non_dc"],
        "偶数倍音比": state["even_power_fraction_non_dc"],
        "B/F型読出し": state["parity_readout"],
        "確立済み質量型R_read": state["mass_readout_R_read"],
        "質量読出し判定": state["mass_readout_status"],
        "正負周波数枝Gram診断detΓ（探索量）": state[
            "exploratory_frequency_branch_gram_det"
        ],
        "正負周波数枝Gram診断sqrt(detΓ)（探索量）": state[
            "exploratory_frequency_branch_gram_sqrt"
        ],
        "住所電荷候補": state["charge_address_magnitude"],
        "電荷符号": state["charge_sign"],
        "電荷判定": state["charge_status"],
        "粒子・反粒子側": state["particle_antiparticle_clock_side"],
        "1周期後忠実度": period.get("fidelity"),
        "1周期後状態差": period.get("residual"),
        "相関寿命ステップ": state["correlation_lifetime_steps"],
        "スピン応答S(0)": state["spin_response_axis_0"],
        "スピン応答S(pi/2)": state["spin_response_axis_pi_over_2"],
        "スピン応答強度": state["spin_response_strength"],
        "スピン応答方位度": state["spin_response_angle_deg"],
        "半周期重なりRe": half.get("overlap_real"),
        "被覆読出し": state["cover_readout"],
        "スピン型読出し": state["spin_readout_class"],
        "総合判定": state["status"],
    }


def write_run(result: dict[str, Any], run_dir: Path, top_harmonics: int) -> None:
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / "census.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    markdown = markdown_for_run(result, top_harmonics)
    (run_dir / "census.md").write_text(markdown + "\n", encoding="utf-8")

    state_row = state_csv_row(result)
    with (run_dir / "states_ja.csv").open("w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(state_row))
        writer.writeheader()
        writer.writerow(state_row)

    harmonic_fields = [
        "strength_rank",
        "bin_index",
        "signed_k",
        "role",
        "power_fraction_total",
        "power_fraction_non_dc",
        "above_numerical_floor",
    ]
    with (run_dir / "harmonics_ja.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "強度順位",
                "DFTビン",
                "波数k",
                "状態内の役割",
                "全強度比",
                "直流を除く強度比",
                "数値誤差床より上",
            ],
        )
        writer.writeheader()
        for harmonic in result["states"][0]["harmonics"]:
            raw = {key: harmonic[key] for key in harmonic_fields}
            writer.writerow(
                {
                    "強度順位": raw["strength_rank"],
                    "DFTビン": raw["bin_index"],
                    "波数k": raw["signed_k"],
                    "状態内の役割": raw["role"],
                    "全強度比": raw["power_fraction_total"],
                    "直流を除く強度比": raw["power_fraction_non_dc"],
                    "数値誤差床より上": "はい" if raw["above_numerical_floor"] else "いいえ",
                }
            )


def combined_markdown(results: list[dict[str, Any]]) -> str:
    lines = [
        "# N=5・N=40 粒子状態一覧表 v2",
        "",
        "## 結論",
        "",
        "現行の白色零閉鎖生成器は、N=5・N=40のどちらでも、多数の独立粒子状態を生成していない。",
        "得られたのは、全関係成分を必要とする一つの零閉鎖倍音束であり、基底波以外の高調波は含むが、",
        "全体系とは別に閉じる部分状態も、有限位数回帰を完走した粒子状態も0である。",
        "",
        "| N体数 | M | 零閉鎖波状態 | 内部独立状態 | 粒子状態 | 基底波のみ | 倍音束 | 粒子所属成分 | 海・未束縛成分 |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for result in results:
        summary = result["state_summary"]
        lines.append(
            "| "
            + " | ".join(
                [
                    str(result["n_body"]),
                    str(result["m_relations"]),
                    str(summary["zero_closed_wave_state_count"]),
                    str(summary["independently_closed_internal_state_count"]),
                    str(summary["finite_order_particle_state_count"]),
                    str(summary["base_only_state_count"]),
                    str(summary["harmonic_bundle_state_count"]),
                    str(summary["particle_assigned_relation_components"]),
                    str(summary["sea_or_unbound_relation_components"]),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## 波形と粒子読出し量",
            "",
            "| N体数 | 波形構成 | 基底波 | 同方向高調波数 | 逆回転枝数 | 基底波強度 | 実効倍音数 | 99%成分数 | 奇数倍音比 | 確立済み質量型 R_read | 正負枝Gram診断（探索量） | 住所電荷候補 | 1周期忠実度 | 半周期重なりRe | スピン応答ρ | 総合判定 |",
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        ]
    )
    for result in results:
        state = result["states"][0]
        lines.append(
            "| "
            + " | ".join(
                [
                    str(result["n_body"]),
                    state["waveform_composition"],
                    f"k={state['base_k']:+d} ({state['base_address']})",
                    str(state["same_direction_harmonic_count"]),
                    str(state["reverse_branch_count"]),
                    f"{100*state['base_power_fraction_non_dc']:.3f}%",
                    fmt(state["effective_harmonic_count_non_dc"], 6),
                    str(state["modes_for_99_percent_non_dc"]),
                    f"{100*state['odd_power_fraction_non_dc']:.3f}%",
                    fmt(state["mass_readout_R_read"], 8),
                    fmt(state["exploratory_frequency_branch_gram_det"], 8),
                    fmt(state["charge_address_magnitude"], 9),
                    fmt(state["period_return"].get("fidelity"), 7),
                    fmt(state["half_period_return"].get("overlap_real"), 7),
                    fmt(state["spin_response_strength"], 7),
                    state["status"],
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "確立済み質量型 R_read は、必要なA/B/C多ゲージ観測契約が入力にないため計算不能である。",
            "正負枝Gram値は質量ではなく探索診断として分離した。住所電荷欄は、基底住所から計算できるモデル内の候補値である。両条件とも有限位数回帰を",
            "完走していないため、電荷として確定していない。詳しい根拠・倍音内訳・外部軸スピン応答は",
            "各Nの `census.md` に分けて表示する。",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    if args.output.exists():
        raise SystemExit(f"出力先が既に存在します。上書きしません: {args.output}")
    if args.top_harmonics < 1:
        raise SystemExit("--top-harmonics は1以上が必要です")
    args.output.mkdir(parents=True, exist_ok=False)

    results = []
    seen_n = set()
    for input_dir in args.input:
        result = analyze_run(input_dir.resolve(), args.top_harmonics)
        n_body = result["n_body"]
        if n_body in seen_n:
            raise SystemExit(f"同じNが複数指定されています: N={n_body}")
        seen_n.add(n_body)
        results.append(result)
        write_run(result, args.output / f"N{n_body}", args.top_harmonics)
        summary = result["state_summary"]
        print(
            f"N={n_body} M={result['m_relations']} "
            f"wave_states={summary['zero_closed_wave_state_count']} "
            f"internal_states={summary['independently_closed_internal_state_count']} "
            f"particle_states={summary['finite_order_particle_state_count']}"
        )

    results.sort(key=lambda item: item["n_body"])
    (args.output / "summary.md").write_text(
        combined_markdown(results) + "\n", encoding="utf-8"
    )
    combined = {
        "schema": CENSUS_VERSION,
        "reader": Path(__file__).name,
        "reader_sha256": sha256_file(Path(__file__)),
        "results": results,
    }
    (args.output / "summary.json").write_text(
        json.dumps(combined, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    rows = [state_csv_row(result) for result in results]
    with (args.output / "summary_ja.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"saved: {args.output.resolve()}")


if __name__ == "__main__":
    main()
