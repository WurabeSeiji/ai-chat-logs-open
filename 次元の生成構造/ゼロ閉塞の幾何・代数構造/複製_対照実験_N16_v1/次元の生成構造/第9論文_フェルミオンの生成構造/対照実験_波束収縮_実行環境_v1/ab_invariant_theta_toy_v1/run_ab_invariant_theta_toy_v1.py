#!/usr/bin/env python3
"""AB 合成波の回転不変量から散乱角を生成する派生トイモデル v1。

外部から R または theta を受け取らず、現在の AB 複素波配列だけから
回転角を毎衝突時に再計算する。元のミラーコードは変更せず、状態生成と
既存指標の読出しにだけ利用する。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
ENV_DIR = HERE.parent
SYSTEM_A_PATH = ENV_DIR / "20260715" / "run_system_A_localization_exchange_R_sweep_preliminary_v1.py"
ENGINE_PATH = ENV_DIR / "20260713" / "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py"
DEFAULT_OUTPUT_DIR = HERE / "result_v1"
NUMERICAL_TOLERANCE = 1.0e-10


def load_local_system_a() -> Any:
    """隔離環境内の系統 A コピーをロードする。リポジトリ外には依存しない。"""

    spec = importlib.util.spec_from_file_location("local_system_a_for_ab_theta_toy_v1", SYSTEM_A_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load local System A module: {SYSTEM_A_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


base = load_local_system_a()
src = base.src


@dataclass(frozen=True)
class ThetaReadout:
    theta: float
    reflection_rate: float
    transmission_rate: float
    fermionic_relation_power: float
    bosonic_relation_power: float
    total_pair_power: float


@dataclass(frozen=True)
class ToyCase:
    name: str
    interpretation: str
    packet_a: tuple[int, ...]
    packet_b: tuple[int, ...]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_cases(high_n: int) -> tuple[ToyCase, ...]:
    if high_n < 3 or high_n % 2 == 0:
        raise ValueError("--high-n must be an odd integer of at least 3")
    odd_packet = tuple(range(1, high_n + 1, 2))
    even_packet = (1, *tuple(range(2, high_n, 2)))
    return (
        ToyCase(
            name="fundamental_control",
            interpretation="A と B がともに基本波だけのボゾン的対照",
            packet_a=(1,),
            packet_b=(1,),
        ),
        ToyCase(
            name=f"even_boson_control_B{high_n - 1}",
            interpretation="B に基本波と等振幅の偶数倍音を与えたボゾン的対照",
            packet_a=(1,),
            packet_b=even_packet,
        ),
        ToyCase(
            name=f"odd_fermion_candidate_B{high_n}",
            interpretation="B に基本波から high_n までの等振幅奇数倍音を与えたフェルミオン候補",
            packet_a=(1,),
            packet_b=odd_packet,
        ),
    )


def to_base_case(case: ToyCase) -> Any:
    return base.explicit_packet_case(
        mode=case.name,
        packet_a=case.packet_a,
        packet_b=case.packet_b,
    )


def combined_chi_power(a: np.ndarray, b: np.ndarray, source_params: Any) -> tuple[np.ndarray, np.ndarray]:
    """符号付き chi 周波数ごとの |A_k|^2 + |B_k|^2 を返す。"""

    shape = (source_params.chi_grid_n, source_params.eta_grid_n)
    a_fft = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
    b_fft = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
    pair_power = np.sum(np.abs(a_fft) ** 2 + np.abs(b_fft) ** 2, axis=1)
    frequencies = np.rint(
        np.fft.fftfreq(source_params.chi_grid_n, d=1.0 / source_params.chi_grid_n)
    ).astype(int)
    return frequencies, pair_power


def theta_from_ab(a: np.ndarray, b: np.ndarray, source_params: Any) -> ThetaReadout:
    """AB 合成スペクトルの保存セクターだけから theta を生成する。

    現行状態には q_A=+1, q_B=-1 の搬送波が入るため、内在的な奇数倍音は
    生の FFT では偶数ビンへ移る。基本波が占める |k| <= 2 を除き、
    偶数かつ |k| >= 4 の AB 合成パワーをフェルミオン関係量 P_f とする。
    残りを P_b とし、theta = atan2(sqrt(P_f), sqrt(P_b)) と定める。
    """

    frequencies, pair_power = combined_chi_power(a, b, source_params)
    abs_frequency = np.abs(frequencies)
    fermionic_mask = (abs_frequency >= 4) & ((abs_frequency % 2) == 0)
    total_power = float(np.sum(pair_power))
    fermionic_power = float(np.sum(pair_power[fermionic_mask]))
    numerical_floor = 1024.0 * np.finfo(float).eps * max(total_power, 1.0)
    if abs(fermionic_power) <= numerical_floor:
        fermionic_power = 0.0
    fermionic_power = min(max(fermionic_power, 0.0), total_power)
    bosonic_power = max(total_power - fermionic_power, 0.0)
    if total_power <= 0.0:
        raise ValueError("AB pair has zero total power")
    theta = math.atan2(math.sqrt(fermionic_power), math.sqrt(bosonic_power))
    reflection_rate = math.sin(theta) ** 2
    return ThetaReadout(
        theta=theta,
        reflection_rate=reflection_rate,
        transmission_rate=math.cos(theta) ** 2,
        fermionic_relation_power=fermionic_power,
        bosonic_relation_power=bosonic_power,
        total_pair_power=total_power,
    )


def rotate_ab(a: np.ndarray, b: np.ndarray, theta: float) -> tuple[np.ndarray, np.ndarray]:
    """実直交 2x2 回転を AB チャネルへ作用させる。

    個別チャネルの再正規化は行わない。これにより AB 合成保存量を保つ。
    """

    cosine = math.cos(theta)
    sine = math.sin(theta)
    return cosine * a - sine * b, sine * a + cosine * b


def pair_hermitian_norm(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.vdot(a, a).real + np.vdot(b, b).real)


def pair_zero_closure(a: np.ndarray, b: np.ndarray) -> complex:
    """非共役二乗和 sum(A_i^2) + sum(B_i^2) の現在値。"""

    return complex(np.sum(a * a) + np.sum(b * b))


def state_metrics(vector: np.ndarray, metric_context: Any) -> dict[str, float]:
    distribution = metric_context.harmonic_distribution(vector)
    n_eff, n_eff_2 = src.effective_n(distribution)
    return {
        "L": src.localization(vector),
        "N_eff": n_eff,
        "N_eff_2": n_eff_2,
    }


def run_case(
    case: ToyCase,
    source_params: Any,
    metric_context: Any,
    max_collision: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    base_case = to_base_case(case)
    a = base.make_case_state(source_params, base_case, "A", hair_enabled=True)
    b = base.make_case_state(source_params, base_case, "B", hair_enabled=True)
    initial_a = a.copy()
    initial_b = b.copy()
    initial_h_a = metric_context.harmonic_distribution(initial_a)
    initial_h_b = metric_context.harmonic_distribution(initial_b)
    _, initial_spectrum = combined_chi_power(a, b, source_params)
    initial_norm = pair_hermitian_norm(a, b)
    initial_closure = pair_zero_closure(a, b)
    initial_readout = theta_from_ab(a, b, source_params)

    rows: list[dict[str, Any]] = []
    for collision in range(max_collision + 1):
        readout = theta_from_ab(a, b, source_params)
        _, current_spectrum = combined_chi_power(a, b, source_params)
        current_norm = pair_hermitian_norm(a, b)
        current_closure = pair_zero_closure(a, b)
        metrics_a = state_metrics(a, metric_context)
        metrics_b = state_metrics(b, metric_context)
        h_a = metric_context.harmonic_distribution(a)
        h_b = metric_context.harmonic_distribution(b)
        rows.append(
            {
                "case": case.name,
                "collision": collision,
                "theta_generated": readout.theta,
                "R_generated": readout.reflection_rate,
                "T_generated": readout.transmission_rate,
                "P_f_relation": readout.fermionic_relation_power,
                "P_b_relation": readout.bosonic_relation_power,
                "pair_power": readout.total_pair_power,
                "theta_drift": abs(readout.theta - initial_readout.theta),
                "pair_norm": current_norm,
                "pair_norm_drift": abs(current_norm - initial_norm),
                "closure_real": current_closure.real,
                "closure_imag": current_closure.imag,
                "closure_drift_abs": abs(current_closure - initial_closure),
                "combined_spectrum_max_drift": float(np.max(np.abs(current_spectrum - initial_spectrum))),
                "L_A": metrics_a["L"],
                "L_B": metrics_b["L"],
                "N_eff_A": metrics_a["N_eff"],
                "N_eff_B": metrics_b["N_eff"],
                "origin_A_in_A": src.projection_weight(a, initial_a),
                "origin_B_in_A": src.projection_weight(a, initial_b),
                "origin_A_in_B": src.projection_weight(b, initial_a),
                "origin_B_in_B": src.projection_weight(b, initial_b),
                "spectrum_A_to_A0": base.distribution_similarity(h_a, initial_h_a),
                "spectrum_A_to_B0": base.distribution_similarity(h_a, initial_h_b),
                "spectrum_B_to_A0": base.distribution_similarity(h_b, initial_h_a),
                "spectrum_B_to_B0": base.distribution_similarity(h_b, initial_h_b),
            }
        )
        if collision < max_collision:
            # theta は固定値として保持せず、その時点の AB 波から毎回読み直す。
            a, b = rotate_ab(a, b, readout.theta)

    max_theta_drift = max(float(row["theta_drift"]) for row in rows)
    max_norm_drift = max(float(row["pair_norm_drift"]) for row in rows)
    max_closure_drift = max(float(row["closure_drift_abs"]) for row in rows)
    max_spectrum_drift = max(float(row["combined_spectrum_max_drift"]) for row in rows)
    invariant_pass = (
        max_theta_drift <= NUMERICAL_TOLERANCE
        and max_norm_drift <= NUMERICAL_TOLERANCE
        and max_closure_drift <= NUMERICAL_TOLERANCE
        and max_spectrum_drift <= NUMERICAL_TOLERANCE
    )
    summary = {
        "case": asdict(case),
        "theta_initial": initial_readout.theta,
        "R_generated_initial": initial_readout.reflection_rate,
        "T_generated_initial": initial_readout.transmission_rate,
        "P_f_relation_initial": initial_readout.fermionic_relation_power,
        "P_b_relation_initial": initial_readout.bosonic_relation_power,
        "pair_norm_initial": initial_norm,
        "closure_initial": {
            "real": initial_closure.real,
            "imag": initial_closure.imag,
            "abs": abs(initial_closure),
        },
        "max_theta_drift": max_theta_drift,
        "max_pair_norm_drift": max_norm_drift,
        "max_closure_drift_abs": max_closure_drift,
        "max_combined_spectrum_drift": max_spectrum_drift,
        "invariant_tolerance": NUMERICAL_TOLERANCE,
        "invariant_verdict": "PASS" if invariant_pass else "CHECK",
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("no rows to write")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, result: dict[str, Any]) -> None:
    lines = [
        "# AB 回転不変量 theta 自動生成トイモデル v1 実行報告",
        "",
        f"- 衝突回数: {result['parameters']['max_collision']}",
        f"- high_n: {result['parameters']['high_n']}",
        "- 外部 R/theta 入力: なし",
        "- 更新: 実直交 2x2 AB 回転、個別正規化なし",
        "",
        "| case | theta(0) | R(0) | max theta drift | max closure drift | verdict |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for item in result["case_summaries"]:
        lines.append(
            "| {case} | {theta:.12g} | {rate:.12g} | {theta_drift:.3e} | "
            "{closure_drift:.3e} | {verdict} |".format(
                case=item["case"]["name"],
                theta=item["theta_initial"],
                rate=item["R_generated_initial"],
                theta_drift=item["max_theta_drift"],
                closure_drift=item["max_closure_drift_abs"],
                verdict=item["invariant_verdict"],
            )
        )
    lines.extend(
        [
            "",
            "この判定は数値的不変性の確認であり、ボゾン・フェルミオン対応の物理的実証ではない。",
            "また theta が一定になることは同じ U の反復を与えるが、有限 n での U^n=I を自動的には保証しない。",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate the 2x2 scattering angle only from the current AB waves"
    )
    parser.add_argument("--high-n", type=int, default=63, help="odd harmonic endpoint of the initial B packet")
    parser.add_argument("--max-collision", type=int, default=32)
    parser.add_argument(
        "--case",
        action="append",
        help="case name to run; repeatable. Omit to run all built-in cases.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_collision < 0:
        raise ValueError("--max-collision must be non-negative")
    cases = build_cases(args.high_n)
    if args.case:
        requested = set(args.case)
        known = {case.name for case in cases}
        unknown = sorted(requested - known)
        if unknown:
            raise ValueError(f"unknown case(s): {', '.join(unknown)}; known={', '.join(sorted(known))}")
        cases = tuple(case for case in cases if case.name in requested)

    params = base.Params(
        high_n=args.high_n,
        recursive_collision_count=args.max_collision,
    )
    source_params = base.build_source_params(params)
    metric_context = base.MetricContext(source_params)
    all_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for case in cases:
        rows, summary = run_case(case, source_params, metric_context, args.max_collision)
        all_rows.extend(rows)
        summaries.append(summary)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "experiment": "ab_invariant_theta_toy_v1",
        "status_boundary": "toy model; numerical invariance is not a particle-statistics proof",
        "parameters": {
            "high_n": args.high_n,
            "max_collision": args.max_collision,
            "chi_grid_n": source_params.chi_grid_n,
            "eta_grid_n": source_params.eta_grid_n,
            "q_A": source_params.q_A,
            "q_B": source_params.q_B,
            "p0": source_params.p0,
        },
        "theta_rule": {
            "fermionic_sector": "raw chi FFT bins with even |k| >= 4",
            "carrier_note": "q_A=+1 and q_B=-1 shift intrinsic odd harmonics to raw even bins",
            "formula": "theta=atan2(sqrt(P_f),sqrt(P_b)); R=sin(theta)^2",
            "external_scattering_parameter": False,
        },
        "rotation_rule": "[A';B']=[[cos(theta),-sin(theta)],[sin(theta),cos(theta)]] [A;B]",
        "local_dependencies": {
            "system_A_path": str(SYSTEM_A_PATH.relative_to(ENV_DIR)),
            "system_A_sha256": sha256(SYSTEM_A_PATH),
            "engine_path": str(ENGINE_PATH.relative_to(ENV_DIR)),
            "engine_sha256": sha256(ENGINE_PATH),
        },
        "case_summaries": summaries,
    }
    write_csv(output_dir / "ab_invariant_theta_toy_rows_v1.csv", all_rows)
    (output_dir / "ab_invariant_theta_toy_result_v1.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    write_report(output_dir / "ab_invariant_theta_toy_report_v1.md", result)
    print(f"saved: {output_dir}")


if __name__ == "__main__":
    main()
