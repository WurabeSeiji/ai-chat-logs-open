#!/usr/bin/env python3
"""モードロック探査・予備テスト v1

三段の手順（逆二乗論文 v1→v2 と同じ規律）:

  [零次 no-go（定理＋数値）]
    現行読出し theta_from_ab はビン別パワー |A_k|^2+|B_k|^2 だけの関数であり、
    実直交回転 rotate_ab はこれを各ビンで厳密保存する。ゆえに theta は
    運動の定数で、theta の動力学（したがってモードロック）は存在しない。
    数値検証: 全振幅で 400 衝突のドリフト = 機械精度 0。

  [保存的作業仮説]
    theta を位相感受的な読出しへ差し替える2変種を試す。挿入ではなく、
    N体系列がもともと基礎に置く「関係波 f_AB」からの読出しという
    フレームワーク内在の候補:
      V1 coherent : 対称チャネル X = A+B のマスク内パワー比から theta
      V2 relational: 交差スペクトル A_k B_k^* のマスク内比から theta
    どちらも回転不変性を持たないため theta_j が発展する（円写像型力学）。

  [独立検証＝本テストの測定]
    初期B振幅（旧来の連続ノブ）を掃引し、後期 theta/pi の漸近値を測る。
    モードロックが在るなら「悪魔の階段」: theta/pi が小分母有理数の
    プラトーに離散的に張り付く。無ければ滑らかに変化する。
    判定はプラトー検出（連続する振幅点で theta/pi が一致）と、
    プラトー値の有理数同定で行う。

注意:
    これは「どの有理数か」(23/124) 以前の、「ロック機構が存在するか」の
    存在テストである。プラトーが一つでも出れば、電荷量子化＝Arnold tongue
    台地の読みがモデル内で初めて実体を持つ。
"""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable

import numpy as np


HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"

HIGH_N = 63
COLLISIONS = 400
TAIL = 150
AMPLITUDES = np.geomspace(0.05, 5.0, 61)
NO_GO_AMPLITUDES = (0.1, 0.5, 1.0, 2.0)
PLATEAU_EPS = 1.0e-6
PLATEAU_MIN_POINTS = 3
RATIONAL_MAX_DEN = 200


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ab_theta_toy_for_mode_locking_probe_v1", TOY_RUNNER_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load toy runner: {TOY_RUNNER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


toy = load_toy_module()
base = toy.base
plt = base.plt


def spectra(a: np.ndarray, b: np.ndarray, source_params: Any):
    shape = (source_params.chi_grid_n, source_params.eta_grid_n)
    a_fft = np.fft.fft(a.reshape(shape), axis=0, norm="ortho")
    b_fft = np.fft.fft(b.reshape(shape), axis=0, norm="ortho")
    frequencies = np.rint(
        np.fft.fftfreq(source_params.chi_grid_n, d=1.0 / source_params.chi_grid_n)
    ).astype(int)
    return frequencies, a_fft, b_fft


def canonical_mask(frequencies: np.ndarray) -> np.ndarray:
    f = np.abs(frequencies)
    return (f >= 4) & (f % 2 == 0)


def theta_from_powers(p_f: float, p_b: float) -> float:
    p_f = max(p_f, 0.0)
    p_b = max(p_b, 0.0)
    if p_f + p_b <= 0.0:
        return 0.0
    return math.atan2(math.sqrt(p_f), math.sqrt(p_b))


def theta_coherent(a: np.ndarray, b: np.ndarray, source_params: Any) -> float:
    """V1: 対称チャネル X=A+B のパワーだけで theta（位相感受的）。"""
    frequencies, a_fft, b_fft = spectra(a, b, source_params)
    x_power = np.sum(np.abs(a_fft + b_fft) ** 2, axis=1)
    mask = canonical_mask(frequencies)
    return theta_from_powers(
        float(np.sum(x_power[mask])), float(np.sum(x_power[~mask]))
    )


def theta_relational(a: np.ndarray, b: np.ndarray, source_params: Any) -> float:
    """V2: 関係波（交差スペクトル A_k B_k^*）の大きさで theta。"""
    frequencies, a_fft, b_fft = spectra(a, b, source_params)
    cross = np.abs(np.sum(a_fft * np.conj(b_fft), axis=1))
    mask = canonical_mask(frequencies)
    return theta_from_powers(float(np.sum(cross[mask])), float(np.sum(cross[~mask])))


VARIANTS: dict[str, Callable[[np.ndarray, np.ndarray, Any], float]] = {
    "V1_coherent_X_power": theta_coherent,
    "V2_relational_cross": theta_relational,
}


def make_templates(source_params: Any) -> tuple[np.ndarray, np.ndarray]:
    case = base.explicit_packet_case(
        mode="mode_locking_b63",
        packet_a=(1,),
        packet_b=tuple(range(1, HIGH_N + 1, 2)),
    )
    a = base.make_case_state(source_params, case, "A", hair_enabled=True)
    b = base.make_case_state(source_params, case, "B", hair_enabled=True)
    return a, b


def no_go_check(
    a_template: np.ndarray,
    b_template: np.ndarray,
    source_params: Any,
) -> dict[str, Any]:
    """零次: 現行読出しでは theta が厳密に不変（数値確認）。"""
    max_drift = 0.0
    for amplitude in NO_GO_AMPLITUDES:
        a = a_template.copy()
        b = amplitude * b_template
        theta0 = toy.theta_from_ab(a, b, source_params).theta
        for _ in range(COLLISIONS):
            readout = toy.theta_from_ab(a, b, source_params)
            max_drift = max(max_drift, abs(readout.theta - theta0))
            a, b = toy.rotate_ab(a, b, readout.theta)
    return {
        "collisions": COLLISIONS,
        "amplitudes": list(NO_GO_AMPLITUDES),
        "max_theta_drift": max_drift,
        "verdict": "NO_GO_CONFIRMED" if max_drift <= 1.0e-12 else "UNEXPECTED_DRIFT",
    }


def run_variant(
    label: str,
    theta_fn: Callable[[np.ndarray, np.ndarray, Any], float],
    a_template: np.ndarray,
    b_template: np.ndarray,
    source_params: Any,
) -> tuple[list[dict[str, Any]], np.ndarray, np.ndarray]:
    rows: list[dict[str, Any]] = []
    final_values = np.zeros(len(AMPLITUDES))
    final_stds = np.zeros(len(AMPLITUDES))
    for idx, amplitude in enumerate(AMPLITUDES):
        a = a_template.copy()
        b = float(amplitude) * b_template
        history = np.zeros(COLLISIONS)
        for j in range(COLLISIONS):
            theta = theta_fn(a, b, source_params)
            history[j] = theta / math.pi
            a, b = toy.rotate_ab(a, b, theta)
        tail = history[-TAIL:]
        final_values[idx] = float(np.mean(tail))
        final_stds[idx] = float(np.std(tail))
        rows.append(
            {
                "variant": label,
                "amplitude": float(amplitude),
                "theta_over_pi_tail_mean": final_values[idx],
                "theta_over_pi_tail_std": final_stds[idx],
                "theta_over_pi_initial": float(history[0]),
            }
        )
    return rows, final_values, final_stds


def find_plateaus(values: np.ndarray) -> list[dict[str, Any]]:
    plateaus: list[dict[str, Any]] = []
    start = 0
    for i in range(1, len(values) + 1):
        if i == len(values) or abs(values[i] - values[start]) > PLATEAU_EPS:
            length = i - start
            if length >= PLATEAU_MIN_POINTS:
                level = float(np.mean(values[start:i]))
                frac = Fraction(level).limit_denominator(RATIONAL_MAX_DEN)
                plateaus.append(
                    {
                        "start_index": start,
                        "n_points": length,
                        "amplitude_from": float(AMPLITUDES[start]),
                        "amplitude_to": float(AMPLITUDES[i - 1]),
                        "theta_over_pi": level,
                        "nearest_rational": f"{frac.numerator}/{frac.denominator}",
                        "distance_to_rational": abs(level - float(frac)),
                        "is_rational_lock": abs(level - float(frac)) <= 1.0e-9,
                    }
                )
            start = i
    return plateaus


def main() -> None:
    params = base.Params(high_n=HIGH_N, recursive_collision_count=COLLISIONS)
    source_params = base.build_source_params(params)
    a_template, b_template = make_templates(source_params)

    no_go = no_go_check(a_template, b_template, source_params)
    print(f"zero-order no-go: {no_go['verdict']} (max drift {no_go['max_theta_drift']:.2e})")

    all_rows: list[dict[str, Any]] = []
    variant_summary: dict[str, Any] = {}
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    for ax, (label, theta_fn) in zip(axes, VARIANTS.items()):
        rows, finals, stds = run_variant(
            label, theta_fn, a_template, b_template, source_params
        )
        all_rows.extend(rows)
        plateaus = find_plateaus(finals)
        converged = float(np.median(stds))
        variant_summary[label] = {
            "n_plateaus": len(plateaus),
            "plateaus": plateaus,
            "median_tail_std": converged,
            "range_theta_over_pi": [float(np.min(finals)), float(np.max(finals))],
        }
        ax.plot(AMPLITUDES, finals, ".-", markersize=3, linewidth=0.8)
        for plateau in plateaus:
            ax.axhspan(
                plateau["theta_over_pi"] - 2e-4,
                plateau["theta_over_pi"] + 2e-4,
                color="tab:red",
                alpha=0.25,
            )
        ax.set_xscale("log")
        ax.set_xlabel("initial B amplitude")
        ax.set_ylabel("tail-mean theta/pi")
        ax.set_title(f"{label}: {len(plateaus)} plateau(s)")
        ax.grid(alpha=0.3)
        print(f"{label}: plateaus={len(plateaus)} median_tail_std={converged:.2e}")
        for plateau in plateaus:
            print(
                f"   theta/pi={plateau['theta_over_pi']:.9f}"
                f" ~ {plateau['nearest_rational']}"
                f" (d={plateau['distance_to_rational']:.1e},"
                f" {plateau['n_points']} pts,"
                f" amp {plateau['amplitude_from']:.3f}..{plateau['amplitude_to']:.3f})"
            )
    fig.suptitle(
        "Mode-locking probe: devil's staircase search over initial-B-amplitude sweep"
    )
    figure_names = []
    for ext in ("png", "svg"):
        path = HERE / f"mode_locking_probe_pre_v1.{ext}"
        fig.savefig(path, dpi=160)
        figure_names.append(path.name)
    plt.close(fig)

    csv_path = HERE / "mode_locking_probe_pre_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)

    payload = {
        "experiment": "mode_locking_probe_pre_v1",
        "procedure": [
            "zero-order no-go (theorem + numeric)",
            "conservative working hypothesis: phase-sensitive theta readouts (coherent X power / relational cross spectrum)",
            "independent verification: devil's-staircase search over amplitude sweep",
        ],
        "conditions": {
            "collisions": COLLISIONS,
            "tail_window": TAIL,
            "n_amplitudes": len(AMPLITUDES),
            "amplitude_range": [float(AMPLITUDES[0]), float(AMPLITUDES[-1])],
            "plateau_eps": PLATEAU_EPS,
            "plateau_min_points": PLATEAU_MIN_POINTS,
        },
        "core_runner": {
            "path": TOY_RUNNER_PATH.name,
            "sha256": toy.sha256(TOY_RUNNER_PATH),
        },
        "zero_order_no_go": no_go,
        "variants": variant_summary,
        "figures": figure_names,
    }
    (HERE / "mode_locking_probe_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
