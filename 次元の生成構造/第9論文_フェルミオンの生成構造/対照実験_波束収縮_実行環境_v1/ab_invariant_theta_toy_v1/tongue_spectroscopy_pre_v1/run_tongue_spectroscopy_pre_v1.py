#!/usr/bin/env python3
"""読出し分類＋タング分光・予備テスト v1

Part A: 読出しの分類（保存 no-go 級 / 動力学級）
    θ 読出しの候補を系統的に並べ、各読出しが「自分自身の駆動する動力学の下で
    保存されるか」を測る。保存なら no-go 級（ロック不能）、動くなら動力学級。
    V1（X パワー）が動き V2（交差スペクトル）が保存された事実の一般化。

Part B: タング分光（V1 読出し）
    振幅掃引を高分解能化し、悪魔の階段上のタングを列挙。各タングの
    有理数・幅・厳密度を測り、幅の分母依存（Arnold スケーリング）を見積もる。

Part C: (23,124) 近傍ズーム
    前回の精密スキャンで、1/3 タング下端のすぐ外（振幅≈2.5）の回転数が
    0.31437 と、θ*/π = 1/2 − 23/124 = 39/124 = 0.314516... の至近だった。
    この帯域を高分解能・長時間でズームし、39/124 タングの存在を探す。
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
TARGET_RHO = Fraction(39, 124)  # 1/2 - 23/124

CLASSIFY_COLLISIONS = 300
CLASSIFY_AMPLITUDES = (0.5, 1.0, 2.0)
CONSERVED_TOLERANCE = 1.0e-12

SCAN_AMPLITUDES = np.geomspace(0.05, 5.0, 161)
SCAN_COLLISIONS = 400
SCAN_TAIL = 150

ZOOM_AMPLITUDES = np.linspace(2.30, 2.62, 49)
ZOOM_COLLISIONS = 900
ZOOM_TAIL = 300

PLATEAU_EPS = 1.0e-6
PLATEAU_MIN_POINTS = 2
RATIONAL_MAX_DEN = 130


def load_toy_module() -> Any:
    spec = importlib.util.spec_from_file_location(
        "ab_theta_toy_for_tongue_spectroscopy_v1", TOY_RUNNER_PATH
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


def ratio_readout(power_fn: Callable[..., np.ndarray]) -> Callable:
    def readout(a: np.ndarray, b: np.ndarray, source_params: Any) -> float:
        frequencies, a_fft, b_fft = spectra(a, b, source_params)
        power = power_fn(a_fft, b_fft)
        mask = canonical_mask(frequencies)
        return theta_from_powers(
            float(np.sum(power[mask])), float(np.sum(power[~mask]))
        )
    return readout


READOUTS: dict[str, Callable] = {
    "R1_diag_powers": ratio_readout(
        lambda af, bf: np.sum(np.abs(af) ** 2 + np.abs(bf) ** 2, axis=1)
    ),
    "R2_X_power_V1": ratio_readout(
        lambda af, bf: np.sum(np.abs(af + bf) ** 2, axis=1)
    ),
    "R3_Y_power": ratio_readout(
        lambda af, bf: np.sum(np.abs(af - bf) ** 2, axis=1)
    ),
    "R4_cross_abs_V2": ratio_readout(
        lambda af, bf: np.abs(np.sum(af * np.conj(bf), axis=1))
    ),
    "R5_A_only": ratio_readout(lambda af, bf: np.sum(np.abs(af) ** 2, axis=1)),
    "R6_B_only": ratio_readout(lambda af, bf: np.sum(np.abs(bf) ** 2, axis=1)),
    "R7_cross_real_abs": ratio_readout(
        lambda af, bf: np.abs(np.real(np.sum(af * np.conj(bf), axis=1)))
    ),
}


def make_templates(source_params: Any) -> tuple[np.ndarray, np.ndarray]:
    case = base.explicit_packet_case(
        mode="tongue_spectroscopy_b63",
        packet_a=(1,),
        packet_b=tuple(range(1, HIGH_N + 1, 2)),
    )
    a = base.make_case_state(source_params, case, "A", hair_enabled=True)
    b = base.make_case_state(source_params, case, "B", hair_enabled=True)
    return a, b


def classify_readout(
    theta_fn: Callable,
    a_template: np.ndarray,
    b_template: np.ndarray,
    source_params: Any,
) -> dict[str, Any]:
    max_drift = 0.0
    for amplitude in CLASSIFY_AMPLITUDES:
        a = a_template.copy()
        b = amplitude * b_template
        theta0 = theta_fn(a, b, source_params)
        for _ in range(CLASSIFY_COLLISIONS):
            theta = theta_fn(a, b, source_params)
            max_drift = max(max_drift, abs(theta - theta0))
            a, b = toy.rotate_ab(a, b, theta)
    return {
        "max_self_drift": max_drift,
        "class": "CONSERVED_NO_GO" if max_drift <= CONSERVED_TOLERANCE else "DYNAMICAL",
    }


def rotation_scan(
    theta_fn: Callable,
    amplitudes: np.ndarray,
    collisions: int,
    tail: int,
    a_template: np.ndarray,
    b_template: np.ndarray,
    source_params: Any,
) -> tuple[np.ndarray, np.ndarray]:
    means = np.zeros(len(amplitudes))
    stds = np.zeros(len(amplitudes))
    for idx, amplitude in enumerate(amplitudes):
        a = a_template.copy()
        b = float(amplitude) * b_template
        history = np.zeros(collisions)
        for j in range(collisions):
            theta = theta_fn(a, b, source_params)
            history[j] = theta / math.pi
            a, b = toy.rotate_ab(a, b, theta)
        means[idx] = float(np.mean(history[-tail:]))
        stds[idx] = float(np.std(history[-tail:]))
    return means, stds


def find_plateaus(amplitudes: np.ndarray, values: np.ndarray) -> list[dict[str, Any]]:
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
                        "n_points": length,
                        "amp_from": float(amplitudes[start]),
                        "amp_to": float(amplitudes[i - 1]),
                        "log_width": float(
                            math.log(amplitudes[i - 1] / amplitudes[start])
                        )
                        if amplitudes[i - 1] > amplitudes[start]
                        else 0.0,
                        "rho": level,
                        "nearest_rational": f"{frac.numerator}/{frac.denominator}",
                        "denominator": frac.denominator,
                        "distance": abs(level - float(frac)),
                        "exact": abs(level - float(frac)) <= 1.0e-9,
                    }
                )
            start = i
    return plateaus


def main() -> None:
    params = base.Params(high_n=HIGH_N, recursive_collision_count=SCAN_COLLISIONS)
    source_params = base.build_source_params(params)
    a_template, b_template = make_templates(source_params)

    # ---- Part A: 読出し分類 ----
    classification: dict[str, Any] = {}
    print("Part A: readout classification")
    for name, theta_fn in READOUTS.items():
        result = classify_readout(theta_fn, a_template, b_template, source_params)
        classification[name] = result
        print(f"  {name:<20} {result['class']:<16} drift={result['max_self_drift']:.2e}")

    # ---- Part B: タング分光（V1 = R2）----
    v1 = READOUTS["R2_X_power_V1"]
    means, stds = rotation_scan(
        v1, SCAN_AMPLITUDES, SCAN_COLLISIONS, SCAN_TAIL,
        a_template, b_template, source_params,
    )
    plateaus = find_plateaus(SCAN_AMPLITUDES, means)
    exact_tongues = [p for p in plateaus if p["exact"]]
    print(f"Part B: {len(plateaus)} plateaus, {len(exact_tongues)} exact rational tongues")
    for p in exact_tongues:
        print(
            f"  rho={p['rho']:.9f} = {p['nearest_rational']}"
            f" width(log)={p['log_width']:.4f} pts={p['n_points']}"
            f" amp {p['amp_from']:.3f}..{p['amp_to']:.3f}"
        )

    # ---- Part C: 39/124 近傍ズーム ----
    zoom_means, zoom_stds = rotation_scan(
        v1, ZOOM_AMPLITUDES, ZOOM_COLLISIONS, ZOOM_TAIL,
        a_template, b_template, source_params,
    )
    target = float(TARGET_RHO)
    zoom_rows = []
    best_idx = int(np.argmin(np.abs(zoom_means - target)))
    for idx, amp in enumerate(ZOOM_AMPLITUDES):
        frac = Fraction(zoom_means[idx]).limit_denominator(RATIONAL_MAX_DEN)
        zoom_rows.append(
            {
                "amplitude": float(amp),
                "rho": float(zoom_means[idx]),
                "tail_std": float(zoom_stds[idx]),
                "dist_to_39_124": float(abs(zoom_means[idx] - target)),
                "nearest_rational": f"{frac.numerator}/{frac.denominator}",
                "dist_to_nearest": float(abs(zoom_means[idx] - float(frac))),
            }
        )
    zoom_plateaus = find_plateaus(ZOOM_AMPLITUDES, zoom_means)
    print(f"Part C: zoom around 39/124={target:.9f}")
    print(
        f"  closest point: amp={ZOOM_AMPLITUDES[best_idx]:.4f}"
        f" rho={zoom_means[best_idx]:.9f}"
        f" dist={abs(zoom_means[best_idx]-target):.2e}"
    )
    for p in zoom_plateaus:
        print(
            f"  zoom plateau rho={p['rho']:.9f} ~ {p['nearest_rational']}"
            f" (d={p['distance']:.1e}, pts={p['n_points']},"
            f" amp {p['amp_from']:.4f}..{p['amp_to']:.4f})"
        )

    # ---- 図 ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    ax1.plot(SCAN_AMPLITUDES, means, ".-", markersize=2.5, linewidth=0.7)
    for p in exact_tongues:
        ax1.axhline(p["rho"], color="tab:red", linewidth=0.5, alpha=0.5)
    ax1.set_xscale("log")
    ax1.set_xlabel("initial B amplitude")
    ax1.set_ylabel("tail-mean theta/pi")
    ax1.set_title(f"Devil's staircase (V1), {len(exact_tongues)} exact tongues")
    ax1.grid(alpha=0.3)
    ax2.plot(ZOOM_AMPLITUDES, zoom_means, ".-", markersize=3, linewidth=0.8)
    ax2.axhline(target, color="tab:red", linestyle=":", label="39/124 (target)")
    ax2.axhline(1.0 / 3.0, color="tab:orange", linestyle=":", label="1/3")
    ax2.set_xlabel("initial B amplitude")
    ax2.set_ylabel("tail-mean theta/pi")
    ax2.set_title("Zoom near the 1/3 tongue lower edge")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    figure_names = []
    for ext in ("png", "svg"):
        path = HERE / f"tongue_spectroscopy_pre_v1.{ext}"
        fig.savefig(path, dpi=160)
        figure_names.append(path.name)
    plt.close(fig)

    csv_path = HERE / "tongue_spectroscopy_zoom_rows_v1.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(zoom_rows[0]))
        writer.writeheader()
        writer.writerows(zoom_rows)

    payload = {
        "experiment": "tongue_spectroscopy_pre_v1",
        "core_runner": {
            "path": TOY_RUNNER_PATH.name,
            "sha256": toy.sha256(TOY_RUNNER_PATH),
        },
        "part_A_readout_classification": classification,
        "part_B_staircase": {
            "n_amplitudes": len(SCAN_AMPLITUDES),
            "collisions": SCAN_COLLISIONS,
            "plateaus": plateaus,
            "exact_tongues": exact_tongues,
        },
        "part_C_zoom_39_124": {
            "target_rho": target,
            "amplitude_range": [float(ZOOM_AMPLITUDES[0]), float(ZOOM_AMPLITUDES[-1])],
            "collisions": ZOOM_COLLISIONS,
            "closest_point": zoom_rows[best_idx],
            "zoom_plateaus": zoom_plateaus,
            "rows": zoom_rows,
        },
        "figures": figure_names,
    }
    (HERE / "tongue_spectroscopy_pre_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
