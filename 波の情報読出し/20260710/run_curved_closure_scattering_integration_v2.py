from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np

from run_curved_closure_stationary_wave_broad_sweep_v2 import (
    base_closed_pairs,
    beta_for_correction,
    conformal_factor,
    correction_names,
    evaluate_case,
    mode_coordinate,
    phase_model_names,
    phase_shape,
)


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "curved_closure_scattering_integration_result_v2"
OUT_DIR.mkdir(exist_ok=True)

MPL_DIR = OUT_DIR / ".matplotlib"
MPL_DIR.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class Params:
    mode_count: int = 50
    delta_max: float = 1.2
    grid_n: int = 2048
    domain_length: float = 200.0
    packet_center: float = -40.0
    packet_k0: float = 5.0
    packet_sigma: float = 2.0
    window_inner: float = 25.0
    window_outer: float = 35.0
    scattering_error_tol: float = 1e-8
    dynamic_expected_tol: float = 1e-8
    uncorrected_min_transmission: float = 5e-2


def normalize(psi: np.ndarray, dx: float) -> np.ndarray:
    norm = math.sqrt(float(np.sum(np.abs(psi) ** 2) * dx))
    return psi / norm


def gaussian_packet(rho: np.ndarray, center: float, k0: float, sigma: float) -> np.ndarray:
    envelope = np.exp(-((rho - center) ** 2) / (4.0 * sigma**2))
    return envelope * np.exp(1j * k0 * rho)


def evolve_free(psi0: np.ndarray, t: float, dx: float) -> np.ndarray:
    n = psi0.size
    k = 2.0 * math.pi * np.fft.fftfreq(n, d=dx)
    ordered = np.fft.ifftshift(psi0)
    evolved = np.fft.ifft(np.fft.fft(ordered) * np.exp(-0.5j * (k**2) * t))
    return np.fft.fftshift(evolved)


def parity_partner(rho: np.ndarray, psi: np.ndarray) -> np.ndarray:
    partner_rho = -rho
    real = np.interp(partner_rho, rho, psi.real, left=0.0, right=0.0)
    imag = np.interp(partner_rho, rho, psi.imag, left=0.0, right=0.0)
    return real + 1j * imag


def interaction_window(rho: np.ndarray, inner: float, outer: float) -> np.ndarray:
    distance = np.abs(rho)
    window = np.zeros_like(rho)
    window[distance <= inner] = 1.0
    edge = (distance > inner) & (distance < outer)
    window[edge] = 0.5 * (1.0 + np.cos(math.pi * (distance[edge] - inner) / (outer - inner)))
    return window


def apply_even_odd_phase(
    rho: np.ndarray,
    psi: np.ndarray,
    delta_f: float,
    inner: float,
    outer: float,
) -> np.ndarray:
    partner = parity_partner(rho, psi)
    local_delta = delta_f * interaction_window(rho, inner, outer)
    t = np.exp(0.5j * local_delta) * np.cos(0.5 * local_delta)
    r = -1j * np.exp(0.5j * local_delta) * np.sin(0.5 * local_delta)
    return t * psi + r * partner


def side_probabilities(rho: np.ndarray, dx: float, psi: np.ndarray) -> Dict[str, float]:
    density = np.abs(psi) ** 2
    left = rho < 0.0
    right = rho > 0.0
    left_prob = float(np.sum(density[left]) * dx)
    right_prob = float(np.sum(density[right]) * dx)
    return {
        "reflection_rate": left_prob,
        "transmission_rate": right_prob,
        "norm_total": left_prob + right_prob,
    }


def build_collision_state(params: Params) -> tuple[np.ndarray, float, float, np.ndarray, float]:
    dx = params.domain_length / params.grid_n
    rho = (np.arange(params.grid_n) - params.grid_n // 2) * dx
    hit_time = -params.packet_center / params.packet_k0
    initial = gaussian_packet(rho, params.packet_center, params.packet_k0, params.packet_sigma)
    initial = normalize(initial, dx)
    collision = evolve_free(initial, hit_time, dx)
    return rho, dx, hit_time, collision, hit_time


def weighted_dynamic_rates(
    params: Params,
    rho: np.ndarray,
    dx: float,
    hit_time: float,
    collision: np.ndarray,
    residual_delta: np.ndarray,
    weights: np.ndarray,
) -> Dict[str, float]:
    weights = weights / float(np.sum(weights))
    reflection = 0.0
    transmission = 0.0
    norm_total = 0.0
    for delta, weight in zip(residual_delta, weights):
        after = apply_even_odd_phase(
            rho,
            collision,
            math.pi + float(delta),
            params.window_inner,
            params.window_outer,
        )
        final = evolve_free(after, hit_time, dx)
        metrics = side_probabilities(rho, dx, final)
        reflection += float(weight) * metrics["reflection_rate"]
        transmission += float(weight) * metrics["transmission_rate"]
        norm_total += float(weight) * metrics["norm_total"]
    return {
        "dynamic_reflection_rate": float(reflection),
        "dynamic_transmission_rate": float(transmission),
        "dynamic_norm_total": float(norm_total),
    }


def expected_rates(residual_delta: np.ndarray, weights: np.ndarray) -> Dict[str, float]:
    weights = weights / float(np.sum(weights))
    transmission = np.sin(residual_delta / 2.0) ** 2
    reflection = np.cos(residual_delta / 2.0) ** 2
    return {
        "expected_reflection_rate": float(np.sum(weights * reflection)),
        "expected_transmission_rate": float(np.sum(weights * transmission)),
    }


def case_residual_and_weights(params: Params, phase_model: str, correction: str) -> tuple[np.ndarray, np.ndarray]:
    u = mode_coordinate(params)
    x0, y0 = base_closed_pairs(params)
    g = conformal_factor(params.delta_max, u, params)
    weights = np.abs(g * x0) ** 2 + np.abs(g * y0) ** 2
    delta = params.delta_max * phase_shape(phase_model, u)
    beta = beta_for_correction(delta, correction, u, weights)
    return delta + beta, weights


def run_cases(params: Params) -> List[Dict[str, float | str]]:
    rho, dx, hit_time, collision, final_time = build_collision_state(params)
    rows: List[Dict[str, float | str]] = []
    for phase_model in phase_model_names():
        for correction in correction_names():
            residual_delta, weights = case_residual_and_weights(params, phase_model, correction)
            dynamic = weighted_dynamic_rates(params, rho, dx, hit_time, collision, residual_delta, weights)
            expected = expected_rates(residual_delta, weights)
            closure_case = evaluate_case(params, phase_model, correction, params.delta_max)
            rows.append(
                {
                    "phase_model": phase_model,
                    "correction": correction,
                    "delta_max": params.delta_max,
                    "closure_pair_rms": float(closure_case["closure_pair_rms"]),
                    "expected_reflection_rate": expected["expected_reflection_rate"],
                    "expected_transmission_rate": expected["expected_transmission_rate"],
                    **dynamic,
                    "reflection_abs_error": abs(
                        dynamic["dynamic_reflection_rate"] - expected["expected_reflection_rate"]
                    ),
                    "transmission_abs_error": abs(
                        dynamic["dynamic_transmission_rate"] - expected["expected_transmission_rate"]
                    ),
                    "norm_abs_error": abs(dynamic["dynamic_norm_total"] - 1.0),
                    "hit_time": hit_time,
                    "final_time": 2.0 * final_time,
                }
            )
    return rows


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def aggregate_by_correction(rows: List[Dict[str, float | str]]) -> Dict[str, Dict[str, float]]:
    aggregate: Dict[str, Dict[str, float]] = {}
    for correction in correction_names():
        selected = [row for row in rows if row["correction"] == correction]
        aggregate[correction] = {
            "max_dynamic_transmission_rate": float(
                max(float(row["dynamic_transmission_rate"]) for row in selected)
            ),
            "max_expected_transmission_rate": float(
                max(float(row["expected_transmission_rate"]) for row in selected)
            ),
            "max_closure_pair_rms": float(max(float(row["closure_pair_rms"]) for row in selected)),
            "max_transmission_abs_error": float(
                max(float(row["transmission_abs_error"]) for row in selected)
            ),
            "max_norm_abs_error": float(max(float(row["norm_abs_error"]) for row in selected)),
        }
    return aggregate


def write_plot(rows: List[Dict[str, float | str]], path: Path) -> None:
    corrections = correction_names()
    xs = np.arange(len(corrections))
    max_dynamic = []
    max_expected = []
    max_closure = []
    for correction in corrections:
        selected = [row for row in rows if row["correction"] == correction]
        max_dynamic.append(max(float(row["dynamic_transmission_rate"]) for row in selected))
        max_expected.append(max(float(row["expected_transmission_rate"]) for row in selected))
        max_closure.append(max(float(row["closure_pair_rms"]) for row in selected))

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axes[0].plot(xs, max_dynamic, marker="o", label="dynamic T leakage")
    axes[0].plot(xs, max_expected, linestyle="--", marker="x", label="expected T leakage")
    axes[0].set_ylabel("max transmission leakage")
    axes[0].set_yscale("symlog", linthresh=1e-18)
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[1].plot(xs, max_closure, marker="o", color="tab:red")
    axes[1].set_ylabel("max closure pair RMS")
    axes[1].set_yscale("symlog", linthresh=1e-18)
    axes[1].set_xticks(xs, corrections, rotation=35, ha="right")
    axes[1].grid(True, alpha=0.3)
    fig.suptitle("Curvature residuals in single-sided exchange scattering")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def build_report(result: Dict[str, object]) -> str:
    verdict = result["verdict"]
    summary = result["summary"]
    outputs = result["outputs"]
    aggregate = summary["aggregate_by_correction"]
    rows = "\n".join(
        f"| {name} | `{data['max_dynamic_transmission_rate']:.16e}` | `{data['max_closure_pair_rms']:.16e}` | `{data['max_transmission_abs_error']:.16e}` |"
        for name, data in aggregate.items()
    )
    return f"""# 曲率付き閉鎖定常波 片側入射散乱統合検証 v2

## 目的

閉鎖定常波レベルで得た曲率相対位相の残差を、片側入射の局所交換干渉写像へ戻し、無補正では通過漏れが出ること、完全な内部位相再選別では完全反射が回復することを確認した。

## 判定

| 項目 | 結果 |
|---|---:|
| 無補正で動的通過漏れを検出 | `{str(verdict["uncorrected_dynamic_leakage_detected"]).lower()}` |
| full 補正で動的完全反射が回復 | `{str(verdict["full_correction_dynamic_reflection_recovered"]).lower()}` |
| 動的散乱が二チャンネル期待式と一致 | `{str(verdict["dynamic_matches_expected_rates"]).lower()}` |
| ノルム保存 | `{str(verdict["dynamic_norm_preserved"]).lower()}` |
| 統合検証の最小判定 | `{str(verdict["integration_valid_minimal"]).lower()}` |

## 補正別集計

| 補正 | 最大動的通過漏れ | 最大閉鎖ペア RMS | 最大期待式誤差 |
|---|---:|---:|---:|
{rows}

## 読み

曲率相対位相の残留量を片側入射の交換干渉写像へ入れると、無補正および限定補正では通過漏れが残った。full 補正では残留位相が消え、動的散乱でも `R=1,T=0` が回復した。

したがって、閉鎖定常波の内部位相再選別は、閉鎖残差だけでなく、完全弾性反射の読出しも回復する。

## 図

![integration]({outputs["plot"]})

## 出力

| 種類 | ファイル |
|---|---|
| JSON | [{Path(outputs["json"]).name}]({outputs["json"]}) |
| CSV | [{Path(outputs["csv"]).name}]({outputs["csv"]}) |
| 図 | [{Path(outputs["plot"]).name}]({outputs["plot"]}) |
"""


def run() -> Dict[str, object]:
    params = Params()
    rows = run_cases(params)
    aggregate = aggregate_by_correction(rows)
    max_transmission_error = max(float(row["transmission_abs_error"]) for row in rows)
    max_reflection_error = max(float(row["reflection_abs_error"]) for row in rows)
    max_norm_error = max(float(row["norm_abs_error"]) for row in rows)

    verdict = {
        "uncorrected_dynamic_leakage_detected": bool(
            aggregate["none"]["max_dynamic_transmission_rate"] > params.uncorrected_min_transmission
        ),
        "full_correction_dynamic_reflection_recovered": bool(
            aggregate["full"]["max_dynamic_transmission_rate"] < params.scattering_error_tol
        ),
        "dynamic_matches_expected_rates": bool(
            max(max_transmission_error, max_reflection_error) < params.dynamic_expected_tol
        ),
        "dynamic_norm_preserved": bool(max_norm_error < params.scattering_error_tol),
    }
    verdict["integration_valid_minimal"] = bool(all(verdict.values()))

    outputs = {
        "json": "curved_closure_scattering_integration_result_v2.json",
        "csv": "curved_closure_scattering_integration_v2.csv",
        "plot": "curved_closure_scattering_integration_v2.png",
        "report": "curved_closure_scattering_integration_report_v2.md",
    }
    result: Dict[str, object] = {
        "experiment": "curved_closure_scattering_integration_v2",
        "params": asdict(params),
        "summary": {
            "aggregate_by_correction": aggregate,
            "max_transmission_abs_error": float(max_transmission_error),
            "max_reflection_abs_error": float(max_reflection_error),
            "max_norm_abs_error": float(max_norm_error),
        },
        "verdict": verdict,
        "output_dir": str(OUT_DIR.relative_to(BASE_DIR)),
        "outputs": outputs,
    }

    write_csv(OUT_DIR / outputs["csv"], rows)
    write_plot(rows, OUT_DIR / outputs["plot"])
    (OUT_DIR / outputs["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT_DIR / outputs["report"]).write_text(build_report(result), encoding="utf-8")
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["verdict"], ensure_ascii=False, indent=2))
