from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Dict, List

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "curved_closure_stationary_wave_broad_sweep_result_v1"
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
    delta_max_sweep_n: int = 31
    delta_max: float = 1.2
    closure_tol: float = 1e-16
    full_tol: float = 1e-16
    nontrivial_leakage_tol: float = 1e-4
    nontrivial_closure_tol: float = 1e-4


def odd_modes(params: Params) -> np.ndarray:
    return 2.0 * np.arange(params.mode_count, dtype=float) + 1.0


def mode_coordinate(params: Params) -> np.ndarray:
    h = odd_modes(params)
    return h / float(np.max(h))


def base_closed_pairs(params: Params) -> tuple[np.ndarray, np.ndarray]:
    amplitude = 1.0 / math.sqrt(2.0 * params.mode_count)
    x = amplitude * np.ones(params.mode_count, dtype=complex)
    return x, 1j * x


def normalize_shape(shape: np.ndarray) -> np.ndarray:
    peak = float(np.max(np.abs(shape)))
    return shape / peak if peak else shape


def phase_shape(name: str, u: np.ndarray) -> np.ndarray:
    if name == "linear":
        shape = u
    elif name == "quadratic_area":
        shape = u**2
    elif name == "cubic_high":
        shape = u**3
    elif name == "quartic_edge":
        shape = u**4
    elif name == "alternating_linear":
        sign = np.where(np.arange(u.size) % 2 == 0, 1.0, -1.0)
        shape = sign * u
    elif name == "sinusoidal_loop":
        shape = np.sin(math.pi * u)
    elif name == "mixed_smooth":
        shape = 0.55 * u + 0.35 * u**2 + 0.10 * np.sin(2.0 * math.pi * u)
    elif name == "rippled_random_like":
        shape = 0.52 * u + 0.21 * np.cos(3.0 * math.pi * u) - 0.17 * np.sin(5.0 * math.pi * u)
    else:
        raise ValueError(f"unknown phase shape: {name}")
    return normalize_shape(shape)


def phase_model_names() -> List[str]:
    return [
        "linear",
        "quadratic_area",
        "cubic_high",
        "quartic_edge",
        "alternating_linear",
        "sinusoidal_loop",
        "mixed_smooth",
        "rippled_random_like",
    ]


def design_matrix(correction: str, u: np.ndarray) -> np.ndarray | None:
    if correction == "none":
        return np.zeros((u.size, 0), dtype=float)
    if correction == "constant":
        return np.column_stack([np.ones_like(u)])
    if correction == "linear":
        return np.column_stack([u])
    if correction == "affine":
        return np.column_stack([np.ones_like(u), u])
    if correction == "quadratic":
        return np.column_stack([np.ones_like(u), u, u**2])
    if correction == "cubic":
        return np.column_stack([np.ones_like(u), u, u**2, u**3])
    if correction == "full":
        return None
    raise ValueError(f"unknown correction: {correction}")


def correction_names() -> List[str]:
    return ["none", "constant", "linear", "affine", "quadratic", "cubic", "full"]


def beta_for_correction(delta: np.ndarray, correction: str, u: np.ndarray, weights: np.ndarray) -> np.ndarray:
    if correction == "none":
        return np.zeros_like(delta)
    if correction == "full":
        return -delta

    basis = design_matrix(correction, u)
    assert basis is not None
    sw = np.sqrt(weights / float(np.sum(weights)))
    a = basis * sw[:, None]
    b = delta * sw
    coef, *_ = np.linalg.lstsq(a, b, rcond=None)
    return -(basis @ coef)


def conformal_factor(delta_max: float, u: np.ndarray, params: Params) -> np.ndarray:
    strength = delta_max / params.delta_max if params.delta_max else 0.0
    area_weight = 1.0 + 0.09 * strength * u**2
    common_phase = 0.25 * delta_max * u
    return np.sqrt(area_weight) * np.exp(1j * common_phase)


def closure_metrics(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    terms = x**2 + y**2
    denom = float(np.sum(np.abs(x) ** 2 + np.abs(y) ** 2))
    return {
        "closure_global_abs": float(abs(complex(np.sum(terms))) / denom if denom else 0.0),
        "closure_pair_rms": float(math.sqrt(float(np.mean(np.abs(terms) ** 2))) / denom if denom else 0.0),
        "closure_pair_max": float(np.max(np.abs(terms)) / denom if denom else 0.0),
    }


def exchange_metrics(effective_delta: np.ndarray, weights: np.ndarray) -> Dict[str, float]:
    weights = weights / float(np.sum(weights))
    transmission = np.sin(effective_delta / 2.0) ** 2
    reflection = np.cos(effective_delta / 2.0) ** 2
    return {
        "transmission_leakage": float(np.sum(weights * transmission)),
        "reflection_rate": float(np.sum(weights * reflection)),
        "max_remaining_delta": float(np.max(np.abs(effective_delta))),
        "rms_remaining_delta": float(math.sqrt(float(np.sum(weights * effective_delta**2)))),
    }


def evaluate_case(
    params: Params,
    phase_model: str,
    correction: str,
    delta_max: float,
) -> Dict[str, float | str]:
    u = mode_coordinate(params)
    x0, y0 = base_closed_pairs(params)
    g = conformal_factor(delta_max, u, params)
    weights = np.abs(g * x0) ** 2 + np.abs(g * y0) ** 2
    delta = delta_max * phase_shape(phase_model, u)
    beta = beta_for_correction(delta, correction, u, weights)
    effective_delta = delta + beta
    x = g * x0
    y = g * y0 * np.exp(1j * effective_delta)
    return {
        "phase_model": phase_model,
        "correction": correction,
        "delta_max": float(delta_max),
        **closure_metrics(x, y),
        **exchange_metrics(effective_delta, weights),
    }


def conformal_control(params: Params, delta_max: float) -> Dict[str, float]:
    u = mode_coordinate(params)
    x0, y0 = base_closed_pairs(params)
    g = conformal_factor(delta_max, u, params)
    metrics = closure_metrics(g * x0, g * y0)
    return {
        "delta_max": float(delta_max),
        **metrics,
    }


def run_sweep(params: Params) -> tuple[List[Dict[str, float | str]], List[Dict[str, float]]]:
    rows: List[Dict[str, float | str]] = []
    control_rows: List[Dict[str, float]] = []
    for delta_max in np.linspace(0.0, params.delta_max, params.delta_max_sweep_n):
        control_rows.append(conformal_control(params, float(delta_max)))
        for phase_model in phase_model_names():
            for correction in correction_names():
                rows.append(evaluate_case(params, phase_model, correction, float(delta_max)))
    return rows, control_rows


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def max_delta_rows(rows: List[Dict[str, float | str]], params: Params) -> List[Dict[str, float | str]]:
    return [row for row in rows if abs(float(row["delta_max"]) - params.delta_max) < 1e-14]


def aggregate_by_correction(rows: List[Dict[str, float | str]], params: Params) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {}
    target = max_delta_rows(rows, params)
    for correction in correction_names():
        corr_rows = [row for row in target if row["correction"] == correction]
        out[correction] = {
            "max_closure_pair_rms": float(max(float(row["closure_pair_rms"]) for row in corr_rows)),
            "max_transmission_leakage": float(max(float(row["transmission_leakage"]) for row in corr_rows)),
            "mean_transmission_leakage": float(np.mean([float(row["transmission_leakage"]) for row in corr_rows])),
            "max_remaining_delta": float(max(float(row["max_remaining_delta"]) for row in corr_rows)),
        }
    return out


def model_correction_matrix(
    rows: List[Dict[str, float | str]],
    params: Params,
    value_key: str,
) -> np.ndarray:
    target = max_delta_rows(rows, params)
    matrix = np.zeros((len(phase_model_names()), len(correction_names())), dtype=float)
    for i, model in enumerate(phase_model_names()):
        for j, correction in enumerate(correction_names()):
            row = next(
                row for row in target if row["phase_model"] == model and row["correction"] == correction
            )
            matrix[i, j] = float(row[value_key])
    return matrix


def write_heatmap(
    matrix: np.ndarray,
    title: str,
    colorbar_label: str,
    path: Path,
) -> None:
    plot_matrix = np.log10(np.maximum(matrix, 1e-18))
    fig, ax = plt.subplots(figsize=(9, 5))
    im = ax.imshow(plot_matrix, aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(correction_names())), correction_names(), rotation=35, ha="right")
    ax.set_yticks(np.arange(len(phase_model_names())), phase_model_names())
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(colorbar_label)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_aggregate_plot(rows: List[Dict[str, float | str]], path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for correction in correction_names():
        xs: List[float] = []
        max_closure: List[float] = []
        max_leakage: List[float] = []
        for delta_max in sorted({float(row["delta_max"]) for row in rows}):
            selected = [
                row
                for row in rows
                if row["correction"] == correction and abs(float(row["delta_max"]) - delta_max) < 1e-14
            ]
            xs.append(delta_max)
            max_closure.append(max(float(row["closure_pair_rms"]) for row in selected))
            max_leakage.append(max(float(row["transmission_leakage"]) for row in selected))
        axes[0].plot(xs, max_closure, marker="o", markersize=3, label=correction)
        axes[1].plot(xs, max_leakage, marker="o", markersize=3, label=correction)

    axes[0].set_ylabel("max closure pair RMS")
    axes[0].set_yscale("symlog", linthresh=1e-18)
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(ncol=4, fontsize=8)
    axes[1].set_xlabel("max curvature relative phase")
    axes[1].set_ylabel("max T leakage")
    axes[1].set_yscale("symlog", linthresh=1e-18)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(ncol=4, fontsize=8)
    fig.suptitle("Correction capacity across curvature phase models")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def build_report(result: Dict[str, object]) -> str:
    summary = result["summary"]
    verdict = result["verdict"]
    outputs = result["outputs"]
    aggregate = summary["aggregate_at_max"]

    aggregate_rows = "\n".join(
        f"| {name} | `{data['max_closure_pair_rms']:.16e}` | `{data['max_transmission_leakage']:.16e}` | `{data['max_remaining_delta']:.16e}` |"
        for name, data in aggregate.items()
    )

    return f"""# 曲率付き閉鎖定常波 広範囲検証 v1

## 目的

最小実験 v1 では、曲率相対位相が固定波に入ると閉鎖残差と通過漏れが出ること、また完全な内部位相再選別でそれらが消えることを確認した。

本検証では、曲率相対位相モデルを複数に増やし、内部補正の自由度を `none`, `constant`, `linear`, `affine`, `quadratic`, `cubic`, `full` に制限して、どの範囲で閉鎖定常波が回復するかを調べた。

## 判定

| 項目 | 結果 |
|---|---:|
| 共通因子型曲率は全掃引で閉鎖保存 | `{str(verdict["conformal_controls_preserve_closure"]).lower()}` |
| 無補正では非自明な曲率漏れを検出 | `{str(verdict["uncorrected_leakage_detected"]).lower()}` |
| 完全補正は全位相モデルで閉鎖回復 | `{str(verdict["full_correction_recovers_closure_all_models"]).lower()}` |
| 完全補正は全位相モデルで完全反射回復 | `{str(verdict["full_correction_recovers_reflection_all_models"]).lower()}` |
| 限定補正はモデル依存の残差を残す | `{str(verdict["limited_correction_is_model_dependent"]).lower()}` |
| 広範囲検証の最小判定 | `{str(verdict["broad_validation_supported"]).lower()}` |

## 最大曲率相対位相での補正別集計

| 補正 | 最大閉鎖ペア RMS | 最大通過漏れ | 最大残留位相 |
|---|---:|---:|---:|
{aggregate_rows}

## 読み

固定波に曲率相対位相を入れると、すべての位相モデルで閉鎖残差と通過漏れが出る。共通因子型の曲率作用は閉鎖を破らない。完全補正は全モデルで閉鎖と完全反射を回復する。一方、定数・線形・二次・三次の限定補正では、補正基底に含まれない位相モデルに残差が残る。

したがって、本検証の範囲では、曲率効果が観測読出しから消えるためには、曲率相対位相を閉鎖定常波の内部位相再選別として吸収するだけの自由度が必要である。

## 図

![aggregate]({outputs["aggregate_plot"]})

![closure heatmap]({outputs["closure_heatmap"]})

![leakage heatmap]({outputs["leakage_heatmap"]})

## 出力

| 種類 | ファイル |
|---|---|
| JSON | [{Path(outputs["json"]).name}]({outputs["json"]}) |
| sweep CSV | [{Path(outputs["sweep_csv"]).name}]({outputs["sweep_csv"]}) |
| conformal control CSV | [{Path(outputs["control_csv"]).name}]({outputs["control_csv"]}) |
| aggregate 図 | [{Path(outputs["aggregate_plot"]).name}]({outputs["aggregate_plot"]}) |
| closure heatmap | [{Path(outputs["closure_heatmap"]).name}]({outputs["closure_heatmap"]}) |
| leakage heatmap | [{Path(outputs["leakage_heatmap"]).name}]({outputs["leakage_heatmap"]}) |
"""


def run() -> Dict[str, object]:
    params = Params()
    rows, control_rows = run_sweep(params)
    aggregate = aggregate_by_correction(rows, params)
    control_max = max(float(row["closure_pair_rms"]) for row in control_rows)

    verdict = {
        "conformal_controls_preserve_closure": bool(control_max < params.closure_tol),
        "uncorrected_leakage_detected": bool(
            aggregate["none"]["max_closure_pair_rms"] > params.nontrivial_closure_tol
            and aggregate["none"]["max_transmission_leakage"] > params.nontrivial_leakage_tol
        ),
        "full_correction_recovers_closure_all_models": bool(
            aggregate["full"]["max_closure_pair_rms"] < params.full_tol
        ),
        "full_correction_recovers_reflection_all_models": bool(
            aggregate["full"]["max_transmission_leakage"] < params.full_tol
        ),
        "limited_correction_is_model_dependent": bool(
            aggregate["cubic"]["max_transmission_leakage"] > params.nontrivial_leakage_tol
            and aggregate["quadratic"]["max_transmission_leakage"] > params.nontrivial_leakage_tol
            and aggregate["linear"]["max_transmission_leakage"] > params.nontrivial_leakage_tol
        ),
    }
    verdict["broad_validation_supported"] = bool(all(verdict.values()))

    outputs = {
        "json": "curved_closure_stationary_wave_broad_sweep_result_v1.json",
        "sweep_csv": "curved_closure_stationary_wave_broad_sweep_v1.csv",
        "control_csv": "curved_closure_stationary_wave_broad_conformal_control_v1.csv",
        "aggregate_plot": "curved_closure_stationary_wave_broad_aggregate_v1.png",
        "closure_heatmap": "curved_closure_stationary_wave_broad_closure_heatmap_v1.png",
        "leakage_heatmap": "curved_closure_stationary_wave_broad_leakage_heatmap_v1.png",
        "report": "curved_closure_stationary_wave_broad_report_v1.md",
    }

    summary = {
        "phase_models": phase_model_names(),
        "corrections": correction_names(),
        "aggregate_at_max": aggregate,
        "conformal_control_max_closure_pair_rms": float(control_max),
    }

    result: Dict[str, object] = {
        "experiment": "curved_closure_stationary_wave_broad_sweep_v1",
        "params": asdict(params),
        "summary": summary,
        "verdict": verdict,
        "output_dir": str(OUT_DIR.relative_to(BASE_DIR)),
        "outputs": outputs,
    }

    write_csv(OUT_DIR / outputs["sweep_csv"], rows)
    write_csv(OUT_DIR / outputs["control_csv"], control_rows)
    write_aggregate_plot(rows, OUT_DIR / outputs["aggregate_plot"])
    write_heatmap(
        model_correction_matrix(rows, params, "closure_pair_rms"),
        "log10 closure pair RMS at max curvature phase",
        "log10 closure pair RMS",
        OUT_DIR / outputs["closure_heatmap"],
    )
    write_heatmap(
        model_correction_matrix(rows, params, "transmission_leakage"),
        "log10 transmission leakage at max curvature phase",
        "log10 T leakage",
        OUT_DIR / outputs["leakage_heatmap"],
    )
    (OUT_DIR / outputs["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT_DIR / outputs["report"]).write_text(build_report(result), encoding="utf-8")
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["verdict"], ensure_ascii=False, indent=2))
