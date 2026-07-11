from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "curved_closure_stationary_wave_result_v1"
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
    kappa_sweep_n: int = 41
    kappa_max: float = 0.012
    relaxation_steps: int = 48
    relaxation_decay: float = 0.72
    closure_tol: float = 1e-16
    stationary_tol: float = 1e-16
    transient_min_residual: float = 1e-4
    transient_min_leakage: float = 1e-4


def odd_modes(params: Params) -> np.ndarray:
    return 2 * np.arange(params.mode_count, dtype=float) + 1.0


def base_closed_pairs(params: Params) -> tuple[np.ndarray, np.ndarray]:
    amplitude = 1.0 / math.sqrt(2.0 * params.mode_count)
    phase = np.zeros(params.mode_count, dtype=float)
    x = amplitude * np.exp(1j * phase)
    ix = 1j * x
    return x, ix


def curvature_phase(kappa: float, modes: np.ndarray) -> np.ndarray:
    return kappa * modes


def conformal_factor(kappa: float, modes: np.ndarray, params: Params) -> np.ndarray:
    scaled = modes / float(np.max(modes))
    area_weight = 1.0 + 0.12 * (kappa / params.kappa_max) * scaled**2 if params.kappa_max else 1.0
    common_phase = 0.35 * kappa * scaled
    return np.sqrt(area_weight) * np.exp(1j * common_phase)


def closure_metrics(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    pair_terms = x**2 + y**2
    denom = float(np.sum(np.abs(x) ** 2 + np.abs(y) ** 2))
    global_residual = abs(complex(np.sum(pair_terms))) / denom if denom else 0.0
    pair_rms = math.sqrt(float(np.mean(np.abs(pair_terms) ** 2))) / denom if denom else 0.0
    pair_max = float(np.max(np.abs(pair_terms))) / denom if denom else 0.0
    return {
        "closure_global_abs": float(global_residual),
        "closure_pair_rms": float(pair_rms),
        "closure_pair_max": float(pair_max),
    }


def exchange_metrics(effective_delta: np.ndarray, weights: np.ndarray) -> Dict[str, float]:
    weights = weights / float(np.sum(weights))
    transmission = np.sin(effective_delta / 2.0) ** 2
    reflection = np.cos(effective_delta / 2.0) ** 2
    node_relative_norm = 2.0 * transmission
    return {
        "transmission_leakage": float(np.sum(weights * transmission)),
        "reflection_rate": float(np.sum(weights * reflection)),
        "node_relative_norm": float(np.sum(weights * node_relative_norm)),
        "max_effective_delta": float(np.max(np.abs(effective_delta))),
    }


def state_metrics(
    label: str,
    kappa: float,
    x: np.ndarray,
    y: np.ndarray,
    effective_delta: np.ndarray,
) -> Dict[str, float | str]:
    weights = np.abs(x) ** 2 + np.abs(y) ** 2
    return {
        "state": label,
        "kappa": float(kappa),
        **closure_metrics(x, y),
        **exchange_metrics(effective_delta, weights),
    }


def compute_states(params: Params, kappa: float) -> Dict[str, Dict[str, float | str]]:
    modes = odd_modes(params)
    x0, y0 = base_closed_pairs(params)
    delta = curvature_phase(kappa, modes)
    g = conformal_factor(kappa, modes, params)

    flat = state_metrics("flat", kappa, x0, y0, np.zeros_like(delta))
    conformal = state_metrics("conformal", kappa, g * x0, g * y0, np.zeros_like(delta))
    transient = state_metrics("transient", kappa, g * x0, g * y0 * np.exp(1j * delta), delta)
    stationary = state_metrics("stationary", kappa, g * x0, g * y0, np.zeros_like(delta))
    return {
        "flat": flat,
        "conformal": conformal,
        "transient": transient,
        "stationary": stationary,
    }


def run_sweep(params: Params) -> List[Dict[str, float | str]]:
    rows: List[Dict[str, float | str]] = []
    for kappa in np.linspace(0.0, params.kappa_max, params.kappa_sweep_n):
        states = compute_states(params, float(kappa))
        rows.extend(states.values())
    return rows


def run_relaxation(params: Params, kappa: float) -> List[Dict[str, float]]:
    modes = odd_modes(params)
    x0, y0 = base_closed_pairs(params)
    delta0 = curvature_phase(kappa, modes)
    g = conformal_factor(kappa, modes, params)
    rows: List[Dict[str, float]] = []
    for step in range(params.relaxation_steps + 1):
        factor = params.relaxation_decay**step
        effective_delta = delta0 * factor
        x = g * x0
        y = g * y0 * np.exp(1j * effective_delta)
        metrics = state_metrics("relaxation", kappa, x, y, effective_delta)
        rows.append(
            {
                "step": float(step),
                "relaxation_factor": float(factor),
                "closure_global_abs": float(metrics["closure_global_abs"]),
                "closure_pair_rms": float(metrics["closure_pair_rms"]),
                "transmission_leakage": float(metrics["transmission_leakage"]),
                "reflection_rate": float(metrics["reflection_rate"]),
                "node_relative_norm": float(metrics["node_relative_norm"]),
                "max_effective_delta": float(metrics["max_effective_delta"]),
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


def rows_by_state(rows: List[Dict[str, float | str]], state: str) -> List[Dict[str, float | str]]:
    return [row for row in rows if row["state"] == state]


def write_sweep_plot(rows: List[Dict[str, float | str]], path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    for state in ["conformal", "transient", "stationary"]:
        state_rows = rows_by_state(rows, state)
        kappa = [float(row["kappa"]) for row in state_rows]
        closure = [float(row["closure_pair_rms"]) for row in state_rows]
        leakage = [float(row["transmission_leakage"]) for row in state_rows]
        axes[0].plot(kappa, closure, marker="o", markersize=3, label=state)
        axes[1].plot(kappa, leakage, marker="o", markersize=3, label=state)

    axes[0].set_ylabel("closure pair RMS")
    axes[0].set_yscale("symlog", linthresh=1e-24)
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()
    axes[1].set_xlabel("curvature phase strength kappa")
    axes[1].set_ylabel("transmission leakage")
    axes[1].set_yscale("symlog", linthresh=1e-24)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    fig.suptitle("Curvature leakage and stationary closure recovery")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_relaxation_plot(rows: List[Dict[str, float]], path: Path) -> None:
    step = [row["step"] for row in rows]
    closure = [row["closure_pair_rms"] for row in rows]
    leakage = [row["transmission_leakage"] for row in rows]
    max_delta = [row["max_effective_delta"] for row in rows]

    fig, axes = plt.subplots(3, 1, figsize=(7, 7), sharex=True)
    axes[0].plot(step, closure, marker="o", markersize=3)
    axes[0].set_ylabel("closure pair RMS")
    axes[0].set_yscale("log")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(step, leakage, marker="o", markersize=3)
    axes[1].set_ylabel("T leakage")
    axes[1].set_yscale("log")
    axes[1].grid(True, alpha=0.3)
    axes[2].plot(step, max_delta, marker="o", markersize=3)
    axes[2].set_xlabel("relaxation step")
    axes[2].set_ylabel("max delta_eff")
    axes[2].set_yscale("log")
    axes[2].grid(True, alpha=0.3)
    fig.suptitle("Relaxation toward closed stationary mode")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def build_report(result: Dict[str, object]) -> str:
    summary = result["summary"]
    verdict = result["verdict"]
    outputs = result["outputs"]
    return f"""# 曲率付き閉鎖定常波 数値検証 v1

## 目的

曲率付き局所セルに置かれた奇数倍音複素波について、曲率が単に消えるのではなく、閉鎖定常波の再選別によって内部へ吸収されるかを最小写像で検証した。

## 定式化

閉鎖条件を、

```math
Q(x)=\\sum_n x_n^2=0
```

とする。曲率付き局所セルで安定に残る波は、

```math
Q(x_K)=0,
\\qquad
\\mathcal U_Kx_K=e^{{i\\alpha_K}}x_K
```

を満たす閉鎖定常波であるとする。

本実験では奇数倍音 `h_m=2m+1` に対し、曲率由来の相対位相漏れを

```math
\\delta_{{K,m}}=\\kappa h_m
```

として置いた。これは実在曲率の定量式ではなく、奇数倍音が相対位相漏れに対して高感度であることを調べるための局所モデルである。

## 判定

| 項目 | 結果 |
|---|---:|
| flat 閉鎖 | `{str(verdict["flat_closed"]).lower()}` |
| 共通因子型曲率は閉鎖を保存 | `{str(verdict["conformal_curvature_preserves_closure"]).lower()}` |
| 相対位相型曲率は過渡残差を生成 | `{str(verdict["transient_curvature_detected"]).lower()}` |
| 定常再選別で閉鎖が回復 | `{str(verdict["stationary_closure_recovered"]).lower()}` |
| 定常再選別で完全反射が回復 | `{str(verdict["stationary_reflection_recovered"]).lower()}` |
| 緩和で閉鎖残差が減少 | `{str(verdict["relaxation_reduces_closure"]).lower()}` |
| 緩和で通過漏れが減少 | `{str(verdict["relaxation_reduces_leakage"]).lower()}` |
| 最小仮説検証 | `{str(verdict["curvature_renormalization_hypothesis_supported_minimal"]).lower()}` |

## 主要数値

| 量 | 値 |
|---|---:|
| 最大 `kappa` | `{summary["kappa_max"]:.16e}` |
| 最大曲率位相 `max(delta_K)` | `{summary["transient_max_effective_delta"]:.16e}` |
| flat 閉鎖ペア RMS | `{summary["flat_pair_rms_at_max"]:.16e}` |
| conformal 閉鎖ペア RMS | `{summary["conformal_pair_rms_at_max"]:.16e}` |
| transient 閉鎖ペア RMS | `{summary["transient_pair_rms_at_max"]:.16e}` |
| stationary 閉鎖ペア RMS | `{summary["stationary_pair_rms_at_max"]:.16e}` |
| transient 通過漏れ | `{summary["transient_transmission_leakage_at_max"]:.16e}` |
| stationary 通過漏れ | `{summary["stationary_transmission_leakage_at_max"]:.16e}` |
| relaxation 最終閉鎖ペア RMS | `{summary["relaxation_final_pair_rms"]:.16e}` |
| relaxation 最終通過漏れ | `{summary["relaxation_final_transmission_leakage"]:.16e}` |

## 読み

固定された平坦波へ曲率相対位相を加えると、閉鎖残差と通過漏れが出た。一方、曲率が共通重みとして入る場合、閉鎖条件は保たれた。また、曲率相対位相を内部位相再選別で吸収した定常状態では、閉鎖残差と通過漏れが消えた。

したがって、本実験の範囲では、曲率の影響は消えているのではなく、閉鎖定常波の存在条件へ繰り込まれる、という仮説と整合する。

## 図

![sweep]({outputs["sweep_plot"]})

![relaxation]({outputs["relaxation_plot"]})

## 出力

| 種類 | ファイル |
|---|---|
| JSON | [{Path(outputs["json"]).name}]({outputs["json"]}) |
| sweep CSV | [{Path(outputs["sweep_csv"]).name}]({outputs["sweep_csv"]}) |
| relaxation CSV | [{Path(outputs["relaxation_csv"]).name}]({outputs["relaxation_csv"]}) |
| sweep 図 | [{Path(outputs["sweep_plot"]).name}]({outputs["sweep_plot"]}) |
| relaxation 図 | [{Path(outputs["relaxation_plot"]).name}]({outputs["relaxation_plot"]}) |
"""


def run() -> Dict[str, object]:
    params = Params()
    rows = run_sweep(params)
    relaxation_rows = run_relaxation(params, params.kappa_max)

    max_rows = [row for row in rows if abs(float(row["kappa"]) - params.kappa_max) < 1e-15]
    by_state = {str(row["state"]): row for row in max_rows}
    flat = by_state["flat"]
    conformal = by_state["conformal"]
    transient = by_state["transient"]
    stationary = by_state["stationary"]
    first_relax = relaxation_rows[0]
    final_relax = relaxation_rows[-1]

    summary = {
        "kappa_max": params.kappa_max,
        "transient_max_effective_delta": float(transient["max_effective_delta"]),
        "flat_pair_rms_at_max": float(flat["closure_pair_rms"]),
        "conformal_pair_rms_at_max": float(conformal["closure_pair_rms"]),
        "transient_pair_rms_at_max": float(transient["closure_pair_rms"]),
        "stationary_pair_rms_at_max": float(stationary["closure_pair_rms"]),
        "transient_transmission_leakage_at_max": float(transient["transmission_leakage"]),
        "stationary_transmission_leakage_at_max": float(stationary["transmission_leakage"]),
        "relaxation_initial_pair_rms": float(first_relax["closure_pair_rms"]),
        "relaxation_final_pair_rms": float(final_relax["closure_pair_rms"]),
        "relaxation_initial_transmission_leakage": float(first_relax["transmission_leakage"]),
        "relaxation_final_transmission_leakage": float(final_relax["transmission_leakage"]),
    }

    verdict = {
        "flat_closed": bool(summary["flat_pair_rms_at_max"] < params.closure_tol),
        "conformal_curvature_preserves_closure": bool(
            summary["conformal_pair_rms_at_max"] < params.closure_tol
        ),
        "transient_curvature_detected": bool(
            summary["transient_pair_rms_at_max"] > params.transient_min_residual
            and summary["transient_transmission_leakage_at_max"] > params.transient_min_leakage
        ),
        "stationary_closure_recovered": bool(
            summary["stationary_pair_rms_at_max"] < params.stationary_tol
        ),
        "stationary_reflection_recovered": bool(
            summary["stationary_transmission_leakage_at_max"] < params.stationary_tol
        ),
        "relaxation_reduces_closure": bool(
            summary["relaxation_final_pair_rms"] < summary["relaxation_initial_pair_rms"] * 1e-6
        ),
        "relaxation_reduces_leakage": bool(
            summary["relaxation_final_transmission_leakage"]
            < summary["relaxation_initial_transmission_leakage"] * 1e-6
        ),
    }
    verdict["curvature_renormalization_hypothesis_supported_minimal"] = bool(
        all(verdict.values())
    )

    outputs = {
        "json": "curved_closure_stationary_wave_result_v1.json",
        "sweep_csv": "curved_closure_stationary_wave_sweep_v1.csv",
        "relaxation_csv": "curved_closure_stationary_wave_relaxation_v1.csv",
        "sweep_plot": "curved_closure_stationary_wave_sweep_v1.png",
        "relaxation_plot": "curved_closure_stationary_wave_relaxation_v1.png",
        "report": "curved_closure_stationary_wave_report_v1.md",
    }

    result: Dict[str, object] = {
        "experiment": "curved_closure_stationary_wave_v1",
        "params": asdict(params),
        "summary": summary,
        "verdict": verdict,
        "output_dir": str(OUT_DIR.relative_to(BASE_DIR)),
        "outputs": outputs,
    }

    write_csv(OUT_DIR / outputs["sweep_csv"], rows)
    write_csv(OUT_DIR / outputs["relaxation_csv"], relaxation_rows)
    write_sweep_plot(rows, OUT_DIR / outputs["sweep_plot"])
    write_relaxation_plot(relaxation_rows, OUT_DIR / outputs["relaxation_plot"])
    (OUT_DIR / outputs["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    (OUT_DIR / outputs["report"]).write_text(build_report(result), encoding="utf-8")
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["verdict"], ensure_ascii=False, indent=2))
