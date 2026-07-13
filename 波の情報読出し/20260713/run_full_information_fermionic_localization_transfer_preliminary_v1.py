from __future__ import annotations

import csv
import json
import math
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np


BASE_DIR = Path(__file__).resolve().parent
OUT_DIR = BASE_DIR / "full_information_fermionic_localization_transfer_preliminary_result_v1"
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
    chi_grid_n: int = 256
    eta_grid_n: int = 16
    high_n: int = 63
    chi_center: float = 0.0
    p0: float = 1.0
    q_A: float = 1.0
    q_B: float = -1.0
    A_A: float = 1.0
    A_B: float = 1.0
    m_A: int = 1
    m_B: int = 2
    delta_f_fermion: float = math.pi
    delta_f_boson: float = 0.0
    p_tol: float = 1.0e-2
    mode_tol: float = 1.0e-2
    copy_distance_tol: float = 1.0e-2
    norm_tol: float = 1.0e-10


def odd_harmonic_kernel(u: np.ndarray, nh: int) -> np.ndarray:
    numerator = np.sin((nh + 1) * u)
    denominator = (nh + 1) * np.sin(u)
    out = np.empty_like(u, dtype=float)
    regular = np.abs(np.sin(u)) > 1e-12
    out[regular] = numerator[regular] / denominator[regular]
    if np.any(~regular):
        k = np.rint(u[~regular] / math.pi).astype(int)
        out[~regular] = np.where(k % 2 == 0, 1.0, -1.0)
    return out


def normalize(v: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(v)
    if norm <= 0.0:
        raise ValueError("zero norm state")
    return v / norm


def make_grids(params: Params) -> Tuple[np.ndarray, np.ndarray]:
    chi = np.linspace(-math.pi, math.pi, params.chi_grid_n, endpoint=False)
    eta = np.linspace(-math.pi, math.pi, params.eta_grid_n, endpoint=False)
    return chi, eta


def make_state(
    params: Params,
    n_chi: int,
    q: float,
    m: int,
    hair_enabled: bool,
    amplitude: float,
) -> np.ndarray:
    chi, eta = make_grids(params)
    chi_part = odd_harmonic_kernel(chi - params.chi_center, n_chi)
    eta_phase = np.exp(1j * m * eta) if hair_enabled else np.ones_like(eta, dtype=complex)
    phase_chi = np.exp(1j * q * params.p0 * (chi - params.chi_center))
    psi = (chi_part * phase_chi)[:, None] * eta_phase[None, :]
    return amplitude * normalize(psi.reshape(-1))


def reshape_state(params: Params, v: np.ndarray) -> np.ndarray:
    return v.reshape(params.chi_grid_n, params.eta_grid_n)


def inner(v: np.ndarray, w: np.ndarray) -> complex:
    return complex(np.vdot(v, w))


def spectral_p_chi(params: Params, w: np.ndarray) -> np.ndarray:
    arr = reshape_state(params, w)
    freqs = np.fft.fftfreq(params.chi_grid_n, d=1.0 / params.chi_grid_n)
    transformed = np.fft.fft(arr, axis=0, norm="ortho")
    p_arr = np.fft.ifft(freqs[:, None] * transformed, axis=0, norm="ortho")
    return p_arr.reshape(-1)


def op_inner_p_chi(params: Params, v: np.ndarray, w: np.ndarray) -> complex:
    return inner(v, spectral_p_chi(params, w))


def chi_frequency_components(params: Params, v: np.ndarray) -> np.ndarray:
    arr = reshape_state(params, v)
    return np.fft.fft(arr, axis=0, norm="ortho")


def eta_frequency_components(params: Params, v: np.ndarray) -> np.ndarray:
    arr = reshape_state(params, v)
    return np.fft.fft(arr, axis=1, norm="ortho")


def frequency_index(freqs: np.ndarray, target: int) -> List[int]:
    return [int(i) for i, freq in enumerate(freqs) if int(round(freq)) == target]


def projector_inner_chi_abs(params: Params, v: np.ndarray, w: np.ndarray, n_abs: int) -> complex:
    fv = chi_frequency_components(params, v)
    fw = chi_frequency_components(params, w)
    freqs = np.fft.fftfreq(params.chi_grid_n, d=1.0 / params.chi_grid_n)
    indices = frequency_index(freqs, n_abs)
    if n_abs != 0:
        indices += frequency_index(freqs, -n_abs)
    if not indices:
        return 0.0 + 0.0j
    total = 0.0 + 0.0j
    for idx in indices:
        total += np.vdot(fv[idx, :], fw[idx, :])
    return complex(total)


def projector_inner_eta(params: Params, v: np.ndarray, w: np.ndarray, m: int) -> complex:
    fv = eta_frequency_components(params, v)
    fw = eta_frequency_components(params, w)
    freqs = np.fft.fftfreq(params.eta_grid_n, d=1.0 / params.eta_grid_n)
    indices = frequency_index(freqs, m)
    if not indices:
        return 0.0 + 0.0j
    total = 0.0 + 0.0j
    for idx in indices:
        total += np.vdot(fv[:, idx], fw[:, idx])
    return complex(total)


def rho_components(alpha: np.ndarray, beta: np.ndarray, delta_f: float) -> Dict[str, Any]:
    c = np.exp(1j * delta_f)
    a2 = float(np.vdot(alpha, alpha).real)
    b2 = float(np.vdot(beta, beta).real)
    s = inner(alpha, beta)
    norm = float((2.0 * a2 * b2 + 2.0 * math.cos(delta_f) * abs(s) ** 2).real)
    if norm <= 1.0e-30:
        raise ValueError("exchange state norm collapsed")
    return {"c": c, "a2": a2, "b2": b2, "s": s, "norm": norm}


def rho_apply(alpha: np.ndarray, beta: np.ndarray, comp: Dict[str, Any], x: np.ndarray) -> np.ndarray:
    c = comp["c"]
    s = comp["s"]
    a2 = comp["a2"]
    b2 = comp["b2"]
    norm = comp["norm"]
    ax = inner(alpha, x)
    bx = inner(beta, x)
    return (
        b2 * alpha * ax
        + a2 * beta * bx
        + np.conj(c) * s * alpha * bx
        + c * np.conj(s) * beta * ax
    ) / norm


def reduced_eigenvalues(alpha: np.ndarray, beta: np.ndarray, comp: Dict[str, Any]) -> np.ndarray:
    q1 = normalize(alpha)
    beta_residual = beta - q1 * inner(q1, beta)
    if np.linalg.norm(beta_residual) <= 1.0e-12:
        basis = [q1]
    else:
        basis = [q1, normalize(beta_residual)]
    matrix = np.zeros((len(basis), len(basis)), dtype=complex)
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            matrix[i, j] = inner(bi, rho_apply(alpha, beta, comp, bj))
    eig = np.linalg.eigvalsh(matrix)
    eig = np.maximum(eig.real, 0.0)
    return eig


def rho_trace_square(alpha: np.ndarray, beta: np.ndarray, comp: Dict[str, Any]) -> float:
    eig = reduced_eigenvalues(alpha, beta, comp)
    return float(np.sum(eig**2))


def rho_expectation(
    alpha: np.ndarray,
    beta: np.ndarray,
    comp: Dict[str, Any],
    bilinear,
) -> complex:
    c = comp["c"]
    s = comp["s"]
    a2 = comp["a2"]
    b2 = comp["b2"]
    norm = comp["norm"]
    return (
        b2 * bilinear(alpha, alpha)
        + a2 * bilinear(beta, beta)
        + np.conj(c) * s * bilinear(beta, alpha)
        + c * np.conj(s) * bilinear(alpha, beta)
    ) / norm


def density_diagonal(alpha: np.ndarray, beta: np.ndarray, comp: Dict[str, Any]) -> np.ndarray:
    c = comp["c"]
    s = comp["s"]
    a2 = comp["a2"]
    b2 = comp["b2"]
    norm = comp["norm"]
    diag = (
        b2 * np.abs(alpha) ** 2
        + a2 * np.abs(beta) ** 2
        + np.conj(c) * s * alpha * np.conj(beta)
        + c * np.conj(s) * beta * np.conj(alpha)
    ) / norm
    return np.maximum(diag.real, 0.0)


def harmonic_distribution(params: Params, alpha: np.ndarray, beta: np.ndarray, comp: Dict[str, Any]) -> Dict[int, float]:
    max_n = min(params.chi_grid_n // 2, params.high_n + 2)
    raw: Dict[int, float] = {}
    total = 0.0
    for n_abs in range(max_n + 1):
        value = rho_expectation(
            alpha,
            beta,
            comp,
            lambda v, w, n=n_abs: projector_inner_chi_abs(params, v, w, n),
        )
        amount = max(float(value.real), 0.0)
        if amount > 1.0e-14:
            raw[n_abs] = amount
            total += amount
    if total <= 0.0:
        return {0: 1.0}
    return {k: v / total for k, v in raw.items()}


def eta_distribution(
    params: Params,
    alpha: np.ndarray,
    beta: np.ndarray,
    comp: Dict[str, Any],
    modes: Iterable[int],
) -> Dict[int, float]:
    raw: Dict[int, float] = {}
    total = 0.0
    for mode in modes:
        value = rho_expectation(
            alpha,
            beta,
            comp,
            lambda v, w, m=mode: projector_inner_eta(params, v, w, m),
        )
        amount = max(float(value.real), 0.0)
        raw[mode] = amount
        total += amount
    if total <= 0.0:
        return {mode: 0.0 for mode in modes}
    return {mode: value / total for mode, value in raw.items()}


def effective_n(distribution: Dict[int, float]) -> Tuple[float, float]:
    n_eff = sum(float(n) * float(weight) for n, weight in distribution.items())
    n_eff_2 = math.sqrt(sum(float(n * n) * float(weight) for n, weight in distribution.items()))
    return n_eff, n_eff_2


def l1_distance(a: Dict[int, float], b: Dict[int, float]) -> float:
    keys = set(a) | set(b)
    return float(sum(abs(float(a.get(k, 0.0)) - float(b.get(k, 0.0))) for k in keys))


def pure_harmonic_distribution(params: Params, vector: np.ndarray) -> Dict[int, float]:
    comp = rho_components(vector, vector, 0.0)
    return harmonic_distribution(params, vector, vector, comp)


def pure_localization(vector: np.ndarray) -> float:
    prob = np.abs(vector) ** 2
    return float(np.sum(prob**2))


def pure_expect_p(params: Params, vector: np.ndarray) -> float:
    return float(op_inner_p_chi(params, vector, vector).real / inner(vector, vector).real)


def copy_distance(alpha: np.ndarray, beta: np.ndarray, comp: Dict[str, Any], pure: np.ndarray) -> float:
    tr_rho2 = rho_trace_square(alpha, beta, comp)
    overlap = float(inner(pure, rho_apply(alpha, beta, comp, pure)).real)
    numerator = max(tr_rho2 + 1.0 - 2.0 * overlap, 0.0)
    denominator = max(tr_rho2, 1.0e-30)
    return float(math.sqrt(numerator / denominator))


def full_model_metrics(
    params: Params,
    model: str,
    delta_f: float,
    n_a: int,
    n_b: int,
    hair_enabled: bool,
) -> List[Dict[str, Any]]:
    alpha = make_state(params, n_a, params.q_A, params.m_A, hair_enabled, params.A_A)
    beta = make_state(params, n_b, params.q_B, params.m_B, hair_enabled, params.A_B)
    comp = rho_components(alpha, beta, delta_f)
    p_read = float(rho_expectation(alpha, beta, comp, lambda v, w: op_inner_p_chi(params, v, w)).real)
    diag = density_diagonal(alpha, beta, comp)
    loc = float(np.sum(diag**2))
    eig = reduced_eigenvalues(alpha, beta, comp)
    trace = float(np.sum(eig))
    kappa = float(np.max(eig) / trace) if trace > 0.0 else 0.0
    harmonics = harmonic_distribution(params, alpha, beta, comp)
    n_eff, n_eff_2 = effective_n(harmonics)
    eta = eta_distribution(params, alpha, beta, comp, [params.m_A, params.m_B])
    copy_a = make_state(params, n_a, -params.q_A, params.m_A, hair_enabled, params.A_A)
    copy_b = make_state(params, n_b, -params.q_B, params.m_B, hair_enabled, params.A_B)
    h_a = pure_harmonic_distribution(params, alpha)
    h_b = pure_harmonic_distribution(params, beta)
    return [
        {
            "stage": "stage0_old_condition_reproduction",
            "model": model,
            "slot": "slot1",
            "delta_f": delta_f,
            "N_A": n_a,
            "N_B": n_b,
            "hair_enabled": hair_enabled,
            "R": trace,
            "p_chi": p_read,
            "p_target": -params.q_A,
            "p_abs_error": abs(p_read + params.q_A),
            "kappa": kappa,
            "L": loc,
            "N_eff": n_eff,
            "N_eff_2": n_eff_2,
            "P_m_A": eta.get(params.m_A, 0.0),
            "P_m_B": eta.get(params.m_B, 0.0),
            "target_mode_probability": eta.get(params.m_A, 0.0),
            "H_distance_to_initial_A": l1_distance(harmonics, h_a),
            "H_distance_to_initial_B": l1_distance(harmonics, h_b),
            "copy_distance_d": copy_distance(alpha, beta, comp, copy_a),
        },
        {
            "stage": "stage0_old_condition_reproduction",
            "model": model,
            "slot": "slot2",
            "delta_f": delta_f,
            "N_A": n_a,
            "N_B": n_b,
            "hair_enabled": hair_enabled,
            "R": trace,
            "p_chi": p_read,
            "p_target": -params.q_B,
            "p_abs_error": abs(p_read + params.q_B),
            "kappa": kappa,
            "L": loc,
            "N_eff": n_eff,
            "N_eff_2": n_eff_2,
            "P_m_A": eta.get(params.m_A, 0.0),
            "P_m_B": eta.get(params.m_B, 0.0),
            "target_mode_probability": eta.get(params.m_B, 0.0),
            "H_distance_to_initial_A": l1_distance(harmonics, h_a),
            "H_distance_to_initial_B": l1_distance(harmonics, h_b),
            "copy_distance_d": copy_distance(alpha, beta, comp, copy_b),
        },
    ]


def pure_model_metrics(
    params: Params,
    model: str,
    n_a: int,
    n_b: int,
    hair_enabled: bool,
    q_a: float,
    q_b: float,
) -> List[Dict[str, Any]]:
    vec_a = make_state(params, n_a, q_a, params.m_A, hair_enabled, params.A_A)
    vec_b = make_state(params, n_b, q_b, params.m_B, hair_enabled, params.A_B)
    rows: List[Dict[str, Any]] = []
    for slot, vec, n, target_m, target_p in [
        ("slot1", vec_a, n_a, params.m_A, -params.q_A),
        ("slot2", vec_b, n_b, params.m_B, -params.q_B),
    ]:
        h = pure_harmonic_distribution(params, vec)
        n_eff, n_eff_2 = effective_n(h)
        eta = {params.m_A: 1.0 if target_m == params.m_A and hair_enabled else 0.0,
               params.m_B: 1.0 if target_m == params.m_B and hair_enabled else 0.0}
        p = pure_expect_p(params, vec)
        rows.append(
            {
                "stage": "stage0_old_condition_reproduction",
                "model": model,
                "slot": slot,
                "delta_f": float("nan"),
                "N_A": n_a,
                "N_B": n_b,
                "hair_enabled": hair_enabled,
                "R": 1.0,
                "p_chi": p,
                "p_target": target_p,
                "p_abs_error": abs(p - target_p),
                "kappa": 1.0,
                "L": pure_localization(vec),
                "N_eff": n_eff,
                "N_eff_2": n_eff_2,
                "P_m_A": eta[params.m_A],
                "P_m_B": eta[params.m_B],
                "target_mode_probability": eta[target_m],
                "H_distance_to_initial_A": 0.0 if slot == "slot1" else float("nan"),
                "H_distance_to_initial_B": 0.0 if slot == "slot2" else float("nan"),
                "copy_distance_d": 0.0 if model == "copy_reflection" else float("nan"),
            }
        )
    return rows


def stage0_rows(params: Params) -> List[Dict[str, Any]]:
    n = params.high_n
    rows: List[Dict[str, Any]] = []
    rows.extend(full_model_metrics(params, "fermionic_full", params.delta_f_fermion, n, n, True))
    rows.extend(full_model_metrics(params, "bosonic_full", params.delta_f_boson, n, n, True))
    rows.extend(pure_model_metrics(params, "copy_reflection", n, n, True, -params.q_A, -params.q_B))
    rows.extend(pure_model_metrics(params, "simple_reflection", n, n, True, -params.q_A, -params.q_B))
    rows.extend(pure_model_metrics(params, "copy_transmission", n, n, True, params.q_A, params.q_B))
    return rows


def compute_stage0_verdict(params: Params, rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    selected = [row for row in rows if row["model"] == "fermionic_full"]
    p_ok = all(float(row["p_abs_error"]) <= params.p_tol for row in selected)
    mode_ok = all(float(row["target_mode_probability"]) >= 1.0 - params.mode_tol for row in selected)
    copy_close = all(float(row["copy_distance_d"]) <= params.copy_distance_tol for row in selected)
    norm_ok = all(abs(float(row["R"]) - 1.0) <= params.norm_tol for row in selected)
    reproduced = bool(p_ok and mode_ok and copy_close and norm_ok)
    return {
        "stage0_reproduced": reproduced,
        "stage0_p_reflection_ok": bool(p_ok),
        "stage0_label_preserved_ok": bool(mode_ok),
        "stage0_copy_distance_small": bool(copy_close),
        "stage0_norm_ok": bool(norm_ok),
        "later_stages_executed": False,
        "stop_reason": None
        if reproduced
        else "Stage 0 failed: full-information exchange did not reproduce the old copy-preserving reflection condition.",
    }


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def serialise_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for row in rows:
        item: Dict[str, Any] = {}
        for key, value in row.items():
            if isinstance(value, (np.integer, np.floating)):
                item[key] = float(value)
            elif isinstance(value, float) and math.isnan(value):
                item[key] = None
            else:
                item[key] = value
        out.append(item)
    return out


def make_plots(rows: List[Dict[str, Any]]) -> Dict[str, str]:
    models = sorted({str(row["model"]) for row in rows})
    slot1 = [row for row in rows if row["slot"] == "slot1"]
    p_by_model = [next(row for row in slot1 if row["model"] == model)["p_chi"] for model in models]
    mode_by_model = [
        next(row for row in slot1 if row["model"] == model)["target_mode_probability"] for model in models
    ]
    copy_by_model = [
        next(row for row in slot1 if row["model"] == model)["copy_distance_d"] for model in models
    ]

    fig, axes = plt.subplots(3, 1, figsize=(9, 9), constrained_layout=True)
    axes[0].bar(models, p_by_model)
    axes[0].axhline(-1.0, color="black", linestyle="--", linewidth=1)
    axes[0].set_ylabel("slot1 p_chi")
    axes[0].set_title("Stage 0 direction readout")
    axes[1].bar(models, mode_by_model)
    axes[1].axhline(1.0, color="black", linestyle="--", linewidth=1)
    axes[1].set_ylabel("target mode probability")
    axes[2].bar(models, [0.0 if isinstance(v, float) and math.isnan(v) else v for v in copy_by_model])
    axes[2].set_ylabel("copy distance d")
    axes[2].tick_params(axis="x", rotation=25)
    path = OUT_DIR / "full_information_stage0_old_condition_diagnostics_v1.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return {"stage0_plot": path.name}


def build_report(result: Dict[str, Any]) -> str:
    verdict = result["verdict"]
    rows = result["stage0_rows"]
    outputs = result["outputs"]
    table_lines = []
    for row in rows:
        table_lines.append(
            f"| {row['model']} | {row['slot']} | {row['p_chi']:.6e} | "
            f"{row['p_target']:.6e} | {row['target_mode_probability']:.6e} | "
            f"{row['kappa']:.6e} | {row['L']:.6e} | {row['copy_distance_d'] if row['copy_distance_d'] is not None else 'nan'} |"
        )

    return f"""# 全情報交換干渉フェルミオン的衝突 予備実験検証メモ v1

## 目的

全情報交換干渉写像を用いた場合に、過去の完全弾性反射系列の旧条件を再現できるかを Stage 0 として検査した。

本実験では、保存コピー近似を主実験に使わず、識別振動、代表振幅、倍音構造、局在カーネルを含む一体波を交換干渉へ入れ、縮約密度から読出し量を再構成した。

## 判定

| 項目 | 結果 |
|---|---:|
| Stage 0 再現 | `{str(verdict['stage0_reproduced']).lower()}` |
| p 反転 | `{str(verdict['stage0_p_reflection_ok']).lower()}` |
| 識別保存 | `{str(verdict['stage0_label_preserved_ok']).lower()}` |
| 保存コピー近似との差が小さい | `{str(verdict['stage0_copy_distance_small']).lower()}` |
| ノルム | `{str(verdict['stage0_norm_ok']).lower()}` |
| 後続 Stage 実行 | `{str(verdict['later_stages_executed']).lower()}` |

## Stage 0 結果

| model | slot | p_chi | p_target | target_mode_probability | kappa | L | copy_distance_d |
|---|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(table_lines)}

![stage0 diagnostics]({outputs['stage0_plot']})

## 解釈

全情報交換干渉をそのまま適用すると、旧モデルが保存コピーしていた識別振動および個体別出力は保存されなかった。

これは過去実験の否定ではない。

過去実験は、方向読出し反転、p/E/R 保存、多ゲージ読出しを検査する有効近似であった。

一方、本実験の目的である局在性および倍音移乗を調べるには、保存コピー近似を使えない。

今回の Stage 0 では、現在の静的な全情報交換干渉縮約だけでは、旧条件を再現できなかった。

仕様書に従い、Stage 1 以降の低局在性底探索、観測停止対照、名前の毛除去、非対称次数実験には進まない。

## 次の課題

次に必要なのは、静的な二体交換合成だけでなく、過去の反射実験で使われた片側入射、局所相互作用窓、偶奇チャンネル分解、自由伝播後読出しを、全情報状態に拡張した動的全情報交換干渉写像である。

これを作らない限り、局在性移乗実験は旧反射系列と比較可能な形にならない。

## 出力

| 種別 | ファイル |
|---|---|
| JSON | `{outputs['json']}` |
| CSV | `{outputs['stage0_csv']}` |
| 図 | `{outputs['stage0_plot']}` |
| report | `{outputs['report']}` |
"""


def run() -> Dict[str, Any]:
    params = Params()
    rows = stage0_rows(params)
    verdict = compute_stage0_verdict(params, rows)
    outputs = {
        "json": "full_information_fermionic_localization_transfer_preliminary_result_v1.json",
        "stage0_csv": "full_information_stage0_old_condition_rows_v1.csv",
        "stage0_plot": "full_information_stage0_old_condition_diagnostics_v1.png",
        "report": "full_information_fermionic_localization_transfer_preliminary_report_v1.md",
    }
    plot_outputs = make_plots(rows)
    outputs.update(plot_outputs)
    result: Dict[str, Any] = {
        "experiment": "full_information_fermionic_localization_transfer_preliminary_v1",
        "params": asdict(params),
        "verdict": verdict,
        "stage0_rows": serialise_rows(rows),
        "note": "Stage 0 is executed first. Later stages are intentionally not executed when Stage 0 fails.",
        "outputs": outputs,
    }
    write_csv(OUT_DIR / outputs["stage0_csv"], rows)
    (OUT_DIR / outputs["json"]).write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    report = build_report(result)
    (OUT_DIR / outputs["report"]).write_text(report, encoding="utf-8")
    (BASE_DIR / "全情報交換干渉フェルミオン的衝突における低局在性・倍音移乗予備実験検証メモ_v1.md").write_text(
        report,
        encoding="utf-8",
    )
    return result


if __name__ == "__main__":
    data = run()
    print(json.dumps(data["verdict"], ensure_ascii=False, indent=2))

