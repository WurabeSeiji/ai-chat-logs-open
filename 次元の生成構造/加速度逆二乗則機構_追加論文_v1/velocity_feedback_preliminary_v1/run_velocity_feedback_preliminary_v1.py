"""AB二体閉鎖位相系・速度フィードバック（v=aτ）予備実験 v1

第1論文（逆二乗則 v1）の零次構造を対照とし、読出し加速度が
キャリア角速度へ蓄積する一階フィードバックを検証する。

零次（v1対照, kappa=0）:
    omega_s = omega_0 固定、偏差は Δ²χ = -omega_d² χ の調和振動。
一階（フィードバック, kappa≠0）:
    a_s      = -4 sin²(omega_s/2) · χ_s        （符号付き読出し加速度）
    χ_{s+1}  = 2χ_s - χ_{s-1} + a_s            （偏差の二階更新）
    omega_{s+1} = omega_s + kappa · a_s        （v = ∫a dτ の蓄積）
    Phi_{s+1}   = Phi_s + omega_{s+1}          （キャリア進行）

閉鎖表現: x1 = ρ e^{i(Phi+χ/2)}, x2 = ρ e^{i(Phi-χ/2+π/2)}
    厳密二波閉鎖残差 |x1²+x2²| = 2ρ²|sin χ|（χ=0で厳密閉鎖）。

検査項目:
  J1 運動学的整合   omega_s - omega_0 = kappa Σ a  （機械精度）
  J2 読出し成立     窓ごとの回帰 Δ²χ vs χ の傾き = -4sin²(ω̄/2)
  J3 有界性         omega・χ包絡・閉鎖残差包絡の長時間有界性
  J4 閉鎖残差       対照とフィードバックの残差包絡の比較
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent
RESULT_DIR = BASE_DIR / "velocity_feedback_preliminary_result_v1"

PERIOD_STEPS = 96
OMEGA_0 = 2.0 * math.pi / PERIOD_STEPS
STEP_COUNT = 96 * 200
INITIAL_DEVIATIONS_DEG = (2.0, 5.0, 10.0, 20.0)
KAPPA_VALUES = (0.0, 0.05, 0.2, 0.5, -0.2)
RHO = 1.0
WINDOW = PERIOD_STEPS


def simulate(chi0_deg: float, kappa: float) -> dict:
    chi0 = math.radians(chi0_deg)
    omega = np.empty(STEP_COUNT + 1)
    chi = np.empty(STEP_COUNT + 1)
    accel = np.empty(STEP_COUNT + 1)
    omega[0] = OMEGA_0
    chi[0] = chi0
    accel[0] = -4.0 * math.sin(omega[0] / 2.0) ** 2 * chi[0]
    # 初速0のleapfrog初期化
    chi[1] = chi[0] + 0.5 * accel[0]
    omega[1] = omega[0] + kappa * accel[0]
    for s in range(1, STEP_COUNT):
        accel[s] = -4.0 * math.sin(omega[s] / 2.0) ** 2 * chi[s]
        chi[s + 1] = 2.0 * chi[s] - chi[s - 1] + accel[s]
        omega[s + 1] = omega[s] + kappa * accel[s]
    accel[STEP_COUNT] = (
        -4.0 * math.sin(omega[STEP_COUNT] / 2.0) ** 2 * chi[STEP_COUNT]
    )

    # J1 運動学的整合: omega - omega0 = kappa * cumsum(a)
    cumulative_accel = np.concatenate([[0.0], np.cumsum(accel[:-1])])
    kinematic_error = float(
        np.max(np.abs((omega - OMEGA_0) - kappa * cumulative_accel))
    )

    # J2 読出し成立: 窓ごとの回帰傾き vs 理論係数
    second_diff = chi[2:] - 2.0 * chi[1:-1] + chi[:-2]
    chi_mid = chi[1:-1]
    omega_mid = omega[1:-1]
    slope_rel_errors = []
    for start in range(0, len(chi_mid) - WINDOW, WINDOW):
        cw = chi_mid[start:start + WINDOW]
        dw = second_diff[start:start + WINDOW]
        ow = omega_mid[start:start + WINDOW]
        denom = float(np.dot(cw, cw))
        if denom < 1.0e-24:
            continue
        slope = float(np.dot(cw, dw)) / denom
        theory = -4.0 * math.sin(float(np.mean(ow)) / 2.0) ** 2
        slope_rel_errors.append(abs(slope - theory) / abs(theory))
    max_slope_rel_error = float(np.max(slope_rel_errors))

    # J3/J4 閉鎖残差と包絡
    closure_residual = 2.0 * RHO * RHO * np.abs(np.sin(chi))
    half = STEP_COUNT // 2
    residual_env_first = float(np.max(closure_residual[:half]))
    residual_env_second = float(np.max(closure_residual[half:]))
    chi_env_first = float(np.max(np.abs(chi[:half])))
    chi_env_second = float(np.max(np.abs(chi[half:])))
    omega_min = float(np.min(omega))
    omega_max = float(np.max(omega))

    return {
        "chi0_deg": chi0_deg,
        "kappa": kappa,
        "kinematic_error": kinematic_error,
        "max_slope_rel_error": max_slope_rel_error,
        "omega_min": omega_min,
        "omega_max": omega_max,
        "omega_drift_ratio": (omega_max - omega_min) / OMEGA_0,
        "chi_envelope_first_half": chi_env_first,
        "chi_envelope_second_half": chi_env_second,
        "chi_envelope_growth": chi_env_second / chi_env_first,
        "closure_residual_max_first_half": residual_env_first,
        "closure_residual_max_second_half": residual_env_second,
        "closure_residual_growth": residual_env_second / residual_env_first,
        "series": {
            "omega": omega,
            "chi": chi,
            "accel": accel,
            "cumulative_accel": cumulative_accel,
            "closure_residual": closure_residual,
        },
    }


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    kept_series = {}
    for chi0_deg in INITIAL_DEVIATIONS_DEG:
        for kappa in KAPPA_VALUES:
            result = simulate(chi0_deg, kappa)
            series = result.pop("series")
            rows.append(result)
            if chi0_deg == 10.0 and kappa in (0.0, 0.2):
                kept_series[f"chi10_kappa{kappa}"] = series
    with (RESULT_DIR / "velocity_feedback_trials_v1.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    np.savez_compressed(
        RESULT_DIR / "velocity_feedback_selected_series_v1.npz",
        **{
            f"{name}_{key}": value
            for name, series in kept_series.items()
            for key, value in series.items()
        },
    )
    payload = {
        "parameters": {
            "period_steps": PERIOD_STEPS,
            "omega_0": OMEGA_0,
            "step_count": STEP_COUNT,
            "initial_deviations_deg": INITIAL_DEVIATIONS_DEG,
            "kappa_values": KAPPA_VALUES,
            "rho": RHO,
            "window": WINDOW,
        },
        "trials": rows,
        "max_kinematic_error": max(row["kinematic_error"] for row in rows),
        "max_slope_rel_error_all": max(row["max_slope_rel_error"] for row in rows),
        "max_closure_growth": max(row["closure_residual_growth"] for row in rows),
        "max_omega_drift_ratio": max(row["omega_drift_ratio"] for row in rows),
    }
    with (RESULT_DIR / "velocity_feedback_preliminary_result_v1.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print(json.dumps(
        {k: payload[k] for k in [
            "max_kinematic_error", "max_slope_rel_error_all",
            "max_closure_growth", "max_omega_drift_ratio"]},
        indent=2))
    for row in rows:
        print(f"chi0={row['chi0_deg']:>4}deg kappa={row['kappa']:>5} "
              f"kin_err={row['kinematic_error']:.2e} "
              f"slope_err={row['max_slope_rel_error']:.2e} "
              f"omega_drift={row['omega_drift_ratio']:.3e} "
              f"chi_growth={row['chi_envelope_growth']:.6f} "
              f"closure_growth={row['closure_residual_growth']:.6f}")


if __name__ == "__main__":
    main()
