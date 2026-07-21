"""速度フィードバック下の逆二乗則再検証と包絡減少の判別実験 v1

実験K: 倍音 n=1..8 で調和閉鎖 omega_n = n·omega_1 を初期値に取り、
       フィードバック動力学の実効係数から距離指数を測る。
       対照 kappa=0 では代数的に -2。フィードバック下での指数を測定。

実験L: 包絡減少（kappa=0.5, chi0=20°で約13%/200周期）の判別。
       同じ連立系 χ̈ = -4sin²(ω/2)χ, ω̇ = κ·χ̈ を RK4 で刻み幅
       dt = 1, 0.1, 0.01 と細分して積分し、包絡減少率が dt に依存する
       （離散化artifact）か、dt→0 で残る（実力学）かを判別する。
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent
RESULT_DIR = BASE_DIR / "velocity_feedback_preliminary_result_v1"

PERIOD_STEPS = 96
OMEGA_1 = 2.0 * math.pi / PERIOD_STEPS
STEP_COUNT = 96 * 200
CHI0 = math.radians(10.0)


def run_map(omega0: float, kappa: float, chi0: float, steps: int):
    omega = np.empty(steps + 1)
    chi = np.empty(steps + 1)
    omega[0] = omega0
    chi[0] = chi0
    a0 = -4.0 * math.sin(omega[0] / 2.0) ** 2 * chi[0]
    chi[1] = chi[0] + 0.5 * a0
    omega[1] = omega[0] + kappa * a0
    for s in range(1, steps):
        a = -4.0 * math.sin(omega[s] / 2.0) ** 2 * chi[s]
        chi[s + 1] = 2.0 * chi[s] - chi[s - 1] + a
        omega[s + 1] = omega[s] + kappa * a
    return omega, chi


def run_map_integrated(omega0: float, kappa: float, chi0: float, steps: int):
    """積分形: omega_s = omega0 + kappa*(chi_s - chi_{s-1})。
    連続系の厳密積分 omega = omega0 + kappa*chi_dot の離散版。"""
    chi = np.empty(steps + 1)
    omega = np.empty(steps + 1)
    chi[0] = chi0
    omega[0] = omega0
    a0 = -4.0 * math.sin(omega0 / 2.0) ** 2 * chi0
    chi[1] = chi0 + 0.5 * a0
    for s in range(1, steps):
        omega[s] = omega0 + kappa * (chi[s] - chi[s - 1])
        a = -4.0 * math.sin(omega[s] / 2.0) ** 2 * chi[s]
        chi[s + 1] = 2.0 * chi[s] - chi[s - 1] + a
    omega[steps] = omega0 + kappa * (chi[steps] - chi[steps - 1])
    return omega, chi


# ---------------------------------------------------------------- 実験K
def experiment_k() -> dict:
    rows = []
    for kappa, runner, scheme in (
        (0.0, run_map, "euler"),
        (0.2, run_map, "euler"),
        (0.2, run_map_integrated, "integrated"),
    ):
        log_dtheta = []
        log_alpha = []
        for n in range(1, 9):
            omega0 = n * OMEGA_1
            omega, chi = runner(omega0, kappa, CHI0, STEP_COUNT)
            second = chi[2:] - 2.0 * chi[1:-1] + chi[:-2]
            chi_mid = chi[1:-1]
            slope = float(np.dot(chi_mid, second) / np.dot(chi_mid, chi_mid))
            omega_eff = 2.0 * math.asin(min(1.0, math.sqrt(-slope) / 2.0))
            alpha_eff = omega_eff ** 2  # R=1
            dtheta = 2.0 * math.pi / n
            log_dtheta.append(math.log(dtheta))
            log_alpha.append(math.log(alpha_eff))
            rows.append({
                "kappa": kappa, "scheme": scheme, "n": n, "delta_theta": dtheta,
                "omega_eff": omega_eff, "alpha_eff": alpha_eff,
            })
        x = np.array(log_dtheta)
        y = np.array(log_alpha)
        exponent = float(np.polyfit(x, y, 1)[0])
        for row in rows:
            if row["kappa"] == kappa and row["scheme"] == scheme and "exponent" not in row:
                row["exponent"] = exponent
    exponents = {
        f"kappa={row['kappa']},{row['scheme']}": row["exponent"]
        for row in rows if "exponent" in row
    }
    return {"rows": rows, "exponents": exponents}


# ---------------------------------------------------------------- 実験L
def rk4_envelope(kappa: float, chi0: float, dt: float, total_tau: float) -> dict:
    def force(omega, chi):
        return -4.0 * math.sin(omega / 2.0) ** 2 * chi

    def deriv(state):
        chi, u, omega = state
        a = force(omega, chi)
        return np.array([u, a, kappa * a])

    steps = int(round(total_tau / dt))
    state = np.array([chi0, 0.0, OMEGA_1])
    env_first = 0.0
    env_second = 0.0
    half = steps // 2
    for i in range(steps):
        k1 = deriv(state)
        k2 = deriv(state + 0.5 * dt * k1)
        k3 = deriv(state + 0.5 * dt * k2)
        k4 = deriv(state + dt * k3)
        state = state + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        magnitude = abs(state[0])
        if i < half:
            env_first = max(env_first, magnitude)
        else:
            env_second = max(env_second, magnitude)
    return {
        "dt": dt,
        "chi_envelope_first_half": env_first,
        "chi_envelope_second_half": env_second,
        "envelope_growth": env_second / env_first,
    }


def experiment_l() -> dict:
    kappa = 0.5
    chi0 = math.radians(20.0)
    total_tau = float(STEP_COUNT)
    # 離散写像（本実験の動力学）
    omega, chi = run_map(OMEGA_1, kappa, chi0, STEP_COUNT)
    half = STEP_COUNT // 2
    map_growth = float(
        np.max(np.abs(chi[half:])) / np.max(np.abs(chi[:half]))
    )
    rows = [{"dt": 1.0, "scheme": "discrete map (leapfrog + Euler omega)",
             "envelope_growth": map_growth}]
    omega_i, chi_i = run_map_integrated(OMEGA_1, kappa, chi0, STEP_COUNT)
    integrated_growth = float(
        np.max(np.abs(chi_i[half:])) / np.max(np.abs(chi_i[:half]))
    )
    rows.append({"dt": 1.0, "scheme": "discrete map (integrated omega form)",
                 "envelope_growth": integrated_growth})
    for dt in (1.0, 0.1, 0.01):
        result = rk4_envelope(kappa, chi0, dt, total_tau)
        result["scheme"] = "RK4 continuous"
        rows.append(result)
    return {"rows": rows}


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    result_k = experiment_k()
    result_l = experiment_l()
    payload = {"experiment_k_inverse_square": result_k,
               "experiment_l_discrimination": result_l}
    with (RESULT_DIR / "feedback_inverse_square_discrimination_v1.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    print("=== K: distance exponent ===")
    for kappa, exponent in result_k["exponents"].items():
        print(f"kappa={kappa}: exponent={exponent:.9f}")
    print("=== L: envelope growth (kappa=0.5, chi0=20deg, 200 cycles) ===")
    for row in result_l["rows"]:
        print(f"{row['scheme']:>36}  dt={row['dt']:<5} growth={row['envelope_growth']:.6f}")


if __name__ == "__main__":
    main()
