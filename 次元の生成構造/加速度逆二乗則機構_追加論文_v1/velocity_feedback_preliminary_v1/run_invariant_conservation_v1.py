"""実験M: 速度フィードバック連続系の保存量 𝓘 の検証 v1

連続系  χ̈ = -g(ω)χ,  ω̇ = κχ̈,  g(ω) = 4 sin²(ω/2)  に対し、

  第一積分   J1 = ω - ω₀ - κχ̇
  保存量     𝓘  = χ²/2 + (1/κ²)∫_{ω₀}^{ω} (s-ω₀)/g(s) ds

の保存を検証する。積分項は閉形式の原始関数

  F(s) = ln sin(s/2) - ((s-ω₀)/2)·cot(s/2)

で厳密に評価する（dF/ds = (s-ω₀)/(4sin²(s/2)) を解析的に確認済み）。

検査:
  M1  RK4連続積分で J1・𝓘 が保存され、ドリフトが dt→0 で消えること
  M2  離散写像（Euler形・積分形）の 𝓘 ドリフトが包絡減少と対応すること
      （離散化アーティファクトの保存量による定量化）
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent
RESULT_DIR = BASE_DIR / "velocity_feedback_preliminary_result_v1"

PERIOD_STEPS = 96
OMEGA_0 = 2.0 * math.pi / PERIOD_STEPS
TOTAL_TAU = 96.0 * 200.0


def g_of(omega: float) -> float:
    return 4.0 * math.sin(omega / 2.0) ** 2


def antiderivative(s: float, omega0: float) -> float:
    """F(s) = ln sin(s/2) - ((s-omega0)/2) cot(s/2)"""
    half = s / 2.0
    return math.log(math.sin(half)) - (s - omega0) / 2.0 / math.tan(half)


def invariant(chi: float, chi_dot: float, omega: float, kappa: float) -> float:
    integral = antiderivative(omega, OMEGA_0) - antiderivative(OMEGA_0, OMEGA_0)
    return chi * chi / 2.0 + integral / (kappa * kappa)


def rk4_run(kappa: float, chi0: float, dt: float):
    def deriv(state):
        chi, u, omega = state
        a = -g_of(omega) * chi
        return np.array([u, a, kappa * a])

    steps = int(round(TOTAL_TAU / dt))
    state = np.array([chi0, 0.0, OMEGA_0])
    j1_0 = state[2] - OMEGA_0 - kappa * state[1]
    inv_0 = invariant(state[0], state[1], state[2], kappa)
    max_j1 = 0.0
    max_inv = 0.0
    sample = max(1, steps // 4000)
    for i in range(steps):
        k1 = deriv(state)
        k2 = deriv(state + 0.5 * dt * k1)
        k3 = deriv(state + 0.5 * dt * k2)
        k4 = deriv(state + dt * k3)
        state = state + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if i % sample == 0 or i == steps - 1:
            j1 = state[2] - OMEGA_0 - kappa * state[1]
            inv = invariant(state[0], state[1], state[2], kappa)
            max_j1 = max(max_j1, abs(j1 - j1_0))
            max_inv = max(max_inv, abs(inv - inv_0))
    return {"dt": dt, "J1_drift": max_j1, "I_drift": max_inv, "I0": inv_0}


def discrete_map_run(kappa: float, chi0: float, scheme: str):
    steps = int(TOTAL_TAU)
    chi = np.empty(steps + 1)
    omega = np.empty(steps + 1)
    chi[0] = chi0
    omega[0] = OMEGA_0
    a0 = -g_of(omega[0]) * chi[0]
    chi[1] = chi[0] + 0.5 * a0
    if scheme == "euler":
        omega[1] = omega[0] + kappa * a0
        for s in range(1, steps):
            a = -g_of(omega[s]) * chi[s]
            chi[s + 1] = 2.0 * chi[s] - chi[s - 1] + a
            omega[s + 1] = omega[s] + kappa * a
    else:  # integrated
        for s in range(1, steps):
            omega[s] = OMEGA_0 + kappa * (chi[s] - chi[s - 1])
            a = -g_of(omega[s]) * chi[s]
            chi[s + 1] = 2.0 * chi[s] - chi[s - 1] + a
        omega[steps] = OMEGA_0 + kappa * (chi[steps] - chi[steps - 1])
    # 中心差分による χ̇ 推定で 𝓘 を評価
    chi_dot = (chi[2:] - chi[:-2]) / 2.0
    inv_series = np.array([
        invariant(chi[s], chi_dot[s - 1], omega[s], kappa)
        for s in range(1, steps, 8)
    ])
    half = len(inv_series) // 2
    env_first = float(np.max(np.abs(chi[: steps // 2])))
    env_second = float(np.max(np.abs(chi[steps // 2:])))
    return {
        "scheme": scheme,
        "I_drift": float(np.max(np.abs(inv_series - inv_series[0]))),
        "I0": float(inv_series[0]),
        "envelope_growth": env_second / env_first,
    }


def main() -> None:
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    rows_rk4 = []
    rows_map = []
    for kappa in (0.05, 0.2, 0.5):
        for chi0_deg in (5.0, 10.0, 20.0):
            chi0 = math.radians(chi0_deg)
            for dt in (1.0, 0.1, 0.01):
                r = rk4_run(kappa, chi0, dt)
                r.update({"kappa": kappa, "chi0_deg": chi0_deg})
                rows_rk4.append(r)
            for scheme in ("euler", "integrated"):
                r = discrete_map_run(kappa, chi0, scheme)
                r.update({"kappa": kappa, "chi0_deg": chi0_deg})
                rows_map.append(r)

    payload = {
        "experiment": "M: invariant conservation",
        "rk4": rows_rk4,
        "discrete_maps": rows_map,
        "max_rk4_J1_drift": max(r["J1_drift"] for r in rows_rk4),
        "max_rk4_I_drift_dt001": max(
            r["I_drift"] for r in rows_rk4 if r["dt"] == 0.01
        ),
    }
    with (RESULT_DIR / "invariant_conservation_v1.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)

    print("=== RK4 (J1 drift / I drift, relative to I0) ===")
    for r in rows_rk4:
        print(f"kappa={r['kappa']:<5} chi0={r['chi0_deg']:>4} dt={r['dt']:<5} "
              f"J1={r['J1_drift']:.2e} I={r['I_drift']:.2e} I/I0={r['I_drift']/r['I0']:.2e}")
    print("=== discrete maps (I drift vs envelope) ===")
    for r in rows_map:
        print(f"kappa={r['kappa']:<5} chi0={r['chi0_deg']:>4} {r['scheme']:>10} "
              f"I_drift/I0={r['I_drift']/r['I0']:.3e} env_growth={r['envelope_growth']:.6f}")


if __name__ == "__main__":
    main()
