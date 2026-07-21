"""速度フィードバック予備実験の論文用図生成 v1"""

from __future__ import annotations

import math
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
RESULT_DIR = BASE_DIR / "velocity_feedback_preliminary_result_v1"
MPL_DIR = RESULT_DIR / ".matplotlib"
MPL_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import run_velocity_feedback_preliminary_v1 as exp

SERIES_BLUE = "#2a78d6"
SERIES_GREEN = "#008300"
SERIES_MAGENTA = "#e87ba4"
NEUTRAL_GRAY = "#52514e"

OMEGA_0 = exp.OMEGA_0
PERIOD = exp.PERIOD_STEPS


def figure_v_equals_a_tau() -> None:
    control = exp.simulate(10.0, 0.0)["series"]
    feedback = exp.simulate(10.0, 0.2)["series"]
    steps = np.arange(len(feedback["omega"]))
    show = 5 * PERIOD

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(steps[:show], control["omega"][:show] / OMEGA_0, color=NEUTRAL_GRAY,
            linewidth=2, label=r"zero order (v1, $\kappa=0$): $\omega$ frozen")
    ax.plot(steps[:show], feedback["omega"][:show] / OMEGA_0, color=SERIES_BLUE,
            linewidth=2, label=r"first order ($\kappa=0.2$): $\omega(\tau)$")
    overlay = 1.0 + 0.2 * feedback["cumulative_accel"][:show] / OMEGA_0
    ax.plot(steps[:show], overlay, color=SERIES_MAGENTA, linewidth=1.2,
            linestyle="--",
            label=r"$1+\kappa\int a\,d\tau\ /\ \omega_0$ (kinematic overlay)")
    ax.set_xlabel(r"step $\tau$")
    ax.set_ylabel(r"carrier angular velocity  $\omega(\tau)/\omega_0$")
    ax.set_title(
        r"$v=\int a\,d\tau$: readout acceleration accumulates into carrier velocity"
    )
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "v_equals_a_tau_v1.png", dpi=180)
    plt.close(fig)


def figure_closure_stability() -> None:
    control = exp.simulate(10.0, 0.0)["series"]
    feedback = exp.simulate(10.0, 0.5)["series"]
    residual_c = control["closure_residual"]
    residual_f = feedback["closure_residual"]
    # 周期ごとの包絡（最大値）
    count = len(residual_c) // PERIOD
    env_c = residual_c[:count * PERIOD].reshape(count, PERIOD).max(axis=1)
    env_f = residual_f[:count * PERIOD].reshape(count, PERIOD).max(axis=1)
    cycles = np.arange(count)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(cycles, env_c, color=NEUTRAL_GRAY, linewidth=2,
            label=r"zero order ($\kappa=0$)")
    ax.plot(cycles, env_f, color=SERIES_BLUE, linewidth=2,
            label=r"first order ($\kappa=0.5$, strongest feedback)")
    ax.set_xlabel("cycle index (200 cycles = 19200 steps)")
    ax.set_ylabel(r"closure-residual envelope  $\max_\mathrm{cycle}|x_1^2+x_2^2|$")
    ax.set_title("Closure residual stays bounded under velocity feedback")
    ax.set_ylim(0, float(max(env_c.max(), env_f.max())) * 1.25)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower left")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "closure_residual_stability_v1.png", dpi=180)
    plt.close(fig)


def figure_readout_validity() -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for kappa, color in ((0.0, NEUTRAL_GRAY), (0.2, SERIES_BLUE),
                         (0.5, SERIES_GREEN)):
        series = exp.simulate(10.0, kappa)["series"]
        chi = series["chi"]
        omega = series["omega"]
        second = chi[2:] - 2.0 * chi[1:-1] + chi[:-2]
        chi_mid = chi[1:-1]
        omega_mid = omega[1:-1]
        ratios = []
        for start in range(0, len(chi_mid) - PERIOD, PERIOD):
            cw = chi_mid[start:start + PERIOD]
            dw = second[start:start + PERIOD]
            ow = omega_mid[start:start + PERIOD]
            denom = float(np.dot(cw, cw))
            if denom < 1.0e-24:
                continue
            slope = float(np.dot(cw, dw)) / denom
            theory = -4.0 * math.sin(float(np.mean(ow)) / 2.0) ** 2
            ratios.append(slope / theory)
        ax.plot(np.arange(len(ratios)), ratios, color=color, linewidth=1.5,
                label=rf"$\kappa={kappa}$")
    ax.axhline(1.0, color=NEUTRAL_GRAY, linestyle="--", linewidth=1.0)
    ax.set_xlabel("window index (window = 96 steps)")
    ax.set_ylabel(r"regression slope / theory  $-4\sin^2(\bar\omega/2)$")
    ax.set_title("Harmonic acceleration readout survives velocity feedback")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "readout_validity_under_feedback_v1.png", dpi=180)
    plt.close(fig)


def main() -> None:
    figure_v_equals_a_tau()
    figure_closure_stability()
    figure_readout_validity()
    for name in ("v_equals_a_tau_v1.png", "closure_residual_stability_v1.png",
                 "readout_validity_under_feedback_v1.png"):
        print(RESULT_DIR / name)


if __name__ == "__main__":
    main()
