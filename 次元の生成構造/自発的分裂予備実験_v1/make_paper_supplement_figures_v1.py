#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
論文用補足図の生成:
  (1) δ=0（厳密周期軌道）からの丸め誤差起源の自発発生曲線
  (2) 成長率の δ 非依存性（summary_v1.json から）
"""

import json
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_spontaneous_splitting_preliminary_v1 as base

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE_DIR, "spontaneous_splitting_result_v1")
STEPS = 4320


def main():
    A = base.line_graph_adjacency(base.N_BODIES)
    m = A.shape[0]

    # ---- (1) δ=0 丸め誤差起源の発生 ----
    rng = np.random.default_rng(20260721)
    Z0, _ = base.prepare_initial_state(rng, A, m, 0.0)
    K_ref = base.sine_generator(np.angle(Z0), A)
    planes_ref, _ = base.plane_decomposition(K_ref)
    rec = base.run_dynamics(Z0, A, planes_ref, STEPS, sequential=True)
    f = base.dormant_fraction(rec)
    tau = np.arange(STEPS + 1)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.semilogy(tau, np.maximum(np.abs(f), 1e-36), color="tab:red")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"dormant fraction $f(\tau)$")
    ax.set_title(
        "Spontaneous onset from floating-point rounding noise alone\n"
        r"(exact periodic parent orbit, $\delta = 0$, no perturbation applied)",
        fontsize=10,
    )
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "zero_delta_machine_noise_onset_v1.png"), dpi=160)
    plt.close(fig)

    # ---- (2) 成長率の δ 非依存性 ----
    with open(os.path.join(RESULT_DIR, "summary_v1.json")) as fh:
        summary = json.load(fh)
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    deltas = summary["params"]["deltas"]
    for delta in deltas:
        rates = [
            r["sequential"]["early_growth_rate_per_step"]
            for r in summary["runs"]
            if r["delta"] == delta
            and r["sequential"]["early_growth_rate_per_step"] is not None
        ]
        ax.scatter([delta] * len(rates), rates, alpha=0.65, color="tab:blue", s=42)
    med_all = np.median(
        [
            r["sequential"]["early_growth_rate_per_step"]
            for r in summary["runs"]
            if r["sequential"]["early_growth_rate_per_step"] is not None
        ]
    )
    ax.axhline(med_all, color="k", ls="--", alpha=0.7,
               label=rf"median $= {med_all:.4f}$ per step")
    ax.set_xscale("log")
    ax.set_xlabel(r"relative seed amplitude $\delta$")
    ax.set_ylabel("early exponential growth rate of $f$ (per step)")
    ax.set_ylim(0.0, 0.08)
    ax.set_title(
        r"Growth rate is independent of seed amplitude over three decades of $\delta$",
        fontsize=10,
    )
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "growth_rate_vs_delta_v1.png"), dpi=160)
    plt.close(fig)

    print("saved: zero_delta_machine_noise_onset_v1.png, growth_rate_vs_delta_v1.png")
    print(f"median growth rate = {med_all:.6f} per step")


if __name__ == "__main__":
    main()
