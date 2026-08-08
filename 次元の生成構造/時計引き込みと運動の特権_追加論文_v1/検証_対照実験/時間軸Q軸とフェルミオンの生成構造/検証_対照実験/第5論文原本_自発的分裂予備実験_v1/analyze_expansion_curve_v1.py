#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
拡大クラスのカーブ形状の定量化（実験O5予備）

N=5,6 の parent 族（拡大クラス）について f(τ), σ2(τ), PR(τ) を密に記録し、
  - 初期指数成長率（f ∈ [1e-5, 1e-2] 窓の log 傾き）
  - 立ち上がり時間（f が飽和値の10%→90%に達する τ）
  - 飽和値と等分配予言 f_eq = 1 - 1/M, PR_eq = M の比較
  - PR の緩和（M - PR(τ) の指数減衰時間）
を測定する。カーブが「指数→ロジスティック型飽和→等分配プラトー」かを判定する。
"""

import json
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import run_spontaneous_splitting_preliminary_v1 as base
import run_outcome_classification_preliminary_v1 as oc

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE_DIR, "outcome_classification_result_v1")
STEPS = 20000
RECORD_EVERY = 2


def analyze(n_bodies, seed):
    A = base.line_graph_adjacency(n_bodies)
    m = A.shape[0]
    rng = np.random.default_rng(30260721 + 1000 * n_bodies + seed)
    Z0, residual = oc.prepare_parent_perturbed(rng, A, m, oc.DELTA)
    f, s2, pr, dev_c = oc.run_and_measure(Z0, A, STEPS, RECORD_EVERY)
    taus = np.arange(len(f)) * RECORD_EVERY

    # 初期指数成長率
    mask = (f > 1e-5) & (f < 1e-2)
    idx = np.where(mask)[0]
    early_rate = None
    if len(idx) >= 5:
        early_rate = float(np.polyfit(taus[idx], np.log(f[idx]), 1)[0])

    # プラトー（後半25%）
    tail = slice(3 * len(f) // 4, len(f))
    f_plateau = float(np.mean(f[tail]))
    f_plateau_std = float(np.std(f[tail]))
    pr_plateau = float(np.mean(pr[tail]))
    s2_plateau = float(np.mean(s2[tail]))
    f_eq = 1.0 - 1.0 / m

    # 立ち上がり 10%→90%（プラトー値基準、初回到達）
    t10 = taus[np.argmax(f >= 0.1 * f_plateau)]
    t90 = taus[np.argmax(f >= 0.9 * f_plateau)]

    # PR 緩和: M - PR の指数減衰時間（PR が M の 10%→ 1% 残差の区間）
    resid_pr = np.maximum(m - pr, 1e-12)
    mask2 = (resid_pr < 0.5 * (m - 1)) & (resid_pr > 1e-3)
    idx2 = np.where(mask2)[0]
    pr_decay_rate = None
    if len(idx2) >= 5:
        pr_decay_rate = float(np.polyfit(taus[idx2], np.log(resid_pr[idx2]), 1)[0])

    return {
        "n_bodies": n_bodies,
        "m": m,
        "seed": seed,
        "parent_residual": residual,
        "early_rate_per_step": early_rate,
        "rise_t10": int(t10),
        "rise_t90": int(t90),
        "f_plateau_mean": f_plateau,
        "f_plateau_std": f_plateau_std,
        "f_equipartition_pred": f_eq,
        "pr_plateau_mean": pr_plateau,
        "pr_equipartition_pred": float(m),
        "sigma2_plateau_mean": s2_plateau,
        "pr_relaxation_rate_per_step": pr_decay_rate,
        "max_dev_closure": dev_c,
    }, (taus, f, s2, pr)


def main():
    results = []
    curves = {}
    for n_bodies in [5, 6]:
        for seed in [0, 1]:
            entry, cv = analyze(n_bodies, seed)
            results.append(entry)
            curves[(n_bodies, seed)] = cv
            print(
                f"N={n_bodies} seed={seed}: early_rate={entry['early_rate_per_step']}"
                f" rise {entry['rise_t10']}→{entry['rise_t90']}"
                f" f_plateau={entry['f_plateau_mean']:.4f}±{entry['f_plateau_std']:.4f}"
                f" (eq pred {entry['f_equipartition_pred']:.4f})"
                f" PR={entry['pr_plateau_mean']:.3f}/(pred {entry['m']})"
                f" σ2={entry['sigma2_plateau_mean']:.3f}"
                f" PR_relax={entry['pr_relaxation_rate_per_step']}"
            )

    with open(os.path.join(RESULT_DIR, "expansion_curve_analysis_v1.json"), "w") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ax = axes[0]
    for (n_bodies, seed), (taus, f, s2, pr) in curves.items():
        if seed != 0:
            continue
        ax.semilogy(taus, np.maximum(f, 1e-8), label=rf"$N={n_bodies}$: $f(\tau)$")
        ax.axhline(1.0 - 1.0 / (n_bodies * (n_bodies - 1) // 2), ls=":", alpha=0.5)
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$f(\tau)$")
    ax.set_title("Expansion curve: exponential onset, logistic-type saturation",
                 fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    for (n_bodies, seed), (taus, f, s2, pr) in curves.items():
        if seed != 0:
            continue
        m = n_bodies * (n_bodies - 1) // 2
        ax.semilogy(taus, np.maximum(m - pr, 1e-12),
                    label=rf"$N={n_bodies}$: $M-\mathrm{{PR}}(\tau)$")
    ax.set_xlabel(r"$\tau$")
    ax.set_ylabel(r"$M-\mathrm{PR}$ (distance to equipartition)")
    ax.set_title("Relaxation toward exact equipartition", fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "expansion_curve_analysis_v1.png"), dpi=160)
    plt.close(fig)
    print("saved: expansion_curve_analysis_v1.{json,png}")


if __name__ == "__main__":
    main()
