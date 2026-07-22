#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分裂後帰結の分類 予備実験 v1（実験O4）

問い: 自発的分裂の開始後、系の帰結は分類できるか。
      候補クラス: 拡大（新回転平面の定着）/ 中間（有界混合振動）/
      回帰（親平面近傍への反復回帰＝縮小エピソード）。

設計:
  - N = 3,4,5,6（M = 3,6,10,15）。拡大には平面数の余地が必要なので N を振る。
  - 初期条件2族:
      parent 族: 自己無撞着固有モード親 + 相対種 δ=1e-3（親平面の直交補空間へ複素種）
      generic 族: 一様ランダム複素状態（親構造なし。別アトラクタの探索）
  - 力学: 逐次再構成正弦生成子（凍結解除）、実直交 Cayley 更新。
  - 分類計量（後期窓 = 後半）:
      f(τ) = 1 - h_親平面/h_total の後期 min/mean/max と傾き
      σ2(τ): 瞬時生成子第2特異方向。σ2 > 0.05（相対閾値）の滞在率と最長連続
      PR(τ): 辺振幅の関与比 (Σ|Z_e|²)²/Σ|Z_e|⁴ ∈ [1, M]（配分の広がり、比率量）
  - 分類規則（予備）:
      expansion    : σ2 活性滞在率 > 0.5（第2平面が定着）
      recurrent    : f_late_min < 0.05（親近傍への深い回帰あり）
      intermediate : それ以外（有界混合振動）

出力: outcome_classification_result_v1/ に JSON・図・分類表。
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
RESULT_DIR = os.path.join(BASE_DIR, "outcome_classification_result_v1")

GAMMA = base.GAMMA
N_LIST = [3, 4, 5, 6]
SEEDS = list(range(6))
DELTA = 1e-3
STEPS = 20000
RECORD_EVERY = 2
SIGMA2_TOL = 0.05
F_RETURN = 0.05


def participation_ratio(Z):
    a2 = np.abs(Z) ** 2
    s2 = float(np.sum(a2))
    s4 = float(np.sum(a2 * a2))
    return s2 * s2 / s4 if s4 > 0 else 1.0


def prepare_parent_perturbed(rng, A, m, delta):
    """自己無撞着固有モード親 + 親平面の直交補空間への零閉鎖種 g=(u+iw)/√2。

    u, w は親平面 {p,q} の実直交補内の実正規直交対。g^T g = 0 厳密、v^T g = 0。"""
    v, residual = base.prepare_initial_state(rng, A, m, 0.0)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    def proj_out(x):
        return x - (x @ p) * p - (x @ q) * q

    u = proj_out(rng.normal(size=m))
    u = u / np.linalg.norm(u)
    w = proj_out(rng.normal(size=m))
    w = w - (w @ u) * u
    nw = np.linalg.norm(w)
    if nw > 1e-8:
        w = w / nw
        g = (u + 1j * w) / np.sqrt(2.0)  # g^T g = 0 厳密
    else:
        # 補空間が1次元（N=3）: 零閉鎖種は存在しない（等方ベクトルに2実次元必要）。
        # 実種 g=u を用い、閉鎖残差 O(δ²) を明示的に許容する。
        g = u.astype(complex)
    Z0 = v + delta * g
    return Z0 / np.linalg.norm(Z0), residual


def run_and_measure(Z0, A, steps, record_every):
    """逐次力学を走らせ、分類計量の時系列（サブサンプル）と保存量偏差を返す。"""
    m = A.shape[0]
    K_ref = base.sine_generator(np.angle(Z0), A)
    planes_ref, _ = base.plane_decomposition(K_ref)
    pl0 = planes_ref[0]
    Z = Z0.copy()
    c0 = Z @ Z
    n_rec = steps // record_every + 1
    f_hist = np.empty(n_rec)
    s2_hist = np.empty(n_rec)
    pr_hist = np.empty(n_rec)
    max_dev_c = 0.0
    r = 0
    for t in range(steps + 1):
        K_now = base.sine_generator(np.angle(Z), A)
        if t % record_every == 0:
            a = pl0["p"] @ Z
            b = pl0["q"] @ Z
            h1 = abs(a) ** 2 + abs(b) ** 2
            h_tot = float(np.real(np.conj(Z) @ Z))
            f_hist[r] = 1.0 - h1 / h_tot
            mu = np.linalg.eigvalsh(1j * K_now)
            s2_hist[r] = float(np.sort(mu)[-2]) if m >= 2 else 0.0
            pr_hist[r] = participation_ratio(Z)
            c = Z @ Z
            max_dev_c = max(max_dev_c, abs(c - c0))
            r += 1
        if t < steps:
            Z = base.cayley(K_now) @ Z
    return f_hist, s2_hist, pr_hist, float(max_dev_c)


def classify(f, s2, family):
    """後期窓の計量から帰結クラスを判定。"""
    n = len(f)
    late = slice(n // 2, n)
    f_late = f[late]
    s2_late = s2[late]
    active = s2_late > SIGMA2_TOL
    active_frac = float(np.mean(active))
    # 最長連続活性（記録点数）
    longest = 0
    cur = 0
    for flag in active:
        cur = cur + 1 if flag else 0
        longest = max(longest, cur)
    taus = np.arange(n // 2, n) * RECORD_EVERY
    slope = float(np.polyfit(taus, f_late, 1)[0]) if len(f_late) > 8 else 0.0
    metrics = {
        "f_late_min": float(np.min(f_late)),
        "f_late_mean": float(np.mean(f_late)),
        "f_late_max": float(np.max(f_late)),
        "f_late_slope_per_step": slope,
        "sigma2_active_frac": active_frac,
        "sigma2_longest_run_steps": int(longest * RECORD_EVERY),
        "sigma2_late_max": float(np.max(s2_late)),
    }
    if active_frac > 0.5:
        cls = "expansion"
    elif metrics["f_late_min"] < F_RETURN:
        cls = "recurrent"
    else:
        cls = "intermediate"
    metrics["class"] = cls
    return metrics


def main():
    os.makedirs(RESULT_DIR, exist_ok=True)
    summary = {
        "params": {
            "n_list": N_LIST,
            "seeds": SEEDS,
            "delta": DELTA,
            "steps": STEPS,
            "record_every": RECORD_EVERY,
            "gamma": GAMMA,
            "sigma2_tol": SIGMA2_TOL,
            "f_return": F_RETURN,
        },
        "runs": [],
    }
    rep_curves = {}

    for n_bodies in N_LIST:
        A = base.line_graph_adjacency(n_bodies)
        m = A.shape[0]
        for family in ["parent", "generic"]:
            for seed in SEEDS:
                rng = np.random.default_rng(30260721 + 1000 * n_bodies + seed)
                if family == "parent":
                    Z0, residual = prepare_parent_perturbed(rng, A, m, DELTA)
                else:
                    X = rng.normal(size=m)
                    Y = rng.normal(size=m)
                    Y = Y - (X @ Y) / (X @ X) * X
                    Y = Y * (np.linalg.norm(X) / np.linalg.norm(Y))
                    Z0 = (X + 1j * Y) / np.linalg.norm(X + 1j * Y)  # Z^T Z = 0 厳密
                    residual = None
                f, s2, pr, dev_c = run_and_measure(Z0, A, STEPS, RECORD_EVERY)
                metrics = classify(f, s2, family)
                entry = {
                    "n_bodies": n_bodies,
                    "m": m,
                    "family": family,
                    "seed": seed,
                    "abs_ztz_initial": float(abs(complex(Z0 @ Z0))),
                    "parent_residual": residual,
                    "f_initial": float(f[0]),
                    "pr_initial": float(pr[0]),
                    "pr_late_mean": float(np.mean(pr[len(pr) // 2 :])),
                    "pr_late_max": float(np.max(pr[len(pr) // 2 :])),
                    "max_dev_closure": dev_c,
                    **metrics,
                }
                summary["runs"].append(entry)
                if seed == 0:
                    rep_curves[(n_bodies, family)] = (f, s2, pr)

    with open(os.path.join(RESULT_DIR, "summary_v1.json"), "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    taus = np.arange(STEPS // RECORD_EVERY + 1) * RECORD_EVERY

    # ---- 図1: parent 族の f(τ)（N別代表） ----
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2), sharex=True)
    for ax, n_bodies in zip(axes.flat, N_LIST):
        f, s2, pr = rep_curves[(n_bodies, "parent")]
        ax.semilogy(taus, np.maximum(f, 1e-18), label=r"$f(\tau)$")
        ax.semilogy(taus, np.maximum(s2, 1e-18), alpha=0.7, label=r"$\sigma_2(\tau)$")
        ax.axhline(SIGMA2_TOL, color="k", ls=":", alpha=0.6)
        ax.axhline(F_RETURN, color="gray", ls="--", alpha=0.6)
        ax.set_title(rf"$N={n_bodies}$ ($M={n_bodies*(n_bodies-1)//2}$), parent family",
                     fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")
    for ax in axes[1]:
        ax.set_xlabel(r"$\tau$")
    fig.suptitle("Post-onset fate: dormant fraction and second singular direction",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "fate_curves_parent_family_v1.png"), dpi=160)
    plt.close(fig)

    # ---- 図2: 分類マップ（f_late_min vs σ2活性率） ----
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    colors = {3: "tab:blue", 4: "tab:orange", 5: "tab:green", 6: "tab:red"}
    markers = {"parent": "o", "generic": "^"}
    for run in summary["runs"]:
        ax.scatter(
            max(run["f_late_min"], 1e-6),
            run["sigma2_active_frac"],
            c=colors[run["n_bodies"]],
            marker=markers[run["family"]],
            alpha=0.75,
            s=46,
        )
    ax.axvline(F_RETURN, color="gray", ls="--", alpha=0.7, label="recurrent boundary")
    ax.axhline(0.5, color="k", ls=":", alpha=0.7, label="expansion boundary")
    ax.set_xscale("log")
    ax.set_xlabel(r"late-window $\min f$ (return depth toward parent)")
    ax.set_ylabel(r"late-window fraction with $\sigma_2 > 0.05$")
    ax.set_title("Outcome-class map (colors: N; circles: parent, triangles: generic)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "outcome_class_map_v1.png"), dpi=160)
    plt.close(fig)

    # ---- 図3: PR の後期平均（配分の広がり） ----
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for family, mk in markers.items():
        xs, ys = [], []
        for run in summary["runs"]:
            if run["family"] == family:
                xs.append(run["m"] + (0.15 if family == "generic" else -0.15))
                ys.append(run["pr_late_mean"])
        ax.scatter(xs, ys, marker=mk, alpha=0.75, label=family)
    ms = [n * (n - 1) // 2 for n in N_LIST]
    ax.plot(ms, ms, "k--", alpha=0.5, label="full spread (PR = M)")
    ax.set_xlabel(r"$M$ (number of relation waves)")
    ax.set_ylabel("late-window mean participation ratio")
    ax.set_title("Spread of amplitude allocation after onset")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "participation_ratio_v1.png"), dpi=160)
    plt.close(fig)

    # ---- コンソール分類表 ----
    print("=== 分裂後帰結の分類 予備実験 v1 ===")
    print(f"steps={STEPS}, delta={DELTA}, seeds={len(SEEDS)}/family")
    dev = max(r["max_dev_closure"] for r in summary["runs"])
    print(f"閉鎖保存の最大偏差（全run）: {dev:.3e}")
    header = f"{'N':>2} {'M':>3} {'family':>8} | {'expansion':>9} {'intermediate':>12} {'recurrent':>9} | " \
             f"{'f_min(med)':>10} {'σ2act(med)':>10} {'PR_late(med)':>12}"
    print(header)
    print("-" * len(header))
    for n_bodies in N_LIST:
        for family in ["parent", "generic"]:
            runs = [
                r
                for r in summary["runs"]
                if r["n_bodies"] == n_bodies and r["family"] == family
            ]
            cnt = {"expansion": 0, "intermediate": 0, "recurrent": 0}
            for r in runs:
                cnt[r["class"]] += 1
            print(
                f"{n_bodies:>2} {runs[0]['m']:>3} {family:>8} |"
                f" {cnt['expansion']:>9} {cnt['intermediate']:>12} {cnt['recurrent']:>9} |"
                f" {np.median([r['f_late_min'] for r in runs]):>10.3e}"
                f" {np.median([r['sigma2_active_frac'] for r in runs]):>10.3f}"
                f" {np.median([r['pr_late_mean'] for r in runs]):>12.2f}"
            )
    res = [
        r["parent_residual"]
        for r in summary["runs"]
        if r["parent_residual"] is not None
    ]
    print(f"親状態の固有モード残差: median={np.median(res):.3e} max={max(res):.3e}")
    print(f"結果出力先: {RESULT_DIR}")


if __name__ == "__main__":
    main()
