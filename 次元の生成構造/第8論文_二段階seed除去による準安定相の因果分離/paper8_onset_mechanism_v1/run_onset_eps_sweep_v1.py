#!/usr/bin/env python3
"""ONS-1: 摂動振幅掃引——潜伏＝線形不安定性の対数則の直接検証

目的（2026-08-05 木原氏合意・開始様式判別論文の機構検証その1）:
    潜伏バーストが相対平衡の線形不安定性であるなら、初期摂動振幅 ε に対して
    crossing 時刻は対数則に従うはずである。

予言（実行前固定・事後変更禁止）:
    f(0) ≈ ε² が f_cross=0.05 まで f ∝ e^{rate_f·t} で成長するなら
        t_cross(ε) ≈ const + (−2 ln ε)/rate_f
    すなわち t_cross 対 (−ln ε) は直線で、傾き = 2/rate_f。
    参照 rate_f はプロファイル測定（paper8_em9r_profile_*_v1.json）の
    control 窓フィット: N=5→0.04935/step、N=40→0.03499/step。
    予言傾き: N=5→40.5 step/e-fold、N=40→57.2 step/e-fold。
    さらに各 run の窓フィット rate は ε に依存しないはずである。

方法:
    control 親（build_init(n, False) の v、無seed）に対し、親平面 {p,q} に
    直交する複素ランダム方向 η（‖η‖=1、rng 93000+j、j=0,1,2。η は ε に
    依らず共通）を ε 倍で加え規格化: Z0=(v+εη)/‖·‖。
    ε ∈ {1e-4,1e-6,1e-8,1e-10,1e-12,1e-14}。f(t) を毎 step 記録し、
    crossing（f>0.05 初回）まで、または T_MAX=6000 まで発展。
    窓フィット: f ∈ [max(100·f(0),1e-20), 1e-2] で ln f 対 t 最小二乗。

再現性: 第7論文力学 read-only import（SHA-256 記録）。乱数は η 用
    rng(93000+j) のみ。全 t_cross・フィット・傾き回帰を JSON 保存、図出力。

使い方: python3 run_onset_eps_sweep_v1.py <N>
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
EPS_LIST = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14]
ETA_SEEDS = [0, 1, 2]
T_MAX = 6000
F_CROSS = 0.05
FIT_HI = 1e-2
RATE_REF = {5: 0.04935428519548658, 40: 0.034990401406084726}  # profile control 窓フィット

spec = importlib.util.spec_from_file_location("abl_ons1", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run_one(n, sys_lr, v, p, q, eta, eps, wp0):
    Z = v + eps * eta
    Z = Z / np.linalg.norm(Z)
    wp = wp0.copy()

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    fs = [fval(Z)]
    t_cross = None
    for t in range(1, T_MAX + 1):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        f = fval(Z)
        fs.append(f)
        if t_cross is None and f > F_CROSS:
            t_cross = t
            break
    fs = np.array(fs)
    f0 = fs[0]
    lo = max(100.0 * f0, 1e-20)
    win = [(t, f) for t, f in enumerate(fs) if lo <= f <= FIT_HI]
    fit = {"n_points": len(win), "rate_per_step": None, "r2": None, "decades": 0.0}
    if len(win) >= 5:
        tt = np.array([w[0] for w in win], float)
        lf = np.log(np.array([w[1] for w in win]))
        A = np.vstack([tt, np.ones_like(tt)]).T
        coef, _, _, _ = np.linalg.lstsq(A, lf, rcond=None)
        pred = A @ coef
        ss_res = float(np.sum((lf - pred) ** 2))
        ss_tot = float(np.sum((lf - lf.mean()) ** 2))
        fit = {"n_points": len(win),
               "rate_per_step": float(coef[0]),
               "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0,
               "decades": float((lf.max() - lf.min()) / np.log(10))}
    return {"eps": eps, "f0": float(f0), "t_cross": t_cross, "fit": fit,
            "f_sampled": [[int(t), float(fs[t])] for t in range(0, len(fs), 10)]}


def main() -> None:
    t0 = time.time()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    print(f"ONS-1 ε掃引 N={n}  ABL {sha256(ABL)[:16]}…")
    sys_lr, v, _, _, _, p, q, _, wp0 = abl.build_init(n, False)
    m = sys_lr.m

    runs = []
    for j in ETA_SEEDS:
        rng = np.random.default_rng(93000 + j)
        xi = rng.normal(size=m) + 1j * rng.normal(size=m)
        eta = xi - p * (p @ xi) - q * (q @ xi)
        eta = eta / np.linalg.norm(eta)
        for eps in EPS_LIST:
            r = run_one(n, sys_lr, v, p, q, eta, eps, wp0)
            r["eta_seed"] = j
            runs.append(r)
            print(f"  j={j} ε=1e{int(np.log10(eps)):+d}: t_cross={r['t_cross']} "
                  f"rate={r['fit']['rate_per_step']} R²={r['fit']['r2']}")

    # 傾き回帰: t_cross 対 −ln ε（全 run プール＋シード別）
    def slope_fit(pairs):
        if len(pairs) < 2:
            return None
        x = np.array([-np.log(e) for e, t in pairs])
        y = np.array([t for e, t in pairs], float)
        A = np.vstack([x, np.ones_like(x)]).T
        coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        pred = A @ coef
        ss_res = float(np.sum((y - pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        return {"slope": float(coef[0]), "intercept": float(coef[1]),
                "r2": 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0,
                "n": len(pairs)}

    pooled = slope_fit([(r["eps"], r["t_cross"]) for r in runs if r["t_cross"]])
    per_seed = {str(j): slope_fit([(r["eps"], r["t_cross"]) for r in runs
                                   if r["t_cross"] and r["eta_seed"] == j])
                for j in ETA_SEEDS}
    pred_slope = 2.0 / RATE_REF[n]
    rates = [r["fit"]["rate_per_step"] for r in runs if r["fit"]["rate_per_step"]]
    out = {"N": n, "imports": {"abl": sha256(ABL)},
           "criteria": {"EPS_LIST": EPS_LIST, "ETA_SEEDS": ETA_SEEDS, "T_MAX": T_MAX,
                         "F_CROSS": F_CROSS, "FIT_HI": FIT_HI,
                         "prediction": {"slope_vs_neg_ln_eps": pred_slope,
                                         "rate_ref": RATE_REF[n]}},
           "runs": runs,
           "slope_pooled": pooled, "slope_per_seed": per_seed,
           "rate_stats": {"mean": float(np.mean(rates)) if rates else None,
                           "std": float(np.std(rates)) if rates else None,
                           "min": float(np.min(rates)) if rates else None,
                           "max": float(np.max(rates)) if rates else None},
           "runtime_sec": time.time() - t0}
    (HERE / f"onset_eps_sweep_N{n:05d}_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")

    fig, ax = plt.subplots(figsize=(7, 5))
    marks = {0: "o", 1: "s", 2: "^"}
    for j in ETA_SEEDS:
        xs = [-np.log(r["eps"]) for r in runs if r["t_cross"] and r["eta_seed"] == j]
        ys = [r["t_cross"] for r in runs if r["t_cross"] and r["eta_seed"] == j]
        ax.plot(xs, ys, marks[j], ms=6, label=f"η seed {j}")
    if pooled:
        xx = np.linspace(min(-np.log(e) for e in EPS_LIST) - 2,
                         max(-np.log(e) for e in EPS_LIST) + 2, 10)
        ax.plot(xx, pooled["slope"] * xx + pooled["intercept"], "k-", lw=1,
                label=f"fit slope={pooled['slope']:.2f}")
        ax.plot(xx, pred_slope * xx + pooled["intercept"]
                + (pooled["slope"] - pred_slope) * float(np.mean(xx)), "r--", lw=1,
                label=f"predicted slope 2/rate_f={pred_slope:.2f}")
    ax.set_xlabel("−ln ε")
    ax.set_ylabel("t_cross (step)")
    ax.set_title(f"N={n} ONS-1: crossing time vs perturbation amplitude")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / f"fig_onset_eps_sweep_N{n:05d}.png", dpi=130)
    plt.close(fig)
    print(f"pooled slope={pooled['slope']:.2f} (R²={pooled['r2']:.5f}) "
          f"predicted={pred_slope:.2f}  rate mean={out['rate_stats']['mean']}")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
