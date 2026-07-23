#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
有理梯子検査（公理1.5候補 U^n=I の最初の審判）

問い: 動力学が選んだ状態（自己無撞着親・拡大定常状態）のσスペクトルは、
      有理比の梯子（共通位数 n: θ_j = 2πm_j/n）に乗っているか。

判定:
  (1) 比 σ_j/σ_1 の有理近似（分母≤20）と誤差——σ2/σ1 ≈ 1/2 の精密検査を含む
  (2) 共通位数検査: 回転角 θ_j = 2 arctan(γ σ_j/σ_1) について、
      err(n) = max_j dist(nθ_j/2π, Z) を n=2..N_MAX で最小化。
      同数のランダム角（帰無分布、モンテカルロ）の best-err と比較する。

対象: parent（自己無撞着親）/ stationary（零閉鎖一般状態を逐次力学で緩和させた終状態）
      / random（無進化の零閉鎖一般状態）——random が帰無側の実測対照。
"""

import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/自発的分裂予備実験_v1")
import run_n_scaling_lowrank_v1 as lr

GAMMA = lr.GAMMA
N_MAX_ORDER = 500
NULL_TRIALS = 200


def thetas_from_sigmas(sig):
    return 2.0 * np.arctan(GAMMA * sig / sig[0])


def ladder_err(thetas, n):
    x = thetas * n / (2.0 * math.pi)
    return float(np.max(np.abs(x - np.round(x))))


def best_order(thetas, n_max=N_MAX_ORDER):
    best_n, best_e = None, np.inf
    for n in range(2, n_max + 1):
        e = ladder_err(thetas, n)
        if e < best_e:
            best_n, best_e = n, e
    return best_n, best_e


def null_best_err(k, theta_max, rng, trials=NULL_TRIALS, n_max=N_MAX_ORDER):
    """同数 k 個の一様ランダム角に同じ検査を適用した best-err の分布。"""
    outs = []
    for _ in range(trials):
        th = rng.uniform(0.0, theta_max, k)
        outs.append(best_order(th, n_max)[1])
    return np.median(outs), np.percentile(outs, 5)


def rational_report(sig, k=8):
    rows = []
    for j in range(1, min(k, len(sig))):
        r = float(sig[j] / sig[0])
        fr = Fraction(r).limit_denominator(20)
        rows.append((j + 1, r, f"{fr.numerator}/{fr.denominator}", abs(r - float(fr))))
    return rows


def stationary_state(n, seed, steps=3000):
    sys_lr = lr.LowRankSystem(n)
    rng = np.random.default_rng(70260722 + seed)
    Z = lr.zero_closure_generic(rng, sys_lr.m)
    wp = rng.normal(size=sys_lr.m)
    for _ in range(steps):
        sys_lr.set_theta(np.angle(Z))
        sig_est, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, sig_est)
    sys_lr.set_theta(np.angle(Z))
    return sys_lr.sigma_spectrum()


def random_state_spectrum(n, seed):
    sys_lr = lr.LowRankSystem(n)
    rng = np.random.default_rng(80260722 + seed)
    Z = lr.zero_closure_generic(rng, sys_lr.m)
    sys_lr.set_theta(np.angle(Z))
    return sys_lr.sigma_spectrum()


def analyze(label, sig, rng_null):
    sig = sig[sig / sig[0] > 1e-6]
    th = thetas_from_sigmas(sig)
    n_best, e_best = best_order(th)
    null_med, null_p5 = null_best_err(len(th), float(th[0]), rng_null)
    print(f"  [{label}] モード数={len(sig)}")
    print(f"    σ_j/σ_1 (j=2..): " + ", ".join(
        f"{r:.6f}≈{fr}(err {e:.1e})" for _, r, fr, e in rational_report(sig, 6)))
    print(f"    共通位数検査: best n={n_best} err={e_best:.4f}"
          f" ｜帰無 best-err: median={null_med:.4f}, 5%点={null_p5:.4f}"
          f" ｜判定: {'梯子の証拠' if e_best < null_p5 * 0.5 else '帰無と区別できず'}")
    return {"label": label, "n_modes": int(len(sig)), "best_n": n_best,
            "best_err": e_best, "null_median": null_med, "null_p5": null_p5,
            "sigma2_over_sigma1": float(sig[1] / sig[0]) if len(sig) > 1 else None}


def main():
    rng_null = np.random.default_rng(12345)
    print("=== 有理梯子検査（U^n=I／公理1.5候補の一次審判） ===")
    print(f"γ = tan(π/144)、位数走査 n ≤ {N_MAX_ORDER}、帰無 {NULL_TRIALS} 試行\n")
    for n in [10, 20, 50, 100]:
        print(f"--- N = {n} (M = {n*(n-1)//2}) ---")
        sys_lr = lr.LowRankSystem(n)
        rng = np.random.default_rng(60260722)
        v, res, sig_parent = lr.make_parent(sys_lr, rng)
        print(f"  親残差 = {res:.2e}")
        analyze("parent", sig_parent, rng_null)
        sig_st = stationary_state(n, 0)
        analyze("stationary", sig_st, rng_null)
        sig_rd = random_state_spectrum(n, 0)
        analyze("random(対照)", sig_rd, rng_null)
        print()


if __name__ == "__main__":
    main()
