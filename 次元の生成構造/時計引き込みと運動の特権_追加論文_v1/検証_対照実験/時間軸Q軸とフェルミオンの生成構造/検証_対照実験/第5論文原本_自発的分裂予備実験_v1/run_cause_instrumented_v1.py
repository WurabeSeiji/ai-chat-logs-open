#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""拡大と停止の「原因」を計装する再走行：f と σ₂/σ₁ を同時記録する。

原本 run_spontaneous_splitting_largeN_v1.py の run() ループを忠実に複製し
（力学は一切変えない：同一種・同一 power 反復 σ₁・同一 Cayley・同一停止）、
毎ステップ（大Nは record_every ごと）に第2特異値比 σ₂/σ₁ を追加記録する。
σ₂/σ₁ は診断専用で、ステップの正規化には使わない（力学不変）。

再現性規約: 走行後、記録した f を、同一パラメータで作った正本
metastable_series_result_v1/fcurve_*.csv と突き合わせる対照テストを行い、
一致しなければ終了コード1で失敗する。

使い方:
    python3 run_cause_instrumented_v1.py 5    --after=20000 --cap=30000 --record-every=1
    python3 run_cause_instrumented_v1.py 40   --after=20000 --cap=30000 --record-every=1
    python3 run_cause_instrumented_v1.py 300  --after=20000 --cap=30000 --record-every=50
"""

import argparse
import csv
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_n_scaling_lowrank_v1 import (  # 原本エンジン（md5固定）
    LowRankSystem, make_parent, zero_closure_kernel_seed,
)

RESULT_DIR = HERE / "cause_instrumented_result_v1"
REFERENCE_DIR = HERE / "metastable_series_result_v1"


def sigma2_over_sigma1(sys_lr):
    sig = sys_lr.sigma_spectrum()
    if len(sig) >= 2 and sig[0] > 0:
        return float(sig[1] / sig[0])
    return 0.0


def run(n, delta, seed, cap, after, record_every, tol=1e-12):
    # ---- run_spontaneous_splitting_largeN_v1.run() の忠実複製（診断行のみ追加）----
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n + seed)  # 原本と同一種
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=tol)
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + delta * g
    Z = Z / np.linalg.norm(Z)

    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    f_hist = []          # 全ステップの f（対照テスト用）
    rows = []            # 診断行（record_every ごと: tau, f, sigma2/sigma1）
    crossed_at = None
    wp = rng.normal(size=sys_lr.m)
    t0 = time.time()
    for t in range(cap + 1):
        Zp = Z - p * (p @ Z) - q * (q @ Z)
        htot = float(np.real(np.conj(Z) @ Z))
        f = float(np.real(np.conj(Zp) @ Zp)) / htot
        f_hist.append(f)
        if crossed_at is None and f > 0.05:
            crossed_at = t
        if crossed_at is not None and t >= crossed_at + after:
            break
        sys_lr.set_theta(np.angle(Z))
        if t % record_every == 0:
            rows.append((t, f, sigma2_over_sigma1(sys_lr)))
        sig_est, wp = sys_lr.sigma_max_power(wp)
        Z = sys_lr.cayley_step(Z, sig_est)
    # ---- 複製ここまで ----

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"N{n:05d}_delta{delta:.0e}_seed{seed}"
    out_csv = RESULT_DIR / f"cause_{tag}.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tau", "f", "sigma2_over_sigma1"])
        w.writerows(rows)

    # ---- 対照テスト: f が正本 metastable fcurve と一致するか ----
    ref_csv = REFERENCE_DIR / f"fcurve_{tag}.csv"
    control = {"reference": str(ref_csv), "checked": False}
    if ref_csv.exists():
        ref = {}
        with open(ref_csv) as fh:
            for r in csv.DictReader(fh):
                ref[int(r["tau"])] = float(r["f"])
        max_dev = 0.0
        for (t, f, _) in rows:
            if t in ref:
                max_dev = max(max_dev, abs(f - ref[t]))
        control = {"reference": str(ref_csv), "checked": True,
                   "max_f_deviation": max_dev, "passed": max_dev < 1e-12}
    summary = {
        "n": n, "m": sys_lr.m, "delta": delta, "seed": seed,
        "cap": cap, "after": after, "record_every": record_every,
        "parent_residual": residual, "crossing_tau": crossed_at,
        "steps_run": len(f_hist) - 1, "runtime_sec": time.time() - t0,
        "control_test": control,
    }
    with open(RESULT_DIR / f"cause_summary_{tag}.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print(json.dumps({"n": n, "crossing_tau": crossed_at,
                      "steps_run": len(f_hist) - 1,
                      "control_test": control}, ensure_ascii=False))
    if control.get("checked") and not control["passed"]:
        raise SystemExit(f"対照テスト失敗: f 最大偏差 {control['max_f_deviation']:.2e}")
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--delta", type=float, default=1e-15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cap", type=int, default=30000)
    ap.add_argument("--after", type=int, default=20000)
    ap.add_argument("--record-every", type=int, default=1)
    ap.add_argument("--tol", type=float, default=1e-12)
    args = ap.parse_args()
    run(args.n, args.delta, args.seed, args.cap, args.after,
        args.record_every, tol=args.tol)


if __name__ == "__main__":
    main()
