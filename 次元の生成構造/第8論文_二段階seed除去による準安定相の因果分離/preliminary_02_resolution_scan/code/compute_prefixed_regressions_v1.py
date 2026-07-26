#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第2予備実験 §11 事前固定回帰。解釈・選別・評価をしてはならない（回帰値の算出のみ）。

固定振幅帯 B1..B6（境界固定, §11.1）で log a_outside(t)=α+γt の通常最小二乗を行う。
有効点数<20 は insufficient_points。区間の追加・結合・移動はしない。
出力: raw/<run_id>/regression_by_fixed_band.csv（各run）と summary/all_fixed_band_regressions.csv（結合）。
γ_local は各run local_growth.csv に既出。
"""
import csv
import json
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
P2 = CODE.parent
RAW = P2 / "raw"
SUM = P2 / "summary"; SUM.mkdir(exist_ok=True)

BANDS = [("B1", 1e-14, 1e-12), ("B2", 1e-12, 1e-10), ("B3", 1e-10, 1e-8),
         ("B4", 1e-8, 1e-6), ("B5", 1e-6, 1e-4), ("B6", 1e-4, 1e-2)]
MIN_POINTS = 20


def ols(t, y):
    t = np.asarray(t, float); y = np.asarray(y, float)
    A = np.column_stack([np.ones_like(t), t])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    alpha, gamma = coef
    yhat = A @ coef
    ss_res = float(np.sum((y - yhat) ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = (1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(ss_res / len(y)))
    return float(alpha), float(gamma), r2, rmse


HDR = ["run_id", "N", "p", "Delta_ref", "Delta_actual", "resolution_operator", "band",
       "a_lo", "a_hi", "n_points", "first_step", "last_step", "alpha", "gamma", "exp_gamma",
       "r_squared", "rmse", "min_a", "max_a", "status"]


def band_rows(cfg, steps, a):
    out = []
    for name, lo, hi in BANDS:
        idx = np.where((a >= lo) & (a < hi) & (a > 0))[0]
        base = [cfg["run_id"], cfg["N"], cfg["p"], cfg["Delta_ref"], cfg["Delta_actual"],
                cfg["resolution_operator"], name, "%.0e" % lo, "%.0e" % hi]
        if idx.size < MIN_POINTS:
            out.append(base + [idx.size, (int(steps[idx[0]]) if idx.size else ""),
                               (int(steps[idx[-1]]) if idx.size else ""),
                               "", "", "", "", "", "", "", "insufficient_points"])
            continue
        tt = steps[idx].astype(float); yy = np.log(a[idx])
        alpha, gamma, r2, rmse = ols(tt, yy)
        out.append(base + [idx.size, int(steps[idx[0]]), int(steps[idx[-1]]),
                           "%.12e" % alpha, "%.12e" % gamma, "%.12e" % np.exp(gamma),
                           "%.12e" % r2, "%.12e" % rmse, "%.12e" % a[idx].min(), "%.12e" % a[idx].max(), "ok"])
    return out


def main():
    runs = sorted([d for d in RAW.iterdir() if d.is_dir() and (d / "timeseries.csv").exists()])
    combined = SUM / "all_fixed_band_regressions.csv"
    with open(combined, "w", newline="") as fh:
        wc = csv.writer(fh); wc.writerow(HDR)
        for d in runs:
            cfg = json.load(open(d / "run_config.json"))
            rows = list(csv.DictReader(open(d / "timeseries.csv")))
            steps = np.array([int(r["step"]) for r in rows])
            a = np.array([float(r["a_outside"]) for r in rows])
            brows = band_rows(cfg, steps, a)
            with open(d / "regression_by_fixed_band.csv", "w", newline="") as pf:
                wp = csv.writer(pf); wp.writerow(HDR); wp.writerows(brows)
            wc.writerows(brows)
    print(f"[regress] {combined}  runs={len(runs)}")


if __name__ == "__main__":
    main()
