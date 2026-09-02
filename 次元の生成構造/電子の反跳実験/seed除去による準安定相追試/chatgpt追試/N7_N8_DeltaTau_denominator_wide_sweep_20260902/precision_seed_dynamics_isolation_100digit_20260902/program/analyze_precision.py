#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Aggregate A/B/C precision-isolation runs.

All threshold logic runs on the log10_Hperp_frac column (strings written at
full precision by the runs; parsed as float64 here, which is exact enough for
threshold crossing and fits since only the exponent scale matters).

Mechanical definitions (recorded, not tuned):
- first_step_above_X: first step with log10(f) > log10(X).
- minimum_Hperp_frac_before_growth: minimum f over steps strictly before
  first_step_above_1e-6 (whole run if that threshold is never crossed).
- initial effective seed floor (fig P7, summary column): f at step 1,
  i.e. the perpendicular fraction expressed after exactly one dynamics step.
- common growth fit window: 1e-20 <= f <= 1e-6 (A/B/C comparable).
- deep fit window (C only): 1e-160 <= f <= 1e-40.
Fit: least squares ln f = a + gamma*x for x in {step, tau}; NA if n < 5.
"""
from __future__ import annotations
import csv, json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parents[1]
CONDS = ["A_IC64_DYN64", "B_IC64_DYN100", "C_IC100_DYN100"]
THRESH = ["1e-180", "1e-150", "1e-120", "1e-90", "1e-60", "1e-30",
          "1e-12", "1e-6", "1e-3", "0.05"]


def load(N, cond):
    p = HERE / "data" / f"N{N}_D{N}" / cond / "timeseries.csv"
    if not p.exists():
        return None
    with open(p, newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        hdr = next(r)
        rows = list(r)
    idx = {h: i for i, h in enumerate(hdr)}
    step = np.array([int(x[idx["step"]]) for x in rows])
    tau = np.array([float(x[idx["tau"]]) for x in rows])
    lg = np.array([(-np.inf if x[idx["log10_Hperp_frac"]] == "-inf"
                    else float(x[idx["log10_Hperp_frac"]])) for x in rows])
    drift = np.array([float(x[idx["H_total_rel_drift"]]) for x in rows])
    clo = np.array([float(x[idx["global_closure"]]) for x in rows])
    init_f_str = rows[0][idx["Hperp_frac"]][:40]
    step1_f_str = rows[1][idx["Hperp_frac"]][:40] if len(rows) > 1 else "NA"
    return {"step": step, "tau": tau, "lg": lg, "drift": drift,
            "closure": clo, "init_f": init_f_str, "step1_f": step1_f_str,
            "n_rows": len(rows)}


def first_above(ts, x):
    lx = math.log10(float(x))
    idx = np.flatnonzero(ts["lg"] > lx)
    return int(ts["step"][idx[0]]) if idx.size else None


def fit(ts, lo, hi):
    llo, lhi = math.log10(lo), math.log10(hi)
    above = np.flatnonzero(ts["lg"] > lhi)
    end = int(above[0]) if above.size else len(ts["lg"])
    sel = np.flatnonzero((ts["lg"][:end] >= llo) & (ts["lg"][:end] <= lhi))
    if sel.size < 5:
        return None
    y = ts["lg"][sel] * math.log(10.0)
    out = {"n": int(sel.size),
           "step_min": int(ts["step"][sel[0]]),
           "step_max": int(ts["step"][sel[-1]]),
           "tau_min": float(ts["tau"][sel[0]]),
           "tau_max": float(ts["tau"][sel[-1]])}
    for name, x in (("step", ts["step"][sel].astype(float)),
                    ("tau", ts["tau"][sel])):
        A = np.vstack([x, np.ones_like(x)]).T
        coef, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
        yh = A @ coef
        sst = float(np.sum((y - y.mean()) ** 2))
        out[f"gamma_{name}"] = float(coef[0])
        out[f"intercept_{name}"] = float(coef[1])
        out[f"R2_{name}"] = (1.0 - float(np.sum((y - yh) ** 2)) / sst
                             if sst > 0 else 0.0)
    return out


def main():
    (HERE / "results").mkdir(exist_ok=True)
    summary, fits, quals = [], [], []
    for N in (8, 7):
        for cond in CONDS:
            ts = load(N, cond)
            if ts is None:
                continue
            row = {"N": N, "D": N, "condition": cond,
                   "steps_recorded": ts["n_rows"] - 1,
                   "initial_Hperp_frac": ts["init_f"],
                   "initial_log10_Hperp_frac": (
                       "-inf" if not np.isfinite(ts["lg"][0])
                       else f"{ts['lg'][0]:.6g}"),
                   "step1_Hperp_frac": ts["step1_f"],
                   "step1_log10_Hperp_frac": (
                       f"{ts['lg'][1]:.6g}" if len(ts["lg"]) > 1 else "NA")}
            s6 = first_above(ts, "1e-6")
            end = (np.flatnonzero(ts["step"] == s6)[0] if s6 is not None
                   else len(ts["lg"]))
            pre = ts["lg"][:end]
            fin = pre[np.isfinite(pre)]
            if fin.size:
                mi = float(fin.min())
                mstep = int(ts["step"][:end][np.isfinite(pre)][fin.argmin()])
                row["min_log10_Hperp_before_growth"] = f"{mi:.6g}"
                row["minimum_Hperp_step"] = mstep
            else:
                row["min_log10_Hperp_before_growth"] = "NA"
                row["minimum_Hperp_step"] = "NA"
            for th in THRESH:
                fa = first_above(ts, th)
                row[f"first_step_above_{th}"] = fa if fa is not None else "NA"
                row[f"first_tau_above_{th}"] = (
                    f"{float(ts['tau'][np.flatnonzero(ts['step']==fa)[0]]):.8g}"
                    if fa is not None else "NA")
            row["Htotal_max_rel_drift"] = f"{float(ts['drift'].max()):.6g}"
            row["closure_final"] = f"{float(ts['closure'][-1]):.8g}"
            summary.append(row)
            f1 = fit(ts, 1e-20, 1e-6)
            if f1:
                fits.append({"N": N, "condition": cond,
                             "window": "1e-20..1e-6", **f1})
            if cond == "C_IC100_DYN100":
                f2 = fit(ts, 1e-160, 1e-40)
                if f2:
                    fits.append({"N": N, "condition": cond,
                                 "window": "1e-160..1e-40", **f2})
        qp = (HERE / "data" / f"N{N}_D{N}" / "C_IC100_DYN100" /
              "ic100_qualification.json")
        if qp.exists():
            q = json.load(open(qp))
            quals.append({k: q[k] for k in
                          ("N", "M_edge", "norm", "mean_amp2",
                           "global_closure_normalized",
                           "local_closure_max_abs", "H_eigen_residual",
                           "mu", "mu_minus_theory_abs",
                           "H_hermiticity_error", "qualification_passed")})

    def w(path, rows):
        if not rows:
            return
        keys = list(rows[0].keys())
        with open(path, "w", newline="", encoding="utf-8") as f:
            wr = csv.DictWriter(f, fieldnames=keys)
            wr.writeheader()
            wr.writerows(rows)
    w(HERE / "results" / "precision_summary.csv", summary)
    w(HERE / "results" / "growth_fits.csv", fits)
    w(HERE / "results" / "qualification_ic100.csv", quals)
    print(f"summary rows: {len(summary)}, fits: {len(fits)}, "
          f"quals: {len(quals)}")


if __name__ == "__main__":
    main()
