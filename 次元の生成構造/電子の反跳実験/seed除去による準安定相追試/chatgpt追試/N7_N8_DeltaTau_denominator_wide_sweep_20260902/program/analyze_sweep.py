#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Aggregate D-sweep runs into stage summaries + sweep_summary.csv.

Mechanical rules only (no reinterpretation):
- onset(thr): first step with Hperp_frac > thr, thr in {1e-6, 1e-3, 0.05}.
- growth fit: points with 1e-12 <= f <= 1e-4 and index < first index f > 1e-4;
  least squares ln f = a + gamma*x for x = step and x = tau; NA if n < 5.
- saturation: last 20% of recorded rows; NA unless onset(0.05) exists and
  tail start step > onset step.
- anomaly candidates (Stage C trigger), evaluated per N on stage A
  (dense integer D) for onset_tau_0p05 and gamma_tau:
  |value - median(D-2..D+2 neighbors)| / median > 0.20  -> flagged.
  Additionally D=124 is compared against {96,112,128,160} (stage B).
"""
from __future__ import annotations
import csv, glob, hashlib, json, math, os, platform
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = [("1e-6", 1e-6), ("1e-3", 1e-3), ("0p05", 0.05)]


def load_ts(path: Path):
    cols = {}
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        header = next(r)
        data = list(r)
    arr = {h: np.array([row[i] for row in data]) for i, h in enumerate(header)}
    for k in arr:
        if k == "finite":
            arr[k] = arr[k] == "True"
        elif k in ("N", "D", "step"):
            arr[k] = arr[k].astype(int)
        else:
            arr[k] = arr[k].astype(float)
    return arr


def onset(ts, thr):
    idx = np.flatnonzero(ts["Hperp_frac"] > thr)
    if idx.size == 0:
        return None
    i = int(idx[0])
    return {"step": int(ts["step"][i]), "tau": float(ts["tau"][i]),
            "cycles": float(ts["cycles"][i]), "chi": float(ts["chi"][i])}


def growth_fit(ts):
    f = ts["Hperp_frac"]
    above = np.flatnonzero(f > 1e-4)
    end = int(above[0]) if above.size else len(f)
    sel = np.flatnonzero((f[:end] >= 1e-12) & (f[:end] <= 1e-4))
    if sel.size < 5:
        return None
    y = np.log(f[sel])
    out = {}
    for name, x in (("step", ts["step"][sel].astype(float)),
                    ("tau", ts["tau"][sel])):
        A = np.vstack([x, np.ones_like(x)]).T
        coef, res, _, _ = np.linalg.lstsq(A, y, rcond=None)
        yhat = A @ coef
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - float(np.sum((y - yhat) ** 2)) / ss_tot if ss_tot > 0 else 0.0
        out[name] = {"slope": float(coef[0]), "intercept": float(coef[1]),
                     "R2": r2}
    out["n"] = int(sel.size)
    out["fit_step_min"] = int(ts["step"][sel[0]])
    out["fit_step_max"] = int(ts["step"][sel[-1]])
    out["fit_tau_min"] = float(ts["tau"][sel[0]])
    out["fit_tau_max"] = float(ts["tau"][sel[-1]])
    return out


def tail_stats(ts, key, i0):
    v = ts[key][i0:]
    return {"mean": float(v.mean()), "std": float(v.std()),
            "min": float(v.min()), "max": float(v.max()),
            "q05": float(np.quantile(v, 0.05)),
            "q95": float(np.quantile(v, 0.95))}


def summarize_run(path: Path, stage: str, shas, prog_sha, env):
    ts = load_ts(path)
    N = int(ts["N"][0]); D = int(ts["D"][0]); M = N * (N - 1) // 2
    row = {"N": N, "M_edge": M, "D": D, "D_over_N": D / N,
           "r2bar": float(ts["H_total"][0]) / M, "stage": stage,
           "steps_run": int(ts["step"][-1]),
           "tau_max": float(ts["tau"][-1]),
           "initial_Hperp_frac": float(ts["Hperp_frac"][0])}
    ons = {}
    for name, thr in THRESHOLDS:
        o = onset(ts, thr)
        ons[name] = o
        row[f"onset_step_{name}"] = o["step"] if o else "NA"
        row[f"onset_tau_{name}"] = f"{o['tau']:.10g}" if o else "NA"
    o05 = ons["0p05"]
    row["onset_cycles_0p05"] = f"{o05['cycles']:.10g}" if o05 else "NA"
    row["onset_chi_0p05"] = f"{o05['chi']:.10g}" if o05 else "NA"
    fit = growth_fit(ts)
    if fit:
        row["gamma_step"] = f"{fit['step']['slope']:.10g}"
        row["gamma_tau"] = f"{fit['tau']['slope']:.10g}"
        row["fit_R2"] = f"{fit['tau']['R2']:.8g}"
        row["fit_n"] = fit["n"]
        row["fit_step_min"] = fit["fit_step_min"]
        row["fit_step_max"] = fit["fit_step_max"]
        row["fit_tau_min"] = f"{fit['fit_tau_min']:.10g}"
        row["fit_tau_max"] = f"{fit['fit_tau_max']:.10g}"
        pred = fit["step"]["slope"] * D / (2 * math.pi)
        row["gamma_crosscheck_rel"] = (
            f"{abs(fit['tau']['slope'] - pred) / abs(pred):.3g}"
            if pred != 0 else "NA")
    else:
        for k in ("gamma_step", "gamma_tau", "fit_R2", "fit_n",
                  "fit_step_min", "fit_step_max", "fit_tau_min",
                  "fit_tau_max", "gamma_crosscheck_rel"):
            row[k] = "NA"
    n_rows = len(ts["step"])
    i0 = int(math.floor(0.8 * n_rows))
    if o05 and int(ts["step"][i0]) > o05["step"]:
        s = tail_stats(ts, "Hperp_frac", i0)
        for k, v in s.items():
            row[f"sat_{k}_Hperp_frac"] = f"{v:.10g}"
        c = tail_stats(ts, "global_closure", i0)
        row["closure_tail_mean"] = f"{c['mean']:.10g}"
        row["closure_tail_std"] = f"{c['std']:.10g}"
        prm = tail_stats(ts, "PR_over_M", i0)
        row["PR_over_M_tail_mean"] = f"{prm['mean']:.10g}"
        row["PR_over_M_tail_std"] = f"{prm['std']:.10g}"
        row["sat_window_note"] = "ok"
    else:
        for k in ("sat_mean_Hperp_frac", "sat_std_Hperp_frac",
                  "sat_min_Hperp_frac", "sat_max_Hperp_frac",
                  "sat_q05_Hperp_frac", "sat_q95_Hperp_frac",
                  "closure_tail_mean", "closure_tail_std",
                  "PR_over_M_tail_mean", "PR_over_M_tail_std"):
            row[k] = "NA"
        row["sat_window_note"] = ("no_onset" if not o05
                                  else "tail_overlaps_onset_needs_longer_run")
    row["Htotal_max_rel_drift"] = f"{float(ts['H_total_rel_drift'].max()):.6g}"
    row["numpy_version"] = env["numpy"]
    row["python_version"] = env["python"]
    row["platform"] = env["platform"]
    row["input_sha256"] = shas[N]
    row["program_sha256"] = prog_sha
    return row


def write_csv(path: Path, rows):
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def anomalies(rows_by_stage):
    """Mechanical anomaly candidates for Stage C."""
    out = []
    for N in (7, 8):
        # dense integer scan on stage A
        sa = sorted((r for r in rows_by_stage["A"]
                     if r["N"] == N and r["onset_tau_0p05"] != "NA"
                     and 2 <= r["D"] <= 256), key=lambda r: r["D"])
        vals = {r["D"]: float(r["onset_tau_0p05"]) for r in sa}
        for r in sa:
            D = r["D"]
            nb = [vals[d] for d in range(D - 2, D + 3)
                  if d != D and d in vals]
            if len(nb) < 3:
                continue
            med = float(np.median(nb))
            if med > 0 and abs(float(r["onset_tau_0p05"]) - med) / med > 0.20:
                out.append({"N": N, "D": D, "quantity": "onset_tau_0p05",
                            "value": r["onset_tau_0p05"],
                            "neighborhood_median": f"{med:.10g}",
                            "rel_dev": f"{abs(float(r['onset_tau_0p05'])-med)/med:.3g}",
                            "rule": "stageA_pm2_median_20pct"})
        # 124 vs neighbors on stage B
        sb = {r["D"]: r for r in rows_by_stage["B"] if r["N"] == N}
        if 124 in sb:
            for qty in ("onset_tau_0p05", "gamma_tau", "sat_mean_Hperp_frac"):
                v124 = sb[124].get(qty, "NA")
                nbs = [sb[d].get(qty, "NA") for d in (96, 112, 128, 160)
                       if d in sb]
                nbs = [float(x) for x in nbs if x != "NA"]
                if v124 == "NA" or len(nbs) < 3:
                    continue
                v = float(v124)
                lo, hi = min(nbs), max(nbs)
                span = hi - lo
                if v < lo - 0.10 * abs(lo) - span or \
                   v > hi + 0.10 * abs(hi) + span:
                    out.append({"N": N, "D": 124, "quantity": qty,
                                "value": v124,
                                "neighborhood_median":
                                    f"{float(np.median(nbs)):.10g}",
                                "rel_dev": "outside_neighbor_range",
                                "rule": "B_124_vs_96_112_128_160"})
    return out


def main():
    env = {"python": platform.python_version(), "numpy": np.__version__,
           "platform": platform.platform()}
    shas = {}
    for N in (7, 8):
        shas[N] = hashlib.sha256(
            (ROOT / "data" / f"N{N}" / "parent_v.npz").read_bytes()).hexdigest()
    prog_sha = hashlib.sha256(
        (ROOT / "program" / "run_sweep.py").read_bytes()).hexdigest()

    rows_by_stage = {"A": [], "B": []}
    for stage in ("A", "B"):
        for p in sorted(ROOT.glob(f"data/N*/D*/timeseries_stage{stage}.csv")):
            rows_by_stage[stage].append(
                summarize_run(p, stage, shas, prog_sha, env))
        rows_by_stage[stage].sort(key=lambda r: (r["N"], r["D"]))
        write_csv(ROOT / "results" / f"stage{stage}_summary.csv",
                  rows_by_stage[stage])
    write_csv(ROOT / "results" / "sweep_summary.csv",
              rows_by_stage["A"] + rows_by_stage["B"])

    an = anomalies(rows_by_stage)
    path = ROOT / "results" / "anomaly_followups.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["N", "D", "quantity", "value",
                                          "neighborhood_median", "rel_dev",
                                          "rule"])
        w.writeheader()
        w.writerows(an)
    print(f"stageA runs: {len(rows_by_stage['A'])}, "
          f"stageB runs: {len(rows_by_stage['B'])}, anomalies: {len(an)}")


if __name__ == "__main__":
    main()
