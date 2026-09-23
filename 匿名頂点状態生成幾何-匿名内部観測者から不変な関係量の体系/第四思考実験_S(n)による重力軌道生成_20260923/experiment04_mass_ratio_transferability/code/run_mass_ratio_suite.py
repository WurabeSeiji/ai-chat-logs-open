#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 4 Experiment 04: mass-ratio transferability of the same S(n) rule.

Fixed orbital initial condition:
    p0 = 60, e0 = 0.45
Mass ratios:
    q = m_A/m_B in {1, 2, 4, 10}
where q=1 is the reference case already used in Experiment 03 and q=2,4,10
are the three additional tests requested here.

The S(n) functional rule is unchanged.  Only q (therefore nu) changes.
The local 2.5PN radiation-reaction recurrence uses the same formula and the
corresponding nu as a state/input parameter; no case-specific fitting occurs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

import experiment03_core as e3

P0 = 60.0
E0 = 0.45
Q_VALUES = (1.0, 2.0, 4.0, 10.0)
DEFAULT_ORDER = 12
DEFAULT_STEPS = 4000
DEFAULT_RADIAL_ORBITS = 3.0


def q_id(q: float) -> str:
    return f"Q{str(q).replace('.', 'p')}"


def worker(q: float, base: str, order: int, steps: int, radial_orbits: float):
    # Each ProcessPool worker has its own module state, so setting the global
    # mass ratio here cannot leak into another case.
    e3.Q_MASS = float(q)
    case = {"case_id": q_id(q), "p0": P0, "e0": E0}
    m = e3.run_case(case, base, order, steps, radial_orbits)
    m["is_additional_case"] = bool(q != 1.0)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    ap.add_argument("--order", type=int, default=DEFAULT_ORDER)
    ap.add_argument("--steps-per-radial-orbit", type=int, default=DEFAULT_STEPS)
    ap.add_argument("--radial-orbits", type=float, default=DEFAULT_RADIAL_ORBITS)
    ap.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1))
    args = ap.parse_args()

    base = Path(args.out) if args.out else Path(__file__).resolve().parent.parent
    (base / "raw").mkdir(parents=True, exist_ok=True)
    (base / "summary").mkdir(parents=True, exist_ok=True)

    cases = []
    for iq, q in enumerate(Q_VALUES):
        ma = q / (1.0 + q)
        mb = 1.0 / (1.0 + q)
        nu = ma * mb
        cases.append({
            "case_id": q_id(q),
            "q_index": iq,
            "q_mass_ratio": q,
            "nu": nu,
            "p0": P0,
            "e0": E0,
            "is_additional_case": q != 1.0,
        })

    with open(base / "summary" / "cases.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "04 mass-ratio transferability",
            "purpose": "test the same S(n) rule under changed binary mass ratio",
            "same_Sn_rule_all_cases": True,
            "no_case_specific_fitting": True,
            "p0": P0,
            "e0": E0,
            "q_values": list(Q_VALUES),
            "additional_q_values": [2.0, 4.0, 10.0],
            "order_N": args.order,
            "steps_per_radial_orbit": args.steps_per_radial_orbit,
            "radial_orbits": args.radial_orbits,
            "cases": cases,
        }, f, ensure_ascii=False, indent=2)

    metrics = []
    if args.jobs <= 1:
        for c in cases:
            print(f"running {c['case_id']} q={c['q_mass_ratio']:g} ...", flush=True)
            metrics.append(worker(c["q_mass_ratio"], str(base), args.order,
                                  args.steps_per_radial_orbit, args.radial_orbits))
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {
                ex.submit(worker, c["q_mass_ratio"], str(base), args.order,
                          args.steps_per_radial_orbit, args.radial_orbits): c
                for c in cases
            }
            for fut in as_completed(futs):
                c = futs[fut]
                m = fut.result()
                print(f"done {c['case_id']}: q={c['q_mass_ratio']:g}, nu={m['nu']:.8f}, "
                      f"RMS/a={m['radial_rms_over_a']:.6g}, "
                      f"dprec={m['precession_error_rad']:.6g}", flush=True)
                metrics.append(m)

    metrics.sort(key=lambda m: Q_VALUES.index(float(m["q_mass_ratio"])))
    fields = list(metrics[0].keys())
    with open(base / "summary" / "case_metrics.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(metrics)

    with open(base / "summary" / "case_metrics.json", "w", encoding="utf-8") as f:
        json.dump({
            "experiment": "04 mass-ratio transferability",
            "interpretation_scope": "representation/transferability test; not a derivation of GR from S(n)",
            "same_Sn_rule_all_cases": True,
            "no_case_specific_fitting": True,
            "p0": P0,
            "e0": E0,
            "q_values": list(Q_VALUES),
            "metrics": metrics,
        }, f, ensure_ascii=False, indent=2)

    # raw manifest for later analysis
    rows = []
    for c in cases:
        cid = c["case_id"]
        for kind, fn in [
            ("pn_integrator", f"pn_integrator_{cid}.csv"),
            ("pn_plot", f"pn_plot_{cid}.csv"),
            ("sn", f"sn_{cid}.csv"),
        ]:
            path = base / "raw" / fn
            with open(path, encoding="utf-8", newline="") as fh:
                nrows = sum(1 for _ in fh) - 1
            rows.append({"case_id": cid, "q_mass_ratio": c["q_mass_ratio"],
                         "kind": kind, "file": f"raw/{fn}", "data_rows": nrows})
    with open(base / "summary" / "raw_manifest.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    rms = np.array([m["radial_rms_over_a"] for m in metrics], dtype=float)
    dprec = np.array([m["precession_error_rad"] for m in metrics], dtype=float)
    shrink = np.array([m["sn_peri_shrink_fraction"] - m["pn_peri_shrink_fraction"] for m in metrics], dtype=float)
    summary = {
        "mean_radial_rms_over_a": float(np.mean(rms)),
        "max_radial_rms_over_a": float(np.max(rms)),
        "mean_abs_precession_error_rad": float(np.mean(np.abs(dprec))),
        "max_abs_precession_error_rad": float(np.max(np.abs(dprec))),
        "mean_abs_peri_shrink_fraction_error": float(np.mean(np.abs(shrink))),
        "max_abs_peri_shrink_fraction_error": float(np.max(np.abs(shrink))),
        "max_det_S_minus_1": float(max(m["max_det_S_minus_1"] for m in metrics)),
    }
    with open(base / "summary" / "manifest.json", "w", encoding="utf-8") as f:
        json.dump({"experiment": "04 mass-ratio transferability", "summary": summary,
                   "files": {"metrics": "summary/case_metrics.csv", "raw_manifest": "summary/raw_manifest.csv"}},
                  f, ensure_ascii=False, indent=2)
    print(json.dumps({"summary": summary, "metrics": metrics}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
