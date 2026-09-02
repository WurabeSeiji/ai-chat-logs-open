#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Condition C: IC100 + Dynamics100.

The dynamics implementation is imported unchanged from
run_mp100_same_ic64.py, so B -> C changes ONLY the initial state
(instruction section 0). Auto-extension (section 5): if no primary onset
(Hperp/H > 0.05) by step 2000, extend 4000 -> 8000 via checkpoints.
"""
from __future__ import annotations
import argparse, csv, json, math, sys
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_mp100_same_ic64 import (run, selftest, json_to_mpc)  # noqa: E402

mp.mp.dps = 100
LOG10_ONSET = math.log10(0.05)


def has_onset(outdir):
    with open(outdir / "timeseries.csv", newline="") as f:
        r = csv.reader(f)
        hdr = next(r)
        i = hdr.index("log10_Hperp_frac")
        for row in r:
            if row[i] != "-inf" and float(row[i]) > LOG10_ONSET:
                return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True, choices=[7, 8])
    args = ap.parse_args()
    N = args.n
    D = N
    out = HERE / "data" / f"N{N}_D{D}" / "C_IC100_DYN100"
    qual = json.load(open(out / "ic100_qualification.json"))
    if not qual.get("qualification_passed"):
        raise SystemExit("IC100 not qualified; refusing to run")
    st = selftest()
    with open(out / "precision_selftest.json", "w", encoding="utf-8") as f:
        json.dump(st, f, indent=1)
    ic = json.load(open(out / "ic100_state.json"))
    z0 = [json_to_mpc(t) for t in ic["z"]]

    targets = [2000, 4000, 8000]
    done_steps = None
    for tgt in targets:
        run(N, D, z0, out, tgt, "C_IC100_DYN100")
        done_steps = tgt
        if has_onset(out):
            break
    with open(out / "run_info.json", "w", encoding="utf-8") as f:
        json.dump({"condition": "C_IC100_DYN100", "N": N, "D": D,
                   "steps": done_steps, "dps": 100,
                   "onset_reached": has_onset(out),
                   "ic": "ic100_state.json (analytic phases, lifted legacy "
                         "norm amplitude)",
                   "pi": "100-digit mpmath", "dynamics": "mpmath.eighe "
                         "imported from run_mp100_same_ic64.py"},
                  f, indent=1)
    print(f"C N={N} COMPLETE steps={done_steps} onset={has_onset(out)}")


if __name__ == "__main__":
    main()
