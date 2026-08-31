#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent hm rank-3 geometry audit.

Purpose
-------
Recompute edge wavelengths directly from raw data/hm_N*/states_treatment.npz,
without using hm_series_k.csv or any pass13-17 derived result as input.

This script intentionally separates:
  A) observation/readout: dominant edge frequencies and normalized wavelengths;
  B) geometry audit: centered Gram matrix PSD/rank diagnostics;
  C) optional reproduction of the old *group-uniform* integer-k ansatz, only as
     a diagnostic baseline.  It is NOT treated as the physical rank-3 search.

The physical target discussed in this audit is rank(B)=3 for N>4.  A full
edge-wise independent integer search is not claimed by this script; that is a
separate combinatorial problem and must not be confused with the baseline
ansatz reproduced here.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import minimize_scalar

DEFAULT_WINDOWS = (4096, 8192, 16384)
PAD = 8
RANK_TOL = 1e-9
PSD_TOL = 1e-9


def edges(n: int):
    return list(itertools.combinations(range(n), 2))


def dominant_omega(sig: np.ndarray, pad: int = PAD) -> float:
    """Dominant angular frequency, matching the pass17 readout method."""
    sig = np.asarray(sig)
    n = len(sig)
    h = np.hanning(n)
    F = np.fft.fft(sig * h, n * pad)
    P = np.abs(F) ** 2
    fr = np.fft.fftfreq(n * pad) * 2 * math.pi
    om0 = float(fr[int(np.argmax(P))])
    t = np.arange(n)
    bw = 2 * math.pi / n

    def objective(om: float) -> float:
        return -abs(np.sum(sig * h * np.exp(-1j * om * t)))

    res = minimize_scalar(
        objective,
        bounds=(om0 - bw, om0 + bw),
        method="bounded",
        options={"xatol": 1e-9},
    )
    return float(res.x)


def normalized_wavelengths(Z: np.ndarray, window: int) -> tuple[np.ndarray, np.ndarray]:
    if window > len(Z):
        raise ValueError(f"window={window} exceeds samples={len(Z)}")
    W = Z[-window:]
    nus = np.array([dominant_omega(W[:, e]) for e in range(W.shape[1])], dtype=float)
    lam = 1.0 / np.abs(nus)
    lamn = lam / np.min(lam)
    return nus, lamn


def centered_gram(n: int, lengths: np.ndarray) -> np.ndarray:
    E = edges(n)
    if len(lengths) != len(E):
        raise ValueError("edge count mismatch")
    D2 = np.zeros((n, n), dtype=float)
    for e, (i, j) in enumerate(E):
        D2[i, j] = D2[j, i] = float(lengths[e]) ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * J @ D2 @ J


def gram_diagnostics(B: np.ndarray, rank_tol: float = RANK_TOL, psd_tol: float = PSD_TOL):
    ev = np.linalg.eigvalsh(B)[::-1]
    scale = max(float(np.max(np.abs(ev))), 1e-300)
    evn = ev / scale
    rank = int(np.sum(evn > rank_tol))
    psd = bool(np.min(evn) >= -psd_tol)
    return ev, evn, rank, psd


def cluster_lambdas(lamn: np.ndarray, rel_tol: float = 0.02):
    groups: list[list] = []
    for e, value in enumerate(lamn):
        for g in groups:
            if abs(value - g[0]) / g[0] < rel_tol:
                g[1].append(e)
                g[0] = (g[0] * (len(g[1]) - 1) + value) / len(g[1])
                break
        else:
            groups.append([float(value), [e]])
    groups.sort(key=lambda x: x[0])
    return groups


def reproduce_group_uniform_baseline(n: int, lamn: np.ndarray, cmax: int = 100):
    """Reproduce old 2% group-uniform ansatz, but report rank instead of calling any PSD rank 'success'."""
    groups = cluster_lambdas(lamn, 0.02)
    E = edges(n)
    idx = {}
    for gi, g in enumerate(groups):
        for e in g[1]:
            idx[e] = gi

    if len(groups) == 1:
        kg = (1,)
        lengths = np.array([groups[0][0] / 2.0] * len(E))
        B = centered_gram(n, lengths)
        ev, evn, rank, psd = gram_diagnostics(B)
        return [{"kg": kg, "rank": rank, "psd": psd, "evn": evn, "groups": groups}]

    out = []
    seen = set()
    # Same two paths as pass17, but do not stop at first PSD solution.
    for d in (1, 2, 3):
        cmin = max(1, int(0.5 * max(g[0] for g in groups[1:]) * d))
        for c in range(cmin, cmax + 1):
            kg = (c,) + (d,) * (len(groups) - 1)
            if kg in seen:
                continue
            seen.add(kg)
            lengths = np.array([kg[idx[e]] * groups[idx[e]][0] / 2.0 for e in range(len(E))])
            B = centered_gram(n, lengths)
            ev, evn, rank, psd = gram_diagnostics(B)
            out.append({"kg": kg, "rank": rank, "psd": psd, "evn": evn, "groups": groups})

    # Exhaustive group-wise fallback only when tractable, matching pass17 ranges.
    ng = len(groups)
    if ng <= 8:
        for c in range(1, cmax + 1):
            for ks in itertools.product(range(1, 4), repeat=ng - 1):
                kg = (c,) + ks
                if kg in seen:
                    continue
                seen.add(kg)
                lengths = np.array([kg[idx[e]] * groups[idx[e]][0] / 2.0 for e in range(len(E))])
                B = centered_gram(n, lengths)
                ev, evn, rank, psd = gram_diagnostics(B)
                out.append({"kg": kg, "rank": rank, "psd": psd, "evn": evn, "groups": groups})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, required=True, help="project data directory containing hm_N*/states_treatment.npz")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-min", type=int, default=3)
    ap.add_argument("--n-max", type=int, default=16)
    ap.add_argument("--windows", type=int, nargs="*", default=list(DEFAULT_WINDOWS))
    ap.add_argument("--baseline-cmax", type=int, default=100)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    edge_rows = []
    audit_rows = []

    for n in range(args.n_min, args.n_max + 1):
        fp = args.data_root / f"hm_N{n}" / "states_treatment.npz"
        if not fp.exists():
            audit_rows.append({"N": n, "status": "MISSING_RAW", "samples": "", "window": "", "rank3_count_group_baseline": "", "best_psd_rank": "", "best_psd_mu4_over_mu1": "", "best_psd_kg": ""})
            continue
        Z = np.load(fp)["Z"]
        E = edges(n)
        if Z.ndim != 2 or Z.shape[1] != len(E):
            raise ValueError(f"{fp}: Z shape {Z.shape}, expected second axis {len(E)}")

        valid_windows = [w for w in args.windows if w <= len(Z)]
        if not valid_windows:
            raise ValueError(f"{fp}: no valid windows")

        readouts = {}
        for w in valid_windows:
            nus, lamn = normalized_wavelengths(Z, w)
            readouts[w] = (nus, lamn)
            for e, (i, j) in enumerate(E):
                edge_rows.append({"N": n, "edge_index": e, "i": i, "j": j, "window": w, "omega": nus[e], "lambda_norm": lamn[e]})

        # Geometry baseline uses the shortest requested valid window to reduce phase mixing.
        w0 = min(valid_windows)
        lamn = readouts[w0][1]
        candidates = reproduce_group_uniform_baseline(n, lamn, args.baseline_cmax)
        rank3 = [c for c in candidates if c["psd"] and c["rank"] == 3]
        psd = [c for c in candidates if c["psd"]]
        best = None
        if psd:
            # Prefer smallest rank, then smallest fourth eigenvalue ratio, then lexicographic kg.
            def score(c):
                evn = c["evn"]
                mu4 = float(evn[3]) if len(evn) > 3 else 0.0
                return (c["rank"], abs(mu4), tuple(c["kg"]))
            best = min(psd, key=score)
        audit_rows.append({
            "N": n,
            "status": "OK",
            "samples": len(Z),
            "window": w0,
            "rank3_count_group_baseline": len(rank3),
            "best_psd_rank": "" if best is None else best["rank"],
            "best_psd_mu4_over_mu1": "" if best is None or len(best["evn"]) <= 3 else float(best["evn"][3]),
            "best_psd_kg": "" if best is None else repr(tuple(best["kg"])),
        })

    with (args.out / "edge_wavelengths.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["N", "edge_index", "i", "j", "window", "omega", "lambda_norm"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(edge_rows)

    with (args.out / "rank3_group_baseline_audit.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["N", "status", "samples", "window", "rank3_count_group_baseline", "best_psd_rank", "best_psd_mu4_over_mu1", "best_psd_kg"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(audit_rows)

    print(f"wrote {args.out / 'edge_wavelengths.csv'}")
    print(f"wrote {args.out / 'rank3_group_baseline_audit.csv'}")


if __name__ == "__main__":
    main()
