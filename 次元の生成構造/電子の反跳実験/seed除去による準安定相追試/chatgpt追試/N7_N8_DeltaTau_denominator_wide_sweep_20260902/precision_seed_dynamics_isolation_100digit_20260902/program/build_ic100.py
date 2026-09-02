#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build IC100: the hm parent regenerated at 100 decimal digits.

Phase rules are taken unchanged from the canonical generator
(generate_hm_mp_free_N3_N40_20260901.py):
- N=8 (even): K_N 1-factorization, classes c=0..N-2, theta_c = pi*c/(N-1).
- N=7 (odd):  cyclic distance classes d=1..(N-1)/2, c=d-1, theta_c = pi*c/q.
Edge ordering: [(i,j) for i<j] lexicographic (identical to A/B).

Amplitude convention (instruction section 4.2): N=7,8 belong to the
historical-compatibility regime whose LEGACY_NORM[N] is a frozen binary64
constant with no analytic formula. As instructed, that fact is recorded and
the SAME amplitude scale as A/B is used via exact binary64 lift of
LEGACY_NORM[N]; only the phases are generated analytically at 100 digits.
The legacy operation order v * NORM / ||v|| is reproduced at 100 digits.

Qualification (section 10) is computed at 100 digits and gates usage.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parents[1]
mp.mp.dps = 100

LEGACY_NORM = {7: 1.1766901825131457, 8: 1.2815849983171124}
GATE = mp.mpf("1e-90")


def lift64(x):
    num, den = float(x).as_integer_ratio()
    return mp.mpf(num) / mp.mpf(den)


def mpc_to_json(c):
    tr, ti = c.real._mpf_, c.imag._mpf_
    return [[tr[0], str(tr[1]), tr[2], tr[3]],
            [ti[0], str(ti[1]), ti[2], ti[3]]]


def classes(N):
    E = [(i, j) for i in range(N) for j in range(i + 1, N)]
    if N % 2 == 0:
        n = N - 1
        col = {}
        for rr in range(n):
            col[tuple(sorted((rr, N - 1)))] = rr
            for k in range(1, N // 2):
                col[tuple(sorted(((rr - k) % n, (rr + k) % n)))] = rr
        return E, col, n
    q = (N - 1) // 2
    col = {}
    for d in range(1, q + 1):
        for i in range(N):
            col[tuple(sorted((i, (i + d) % N)))] = d - 1
    return E, col, q


def build(N):
    E, col, nclass = classes(N)
    M = len(E)
    rbar = mp.sqrt(mp.mpf(1) / mp.mpf(15))
    v = []
    for e in E:
        c = col[e]
        theta = mp.pi * c / nclass
        v.append(rbar * mp.mpc(mp.cos(theta), mp.sin(theta)))
    norm = mp.sqrt(mp.re(sum(mp.conj(x) * x for x in v)))
    scale = lift64(LEGACY_NORM[N]) / norm
    v = [x * scale for x in v]
    return E, col, v


def qualify(N, E, v):
    M = len(E)
    h = mp.re(sum(mp.conj(x) * x for x in v))
    norm = mp.sqrt(h)
    r2 = h / M
    gs = sum(x * x for x in v)
    S = [mp.mpc(0) for _ in range(N)]
    for k, (i, j) in enumerate(E):
        S[i] += v[k] ** 2
        S[j] += v[k] ** 2
    loc = max(abs(s) for s in S)
    adj = [[] for _ in range(M)]
    for a in range(M):
        sa = set(E[a])
        for b in range(M):
            if b != a and sa.intersection(E[b]):
                adj[a].append(b)
    Hv = []
    for a in range(M):
        ca = mp.conj(v[a])
        Hv.append(sum(ca * v[b] * v[b] for b in adj[a]))
    mu = mp.re(sum(mp.conj(v[k]) * Hv[k] for k in range(M))) / h
    res = mp.sqrt(mp.re(sum(mp.conj(Hv[k] - mu * v[k]) *
                            (Hv[k] - mu * v[k]) for k in range(M)))) / norm
    herm = mp.mpf(0)  # H_ab = conj(v_a) v_b on symmetric adjacency
    for a in range(M):
        for b in adj[a]:
            herm = max(herm, abs(mp.conj(v[a]) * v[b] -
                                 mp.conj(mp.conj(v[b]) * v[a])))
    mu_theory = -2 * r2
    q = {"N": N, "M_edge": M,
         "norm": mp.nstr(norm, 110),
         "mean_amp2": mp.nstr(r2, 110),
         "global_closure_abs": mp.nstr(abs(gs), 20),
         "global_closure_normalized": mp.nstr(abs(gs) / h, 20),
         "local_closure_max_abs": mp.nstr(loc, 20),
         "mu": mp.nstr(mu, 110),
         "mu_theory_minus2r2": mp.nstr(mu_theory, 110),
         "mu_minus_theory_abs": mp.nstr(abs(mu - mu_theory), 20),
         "H_eigen_residual": mp.nstr(res, 20),
         "H_hermiticity_error": mp.nstr(herm, 20),
         "amplitude_convention":
             "exact binary64 lift of frozen LEGACY_NORM (no analytic formula "
             "exists for it); legacy operation v*NORM/||v|| reproduced at "
             "dps=100; phases analytic at dps=100",
         "legacy_norm_float64": LEGACY_NORM[N]}
    passed = (abs(gs) / h < GATE and loc / r2 < GATE and res < GATE and
              abs(mu - mu_theory) < GATE)
    q["qualification_passed"] = bool(passed)
    return q, passed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True, choices=[7, 8])
    args = ap.parse_args()
    N = args.n
    E, col, v = build(N)
    q, passed = qualify(N, E, v)
    out = HERE / "data" / f"N{N}_D{N}" / "C_IC100_DYN100"
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "ic100_qualification.json", "w", encoding="utf-8") as f:
        json.dump(q, f, indent=1, ensure_ascii=False)
    with open(out / "ic100_state.json", "w", encoding="utf-8") as f:
        json.dump({"N": N, "edges": E,
                   "class": [col[e] for e in E],
                   "z": [mpc_to_json(c) for c in v]}, f)
    with open(out / "ic100_readable.csv", "w", encoding="utf-8") as f:
        f.write("edge_index,i,j,class,re_110digits,im_110digits\n")
        for k, (i, j) in enumerate(E):
            f.write(f"{k},{i},{j},{col[(i, j)]},"
                    f"{mp.nstr(v[k].real, 110)},{mp.nstr(v[k].imag, 110)}\n")
    print(json.dumps({k: q[k] for k in
                      ("N", "global_closure_normalized",
                       "local_closure_max_abs", "H_eigen_residual",
                       "mu_minus_theory_abs", "qualification_passed")},
                     indent=1))
    if not passed:
        raise SystemExit("IC100 qualification FAILED - not usable as seed")


if __name__ == "__main__":
    main()
