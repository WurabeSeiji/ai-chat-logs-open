#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Condition B: IC64 (exact binary64 lift) + Dynamics100 (mpmath dps=100).

IC lift: every real/imag component x of the float64 parent is lifted via
num, den = float(x).as_integer_ratio(); mp.mpf(num)/mp.mpf(den).
The lift is verified (mathematically identical binary64 values).

Dynamics: H_ef = A_ef conj(z_e) z_f, z' = Q exp(-i(2pi/D)E) Q^H z with
mpmath.eighe at mp.dps = 100. pi is 100-digit. No numpy in the state path.

Checkpoint: exact state (_mpf_ tuples) saved every CKPT steps; resumable.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, platform, time
from pathlib import Path
import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parents[1]
SWEEP = HERE.parent
STEPS = 2000
CKPT = 50
mp.mp.dps = 100


def selftest():
    one = mp.mpf(1)
    eff = None
    for k in (80, 90, 95, 100, 105, 110):
        if one + mp.mpf(10) ** (-k) == one:
            eff = k
            break
    Ht = mp.matrix(3, 3)
    ent = [[mp.mpc(2, 0), mp.mpc(1, 2), mp.mpc(0, 0)],
           [mp.mpc(1, -2), mp.mpc(3, 0), mp.mpc(0, 1)],
           [mp.mpc(0, 0), mp.mpc(0, -1), mp.mpc(1, 0)]]
    for i in range(3):
        for j in range(3):
            Ht[i, j] = ent[i][j]
    E, Q = mp.eighe(Ht)
    R = Ht * Q - Q * mp.diag(E)
    resid = max(abs(R[i, j]) for i in range(3) for j in range(3))
    return {
        "precision_digits_requested": 100,
        "first_k_with_1_plus_10^-k_equal_1": eff,
        "1+1e-80 != 1": bool(one + mp.mpf("1e-80") != one),
        "1+1e-110 == 1": bool(one + mp.mpf("1e-110") == one),
        "pi_string_110digits": mp.nstr(mp.pi, 110),
        "eighe_selftest_max_residual": mp.nstr(resid, 8),
        "eigensolver_backend": f"mpmath.eighe {mp.__version__}",
        "python": platform.python_version(),
        "mpmath_version": mp.__version__,
        "platform": platform.platform(),
        "machine": platform.machine()}


def adjacency_sets(N):
    E = [(i, j) for i in range(N) for j in range(i + 1, N)]
    M = len(E)
    adj = [[] for _ in range(M)]
    for a in range(M):
        sa = set(E[a])
        for b in range(M):
            if b != a and sa.intersection(E[b]):
                adj[a].append(b)
    return M, adj


def lift64(x):
    num, den = float(x).as_integer_ratio()
    return mp.mpf(num) / mp.mpf(den)


def mpc_to_json(c):
    tr, ti = c.real._mpf_, c.imag._mpf_
    return [[tr[0], str(tr[1]), tr[2], tr[3]],
            [ti[0], str(ti[1]), ti[2], ti[3]]]


def json_to_mpc(t):
    re = mp.mpf((t[0][0], int(t[0][1]), t[0][2], t[0][3]))
    im = mp.mpf((t[1][0], int(t[1][1]), t[1][2], t[1][3]))
    return mp.mpc(re, im)


def observe(z, p, q, M, h0):
    h = mp.re(sum(mp.conj(z[k]) * z[k] for k in range(M)))
    pz = sum(p[k] * z[k] for k in range(M))
    qz = sum(q[k] * z[k] for k in range(M))
    zp = [z[k] - p[k] * pz - q[k] * qz for k in range(M)]
    hp = mp.re(sum(mp.conj(x) * x for x in zp))
    f = hp / h
    zz = sum(z[k] * z[k] for k in range(M))
    s2 = [mp.re(mp.conj(z[k]) * z[k]) for k in range(M)]
    ssum = sum(s2)
    pr = ssum ** 2 / sum(x * x for x in s2)
    amps = [mp.sqrt(x) for x in s2]
    amean = sum(amps) / M
    astd = mp.sqrt(sum((a - amean) ** 2 for a in amps) / M)
    return {"f": f, "log10f": (mp.log10(f) if f > 0 else None),
            "h": h, "drift": abs(h - h0) / h0, "closure": abs(zz) / h,
            "pr": pr, "amp_min": min(amps), "amp_max": max(amps),
            "amp_std": astd}


def one_step(z, M, adj, dt):
    H = mp.matrix(M, M)
    for a in range(M):
        ca = mp.conj(z[a])
        for b in adj[a]:
            H[a, b] = ca * z[b]
    w, Q = mp.eighe(H)
    y = [sum(mp.conj(Q[k, j]) * z[k] for k in range(M)) for j in range(M)]
    for j in range(M):
        y[j] = mp.exp(mp.mpc(0, -1) * dt * w[j]) * y[j]
    return [sum(Q[k, j] * y[j] for j in range(M)) for k in range(M)]


def make_plane(z_init, M):
    p_raw = [mp.re(z_init[k]) for k in range(M)]
    np_ = mp.sqrt(sum(x * x for x in p_raw))
    p = [x / np_ for x in p_raw]
    q_raw = [mp.im(z_init[k]) for k in range(M)]
    qp = sum(q_raw[k] * p[k] for k in range(M))
    q_raw = [q_raw[k] - qp * p[k] for k in range(M)]
    nq = mp.sqrt(sum(x * x for x in q_raw))
    q = [x / nq for x in q_raw]
    return p, q


def run(N, D, z_init, outdir, steps_target, label):
    M, adj = adjacency_sets(N)
    dt = 2 * mp.pi / D
    ck = outdir / "checkpoint.json"
    rows = []
    if ck.exists():
        c = json.load(open(ck))
        t_start = c["step"]
        z = [json_to_mpc(t) for t in c["z"]]
        if (outdir / "timeseries.csv").exists():
            with open(outdir / "timeseries.csv", newline="") as f_:
                r = csv.reader(f_)
                next(r)
                rows = [row for row in r if int(row[0]) <= t_start]
    else:
        t_start = 0
        z = list(z_init)
    p, q = make_plane(z_init, M)
    h0 = mp.re(sum(mp.conj(z_init[k]) * z_init[k] for k in range(M)))

    def write_all():
        with open(outdir / "timeseries.csv", "w", newline="",
                  encoding="utf-8") as f_:
            w_ = csv.writer(f_)
            w_.writerow(["step", "tau", "Hperp_frac", "log10_Hperp_frac",
                         "H_parallel_frac", "H_total", "H_total_rel_drift",
                         "global_closure", "PR", "PR_over_M",
                         "amp_min", "amp_max", "amp_std", "finite"])
            w_.writerows(rows)

    have = {int(r[0]) for r in rows}
    t0 = time.time()
    for t in range(t_start, steps_target + 1):
        if t not in have:
            o = observe(z, p, q, M, h0)
            s = lambda x: mp.nstr(x, 110)
            rows.append([t, s(dt * t), s(o["f"]),
                         ("-inf" if o["log10f"] is None else s(o["log10f"])),
                         s(1 - o["f"]), s(o["h"]), mp.nstr(o["drift"], 8),
                         s(o["closure"]), s(o["pr"]), s(o["pr"] / M),
                         s(o["amp_min"]), s(o["amp_max"]), s(o["amp_std"]),
                         True])
        if t % CKPT == 0 or t == steps_target:
            write_all()
            with open(ck, "w") as f_:
                json.dump({"step": t, "label": label,
                           "z": [mpc_to_json(c) for c in z]}, f_)
            if t % 200 == 0:
                print(f"{label} N={N} step {t}/{steps_target} "
                      f"({time.time()-t0:.0f}s)", flush=True)
        if t < steps_target:
            z = one_step(z, M, adj, dt)
    write_all()
    with open(ck, "w") as f_:
        json.dump({"step": steps_target, "label": label,
                   "z": [mpc_to_json(c) for c in z]}, f_)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True, choices=[7, 8])
    args = ap.parse_args()
    N = args.n
    D = N
    pf = SWEEP / "data" / f"N{N}" / "parent_v.npz"
    sha = hashlib.sha256(pf.read_bytes()).hexdigest()
    v64 = np.load(pf)["v"].astype(np.complex128)
    M = len(v64)

    out = HERE / "data" / f"N{N}_D{D}" / "B_IC64_DYN100"
    out.mkdir(parents=True, exist_ok=True)

    st = selftest()
    with open(out / "precision_selftest.json", "w", encoding="utf-8") as f_:
        json.dump(st, f_, indent=1)

    z0 = []
    lift_ok = True
    max_back_diff = 0.0
    for c in v64:
        re = lift64(c.real)
        im = lift64(c.imag)
        if float(re) != c.real or float(im) != c.imag:
            lift_ok = False
        if re != mp.mpf(c.real) or im != mp.mpf(c.imag):
            lift_ok = False
        max_back_diff = max(max_back_diff,
                            abs(float(re) - c.real), abs(float(im) - c.imag))
        z0.append(mp.mpc(re, im))
    with open(out / "lift_verification.json", "w", encoding="utf-8") as f_:
        json.dump({"parent_sha256": sha, "components": 2 * M,
                   "method": "as_integer_ratio -> mp.mpf(num)/mp.mpf(den)",
                   "roundtrip_float_equal": lift_ok,
                   "max_roundtrip_abs_diff": max_back_diff,
                   "also_equal_to_mp.mpf(x)_direct": lift_ok}, f_, indent=1)
    if not lift_ok:
        raise SystemExit("exact binary64 lift verification FAILED")

    run(N, D, z0, out, STEPS, "B_IC64_DYN100")
    with open(out / "run_info.json", "w", encoding="utf-8") as f_:
        json.dump({"condition": "B_IC64_DYN100", "N": N, "D": D,
                   "steps": STEPS, "dps": 100, "parent_sha256": sha,
                   "pi": "100-digit mpmath", "dynamics": "mpmath.eighe"},
                  f_, indent=1)
    print(f"B N={N} COMPLETE")


if __name__ == "__main__":
    main()
