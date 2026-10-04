#!/usr/bin/env python3
"""対照：走行で測った近点移動（一周あたり）を、E–M の解析値（Delaunay 平均 ∂⟨δε⟩/∂L × 周期）と比べる。
解析値は ε の 22 項そのものから計算する（走行とは独立な計算。v1〜v7 の数値は使わない）。"""
import csv, gzip, math, sys
import eight_state_exchange as E

def orbit_avg_delta_eps(a, e, sa, sb, nq=256):
    e = abs(e); L = math.sqrt(E.KAPPA * a * (1 - e * e)); m = L / (2 * E.GM)
    tot = 0.0
    for j in range(nq):
        Ea = 2 * math.pi * (j + 0.5) / nq
        r = a * (1 - e * math.cos(Ea)); rdot = math.sqrt(E.KAPPA / a) * e * math.sin(Ea) / (1 - e * math.cos(Ea))
        epsK = 0.5 * (rdot**2 + (2 * E.GM * m / r)**2) - E.KAPPA / r
        tot += (E.eps_full(r, rdot, m, sa, sb) - epsK) * (1 - e * math.cos(Ea))
    return tot / nq

def analytic_precession_per_orbit(a, e, sa=1, sb=-1):
    L = math.sqrt(E.KAPPA * a * (1 - e * e)); h = 1e-3
    f1 = orbit_avg_delta_eps(a, e + h, sa, sb); f2 = orbit_avg_delta_eps(a, abs(e - h), sa, sb)
    if e > h:
        one_over_e_dde = (f1 - f2) / (2 * h * e)
    else:
        fp = orbit_avg_delta_eps(a, 1.5 * h, sa, sb); fm = orbit_avg_delta_eps(a, 0.5 * h, sa, sb)
        one_over_e_dde = ((fp - fm) / h) / h
    rate = one_over_e_dde * (-L / (E.KAPPA * a))
    T = 2 * math.pi * a**1.5 / math.sqrt(E.KAPPA)
    return rate * T

for p in sys.argv[1:]:
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        rows = list(csv.DictReader(fh))
    a = float(rows[0]["a_orb"]); e = float(rows[0]["ecc"]); n = int(rows[-1]["node"]); turns = n * 23 / 124
    dv = 0.0; prev = float(rows[0]["varpi"])
    for r in rows[1:]:
        cur = float(r["varpi"]); dv += (cur - prev + math.pi) % (2 * math.pi) - math.pi; prev = cur
    ana = analytic_precession_per_orbit(a, e)
    note = "  [e < 1e-4: 近点方向 ϖ は数値的に定義できず、測定値は無意味]" if e < 1e-4 else ""
    print(f"{p}: a = {a:.6f}, e = {e:.3e}; measured per orbit = {dv/turns:.6e} rad; analytic (E–M, Delaunay) = {ana:.6e} rad; ratio = {dv/turns/ana if ana else float('nan'):.4f}{note}")
