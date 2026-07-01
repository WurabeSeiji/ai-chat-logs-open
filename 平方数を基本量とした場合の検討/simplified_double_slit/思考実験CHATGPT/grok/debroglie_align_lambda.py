#!/usr/bin/env python3
"""
de Broglie ALIGNMENT-WAVELENGTH solver by the ABSOLUTE central-alignment A(lambda).

Central screen intensity (s=0, centred source y=0) of the N-harmonic wave:
    I(0) = | 2 sum_{n=1,3,...,N} e^{i 2 pi n r_k/lambda} |^2 = 4 K^2 A(lambda) ,
    A(lambda) = | sum_n e^{i 2 pi n r_k/lambda} |^2 / K^2 ,   K = (N+1)/2 .
A is the FRACTION of full alignment: A = 1 only when every odd harmonic shares a
common phase (a true resonance, r_k/lambda = m/2), and A < 1 otherwise.  It is an
ABSOLUTE measure (max value 4K^2 is fixed by geometry & N, not self-normalised),
so a degraded / partially-aligned pattern gives A << 1 and is correctly rejected.

Interference test (h-free; geometry + trial lambda only):
    A(lambda) >= 0.98  <=>  the centre reaches >= 98% of full alignment => INTERFERES.

Search / decision:
    window = lambda0 +/- lambda0^2/(4 r_k)  (half the resonance spacing).
    Scan A across it; count resonance BANDS (connected runs with A >= 0.98):
      exactly 1 band -> golden-section maximise A in it -> lambda' (A(lambda') ~ 1);
      0 or >= 2 bands -> nearest resonance NOT unique (lambda0 near the midpoint
                         between two resonances) -> lambda' = NA (status "NG").
    h enters ONLY via lambda0.  m* = round(2 r_k/lambda') reported POST-HOC.

Standalone:  python3 debroglie_align_lambda.py --L 0.05 --W 5e-6 --V 100 --N 9
Importable:  from debroglie_align_lambda import align_lambda_search, align_lambda
"""
import argparse
import numpy as np

h  = 6.62607015e-34
me = 9.1093837015e-31
e  = 1.602176634e-19

THR = 0.98   # A threshold: centre within 2% of full alignment => interferes


def lam_from_V(V):
    return h / np.sqrt(2.0 * me * e * V)


def rk(L, W):
    return np.sqrt(L**2 + (W / 2.0)**2)


def _A(r, n_list, K, lam):
    """Absolute central alignment  A = |sum_n e^{i 2 pi n r/lam}|^2 / K^2  (= I(0)/(4K^2))."""
    s = np.sum(np.exp(1j * 2.0 * np.pi * n_list * r / lam))
    return (s * np.conj(s)).real / (K * K)


def _golden(fun, a, b, tol, max_iter):
    """Golden-section maximisation of fun on [a,b].  Returns (x*, n_eval, iters)."""
    invphi = (np.sqrt(5.0) - 1.0) / 2.0
    c = b - invphi * (b - a)
    d = a + invphi * (b - a)
    fc, fd = fun(c), fun(d)
    n_eval, iters = 2, 0
    while (b - a) > tol and iters < max_iter:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = fun(c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = fun(d)
        n_eval += 1
        iters += 1
    return 0.5 * (a + b), n_eval, iters


def align_lambda(L, W, lam0):
    """DETERMINISTIC reference (one-shot round) -- verification only, not the method."""
    r = rk(L, W)
    m = int(round(2.0 * r / lam0))
    lam_det = 2.0 * r / m
    return {"lam0": lam0, "m": m, "lam_star": lam_det, "rk": r,
            "rel": lam_det / lam0 - 1.0, "half_int": r / lam_det}


def align_lambda_search(L, W, lam0, N, y=0.0, n_scan=201, thr=THR, tol_rel=1e-13, max_iter=100):
    """Absolute-A search at the ACTUAL source position y.  h enters ONLY via lam0.
    The alignment geometry is the mean source-slit distance r_avg(y):
      r1=sqrt(L^2+(y-W/2)^2), r2=sqrt(L^2+(y+W/2)^2), r_avg=(r1+r2)/2  (=rk at y=0).
    status "OK" or "NG" (NA)."""
    r1 = np.sqrt(L**2 + (y - W / 2.0)**2)
    r2 = np.sqrt(L**2 + (y + W / 2.0)**2)
    r = 0.5 * (r1 + r2)                                   # r_avg(y);  = rk(L,W) at y=0
    n_list = np.arange(1, N + 1, 2)
    K = len(n_list)
    lam_det = 2.0 * r / round(2.0 * r / lam0)             # verification reference
    fourK2 = 4.0 * K * K                                  # full-alignment I(0)

    if N == 1:                                            # single harmonic: always A=1
        return {"lam0": lam0, "lam_star": lam0, "m": int(round(2.0 * r / lam0)),
                "rk": r, "rel": 0.0, "iters": 0, "n_eval": 0, "scan_pts": 0,
                "window_half": 0.0, "N": N, "lam_det": lam_det, "match": abs(lam0/lam_det - 1.0),
                "A0": 1.0, "A_star": 1.0, "I0_star": fourK2, "n_bands": 1,
                "interferes0": True, "status": "OK"}

    Af = lambda lam: _A(r, n_list, K, lam)
    half = lam0 * lam0 / (4.0 * r)                        # half resonance spacing
    xs = np.linspace(lam0 - half, lam0 + half, n_scan)
    AA = np.array([Af(x) for x in xs])
    n_eval = n_scan
    A0 = Af(lam0); n_eval += 1

    # count resonance bands = connected runs with A >= thr
    above = AA >= thr
    bands, i = [], 0
    while i < n_scan:
        if above[i]:
            k = i
            while k < n_scan and above[k]:
                k += 1
            bands.append((i, k))
            i = k
        else:
            i += 1

    base = {"lam0": lam0, "rk": r, "N": N, "lam_det": lam_det, "scan_pts": n_scan,
            "window_half": half, "A0": A0, "interferes0": bool(A0 >= thr),
            "n_bands": len(bands), "n_eval": n_eval}

    if len(bands) != 1:                                  # 0 or >=2 resonances -> NA
        base.update({"lam_star": None, "m": None, "rel": None, "iters": 0,
                     "match": None, "A_star": None, "I0_star": None, "status": "NG"})
        return base

    lo, hi = bands[0]
    aa = xs[max(0, lo - 1)]
    bb = xs[min(n_scan - 1, hi)]
    lam_star, ev, iters = _golden(Af, aa, bb, tol_rel * lam0, max_iter)
    n_eval += ev
    A_star = Af(lam_star); n_eval += 1
    m = int(round(2.0 * r / lam_star))
    base.update({"lam_star": lam_star, "m": m, "rel": lam_star / lam0 - 1.0,
                 "iters": iters, "n_eval": n_eval, "match": abs(lam_star / lam_det - 1.0),
                 "A_star": A_star, "I0_star": fourK2 * A_star, "status": "OK"})
    return base


def align_lambda_from_V(L, W, V, N):
    d = align_lambda_search(L, W, lam_from_V(V), N)
    d["V"] = V
    return d


def _main():
    ap = argparse.ArgumentParser(description="de Broglie aligned-wavelength search (absolute A)")
    ap.add_argument("--L", type=float, required=True)
    ap.add_argument("--W", type=float, required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--V", type=float)
    g.add_argument("--lam0", type=float)
    ap.add_argument("--N", type=int, default=9)
    args = ap.parse_args()

    lam0 = lam_from_V(args.V) if args.V is not None else args.lam0
    d = align_lambda_search(args.L, args.W, lam0, args.N)
    K = (args.N + 1) // 2
    print("r_k               = %.12e m   (geometry, no h)" % d["rk"])
    if args.V is not None:
        print("V                 = %g V" % args.V)
    print("lambda0 (initial) = %.12e m   (%.6f nm)   [h used ONLY here]" % (lam0, lam0 * 1e9))
    print("test: A(lam)=|sum e^{i2pi n r/lam}|^2/K^2 = I(0)/(4K^2);  >= %.2f = interferes  (4K^2=%d)"
          % (THR, 4 * K * K))
    print("search window     = lam0 +/- lam0^2/(4 r_k) = +/- %.3e m  (half resonance spacing)" % d["window_half"])
    print("A(lam0)           = %.6f   (I(0)=%.3f of %d)  ->  %s"
          % (d["A0"], d["A0"] * 4 * K * K, 4 * K * K, "INTERFERES" if d["interferes0"] else "does NOT interfere"))
    print("resonance bands in window (A>=%.2f): %d   scan pts=%d   A-evals=%d"
          % (THR, d["n_bands"], d["scan_pts"], d["n_eval"]))
    if d["status"] == "NG":
        print("RESULT: lambda' = NA  (NG)  -- nearest resonance NOT unique in window")
        print("        (lambda0 near the midpoint between two resonances)")
        return
    print("lambda' (searched)= %.12e m   (%.6f nm)   A(lam')=%.6f  I(0)=%.4f of %d  (iters=%d)"
          % (d["lam_star"], d["lam_star"] * 1e9, d["A_star"], d["I0_star"], 4 * K * K, d["iters"]))
    print("observed order m* = round(2 r_k/lambda') = %d   (OUTPUT, post-hoc)" % d["m"])
    print("Delta_lambda = lambda' - lambda0 = %+.3e m   (rel %+.3e)" % (d["lam_star"] - lam0, d["rel"]))
    print("verification: analytic resonance lam_det = %.12e m,  |lam'/lam_det-1| = %.2e"
          % (d["lam_det"], d["match"]))


if __name__ == "__main__":
    _main()
