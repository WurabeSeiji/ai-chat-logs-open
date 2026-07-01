#!/usr/bin/env python3
"""
Deterministic de Broglie ALIGNMENT-WAVELENGTH solver (single shot, no iteration).

Given the geometry (L, W) and an approximate wavelength lambda0 (or accelerating
voltage V), return the EXACT wavelength lambda* at which the localized odd-harmonic
two-slit wave is ALIGNED, i.e.  r_k/lambda* = m*/2  (half-integer): every odd
harmonic then shares a common phase (2 pi n r_k/lambda* = pi n m*) and a sharp
localized fringe forms.  Alignment condition = Paper 2 (Zenodo 10.5281/zenodo.21035831,
fig_oddharm_interference.py: dev = |ratio - round(2*ratio)/2|).

Closed form (deterministic; NO scan, NO root-finding):
    r_k     = sqrt(L^2 + (W/2)^2)          (centred source, y = 0)
    m*      = round(2 r_k / lambda0)
    lambda* = 2 r_k / m*

The relative adjustment |lambda* - lambda0|/lambda0 <= 1/(2 m*) is tiny
(~1e-9 at the physical scale L=5cm, W=5um): the slightest wavelength offset
destroys the interference, so the alignment condition pins lambda to
machine-level precision.  (Building lambda* = 2 r_k/m* directly also sidesteps
the float64 loss that r_k/lambda ~ 1e9 would otherwise cause.)

Standalone:
    python3 debroglie_align_lambda.py --L 0.05 --W 5e-6 --V 150
    python3 debroglie_align_lambda.py --L 0.05 --W 5e-6 --lam0 1.0e-10 --N 9
Importable:
    from debroglie_align_lambda import align_lambda, align_lambda_from_V
"""
import argparse
import numpy as np

# physical constants (CODATA)
h  = 6.62607015e-34
me = 9.1093837015e-31
e  = 1.602176634e-19


def lam_from_V(V):
    """de Broglie (approximate) wavelength from accelerating voltage."""
    return h / np.sqrt(2.0 * me * e * V)


def rk(L, W):
    """source-to-slit distance for the centred source (y=0)."""
    return np.sqrt(L**2 + (W / 2.0)**2)


def align_lambda(L, W, lam0):
    """Deterministic aligned wavelength for the approximate lam0.

    Returns a dict:
      lam0     : input approximate wavelength [m]
      m        : aligned order  m* = round(2 r_k/lam0)
      lam_star : aligned wavelength  lambda* = 2 r_k/m*  [m]
      rk       : r_k [m]
      rel      : lambda*/lam0 - 1  (relative adjustment)
      half_int : r_k/lambda*  (= m*/2, the half-integer alignment value)
    """
    r = rk(L, W)
    m = int(round(2.0 * r / lam0))
    lam_star = 2.0 * r / m
    return {"lam0": lam0, "m": m, "lam_star": lam_star, "rk": r,
            "rel": lam_star / lam0 - 1.0, "half_int": r / lam_star}


def align_lambda_from_V(L, W, V):
    """Same as align_lambda but starting from an accelerating voltage V."""
    d = align_lambda(L, W, lam_from_V(V))
    d["V"] = V
    return d


def tol_band(N):
    """Paper-2 alignment tolerance band ~ 1/(2N) (reporting only)."""
    return 1.0 / (2.0 * N)


def _main():
    ap = argparse.ArgumentParser(description="deterministic aligned-wavelength solver")
    ap.add_argument("--L", type=float, required=True, help="source-slit distance [m]")
    ap.add_argument("--W", type=float, required=True, help="slit separation [m]")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--V", type=float, help="accelerating voltage [V]  (lam0 = h/sqrt(2 me e V))")
    g.add_argument("--lam0", type=float, help="approximate wavelength lambda0 [m]")
    ap.add_argument("--N", type=int, default=None, help="highest odd harmonic (reporting only)")
    args = ap.parse_args()

    lam0 = lam_from_V(args.V) if args.V is not None else args.lam0
    d = align_lambda(args.L, args.W, lam0)
    print("r_k              = %.12e m" % d["rk"])
    if args.V is not None:
        print("V                = %g V" % args.V)
    print("lambda0 (approx) = %.12e m   (%.6f nm)" % (d["lam0"], d["lam0"] * 1e9))
    print("m* = round(2 r_k/lambda0) = %d" % d["m"])
    print("lambda* (aligned)= %.12e m   (%.6f nm)" % (d["lam_star"], d["lam_star"] * 1e9))
    print("Delta_lambda = lambda* - lambda0 = %+.3e m   (abs adjustment)" % (d["lam_star"] - d["lam0"]))
    print("rel adjust  lambda*/lambda0 - 1 = %+.3e" % d["rel"])
    print("check  r_k/lambda* = %.6f  (half-integer m*/2 = %.1f)" % (d["half_int"], d["m"] / 2.0))
    if args.N is not None:
        # interference tolerance in lambda: dev<1/(2N) in ratio r_k/lambda ->
        #   |d lambda|_full = (1/N) * lambda^2 / r_k
        dlam_tol = d["lam0"]**2 / d["rk"] / args.N
        print("N = %d   tolerance band 1/(2N) = %.3e (ratio units)" % (args.N, tol_band(args.N)))
        print("         interference band FULL width in lambda ~ lambda^2/(N r_k) = %.3e m" % dlam_tol)


if __name__ == "__main__":
    _main()
