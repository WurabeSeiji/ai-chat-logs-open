#!/usr/bin/env python3
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

OUT = Path(__file__).resolve().parent
SEED = 20260914
SAMPLES_PER_N = 2000

def Pk(z, k):
    z = np.asarray(z, dtype=complex)
    H = np.sum(np.abs(z)**2)
    return np.sum(z**k) / (H**(k/2))

def frame(z, K):
    p1 = Pk(z, 1)
    vals = [abs(p1)]
    for k in range(2, K+1):
        pk = Pk(z, k)
        vals.append(abs(pk))
        vals.append(np.angle(pk*np.conj(p1**k)) if abs(p1)>1e-12 and abs(pk)>1e-18 else np.nan)
    return np.asarray(vals, float)

def gauge_fix_and_powers(z, K):
    z = np.asarray(z, dtype=complex)
    H = np.sum(np.abs(z)**2)
    zn = z/np.sqrt(H)
    p1 = np.sum(zn)
    if abs(p1) < 1e-12:
        return None
    zg = np.exp(-1j*np.angle(p1))*zn
    p = [None] + [np.sum(zg**k) for k in range(1, K+1)]
    return zg, p

def newton(p, m):
    e = [1+0j]
    for k in range(1, m+1):
        s = 0j
        for i in range(1, k+1):
            s += ((-1)**(i-1))*e[k-i]*p[i]
        e.append(s/k)
    return e

def roots_from_e(e):
    m = len(e)-1
    coeff = [1+0j] + [((-1)**k)*e[k] for k in range(1, m+1)]
    return np.roots(np.asarray(coeff, complex))

def match_error(a, b):
    C = np.abs(a.reshape(-1,1)-b.reshape(1,-1))
    r, c = linear_sum_assignment(C)
    return float(np.max(C[r,c]))

def kappa(z):
    z = np.asarray(z, complex)
    return abs(np.sum(z*z))/np.sum(np.abs(z)**2)

def n3_checks(z):
    z1, z2 = z
    r = abs(z1)/abs(z2)
    d = np.angle(z1/z2)
    u = r + 1/r
    c = np.cos(d)
    w = (r - 1/r)*np.sin(d)
    S = z1/z2 + z2/z1
    e1 = z1+z2
    J2 = z1*z2/(e1**2) if abs(e1)>1e-12 else np.nan+1j*np.nan
    kap = kappa(z)
    return {
        "u":u, "c":c, "w":w,
        "relation_error":abs(w*w-(u*u-4)*(1-c*c)),
        "S_error":abs(S-(u*c+1j*w)),
        "J2_error":abs(1/J2-(2+S)) if np.isfinite(J2.real) else np.nan,
        "kappa_identity_error":abs(abs(S)-u*kap),
    }

def main():
    rng = np.random.default_rng(SEED)
    rows, n3rows = [], []
    for N in range(2, 7):
        m = N-1
        for trial in range(SAMPLES_PER_N):
            z = rng.normal(size=m) + 1j*rng.normal(size=m)
            if np.sum(np.abs(z)**2) < 1e-14:
                continue

            lam = 10**rng.uniform(-2,2)
            theta = rng.uniform(-np.pi,np.pi)
            perm = rng.permutation(m)
            z2 = lam*np.exp(1j*theta)*z[perm]

            f1, f2 = frame(z,m), frame(z2,m)
            errs = []
            for idx,(a,b) in enumerate(zip(f1,f2)):
                if np.isnan(a) or np.isnan(b):
                    continue
                is_angle = idx>=2 and idx%2==0
                errs.append(abs(np.angle(np.exp(1j*(a-b)))) if is_angle else abs(a-b))

            gf = gauge_fix_and_powers(z,m)
            if gf is None:
                recerr = np.nan
            else:
                zg,p = gf
                roots = roots_from_e(newton(p,m))
                recerr = match_error(roots,zg)

            rows.append({
                "N":N, "m":m, "trial":trial,
                "abs_P1":abs(Pk(z,1)),
                "kappa":kappa(z),
                "frame_invariance_error":max(errs) if errs else np.nan,
                "reconstruction_error_K_Nminus1":recerr
            })

            if N==3:
                d = {"trial":trial}
                d.update(n3_checks(z))
                n3rows.append(d)

    df = pd.DataFrame(rows)
    n3df = pd.DataFrame(n3rows)

    summary = []
    for N,g in df.groupby("N"):
        summary.append({
            "N":N, "m":N-1, "samples":len(g),
            "max_frame_invariance_error":g["frame_invariance_error"].max(),
            "median_frame_invariance_error":g["frame_invariance_error"].median(),
            "max_reconstruction_error_K_Nminus1":g["reconstruction_error_K_Nminus1"].max(),
            "median_reconstruction_error_K_Nminus1":g["reconstruction_error_K_Nminus1"].median(),
            "min_abs_P1":g["abs_P1"].min(),
            "fraction_abs_P1_lt_1e-3":float((g["abs_P1"]<1e-3).mean()),
            "mean_kappa":g["kappa"].mean(),
            "median_kappa":g["kappa"].median(),
        })
    sdf = pd.DataFrame(summary)

    n3summary = {
        "samples":int(len(n3df)),
        "max_relation_error":float(n3df["relation_error"].max()),
        "max_S_error":float(n3df["S_error"].max()),
        "max_J2_error":float(n3df["J2_error"].dropna().max()),
        "max_kappa_identity_error":float(n3df["kappa_identity_error"].max()),
    }

    df.to_csv(OUT/"lowN_v06_samples.csv", index=False)
    sdf.to_csv(OUT/"lowN_v06_summary.csv", index=False)
    n3df.to_csv(OUT/"N3_identity_checks.csv", index=False)
    (OUT/"N3_identity_summary.json").write_text(json.dumps(n3summary, indent=2), encoding="utf-8")

if __name__ == "__main__":
    main()
