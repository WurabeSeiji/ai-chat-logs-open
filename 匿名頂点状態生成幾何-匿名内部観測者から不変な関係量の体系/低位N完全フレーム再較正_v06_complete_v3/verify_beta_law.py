#!/usr/bin/env python3
"""
Independent validation suite for Paper 2.

Reproduces four claims used in the paper:
1) same-seed regeneration of abs_P1 and kappa against the stored 10,000-row sample file,
2) independent implementation / independent-seed validation of invariance and Newton reconstruction,
3) direct algebraic validation of kappa = (lambda1-lambda2)/(lambda1+lambda2),
4) KS tests of kappa^2 ~ Beta(1,(m-1)/2) for m=2,3,4,5,6,39.

Default KS sample size is 1,000,000 per m, matching the reported validation.
"""

from pathlib import Path
import argparse
import json
import math
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import kstest, beta

OUT = Path(__file__).resolve().parent

def kappa(z):
    z = np.asarray(z, dtype=np.complex128)
    return np.abs(np.sum(z*z, axis=-1)) / np.sum(np.abs(z)**2, axis=-1)

def abs_p1(z):
    z = np.asarray(z, dtype=np.complex128)
    H = np.sum(np.abs(z)**2, axis=-1)
    return np.abs(np.sum(z, axis=-1)) / np.sqrt(H)

def pk_scalar(z, k):
    H = np.sum(np.abs(z)**2)
    return np.sum(z**k)/(H**(k/2))

def frame_scalar(z, K):
    p1 = pk_scalar(z, 1)
    vals = [abs(p1)]
    for k in range(2, K+1):
        pk = pk_scalar(z, k)
        vals.extend([abs(pk), np.angle(pk*np.conj(p1**k))])
    return np.asarray(vals, float)

def gauge_fixed_powers(z):
    z = np.asarray(z, dtype=np.complex128)
    H = np.sum(np.abs(z)**2)
    zn = z/np.sqrt(H)
    p1 = np.sum(zn)
    if abs(p1) < 1e-12:
        return None, None
    zg = np.exp(-1j*np.angle(p1))*zn
    m = len(z)
    p = [None] + [np.sum(zg**k) for k in range(1, m+1)]
    return zg, p

def newton_to_e(p, m):
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

def wrapped_angle_error(a,b):
    return abs(np.angle(np.exp(1j*(a-b))))

def same_seed_regeneration():
    stored = pd.read_csv(OUT/"lowN_v06_samples.csv")
    rng = np.random.default_rng(20260914)
    rows = []
    for N in range(2,7):
        m=N-1
        for trial in range(2000):
            z = rng.normal(size=m) + 1j*rng.normal(size=m)
            if np.sum(np.abs(z)**2) < 1e-14:
                continue
            # consume the exact same RNG calls as the primary program
            _lam = 10**rng.uniform(-2,2)
            _theta = rng.uniform(-np.pi,np.pi)
            _perm = rng.permutation(m)
            rows.append((N,m,trial,float(abs_p1(z)),float(kappa(z))))
    regen = pd.DataFrame(rows, columns=["N","m","trial","abs_P1","kappa"])
    merged = stored[["N","m","trial","abs_P1","kappa"]].merge(
        regen, on=["N","m","trial"], suffixes=("_stored","_regen"), validate="one_to_one"
    )
    return {
        "rows": int(len(merged)),
        "max_abs_P1_difference": float(np.max(np.abs(merged.abs_P1_stored-merged.abs_P1_regen))),
        "max_kappa_difference": float(np.max(np.abs(merged.kappa_stored-merged.kappa_regen))),
        "exact_abs_P1_fraction": float(np.mean(merged.abs_P1_stored.values == merged.abs_P1_regen.values)),
        "exact_kappa_fraction": float(np.mean(merged.kappa_stored.values == merged.kappa_regen.values)),
    }

def independent_seed_crosscheck(seed=12345, samples_per_N=2000):
    rng = np.random.default_rng(seed)
    out = []
    for N in range(2,7):
        m=N-1
        max_inv = 0.0
        max_rec = 0.0
        kappas = []
        for _ in range(samples_per_N):
            z = rng.normal(size=m)+1j*rng.normal(size=m)
            c = 10**rng.uniform(-2,2)*np.exp(1j*rng.uniform(-np.pi,np.pi))
            perm = rng.permutation(m)
            z2 = c*z[perm]
            f1, f2 = frame_scalar(z,m), frame_scalar(z2,m)
            for idx,(a,b) in enumerate(zip(f1,f2)):
                err = wrapped_angle_error(a,b) if (idx>=2 and idx%2==0) else abs(a-b)
                max_inv = max(max_inv, float(err))
            zg,p = gauge_fixed_powers(z)
            if zg is not None:
                roots = roots_from_e(newton_to_e(p,m))
                max_rec = max(max_rec, match_error(roots,zg))
            kappas.append(float(kappa(z)))
        out.append({
            "N":N, "m":m, "samples":samples_per_N,
            "max_invariance_error":max_inv,
            "max_reconstruction_error":max_rec,
            "mean_kappa":float(np.mean(kappas)),
            "median_kappa":float(np.median(kappas)),
        })
    return out

def wishart_identity_check(seed=24680, samples=2000):
    rng = np.random.default_rng(seed)
    max_err = 0.0
    for _ in range(samples):
        m = int(rng.integers(2,40))
        z = rng.normal(size=m)+1j*rng.normal(size=m)
        x, y = z.real, z.imag
        W = np.array([[np.dot(x,x),np.dot(x,y)],[np.dot(x,y),np.dot(y,y)]], float)
        ev = np.linalg.eigvalsh(W)
        kw = (ev[1]-ev[0])/(ev[1]+ev[0])
        max_err = max(max_err, abs(float(kappa(z))-float(kw)))
    return {"samples":samples, "max_error":float(max_err)}

def sample_kappas_complex_gaussian(m, n, rng, chunk=20000):
    q = np.empty(n, float)
    pos = 0
    while pos < n:
        b = min(chunk, n-pos)
        z = rng.normal(size=(b,m)) + 1j*rng.normal(size=(b,m))
        q[pos:pos+b] = kappa(z)
        pos += b
    return q

def exact_mean(m):
    return math.sqrt(math.pi)/2 * math.gamma((m+1)/2)/math.gamma((m+2)/2)

def exact_median(m):
    return math.sqrt(1 - 2**(-2/(m-1)))

def beta_law_tests(n=1_000_000, seed=97531):
    results = []
    for m in [2,3,4,5,6,39]:
        rng = np.random.default_rng(seed+m)
        q = sample_kappas_complex_gaussian(m,n,rng)
        y = q*q
        a, b = 1.0, (m-1)/2
        stat,pvalue = kstest(y, beta(a,b).cdf)
        results.append({
            "m":m,
            "samples":n,
            "ks_statistic":float(stat),
            "ks_pvalue":float(pvalue),
            "sample_mean_kappa":float(np.mean(q)),
            "exact_mean_kappa":float(exact_mean(m)),
            "mean_difference":float(np.mean(q)-exact_mean(m)),
            "sample_median_kappa":float(np.median(q)),
            "exact_median_kappa":float(exact_median(m)),
            "median_difference":float(np.median(q)-exact_median(m)),
        })
    return results

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks-samples", type=int, default=1_000_000)
    args = ap.parse_args()

    result = {
        "same_seed_regeneration": same_seed_regeneration(),
        "independent_seed_crosscheck": independent_seed_crosscheck(),
        "wishart_identity_check": wishart_identity_check(),
        "beta_law_tests": beta_law_tests(args.ks_samples),
    }
    (OUT/"beta_law_validation.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    pd.DataFrame(result["independent_seed_crosscheck"]).to_csv(
        OUT/"independent_crosscheck_summary.csv", index=False
    )
    pd.DataFrame(result["beta_law_tests"]).to_csv(
        OUT/"beta_law_ks_summary.csv", index=False
    )
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
