#!/usr/bin/env python3
from pathlib import Path
import math
import numpy as np
import pandas as pd
from scipy.stats import kstest, beta
from scipy.optimize import linear_sum_assignment

OUT=Path(__file__).resolve().parent
SEED=20260914
KAPPA_SAMPLES=20000
FRAME_SAMPLES=300
RECON_SAMPLES=50

def kappa_batch(z):
    return np.abs(np.sum(z*z,axis=1))/np.sum(np.abs(z)**2,axis=1)

def exact_mean(m):
    return math.sqrt(math.pi)/2*math.gamma((m+1)/2)/math.gamma((m+2)/2)

def exact_median(m):
    return 1.0 if m==1 else math.sqrt(1-2**(-2/(m-1)))

def Pk(z,k):
    H=np.sum(np.abs(z)**2); return np.sum(z**k)/(H**(k/2))

def frame(z,K):
    p1=Pk(z,1); vals=[abs(p1)]
    for k in range(2,K+1):
        pk=Pk(z,k)
        vals.extend([abs(pk),np.angle(pk*np.conj(p1**k)) if abs(p1)>1e-14 and abs(pk)>1e-18 else np.nan])
    return np.asarray(vals,float)

def frame_err(z,rng):
    m=len(z); f1=frame(z,m)
    z2=10**rng.uniform(-2,2)*np.exp(1j*rng.uniform(-np.pi,np.pi))*z[rng.permutation(m)]
    f2=frame(z2,m); errs=[]
    for idx,(a,b) in enumerate(zip(f1,f2)):
        if np.isnan(a) or np.isnan(b): continue
        errs.append(abs(np.angle(np.exp(1j*(a-b)))) if (idx>=2 and idx%2==0) else abs(a-b))
    return max(errs) if errs else np.nan

def gauge_fixed_powers(z):
    zn=z/np.sqrt(np.sum(np.abs(z)**2)); p1=np.sum(zn)
    if abs(p1)<1e-12: return None,None
    zg=np.exp(-1j*np.angle(p1))*zn; m=len(z)
    return zg,[None]+[np.sum(zg**k) for k in range(1,m+1)]

def newton(p,m):
    e=[1+0j]
    for k in range(1,m+1):
        s=0j
        for i in range(1,k+1): s+=((-1)**(i-1))*e[k-i]*p[i]
        e.append(s/k)
    return e

def roots_from_e(e):
    m=len(e)-1
    return np.roots(np.asarray([1+0j]+[((-1)**k)*e[k] for k in range(1,m+1)],complex))

def match_error(a,b):
    C=np.abs(a.reshape(-1,1)-b.reshape(1,-1)); r,c=linear_sum_assignment(C)
    return float(np.max(C[r,c]))

def main():
    rng=np.random.default_rng(SEED); summary=[]; recrows=[]
    for N in range(2,41):
        m=N-1
        z=rng.normal(size=(KAPPA_SAMPLES,m))+1j*rng.normal(size=(KAPPA_SAMPLES,m))
        q=kappa_batch(z)
        ks_stat,ks_p=(0.0,1.0) if m==1 else kstest(q*q,beta(1,(m-1)/2).cdf)
        ferr=[]; p1abs=[]
        for _ in range(FRAME_SAMPLES):
            zz=rng.normal(size=m)+1j*rng.normal(size=m)
            ferr.append(frame_err(zz,rng)); p1abs.append(abs(Pk(zz,1)))
        recerrs=[]; fails=0
        for trial in range(RECON_SAMPLES):
            zz=rng.normal(size=m)+1j*rng.normal(size=m); zg,p=gauge_fixed_powers(zz)
            if zg is None: fails+=1; continue
            try: er=match_error(roots_from_e(newton(p,m)),zg)
            except Exception: er=np.nan; fails+=1
            recerrs.append(er); recrows.append({"N":N,"m":m,"trial":trial,"reconstruction_error":er})
        em=exact_mean(m); med=exact_median(m)
        summary.append({
            "N":N,"m":m,"kappa_samples":KAPPA_SAMPLES,
            "sample_mean_kappa":float(np.mean(q)),"exact_mean_kappa":em,"mean_error":float(np.mean(q)-em),
            "sample_median_kappa":float(np.median(q)),"exact_median_kappa":med,"median_error":float(np.median(q)-med),
            "ks_statistic_beta":float(ks_stat),"ks_pvalue_beta":float(ks_p),
            "frame_samples":FRAME_SAMPLES,"max_frame_invariance_error":float(np.nanmax(ferr)),
            "median_frame_invariance_error":float(np.nanmedian(ferr)),
            "min_abs_P1_frame_samples":float(np.min(p1abs)),
            "reconstruction_samples":RECON_SAMPLES,
            "max_reconstruction_error":float(np.nanmax(recerrs)) if recerrs else np.nan,
            "median_reconstruction_error":float(np.nanmedian(recerrs)) if recerrs else np.nan,
            "reconstruction_failures":fails
        })
    pd.DataFrame(summary).to_csv(OUT/"highN_N2_N40_summary.csv",index=False)
    pd.DataFrame(recrows).to_csv(OUT/"highN_reconstruction_conditioning.csv",index=False)

if __name__=="__main__": main()
