#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np, itertools, math
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
TOL=1e-6

def incidence(N,d):
    A=np.zeros((N,len(d)))
    for e,r in d.iterrows():
        A[int(r.i)-1,e]=1
        A[int(r.j)-1,e]=1
    return A

def star_span_test(mask,A,tol=1e-10):
    c=np.linalg.lstsq(A.T,mask.astype(float),rcond=None)[0]
    err=np.linalg.norm(A.T@c-mask)
    return err<tol,err

def closure_residual(q,idx):
    idx=np.asarray(idx,int)
    return abs(q[idx].sum())/np.abs(q[idx]).sum()

rows=[]
for N in range(3,17):
    d=pd.read_csv(ROOT/f"SOURCE_N{N}_step5000_final_edges.csv")
    q=d.z2_re.to_numpy()+1j*d.z2_im.to_numpy()
    A=incidence(N,d)
    M=len(d)
    for k in range(2,min(5,M-1)+1):
        if math.comb(M,k)>2_500_000:
            continue
        for idx in itertools.combinations(range(M),k):
            r=closure_residual(q,idx)
            if r<TOL:
                mask=np.zeros(M,bool); mask[list(idx)]=True
                trivial,err=star_span_test(mask,A)
                edges=",".join(f"{int(d.iloc[e].i)}-{int(d.iloc[e].j)}" for e in idx)
                rows.append([N,M,k,r,trivial,err,edges])

pd.DataFrame(rows,columns=[
"N","M","k","residual","star_span_trivial","star_span_error","edges"
]).to_csv(ROOT/"reproduced_small_subset_classification.csv",index=False)
