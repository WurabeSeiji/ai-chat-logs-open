#!/usr/bin/env python3
"""Meet-in-the-middle k=5/6 subset search used by the closure surveys.

Unlike the orphaned historical helper, all paths and requested k are explicit.
The search rule is the preserved cKDTree nearest-neighbor rule (16 neighbors).
"""
from __future__ import annotations
import argparse,itertools,math,time
from pathlib import Path
import numpy as np,pandas as pd
from scipy.spatial import cKDTree

def comb_array(n,k):
    return np.fromiter((x for c in itertools.combinations(range(n),k) for x in c),dtype=np.int16,count=math.comb(n,k)*k).reshape(-1,k)

def best(df,k,neighbors=16,batch=25000):
    z=df.z2_re.to_numpy()+1j*df.z2_im.to_numpy(); mag=np.abs(z); n=len(z)
    p=k//2; q=k-p; A=comb_array(n,p); B=comb_array(n,q); SB=z[B].sum(axis=1)
    tree=cKDTree(np.c_[SB.real,SB.imag]); br=np.inf; bi=None
    for s in range(0,len(A),batch):
        X=A[s:s+batch]; SA=z[X].sum(axis=1)
        _,ids=tree.query(np.c_[-SA.real,-SA.imag],k=min(neighbors,len(B)))
        if ids.ndim==1: ids=ids[:,None]
        for row,a in enumerate(X):
            aset=set(map(int,a))
            for bj in ids[row]:
                b=B[int(bj)]
                if any(int(x) in aset for x in b): continue
                idx=tuple(sorted(aset|set(map(int,b))))
                if len(idx)!=k: continue
                rr=abs(z[list(idx)].sum())/mag[list(idx)].sum()
                if rr<br: br=float(rr); bi=idx
    return br,bi

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--N',type=int,required=True); ap.add_argument('--input',type=Path,required=True); ap.add_argument('--ks',type=int,nargs='+',default=[5,6]); ap.add_argument('--out',type=Path,required=True); ap.add_argument('--neighbors',type=int,default=16); a=ap.parse_args()
    df=pd.read_csv(a.input); rows=[]
    for k in a.ks:
        t=time.time(); r,idx=best(df,k,a.neighbors); runtime=time.time()-t
        edges=','.join(f'{int(df.iloc[e].i)}-{int(df.iloc[e].j)}' for e in idx)
        rows.append([a.N,len(df),k,r,edges,runtime]); print(a.N,k,r,edges,runtime)
    pd.DataFrame(rows,columns=['N','M','k','best_residual','edges','runtime_s']).to_csv(a.out,index=False)
if __name__=='__main__': main()
