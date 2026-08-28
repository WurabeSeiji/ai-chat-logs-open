#!/usr/bin/env python3
from __future__ import annotations
import itertools,math
import numpy as np,pandas as pd
from scipy.spatial import cKDTree

def comb_array(n,k):
 return np.fromiter((x for c in itertools.combinations(range(n),k) for x in c),dtype=np.int16,count=math.comb(n,k)*k).reshape(-1,k)
def best(df,k,neighbors=16,batch=25000):
 z=df.z2_re.to_numpy()+1j*df.z2_im.to_numpy(); mag=np.abs(z); n=len(z); p=k//2;q=k-p;A=comb_array(n,p);B=comb_array(n,q);SB=z[B].sum(axis=1);tree=cKDTree(np.c_[SB.real,SB.imag]);br=np.inf;bi=None
 for s in range(0,len(A),batch):
  X=A[s:s+batch];SA=z[X].sum(axis=1);_,ids=tree.query(np.c_[-SA.real,-SA.imag],k=min(neighbors,len(B)))
  if ids.ndim==1:ids=ids[:,None]
  for row,a in enumerate(X):
   aset=set(map(int,a))
   for bj in ids[row]:
    b=B[int(bj)]
    if any(int(x) in aset for x in b):continue
    idx=tuple(sorted(aset|set(map(int,b))))
    if len(idx)!=k:continue
    r=float(abs(z[list(idx)].sum())/mag[list(idx)].sum())
    if r<br:br=r;bi=idx
 return br,bi
