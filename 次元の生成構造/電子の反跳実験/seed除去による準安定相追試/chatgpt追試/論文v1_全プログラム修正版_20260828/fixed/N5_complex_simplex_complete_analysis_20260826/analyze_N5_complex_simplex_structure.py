#!/usr/bin/env python3
import pandas as pd, numpy as np, itertools
from pathlib import Path

HERE=Path(__file__).resolve().parent
df=pd.read_csv(HERE/"N5_step5000_final_complex_distance_classes.csv")
N=5
pairs=[(i,j) for i in range(N) for j in range(i+1,N)]
D2=np.zeros((N,N),complex)
for (_,r),(i,j) in zip(df.iterrows(),pairs):
    D2[i,j]=D2[j,i]=r.z2_re+1j*r.z2_im

J=np.eye(N)-np.ones((N,N))/N
B=-0.5*J@D2@J
s=np.linalg.svd(B,compute_uv=False)
print("Centered complex Gram singular values:",s)
print("rank(tol=1e-10):",np.sum(s>1e-10))
print("|z_ij^2| min/max:",np.abs(D2[np.triu_indices(N,1)]).min(),
      np.abs(D2[np.triu_indices(N,1)]).max())

def preserve(p,tol=1e-7):
    return all(abs(D2[p[i],p[j]]-D2[i,j])<tol
               for i in range(N) for j in range(i+1,N))
auts=[p for p in itertools.permutations(range(N)) if preserve(p)]
print("exact complex-distance-preserving vertex permutations:",len(auts))
print(auts)
