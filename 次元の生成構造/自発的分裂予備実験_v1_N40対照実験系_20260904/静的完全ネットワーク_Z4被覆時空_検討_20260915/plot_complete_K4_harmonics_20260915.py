#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
U^4=I, U=i.
State S_k=(U^k, U^(m k)), k=0,1,2,3.
Build the COMPLETE graph K4: 4 vertices, all 6 pair relations.

Edge length is the chord distance in normalized two-wave complex state space C^2:
    d_ij = sqrt(|U^i-U^j|^2 + |U^(mi)-U^(mj)|^2)

This intentionally ignores absolute amplitude and center phase.
The 3D coordinates are obtained by classical MDS from the exact 6 pair distances.
For 4 points this reproduces the metric exactly up to floating-point error.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from pathlib import Path

OUT = Path(__file__).resolve().parent

def states(m, n=4):
    U = np.exp(2j*np.pi/n)
    return np.array([[U**k, U**(m*k)] for k in range(n)], dtype=complex)

def pairwise_dist(S):
    n=len(S); D=np.zeros((n,n))
    for i,j in combinations(range(n),2):
        D[i,j]=D[j,i]=np.linalg.norm(S[i]-S[j])
    return D

def classical_mds(D, ndim=3):
    n=D.shape[0]
    J=np.eye(n)-np.ones((n,n))/n
    B=-0.5*J@(D**2)@J
    vals,vecs=np.linalg.eigh(B)
    idx=np.argsort(vals)[::-1]
    vals,vecs=vals[idx],vecs[:,idx]
    keep=vals>1e-12
    vals,vecs=vals[keep][:ndim],vecs[:,keep][:,:ndim]
    X=vecs*np.sqrt(vals)
    if X.shape[1]<ndim:
        X=np.pad(X,((0,0),(0,ndim-X.shape[1])))
    return X

def cfmt(z):
    a,b=float(np.real(z)),float(np.imag(z))
    if abs(a)<1e-10 and abs(b-1)<1e-10:return "i"
    if abs(a)<1e-10 and abs(b+1)<1e-10:return "-i"
    if abs(b)<1e-10 and abs(a-1)<1e-10:return "1"
    if abs(b)<1e-10 and abs(a+1)<1e-10:return "-1"
    return f"{z:.4g}"

for m in (2,3):
    S=states(m)
    D=pairwise_dist(S)
    X=classical_mds(D,3)
    rows=[]
    for i,j in combinations(range(4),2):
        temporal=((j-i)%4==1) or ((i-j)%4==1)
        rows.append((i,j,D[i,j],D[i,j]**2,temporal))
    df=pd.DataFrame(rows,columns=["i","j","distance","distance_squared","temporal_neighbor"])
    df.to_csv(OUT/f"edges_m{m}.csv",index=False)

    fig=plt.figure(figsize=(9,8))
    ax=fig.add_subplot(111,projection="3d")
    for i,j,d,d2,temporal in rows:
        ls="-" if temporal else "--"
        lw=2.4 if temporal else 1.6
        ax.plot([X[i,0],X[j,0]],[X[i,1],X[j,1]],[X[i,2],X[j,2]],
                linestyle=ls,linewidth=lw)
        mid=(X[i]+X[j])/2
        ax.text(*mid,f"d={d:.6g}",fontsize=9)
    ax.scatter(X[:,0],X[:,1],X[:,2],s=90)
    for k in range(4):
        ax.text(*X[k],f"S{k}\n({cfmt(S[k,0])},{cfmt(S[k,1])})",fontsize=10)
    ax.set_title(f"Complete state network K4: base + harmonic m={m}\nAll 6 relations")
    ax.set_xlabel("MDS-1"); ax.set_ylabel("MDS-2"); ax.set_zlabel("MDS-3")
    ax.set_box_aspect((1,1,1))
    fig.tight_layout()
    fig.savefig(OUT/f"K4_complete_network_m{m}.png",dpi=220,bbox_inches="tight")
    fig.savefig(OUT/f"K4_complete_network_m{m}.svg",bbox_inches="tight")
    plt.close(fig)
