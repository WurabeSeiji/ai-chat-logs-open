#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-shot generator for the complete U^4=I two-wave study.

For m=2 and m=3 it generates:
1) base wave + harmonic wave
2) composite waveform
3) two complex vectors at all four U^4 states
4) representative two-vector complex diagram at theta=pi/4
5) complete K4 state network with all 6 relations
6) CSV edge tables

Definitions
-----------
U = exp(2*pi*i/4) = i
S_k = (U^k, U^(m*k)),  k=0,1,2,3

Edge metric in normalized two-wave complex state space C^2:
d_ij = sqrt(|U^i-U^j|^2 + |U^(m i)-U^(m j)|^2)

The 3D K4 coordinates are reconstructed from all six exact distances by
classical MDS. Thus the graph layout preserves the edge metric instead of
using an arbitrary schematic placement.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations

OUT = Path(__file__).resolve().parent

def states(m,n=4):
    U=np.exp(2j*np.pi/n)
    return np.array([[U**k,U**(m*k)] for k in range(n)],dtype=complex)

def pairwise_dist(S):
    n=len(S); D=np.zeros((n,n))
    for i,j in combinations(range(n),2):
        D[i,j]=D[j,i]=np.linalg.norm(S[i]-S[j])
    return D

def classical_mds(D,ndim=3):
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
    if abs(b)<1e-10 and abs(a-1)<1e-10:return "1"
    if abs(b)<1e-10 and abs(a+1)<1e-10:return "-1"
    if abs(a)<1e-10 and abs(b-1)<1e-10:return "i"
    if abs(a)<1e-10 and abs(b+1)<1e-10:return "-i"
    return f"{z:.4g}"

def edge_df(m):
    S=states(m); rows=[]
    for i,j in combinations(range(4),2):
        d=np.linalg.norm(S[i]-S[j])
        temporal=((j-i)%4==1) or ((i-j)%4==1)
        rows.append(dict(m=m,i=i,j=j,distance=d,distance_squared=d*d,
                         temporal_neighbor=temporal))
    return pd.DataFrame(rows)

def savefig(fig,stem):
    fig.tight_layout()
    fig.savefig(OUT/f"{stem}.png",dpi=220,bbox_inches="tight")
    fig.savefig(OUT/f"{stem}.svg",bbox_inches="tight")
    plt.close(fig)

def plot_waves(m):
    th=np.linspace(-2*np.pi,2*np.pi,2000)
    y1=np.cos(th); ym=np.cos(m*th); yc=y1+ym
    fig,ax=plt.subplots(figsize=(10,5.5))
    ax.plot(th,y1,label=r"base: $\cos\theta$")
    ax.plot(th,ym,label=rf"harmonic: $\cos({m}\theta)$")
    ax.set(xlabel=r"$\theta$",ylabel="amplitude",title=f"Base wave and harmonic m={m}")
    ax.grid(True,alpha=.25); ax.legend()
    savefig(fig,f"waves_base_harmonic_m{m}")
    fig,ax=plt.subplots(figsize=(10,5.5))
    ax.plot(th,yc,label=rf"$\cos\theta+\cos({m}\theta)$")
    ax.axhline(0,linewidth=.8)
    ax.set(xlabel=r"$\theta$",ylabel="composite amplitude",
           title=f"Composite waveform: base + harmonic m={m}")
    ax.grid(True,alpha=.25); ax.legend()
    savefig(fig,f"composite_wave_m{m}")

def plot_vectors(m):
    U=1j
    fig,axs=plt.subplots(2,2,figsize=(9,9))
    for k,ax in enumerate(axs.ravel()):
        z1=U**k; z2=U**(m*k)
        ax.axhline(0,linewidth=.8); ax.axvline(0,linewidth=.8)
        ax.add_patch(plt.Circle((0,0),1,fill=False,linestyle="--",linewidth=1))
        ax.arrow(0,0,z1.real,z1.imag,width=.015,length_includes_head=True)
        ax.arrow(0,0,z2.real,z2.imag,width=.015,length_includes_head=True)
        ax.text(z1.real*1.08,z1.imag*1.08,"base")
        ax.text(z2.real*1.08,z2.imag*1.08,f"m={m}")
        ax.set_xlim(-1.35,1.35); ax.set_ylim(-1.35,1.35); ax.set_aspect("equal")
        ax.set_title(f"k={k}: ({cfmt(z1)}, {cfmt(z2)})")
        ax.set_xlabel("Re"); ax.set_ylabel("Im"); ax.grid(True,alpha=.2)
    fig.suptitle(f"Two complex vectors at all U^4 states: harmonic m={m}",y=.98)
    savefig(fig,f"complex_vectors_all_states_m{m}")

    th=np.pi/4; z1=np.exp(1j*th); z2=np.exp(1j*m*th)
    fig,ax=plt.subplots(figsize=(6.5,6.5))
    ax.axhline(0,linewidth=.8); ax.axvline(0,linewidth=.8)
    ax.add_patch(plt.Circle((0,0),1,fill=False,linestyle="--",linewidth=1))
    ax.arrow(0,0,z1.real,z1.imag,width=.015,length_includes_head=True)
    ax.arrow(0,0,z2.real,z2.imag,width=.015,length_includes_head=True)
    ax.text(z1.real*1.08,z1.imag*1.08,r"$z=e^{i\theta}$")
    ax.text(z2.real*1.08,z2.imag*1.08,rf"$z^{m}=e^{{i{m}\theta}}$")
    ax.set_xlim(-1.35,1.35); ax.set_ylim(-1.35,1.35); ax.set_aspect("equal")
    ax.set_xlabel("Re"); ax.set_ylabel("Im")
    ax.set_title(rf"Two complex vectors, m={m}, $\theta=\pi/4$")
    ax.grid(True,alpha=.2)
    savefig(fig,f"complex_vectors_theta_pi4_m{m}")

def plot_network(m):
    S=states(m); D=pairwise_dist(S); X=classical_mds(D)
    df=edge_df(m)
    fig=plt.figure(figsize=(9,8)); ax=fig.add_subplot(111,projection="3d")
    for _,r in df.iterrows():
        i,j=int(r.i),int(r.j)
        ls="-" if r.temporal_neighbor else "--"
        lw=2.4 if r.temporal_neighbor else 1.6
        ax.plot([X[i,0],X[j,0]],[X[i,1],X[j,1]],[X[i,2],X[j,2]],
                linestyle=ls,linewidth=lw)
        mid=(X[i]+X[j])/2
        ax.text(*mid,f"d={r.distance:.6g}",fontsize=9)
    ax.scatter(X[:,0],X[:,1],X[:,2],s=90)
    for k in range(4):
        ax.text(*X[k],f"S{k}\n({cfmt(S[k,0])},{cfmt(S[k,1])})",fontsize=10)
    ax.set_title(f"Complete state network K4: base + harmonic m={m}\nAll 6 relations")
    ax.set_xlabel("MDS-1"); ax.set_ylabel("MDS-2"); ax.set_zlabel("MDS-3")
    ax.set_box_aspect((1,1,1))
    savefig(fig,f"K4_complete_network_m{m}")

for m in (2,3):
    plot_waves(m)
    plot_vectors(m)
    plot_network(m)
    edge_df(m).to_csv(OUT/f"edges_m{m}.csv",index=False)

pd.concat([edge_df(2),edge_df(3)],ignore_index=True).to_csv(
    OUT/"all_edges_m2_m3.csv",index=False
)
print("generated all PNG/SVG/CSV outputs in",OUT)
