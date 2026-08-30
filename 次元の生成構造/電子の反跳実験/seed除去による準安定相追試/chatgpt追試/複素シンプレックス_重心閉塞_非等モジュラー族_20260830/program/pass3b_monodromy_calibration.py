#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス3b：共回転モノドロミー ρ の較正（走行後に追加。パス3 の閾値 ρ>1+1e-9 が N=9 を誤判定したため）。
既知の床（中立）親 N=4,5,7（手作り親パッケージで 40000 step 床を実測済）と、N=9,10,11 の等モジュラー親で ρ を計算し、
「床でも ρ−1 ~ 1e-4（刻み由来）」「インフレーションは ρ−1 ≳ 3e-3」という帯を確定する。
ANGLE=(2π/L)|μ| と ANGLE²/16（ρ−1 の刻み項の目安）も併記。出力：results/monodromy_calibration.csv"""
import os, csv, sys, math
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
from common import *
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NORM={4:0.816497,5:0.927173,6:1.036862,7:1.210458,8:1.301688,9:1.40,10:1.50,11:1.60}
L=124
def step_map(N,z,A):
    K=K_of(N,z,A); w,V=np.linalg.eigh(1j*K); return V@(np.exp(-1j*(2*math.pi/L)*w)*(V.conj().T@z))
def mono(N,v,A,h=1e-7):
    M=len(v); x=np.r_[v.real,v.imag]; D=np.zeros((2*M,2*M))
    for j in range(2*M):
        e=np.zeros(2*M); e[j]=h; zp=step_map(N,(x+e)[:M]+1j*(x+e)[M:],A); zm=step_map(N,(x-e)[:M]+1j*(x-e)[M:],A); D[:,j]=np.r_[(zp-zm).real,(zp-zm).imag]/(2*h)
    mu=selfconsistency(N,v,A)['mu']; phi=(2*math.pi/L)*mu; c,s=math.cos(phi),math.sin(phi)
    R=np.block([[c*np.eye(M),-s*np.eye(M)],[s*np.eye(M),c*np.eye(M)]]); ev=np.linalg.eigvals(R@D); return float(np.abs(ev).max()),ev,mu
KNOWN={4:'floor (measured 1.5e-23 @40000)',5:'floor (2.5e-19)',6:'inflating (t50=5269)',7:'floor (1.3e-18)',8:'inflating (t50=3791)',9:'floor (this package, 6.4e-18)',10:'not run',11:'not run'}
rows=[]
for N in [4,5,6,7,8,9,10,11]:
    th,cls,q=phases(N); A=adjacency(N); v=state(N,np.ones(q)); v*=NORM[N]/np.linalg.norm(v)
    rho,ev,mu=mono(N,v,A); ang=(2*math.pi/L)*abs(mu)
    rows.append(dict(N=N,q=q,norm=NORM[N],mu=f'{mu:.6f}',rho=f'{rho:.9f}',rho_minus_1=f'{rho-1:.3e}',lambda_f=f'{2*math.log(rho):.6f}',ANGLE=f'{ang:.5f}',ANGLE2_over_16=f'{ang**2/16:.2e}',n_eig_gt_1=int((np.abs(ev)>1+1e-9).sum()),known_behaviour=KNOWN[N],band=('inflating' if rho-1>3e-3 else 'neutral (rounding-floor)')))
    print(f"N={N:2d} q={q:2d} ρ−1={rho-1:.3e} λ_f={2*math.log(rho):.6f} ANGLE²/16={ang**2/16:.1e} |ev|>1: {rows[-1]['n_eig_gt_1']:3d} → {rows[-1]['band']:26s} 既知: {KNOWN[N]}")
with open(os.path.join(ROOT,'results','monodromy_calibration.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("PASS3b OK")
