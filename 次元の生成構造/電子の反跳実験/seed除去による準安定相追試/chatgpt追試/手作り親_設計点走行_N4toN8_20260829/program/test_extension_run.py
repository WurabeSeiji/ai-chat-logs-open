#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 2 自己無撞着状態テスト②（延長走行）：τ=40000 の状態から力学無変更（振幅込み K、
exp((2π/124)K)、種なし・正規化なし）で +160000 step 継続。200 step ごとに残差・μ・f・位相進みを記録。
データは data/N*/second_state/extension_timeseries.csv、状態は 1000 step ごとに間引き保存。"""
import os, math, sys, csv
import numpy as np
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N=int(sys.argv[1]); EXTRA=160000; L=124
def edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
E=edges(N); M=len(E)
A=np.zeros((M,M))
for a in range(M):
    for b in range(M):
        if a!=b and set(E[a])&set(E[b]): A[a,b]=1.0
S=np.load(os.path.join(ROOT,'data',f'N{N}','states_treatment.npz'))['Z']
pz=np.load(os.path.join(ROOT,'data',f'N{N}','parent_v.npz')); v=pz['v']; mu_parent=float(pz['mu'])
p=v.real/np.linalg.norm(v.real); q=v.imag-(v.imag@p)*p; q/=np.linalg.norm(q)
Z=S[40000].copy()
out=os.path.join(ROOT,'data',f'N{N}','second_state'); os.makedirs(out,exist_ok=True)
rows=[]; keep=[]
for t in range(EXTRA+1):
    if t%200==0:
        K=A*np.imag(np.conj(Z)[:,None]*Z[None,:]); hv=1j*(K@Z)
        mu=(np.vdot(Z,hv)/np.vdot(Z,Z)).real
        r=float(np.linalg.norm(hv-mu*Z)/np.linalg.norm(Z))
        Zp=Z-p*(p@Z)-q*(q@Z); f=float(np.vdot(Zp,Zp).real/np.vdot(Z,Z).real)
        rows.append([40000+t,r,mu,mu/mu_parent,f,float(np.vdot(Z,Z).real),float(abs((Z*Z).sum()))])
    if t%1000==0: keep.append(Z.copy())
    if t==EXTRA: break
    K=A*np.imag(np.conj(Z)[:,None]*Z[None,:])
    if np.linalg.norm(K+K.T)>1e-10: raise RuntimeError('K antisymmetry failure')
    w,V=np.linalg.eigh(1j*K)
    Z=V@(np.exp(-1j*(2*math.pi/L)*w)*(V.conj().T@Z))
with open(os.path.join(out,'extension_timeseries.csv'),'w',newline='') as f:
    w2=csv.writer(f); w2.writerow(['step','selfcons_residual','mu','mu_over_mu_parent','Hperp_frac','H_total','abs_ZT_Z']); w2.writerows(rows)
np.savez_compressed(os.path.join(out,'extension_states_every1000.npz'),Z=np.array(keep))
r0=rows[0]; r1=rows[-1]
print(f"EXT N={N} done: 残差 {r0[1]:.3e}→{r1[1]:.3e}  μ/μ_parent {r0[3]:+.5f}→{r1[3]:+.5f}  f {r0[4]:.3f}→{r1[4]:.3f}  H drift={abs(r1[5]-r0[5]):.2e}")
