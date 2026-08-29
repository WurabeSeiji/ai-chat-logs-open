#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第 2 自己無撞着状態テスト①（磨き）：τ=40000 の飽和状態を、振幅込み K の固有モードへの
混合反復（make_parent 段階 3 と同じ操作、β=0.5、スケール保存、正規化なし）で追い込めるかを試す。
ダメ元テスト：収束しない・別の状態へ逃げる可能性あり。結果は成否によらず data/N*/second_state/ に保存。"""
import os, json, sys
import numpy as np
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def adjacency(N):
    E=edges(N); M=len(E); A=np.zeros((M,M))
    for a in range(M):
        for b in range(M):
            if a!=b and set(E[a])&set(E[b]): A[a,b]=1.0
    return A
def resid(A,Z):
    K=A*np.imag(np.conj(Z)[:,None]*Z[None,:]); hv=1j*(K@Z)
    mu=(np.vdot(Z,hv)/np.vdot(Z,Z)).real
    return float(np.linalg.norm(hv-mu*Z)/np.linalg.norm(Z)), float(mu), K
for N in [6,8]:
    A=adjacency(N)
    S=np.load(os.path.join(ROOT,'data',f'N{N}','states_treatment.npz'))['Z']
    pz=np.load(os.path.join(ROOT,'data',f'N{N}','parent_v.npz')); mu_parent=float(pz['mu'])
    Z=S[40000].copy(); Z0=Z.copy()
    out=os.path.join(ROOT,'data',f'N{N}','second_state'); os.makedirs(out,exist_ok=True)
    hist=[]; beta=0.5
    r,mu,K=resid(A,Z); hist.append((0,r,mu))
    for it in range(1,5001):
        w,U=np.linalg.eigh(1j*K)
        ov=np.abs(U.conj().T@Z); j=int(np.argmax(ov))
        u=U[:,j]; ph=np.vdot(u,Z)
        if abs(ph)>0: u=u*(ph/abs(ph))
        u=u*(np.linalg.norm(Z)/np.linalg.norm(u))
        Z=(1.0-beta)*Z+beta*u
        r,mu,K=resid(A,Z)
        if it%50==0 or r<1e-12: hist.append((it,r,mu))
        if r<1e-12: break
    amp2=np.abs(Z)**2
    res_final,mu_final,_=resid(A,Z)
    summ=dict(N=N,iters=it,residual_initial=hist[0][1],residual_final=res_final,
              mu_parent=mu_parent,mu_initial=hist[0][2],mu_final=mu_final,
              mu_ratio_final=mu_final/mu_parent,
              overlap_with_tau40000=float(abs(np.vdot(Z0,Z))/(np.linalg.norm(Z0)*np.linalg.norm(Z))),
              closure=float(abs((Z*Z).sum())/amp2.sum()),
              equimodular_spread=float((amp2.max()-amp2.min())/amp2.mean()),
              nonzero_edges=int((np.abs(Z)>1e-8*np.abs(Z).max()).sum()),M=len(Z),
              norm=float(np.linalg.norm(Z)))
    np.savez_compressed(os.path.join(out,'polished_state.npz'),Z=Z,Z_start=Z0)
    with open(os.path.join(out,'polish_history.csv'),'w') as f:
        f.write('iter,residual,mu\n')
        for a,b,c in hist: f.write(f'{a},{b},{c}\n')
    with open(os.path.join(out,'polish_summary.json'),'w') as f: json.dump(summ,f,indent=2)
    print(f"N={N}: iters={it} 残差 {hist[0][1]:.3e} → {res_final:.3e}  μ: {hist[0][2]:+.5f} → {mu_final:+.5f}  μ/μ_parent={mu_final/mu_parent:+.5f}  元状態との重なり={summ['overlap_with_tau40000']:.4f}  閉塞={summ['closure']:.1e}  |z|²幅={summ['equimodular_spread']:.3f}")
print("POLISH DONE")
