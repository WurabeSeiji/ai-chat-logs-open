#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス7：自己無撞着の厳密な書き換え（本文 §5 の導出）の数値検証と、対称性を一切使わない親の生成。
恒等式：(iKv)_e = ½[ z̄_e (A z²)_e − z_e (A|z|²)_e ]。S_i := Σ_{k≠i} z_ik²（頂点 i の局所和）、W_i := Σ_{k≠i} |z_ik|²（頂点重み）とすると
(A z²)_e = S_i + S_j − 2 z_e²、(A|z|²)_e = W_i + W_j − 2|z_e|² なので、自己無撞着 iKv = μv は各辺 e=(i,j) について
    S_i + S_j = (2μ + W_i + W_j) · z_e² / |z_e|²        …(★)
と同値。系：S ≡ 0（局所閉塞）かつ W_i ≡ W（頂点重み均等）ならば自己無撞着で μ = −W。逆に自己無撞着かつ S≡0 なら W_i ≡ −μ。
(a) 既存の親・族メンバー・プローブ解で (★) の残差を検証。
(b) 乱数初期値から {S_i=0 (2N 実), W_i=W₀ (N 実)} だけを Newton で解き（対称性ゼロ）、自己無撞着残差・閉塞・rank・丸さ・モノドロミー ρ を測る。
    N=5..8 各 5 個。‖v‖ は同 N の親に合わせる。生成した状態は data/random_N{N}_s{k}/parent_v.npz に保存（力学走行はパス4 と同じ run_dynamics.py）。
出力：results/identity_check.csv, results/balanced_random_parents.csv"""
import os, csv, sys, math, json
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
from common import *
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NORM={4:0.816497,5:0.927173,6:1.036862,7:1.210458,8:1.301688,9:1.40}
L=124
def SW(N,v):
    E=edges(N); S=np.zeros(N,complex); W=np.zeros(N)
    for k,(i,j) in enumerate(E): S[i]+=v[k]**2; S[j]+=v[k]**2; W[i]+=abs(v[k])**2; W[j]+=abs(v[k])**2
    return S,W
def star_residual(N,v,mu):
    E=edges(N); S,W=SW(N,v)
    r=[S[i]+S[j]-(2*mu+W[i]+W[j])*v[k]**2/abs(v[k])**2 for k,(i,j) in enumerate(E)]
    return float(np.linalg.norm(r)/np.linalg.norm(v)**2)
def step_map(N,z,A):
    K=K_of(N,z,A); w,V=np.linalg.eigh(1j*K); return V@(np.exp(-1j*(2*math.pi/L)*w)*(V.conj().T@z))
def mono(N,v,A,h=1e-7):
    M=len(v); x=np.r_[v.real,v.imag]; D=np.zeros((2*M,2*M))
    for j in range(2*M):
        e=np.zeros(2*M); e[j]=h; zp=step_map(N,(x+e)[:M]+1j*(x+e)[M:],A); zm=step_map(N,(x-e)[:M]+1j*(x-e)[M:],A); D[:,j]=np.r_[(zp-zm).real,(zp-zm).imag]/(2*h)
    mu=selfconsistency(N,v,A)['mu']; phi=(2*math.pi/L)*mu; c,s=math.cos(phi),math.sin(phi)
    R=np.block([[c*np.eye(M),-s*np.eye(M)],[s*np.eye(M),c*np.eye(M)]]); ev=np.linalg.eigvals(R@D); return float(np.abs(ev).max())
# (a)
rows=[]
for N in range(4,12):
    th,cls,q=phases(N); A=adjacency(N)
    cases=[('equimodular',state(N,np.ones(q)/15))]
    if q>=4: cases.append(('member_k2_eps0.6',state(N,cos_mode(q,2,0.6)/15)))
    for name,v in cases:
        sc=selfconsistency(N,v,A); S,W=SW(N,v)
        rows.append(dict(N=N,case=name,selfcons_residual=f"{sc['residual']:.1e}",star_residual=f"{star_residual(N,v,sc['mu']):.1e}",max_abs_S=f"{abs(S).max():.1e}",W_spread=f"{np.ptp(W):.1e}",W_mean=f"{W.mean():.6f}",minus_mu=f"{-sc['mu']:.6f}"))
        print(f"(a) N={N} {name:18s} 自己無撞着残差={sc['residual']:.1e} (★)残差={rows[-1]['star_residual']} max|S_i|={abs(S).max():.1e} W の幅={np.ptp(W):.1e} W̄={W.mean():.6f} −μ={-sc['mu']:.6f}")
# 対照：自己無撞着でない状態で (★) が破れること
rng=np.random.default_rng(1); N=6; A=adjacency(6); v=rng.standard_normal(15)+1j*rng.standard_normal(15); sc=selfconsistency(6,v,A)
rows.append(dict(N=6,case='random_control',selfcons_residual=f"{sc['residual']:.2e}",star_residual=f"{star_residual(6,v,sc['mu']):.2e}",max_abs_S='',W_spread='',W_mean='',minus_mu=''))
print(f"(a) 対照 乱数状態 N=6: 自己無撞着残差={sc['residual']:.2e} (★)残差={rows[-1]['star_residual']}")
with open(os.path.join(ROOT,'results','identity_check.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
# (b)
def solve_balanced(N,rng,W0,iters=100):
    E=edges(N); M=len(E); v=rng.standard_normal(M)+1j*rng.standard_normal(M); v*=np.sqrt(N*W0/2)/np.linalg.norm(v)
    def G(x):
        vv=x[:M]+1j*x[M:]; S,W=SW(N,vv); return np.r_[S.real,S.imag,W-W0]
    for it in range(iters):
        x=np.r_[v.real,v.imag]; g=G(x)
        if np.linalg.norm(g)<1e-15: break
        J=np.zeros((3*N,2*M)); h=1e-7
        for j in range(2*M):
            e=np.zeros(2*M); e[j]=h; J[:,j]=(G(x+e)-G(x-e))/(2*h)
        dx=-np.linalg.lstsq(J,g,rcond=1e-12)[0]; x=x+dx; v=x[:M]+1j*x[M:]
    return v,float(np.linalg.norm(G(np.r_[v.real,v.imag])))
rows=[]
for N in [5,6,7,8]:
    A=adjacency(N); M=N*(N-1)//2; W0=2*NORM[N]**2/N   # Σ_i W_i = 2‖v‖²
    rng=np.random.default_rng(100+N)
    for s in range(5):
        v,g=solve_balanced(N,rng,W0)
        sc=selfconsistency(N,v,A); gt=gram_takagi(N,v); ax=gt['axes'][:gt['rank']]; rho=mono(N,v,A)
        S,W=SW(N,v); band='inflating' if rho-1>3e-3 else 'neutral'
        t=f'random_N{N}_s{s}'; dd=os.path.join(ROOT,'data',t); os.makedirs(dd,exist_ok=True)
        np.savez_compressed(os.path.join(dd,'parent_v.npz'),v=v,edges=np.array(edges(N)),color=np.zeros(M,int),theta=np.angle(v),design='random_balanced',r=float(np.sqrt(sc['mean_amp2'])),sigma=np.sort(np.linalg.eigvalsh(1j*K_of(N,v,A)))[::-1],mu=sc['mu'],residual=sc['residual'])
        rows.append(dict(tag=t,N=N,seed=s,constraint_residual=f'{g:.1e}',selfcons_residual=f"{sc['residual']:.1e}",mu=f"{sc['mu']:.6f}",minus_W0=f"{-W0:.6f}",global_closure=f"{sc['global_closure']:.1e}",local_closure=f"{sc['local_closure']:.1e}",
                         amp_min=f'{abs(v).min():.4f}',amp_max=f'{abs(v).max():.4f}',amp_std_over_rms=f'{abs(v).std()/np.sqrt((abs(v)**2).mean()):.3f}',phase_classes='none (generic)',rank=gt['rank'],roundness=f'{(ax.max()-ax.min())/ax.max():.3f}',axes=' '.join(f'{x:.4f}' for x in ax),vertex_null=f"{gt['vertex_null_dev']:.1e}",rho=f'{rho:.9f}',rho_minus_1=f'{rho-1:.2e}',lambda_f=f'{2*math.log(rho):.5f}',pred_band=band))
        print(f"(b) {t}: 拘束残差={g:.1e} 自己無撞着残差={sc['residual']:.1e} μ={sc['mu']:+.6f} (−W₀={-W0:+.6f}) 大域={sc['global_closure']:.1e} 局所={sc['local_closure']:.1e} |z|∈[{abs(v).min():.3f},{abs(v).max():.3f}] rank={gt['rank']} 丸さ={(ax.max()-ax.min())/ax.max():.3f} ρ−1={rho-1:.2e} λ_f={2*math.log(rho):.5f} → {band}")
with open(os.path.join(ROOT,'results','balanced_random_parents.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("PASS7 OK")
