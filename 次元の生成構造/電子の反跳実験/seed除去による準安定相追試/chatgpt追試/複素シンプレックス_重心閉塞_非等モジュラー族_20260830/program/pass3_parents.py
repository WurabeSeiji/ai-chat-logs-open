#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス3：力学走行用の親（等モジュラー対照＋非等モジュラー族メンバー）の生成・受け入れ検査・保存と、
共回転モノドロミーによる成長率の事前予測（走行前に results/parents_predictions.csv に固定＝事前登録）。
親：N6: ε=0/0.3/0.6/0.9 (k=2)、N8: ε=0 / 0.6(k=2) / 0.6(k=3)、N9: ε=0 / 0.6(k=2)。
‖v‖ は N=4..8 は手作り親パッケージ（fixed_equimodular 親の実測値）に一致、N=9 は 1.40（系列に沿う外挿。スケールは規約で速さにしか効かない）。
受け入れ検査：残差<1e-12、大域閉塞<1e-13、局所閉塞<1e-13、μ/⟨|z|²⟩=−(N−1) (1e-9)、‖v‖ 一致。
予測：1 step 写像 Φ(z)=exp(−i(2π/L)·iK(z)) z の実ヤコビアン DΦ を中心差分で求め、共回転 G=R(+φ)DΦ（φ=(2π/L)μ）の
最大固有値絶対値 ρ から λ_f=2 ln ρ（H⊥ の per-step 成長率）。ρ≤1+1e-9 なら中立（床）と予測。
対照テスト：ε=0 の N6/N8 で先行パッケージの予測 0.01466/0.01861（実測 0.01431/0.01825）と一致することを確認。"""
import os, csv, json, sys, math
import numpy as np
sys.path.insert(0,os.path.dirname(__file__))
from common import *
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NORM={4:0.816497,5:0.927173,6:1.036862,7:1.210458,8:1.301688,9:1.40}
L=124
PARENTS=[(6,0.0,2),(6,0.3,2),(6,0.6,2),(6,0.9,2),(8,0.0,2),(8,0.6,2),(8,0.6,3),(9,0.0,2),(9,0.6,2)]
def tag(N,eps,k): return f"N{N}_eps{eps:.2f}_k{k}"
def step_map(N,z,A):
    K=K_of(N,z,A); H=1j*K; w,V=np.linalg.eigh(H); return V@(np.exp(-1j*(2*math.pi/L)*w)*(V.conj().T@z))
def monodromy_growth(N,v,A,h=1e-7):
    M=len(v); x=np.r_[v.real,v.imag]; D=np.zeros((2*M,2*M))
    for j in range(2*M):
        e=np.zeros(2*M); e[j]=h
        zp=step_map(N,(x+e)[:M]+1j*(x+e)[M:],A); zm=step_map(N,(x-e)[:M]+1j*(x-e)[M:],A)
        D[:,j]=np.r_[(zp-zm).real,(zp-zm).imag]/(2*h)
    mu=selfconsistency(N,v,A)['mu']; phi=(2*math.pi/L)*mu   # 親は z→e^{−iφ}z
    c,s=math.cos(phi),math.sin(phi); R=np.block([[c*np.eye(M),-s*np.eye(M)],[s*np.eye(M),c*np.eye(M)]])  # e^{+iφ} を掛けて共回転へ
    G=R@D; ev=np.linalg.eigvals(G); rho=float(np.abs(ev).max())
    return rho,2*math.log(rho),int((np.abs(ev)>1+1e-9).sum())
rows=[]
for N,eps,k in PARENTS:
    th,cls,q=phases(N); A=adjacency(N); M=len(cls)
    a=cos_mode(q,k,eps); v=state(N,a); v*=NORM[N]/np.linalg.norm(v)
    sc=selfconsistency(N,v,A); gt=gram_takagi(N,v)
    ok=dict(residual=sc['residual']<1e-12,global_closure=sc['global_closure']<1e-13,local_closure=sc['local_closure']<1e-13,
            mu=abs(sc['mu_over_mean_amp2']+(N-1))<1e-9,norm=abs(np.linalg.norm(v)-NORM[N])<1e-9)
    rho,lam,nunst=monodromy_growth(N,v,A)
    kind='inflating' if rho>1+1e-9 else 'neutral'
    t=tag(N,eps,k); dd=os.path.join(ROOT,'data',t); os.makedirs(dd,exist_ok=True)
    print(f"{t}: M={M} q={q} |z|∈[{abs(v).min():.6f},{abs(v).max():.6f}] 残差={sc['residual']:.1e} μ={sc['mu']:+.6f} μ/⟨|z|²⟩={sc['mu_over_mean_amp2']:+.6f} 大域={sc['global_closure']:.1e} 局所={sc['local_closure']:.1e} rank={gt['rank']} 軸={np.round(gt['axes'][:gt['rank']],5)} | 予測 ρ={rho:.8f} λ_f={lam:.5f} 不安定固有値数={nunst} → {kind} 判定={ok}")
    if not all(ok.values()): raise SystemExit(f"ABORT: {t} acceptance failed {ok}")
    E=edges(N)
    np.savez_compressed(os.path.join(dd,'parent_v.npz'),v=v,edges=np.array(E),color=cls,theta=th,design=f'{"A_matching" if N%2==0 else "C_distance"}_eps{eps}_k{k}',r=float(np.sqrt(sc['mean_amp2'])),sigma=np.sort(np.linalg.eigvalsh(1j*K_of(N,v,A)))[::-1],mu=sc['mu'],residual=sc['residual'],eps=eps,k=k,a_class=a*(NORM[N]**2/ (a[cls].sum())))
    with open(os.path.join(dd,'parent_v.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['edge_index','i','j','class','theta_deg','a_Re','b_Im','abs_z','abs_z_sq'])
        for kk,(e,cc) in enumerate(zip(E,cls)): w.writerow([kk,e[0],e[1],int(cc),np.degrees(th[kk]),v[kk].real,v[kk].imag,abs(v[kk]),abs(v[kk])**2])
    rep=dict(N=N,eps=eps,k=k,M=M,q=q,norm=float(np.linalg.norm(v)),amp_min=float(abs(v).min()),amp_max=float(abs(v).max()),**{kk:(float(vv) if not isinstance(vv,dict) else vv) for kk,vv in sc.items()},rank=gt['rank'],axes=[float(x) for x in gt['axes'][:gt['rank']]],vertex_null=gt['vertex_null_dev'],ok={kk:bool(vv) for kk,vv in ok.items()},pred_rho=rho,pred_lambda_f=lam,pred_n_unstable=nunst,pred_kind=kind)
    with open(os.path.join(dd,'parent_checks.json'),'w') as f: json.dump(rep,f,indent=2,ensure_ascii=False)
    rows.append(dict(tag=t,N=N,eps=eps,k=k,M=M,q=q,norm=rep['norm'],amp_min=rep['amp_min'],amp_max=rep['amp_max'],residual=sc['residual'],mu=sc['mu'],global_closure=sc['global_closure'],local_closure=sc['local_closure'],rank=gt['rank'],
                     axes=' '.join(f'{x:.6f}' for x in gt['axes'][:gt['rank']]),pred_rho=f'{rho:.9f}',pred_lambda_f=f'{lam:.6f}',pred_n_unstable=nunst,pred_kind=kind,
                     pred_t50=f'{(math.log(0.5)+72.6)/lam:.0f}' if lam>1e-9 else ''))
with open(os.path.join(ROOT,'results','parents_predictions.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("PASS3 OK（予測は走行前に固定）")
