#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス1：手作り親の生成・受け入れ検査・保存（N=4..8、乱数ゼロ・決定論）。
設計 A（偶数 N）：円法 1-因子分解、位相 θ_e = color·π/(N−1)。
設計 C（奇数 N）：距離巡回ハミルトン閉路分解、位相 θ_e = (d−1)·π/((N−1)/2)。
振幅：全辺 r = ‖v‖/√M、‖v‖ は同 N の fixed_equimodular 親の実測値に一致させる。
検査（1 つでも不合格なら abort、実験を走らせない）：
  残差 <1e-12 / 閉塞 |Σz²|/H <1e-13 / |z|² 相対幅 <1e-12 / μ/r²=−(N−1) (1e-9) / (σ/r²)² が閉形式 (1e-6)"""
import os, csv, json, sys
import numpy as np
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NORM={4:0.816497,5:0.927173,6:1.036862,7:1.210458,8:1.301688}  # fixed_equimodular 親の ‖v‖
def edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def adjacency(N):
    E=edges(N); M=len(E); A=np.zeros((M,M))
    for a in range(M):
        for b in range(M):
            if a!=b and set(E[a])&set(E[b]): A[a,b]=1.0
    return A
def one_factorization(N):
    n=N-1; col={}
    for r in range(n):
        col[tuple(sorted((r,N-1)))]=r
        for k in range(1,N//2): col[tuple(sorted(((r-k)%n,(r+k)%n)))]=r
    return col
def hamilton_decomposition(N):
    n=(N-1)//2
    assert all(np.gcd(d,N)==1 for d in range(1,n+1)), f"N={N}: distance circulant not Hamiltonian"
    col={}
    for d in range(1,n+1):
        for i in range(N): col[tuple(sorted((i,(i+d)%N)))]=d-1
    return col
def make(N):
    E=edges(N); M=len(E)
    if N%2==0: col=one_factorization(N); th=np.array([col[e]*np.pi/(N-1) for e in E]); typ='A_matching'
    else: col=hamilton_decomposition(N); th=np.array([col[e]*np.pi/((N-1)//2) for e in E]); typ='C_cycle'
    r=NORM[N]/np.sqrt(M)
    return E,np.array([col[e] for e in E]),th,r*np.exp(1j*th),typ,r
def checks(N,v):
    A=adjacency(N); K=A*np.imag(np.conj(v)[:,None]*v[None,:]); hv=1j*(K@v)
    mu=float((np.vdot(v,hv)/np.vdot(v,v)).real); res=float(np.linalg.norm(hv-mu*v)/np.linalg.norm(v))
    r2=float((np.abs(v)**2).mean()); H=float((np.abs(v)**2).sum())
    clo=float(abs((v*v).sum())/H); amp2=np.abs(v)**2; spread=float((amp2.max()-amp2.min())/amp2.mean())
    ev=np.linalg.eigvalsh(1j*K); s=np.sort(ev[ev>1e-9])[::-1]/r2
    if N%2==0: pred=sorted({float((N-1)**2),(N-1)*(N-3)/4.0,(N-1)*(N-4)/4.0}) if N>4 else [9.0]
    else: pred=sorted({float((N-1)**2),(N-1)*(N-1)/8.0}) if N==7 else sorted({16.0,1.0})
    got=sorted(set(np.round(s*s,6)))
    ok_sig=len(got)==len(pred) and max(abs(g-p) for g,p in zip(got,pred))<1e-6
    rep=dict(residual=res,mu=mu,mu_over_r2=mu/r2,closure=clo,spread=spread,r2=r2,
             sigma2_over_r2=got,sigma2_pred=pred,
             ok=dict(residual=res<1e-12,closure=clo<1e-13,spread=spread<1e-12,
                     mu=abs(mu/r2+(N-1))<1e-9,sigma2=ok_sig))
    return rep,s*r2
for N in [4,5,6,7,8]:
    E,c,th,v,typ,r=make(N)
    rep,sig=checks(N,v)
    print(f"N={N} {typ} M={len(E)} r={r:.9f} 残差={rep['residual']:.1e} μ/r²={rep['mu_over_r2']:+.6f} 閉塞={rep['closure']:.1e} (σ/r²)²={rep['sigma2_over_r2']} 判定={rep['ok']}")
    if not all(rep['ok'].values()):
        raise SystemExit(f"ABORT: N={N} acceptance failed: {rep['ok']}")
    dd=os.path.join(ROOT,'data',f'N{N}'); os.makedirs(dd,exist_ok=True)
    np.savez_compressed(os.path.join(dd,'parent_v.npz'),v=v,edges=np.array(E),color=c,theta=th,design=typ,r=r,sigma=sig,mu=rep['mu'],residual=rep['residual'])
    with open(os.path.join(dd,'parent_v.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['edge_index','i','j','color','theta_deg','a_Re','b_Im','abs_z'])
        for k,(e,cc) in enumerate(zip(E,c)): w.writerow([k,e[0],e[1],int(cc),np.degrees(th[k]),v[k].real,v[k].imag,abs(v[k])])
    with open(os.path.join(dd,'parent_checks.json'),'w') as f: json.dump(rep,f,indent=2,default=lambda o: bool(o) if isinstance(o,np.bool_) else float(o))
print("PASS1 OK")
