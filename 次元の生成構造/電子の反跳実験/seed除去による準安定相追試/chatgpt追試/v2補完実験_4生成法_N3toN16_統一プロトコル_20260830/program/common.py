#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""共通モジュール：辺・隣接・位相クラス構成（偶数=1-因子分解／奇数=距離クラス）・K・自己無撞着検査・B/Takagi。
位相配置は 複素シンプレックス基礎_N別全展開_20260830/program/make_Ngeneric.py と同一。
本パッケージの新規点は振幅：辺 e の振幅二乗 a_e = r_e² をクラス別 a_c に置く（等モジュラーは a_c ≡ 一定の特殊例）。"""
import numpy as np

def edges(N): return [(i,j) for i in range(N) for j in range(i+1,N)]

def adjacency(N):
    E=edges(N); M=len(E); A=np.zeros((M,M))
    for a in range(M):
        for b in range(M):
            if a!=b and set(E[a])&set(E[b]): A[a,b]=1.0
    return A

def classes(N):
    """辺 → 位相クラス番号。偶数 N: 円法 1-因子分解（q=N−1 色、各頂点に各色 1 本）。
    奇数 N: 距離クラス（q=(N−1)/2、各頂点に各クラス 2 本）。N=3 は距離クラス 1 個（q=1）なので別扱い（Z3 配置）。"""
    E=edges(N); col={}
    if N%2==0:
        n=N-1
        for rr in range(n):
            col[tuple(sorted((rr,N-1)))]=rr
            for k in range(1,N//2): col[tuple(sorted(((rr-k)%n,(rr+k)%n)))]=rr
        q=n
    else:
        q=(N-1)//2
        for d in range(1,q+1):
            for i in range(N): col[tuple(sorted((i,(i+d)%N)))]=d-1
    return np.array([col[e] for e in E]), q

def phases(N):
    """クラス c の位相 θ_c = π c / q（d² の角度は 2θ_c = 2π c/q、つまり 1 の q 乗根）。N=3 は 0/60/120°。"""
    if N==3: return np.array([0.0,np.pi/3,2*np.pi/3]), np.array([0,1,2]), 3
    cls,q=classes(N); return np.pi*cls/q, cls, q

def state(N,a_class):
    """クラス別振幅二乗 a_class[c] から状態 v を作る。"""
    th,cls,q=phases(N); return np.sqrt(np.asarray(a_class,float)[cls])*np.exp(1j*th)

def K_of(N,v,A=None):
    if A is None: A=adjacency(N)
    return A*np.imag(np.conj(v)[:,None]*v[None,:])

def selfconsistency(N,v,A=None):
    """残差 ‖iKv−μv‖/‖v‖、μ、大域閉塞 |Σz²|/H、局所閉塞 max_i|Σ_{j} z_ij²|/H、H=‖v‖²。"""
    K=K_of(N,v,A); hv=1j*(K@v); H=float(np.vdot(v,v).real)
    mu=float((np.vdot(v,hv)/np.vdot(v,v)).real); res=float(np.linalg.norm(hv-mu*v)/np.sqrt(H))
    d2=v*v; E=edges(N)
    glob=float(abs(d2.sum())/H); loc=max(float(abs(sum(d2[k] for k,(a,b) in enumerate(E) if i in (a,b)))/H) for i in range(N))
    return dict(residual=res,mu=mu,H=H,mean_amp2=H/len(v),mu_over_mean_amp2=mu/(H/len(v)),global_closure=glob,local_closure=loc)

def gram_takagi(N,v):
    """D²→B=−½JD²J→特異値（Takagi 値）σ_k、軸スケール √σ_k、rank、埋め込み座標 x（(x_i−x_j)·(x_i−x_j)=d²_ij）。"""
    E=edges(N); d2=v*v; D=np.zeros((N,N),complex)
    for val,(i,j) in zip(d2,E): D[i,j]=D[j,i]=val
    J=np.eye(N)-np.ones((N,N))/N; B=-0.5*J@D@J
    sv=np.linalg.svd(B,compute_uv=False); rank=int((sv>1e-12*max(sv.max(),1e-300)).sum())
    Mre=np.block([[B.real,B.imag],[B.imag,-B.real]]); w,V=np.linalg.eigh(Mre)
    idx=np.argsort(w)[::-1][:N]; X=[]
    for j in idx:
        if w[j]<=1e-14: continue
        vv=V[:N,j]+1j*V[N:,j]; vv/=np.linalg.norm(vv); X.append(np.sqrt(w[j])*vv)
    X=np.array(X).T if X else np.zeros((N,0),complex)
    err=max(abs(((X[i]-X[j])@(X[i]-X[j]))-D[i,j]) for i,j in E) if X.size else float('nan')
    nulldev=max(abs(X[i]@X[i]) for i in range(N)) if X.size else float('nan')
    return dict(sigma=sv,axes=np.sqrt(np.maximum(sv,0)),rank=rank,X=X,embed_err=float(err),vertex_null_dev=float(nulldev),B=B)

def is_round(axes,rank,tol=1e-8):
    ax=axes[:rank]; return bool(len(ax)>0 and (ax.max()-ax.min())<tol*ax.max())

def family_dim(q):
    """クラス別振幅の自己無撞着条件 Σ_c a_c ω^c = 0（ω=e^{2πi/q}）の実自由度（全体スケールを除く）。
    q=1: 条件は a_0=0 → 非自明解なし（N=3 は別配置）。q=2: ω=−1 で a_0=a_1 → 0。q≥3: 実条件 2 本 → q−3。"""
    if q<=1: return 0
    if q==2: return 0
    return q-3

def cos_mode(q,k,eps):
    """a_c = 1 + ε cos(2πkc/q)。k ≢ 0,±1 (mod q) なら Σ a_c ω^c = 0 を厳密に満たす（三角和の直交性）。"""
    return 1.0+eps*np.cos(2*np.pi*k*np.arange(q)/q)
