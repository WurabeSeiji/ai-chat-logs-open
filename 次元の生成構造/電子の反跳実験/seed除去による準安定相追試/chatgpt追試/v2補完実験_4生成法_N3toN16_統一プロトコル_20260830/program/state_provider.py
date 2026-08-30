#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""非等モジュラー版の各 N の代表状態（決定論・乱数は固定 seed）。
規約：全 N で平均振幅二乗 〈|z|²〉 = r̄² = 1/15（等モジュラー版の r² = 1/15 と同じ ‖v‖²、同じ μ）。全体位相：クラス 0 の最初の組を θ=0°。
(A) クラス重み付き族（q ≥ 4、N=6,8,…,16）：位相配置は等モジュラー版と同一、クラス c の振幅二乗 a_c = r̄²(1 + ε cos(4πc/q))、ε=0.6。
    条件 Σ_c a_c ω^c = 0（ω=e^{2πi/q}）を三角和の直交性で厳密に満たす（k=2 モード、q≥4 で k=2 ≢ 0,±1）。
(B) 多様体上の代表点（q ≤ 3、N=3,4,5,7：クラス重み付き族が点に退化する N）：等モジュラー点から自己無撞着写像
    F(v)=iK(v)v−μ(v)v のヤコビアン核（スケール・全体位相を除く）の固定方向へ δ 動かし、Gauss–Newton で F=0 に戻した状態。
    N=5,7 は局所閉塞 S_i=0 を保つ分枝（{S=0, W=const}）を Newton の拘束に含めて選ぶ。N=3,4 は局所閉塞を保つ非等モジュラー解が存在しないので F=0 のみ。"""
import numpy as np
EPS=0.6; K_MODE=2; RBAR2=1.0/15.0
def edges(N): return [(i,j) for i in range(N) for j in range(i+1,N)]
def adjacency(N):
    E=edges(N); M=len(E); A=np.zeros((M,M))
    for a in range(M):
        for b in range(M):
            if a!=b and set(E[a])&set(E[b]): A[a,b]=1.0
    return A
def classes(N):
    """等モジュラー版（make_Ngeneric.py / make_N3.py / make_N4.py）と同一の位相クラスと位相刻み。"""
    E=edges(N)
    if N==3: return {(0,2):0,(1,2):1,(0,1):2}, 3, 60.0
    if N==4: return {(0,3):0,(1,2):0,(0,2):1,(1,3):1,(0,1):2,(2,3):2}, 3, 60.0
    if N%2==0:
        n=N-1; col={}
        for rr in range(n):
            col[tuple(sorted((rr,N-1)))]=rr
            for k in range(1,N//2): col[tuple(sorted(((rr-k)%n,(rr+k)%n)))]=rr
        return col, n, 180.0/(N-1)
    nn=(N-1)//2; col={}
    for d in range(1,nn+1):
        for i in range(N): col[tuple(sorted((i,(i+d)%N)))]=d-1
    return col, nn, 180.0/nn
def equimodular(N):
    E=edges(N); col,q,step=classes(N)
    th=np.array([np.radians(step*col[e]) for e in E]); return np.sqrt(RBAR2)*np.exp(1j*th)
def K_of(N,v,A): return A*np.imag(np.conj(v)[:,None]*v[None,:])
def mu_res(N,v,A):
    K=K_of(N,v,A); hv=1j*(K@v); mu=(np.vdot(v,hv)/np.vdot(v,v)).real; return float(mu), float(np.linalg.norm(hv-mu*v)/np.linalg.norm(v))
def SW(N,v):
    E=edges(N); S=np.zeros(N,complex); W=np.zeros(N)
    for k,(i,j) in enumerate(E): S[i]+=v[k]**2; S[j]+=v[k]**2; W[i]+=abs(v[k])**2; W[j]+=abs(v[k])**2
    return S,W
def _F(N,x,A,keep_local):
    M=len(x)//2; v=x[:M]+1j*x[M:]; K=K_of(N,v,A); hv=1j*(K@v); mu=(np.vdot(v,hv)/np.vdot(v,v)).real; r=hv-mu*v
    out=np.r_[r.real,r.imag]
    if keep_local:
        S,W=SW(N,v); out=np.r_[out,S.real,S.imag]
    return out
def _jac(N,x,A,keep_local,h=1e-6):
    f0=_F(N,x,A,keep_local); J=np.zeros((len(f0),len(x)))
    for j in range(len(x)):
        e=np.zeros(len(x)); e[j]=h; J[:,j]=(_F(N,x+e,A,keep_local)-_F(N,x-e,A,keep_local))/(2*h)
    return J
def balanced_point(N,delta=0.3,seed=0):
    """N=5,7 用：等モジュラー点から核方向へ δ‖v0‖ 動かし、{S_i=0, W_i=W0}（局所閉塞＋頂点重み均等 ⇒ 自己無撞着、μ=−W0）へ Newton。"""
    E=edges(N); M=len(E); A=adjacency(N); v0=equimodular(N); x0=np.r_[v0.real,v0.imag]; W0=2*np.linalg.norm(v0)**2/N
    J=_jac(N,x0,A,False); U,s,Vt=np.linalg.svd(J); null=int((s<s.max()*1e-7).sum()); Nsp=Vt[-null:].T
    triv=np.c_[x0,np.r_[-v0.imag,v0.real]]; triv,_=np.linalg.qr(triv); P=Nsp-triv@(triv.T@Nsp); Q,_=np.linalg.qr(P); Q=Q[:,:null-2]
    rng=np.random.default_rng(seed); d=Q@rng.standard_normal(Q.shape[1]); d/=np.linalg.norm(d)
    x=x0+delta*np.linalg.norm(x0)*d
    def G(x):
        vv=x[:M]+1j*x[M:]; S,W=SW(N,vv); return np.r_[S.real,S.imag,W-W0]
    for it in range(100):
        g=G(x)
        if np.linalg.norm(g)<1e-15: break
        Jm=np.zeros((3*N,2*M)); h=1e-7
        for j in range(2*M):
            e=np.zeros(2*M); e[j]=h; Jm[:,j]=(G(x+e)-G(x-e))/(2*h)
        x=x-np.linalg.lstsq(Jm,g,rcond=1e-12)[0]
    return x[:M]+1j*x[M:]
def manifold_point(N,delta=0.3,seed=0,keep_local=False):
    """N=3,4 用：等モジュラー点 v0 から核方向へ δ‖v0‖ 動かし F=0 へ Newton。‖v‖ と全体位相を固定。"""
    E=edges(N); M=len(E); A=adjacency(N); v0=equimodular(N); x0=np.r_[v0.real,v0.imag]
    J=_jac(N,x0,A,False); U,s,Vt=np.linalg.svd(J); null=int((s<s.max()*1e-7).sum()); Nsp=Vt[-null:].T
    triv=np.c_[x0,np.r_[-v0.imag,v0.real]]; triv,_=np.linalg.qr(triv); P=Nsp-triv@(triv.T@Nsp); Q,_=np.linalg.qr(P); Q=Q[:,:null-2]
    rng=np.random.default_rng(seed); d=Q@rng.standard_normal(Q.shape[1]); d/=np.linalg.norm(d)
    x=x0+delta*np.linalg.norm(x0)*d; n0=np.linalg.norm(x0)
    for it in range(200):
        r=_F(N,x,A,keep_local)
        if np.linalg.norm(r)<1e-15: break
        Jm=_jac(N,x,A,keep_local); v=x[:M]+1j*x[M:]
        g1=x/np.linalg.norm(x); g2=np.r_[-v0.imag,v0.real]
        Ja=np.vstack([Jm,g1,g2]); ra=np.r_[r,np.linalg.norm(x)-n0,float(np.vdot(v0,v).imag)]
        x=x-np.linalg.lstsq(Ja,ra,rcond=1e-15)[0]
    # 仕上げ：scipy の Levenberg–Marquardt で F=0（＋‖v‖・全体位相の拘束）を機械精度まで磨く
    from scipy.optimize import least_squares
    def resid(x):
        v=x[:M]+1j*x[M:]; return np.r_[_F(N,x,A,False),np.linalg.norm(x)-n0,float(np.vdot(v0,v).imag)]
    x=least_squares(resid,x,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=20000).x
    v=x[:M]+1j*x[M:]
    return v
def state(N):
    """(v, kind, col, q, step)。kind='class' はクラス重み付き族、'manifold' は多様体上の代表点。"""
    E=edges(N); col,q,step=classes(N)
    if q>=4:
        cls=np.array([col[e] for e in E]); a=RBAR2*(1+EPS*np.cos(2*np.pi*K_MODE*cls/q))
        # 各クラス同数なので平均は r̄² のまま。全体位相：クラス 0 は θ=0
        v=np.sqrt(a)*np.exp(1j*np.radians(step*cls)); kind='class'
    else:
        v=balanced_point(N) if N in (5,7) else manifold_point(N); kind='manifold'
        # 全体位相：クラス 0 の最初の組を θ=0 に
        k0=[k for k,e in enumerate(E) if col[e]==0][0]; v=v*np.exp(-1j*np.angle(v[k0]))
    v=v*np.sqrt(len(E)*RBAR2)/np.linalg.norm(v)
    return v,kind,col,q,step
if __name__=='__main__':
    for N in [3,4,5,7]:
        v,kind,col,q,step=state(N); A=adjacency(N); mu,res=mu_res(N,v,A); S,W=SW(N,v); E=edges(N)
        print(f"N={N} {kind}: μ={mu:+.9f} (等モジュラー −(N−1)/15={-(N-1)/15:+.9f}) 残差={res:.1e} |z|∈[{abs(v).min():.4f},{abs(v).max():.4f}] max|S_i|/r̄²={abs(S).max()*15:.2e} W 幅/r̄²={np.ptp(W)*15:.2e} |Σz²|={abs((v*v).sum()):.1e}")
        for k,e in enumerate(E): print(f"   組{e} クラス{col[e]}: |z|={abs(v[k]):.6f} θ={np.degrees(np.angle(v[k]))%360:7.3f}°  (等モジュラー θ={step*col[e]:.0f}°)")
