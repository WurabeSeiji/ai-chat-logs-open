#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス3：図化・τ 判定・分析（読出し専用、力学に無関係。何度でも再実行可）。
- 従来図（H⊥/H・PR/M・amp_std・|ZᵀZ|、局面マーカー付き）
- 平面占有スペクトル（iK(v) の ±σ 固有平面ごとの |射影|²、基底は step 0 の親で固定）
- 5 局面（①初期 ②急拡大の中央 ③緩和の膝の中央 ④再直線化の開始 ⑤最終）で複素平面図＋幾何構造図。
  中立走行（max f < 1e-6）は ①・τ=20000・⑤ の 3 点に縮退。
- 共回転ゲージ：⟨v,Z⟩ が正実になる全体位相を掛けてから描く。
- results/stage_taus.csv, results/analysis_summary.csv を出力。"""
import os, csv, json, sys, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CMAP=plt.get_cmap('twilight')
def edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def adjacency(N):
    E=edges(N); M=len(E); A=np.zeros((M,M))
    for a in range(M):
        for b in range(M):
            if a!=b and set(E[a])&set(E[b]): A[a,b]=1.0
    return A
def gauge(v,Z):
    ph=np.vdot(v,Z)
    return Z*np.exp(-1j*np.angle(ph)) if abs(ph)>0 else Z
def stage_taus(f,steps):
    """局面 τ の判定。戻り値 [(label, tau), ...]"""
    if np.nanmax(f)<1e-6:
        return [('S1_initial',0),('S3_mid',steps//2),('S5_final',steps)],'neutral'
    lf=np.log(np.maximum(f,1e-300)); t=np.arange(len(f))
    on=np.where(f>1e-8)[0]; t_on=int(on[0]) if len(on) else 0
    plateau=float(np.nanmedian(lf[-max(len(f)//8,100):]))
    mid_lf=(lf[t_on]+plateau)/2.0
    seg=np.where((t>=t_on)&(lf>=mid_lf))[0]; t2=int(seg[0]) if len(seg) else t_on
    # 平滑化した傾き
    w=51; k=np.ones(w)/w; sl=np.gradient(np.convolve(lf,k,mode='same'))
    lam=np.nanmax(sl[t_on:min(t_on+8000,len(f)-1)]) if t_on<len(f)-1 else np.nanmax(sl)
    after=np.arange(t2,len(f))
    a50=after[sl[after]<0.5*lam]; tA=int(a50[0]) if len(a50) else t2
    a10=after[(after>=tA)&(sl[after]<0.1*lam)]; tB=int(a10[0]) if len(a10) else tA
    t3=(tA+tB)//2
    tail=float(np.nanmedian(sl[-max(len(f)//8,100):]))
    cand=np.arange(tB,len(f)); good=cand[np.abs(sl[cand]-tail)<=max(abs(tail)*0.5,1e-7)]
    t4=int(good[0]) if len(good) else tB
    return [('S1_initial',0),('S2_growth_mid',t2),('S3_knee_mid',t3),('S4_relinear',t4),('S5_final',steps)],'inflating'
def complex_plane_fig(v,Z,r0,path,title,lim=None):
    """log-polar（平方根ゲージ）複素平面図・2 色重ね：
      青 = 状態 z_e、赤 = 親平面と直交な面に乗る成分 (Z⊥)_e（写像）
    半径 ρ = 1 − sqrt(u/UMAX)、u = −log10(|z|/Rtop)、UMAX=34（倍精度床 1e-32 まで窓内。
    対数軸より深部を強く引き伸ばすゲージ：1e-31→1e-23 の 8 桁移動が半径の ~13% で見える）。
    角度=位相（共回転ゲージ、青赤同一）。橙破線=親振幅 r0。厳密ゼロのみ中心マーカー。"""
    UMAX=34.0
    Rtop=float(np.linalg.norm(v))
    ph_g=np.vdot(v,Z); gph=np.exp(-1j*np.angle(ph_g)) if abs(ph_g)>0 else 1.0
    pv=v.real/np.linalg.norm(v.real); qv=v.imag-(v.imag@pv)*pv; qv/=np.linalg.norm(qv)
    Zp=Z-pv*(pv@Z)-qv*(qv@Z)
    def rad(amp):
        with np.errstate(divide='ignore'):
            u=np.clip(-(np.log10(np.maximum(amp,1e-300))-np.log10(Rtop)),0.0,UMAX)
        return 1.0-np.sqrt(u/UMAX)
    fig,ax=plt.subplots(figsize=(6.8,6.8))
    th=np.linspace(0,2*np.pi,241)
    for k in [2,4,8,16,24,32]:
        rr=1.0-np.sqrt(k/UMAX)
        ax.plot(rr*np.cos(th),rr*np.sin(th),color='gray',lw=0.5,alpha=0.35)
        ax.annotate(f'1e-{k}',(rr*0.707+0.02,rr*0.707+0.02),fontsize=6,color='gray',alpha=0.85)
    rr0=rad(np.array([r0]))[0]
    ax.plot(rr0*np.cos(th),rr0*np.sin(th),'--',color='tab:orange',lw=1.0,alpha=0.85)
    def put(W,color,marker,label,dy):
        Wg=W*gph; amp=np.abs(Wg); ang=np.angle(Wg); rho=rad(amp)
        zero=amp<=0
        ax.scatter(rho[~zero]*np.cos(ang[~zero]),rho[~zero]*np.sin(ang[~zero]),s=55,alpha=0.4,color=color,edgecolors='none',label=label)
        if zero.any():
            ax.scatter([0],[0],marker=marker,s=80,color=color,alpha=0.8)
            ax.annotate(f'=0 ×{int(zero.sum())}',(0.03,dy),fontsize=7,color=color)
    put(Z,'tab:blue','v','state z_e',-0.07)
    put(Zp,'tab:red','^','transverse (Z_perp)_e',-0.14)
    ax.axhline(0,color='k',lw=0.3,alpha=0.3); ax.axvline(0,color='k',lw=0.3,alpha=0.3)
    ax.set_xlim(-1.25,1.25); ax.set_ylim(-1.25,1.25); ax.set_aspect('equal'); ax.set_axis_off()
    ax.legend(loc='lower left',fontsize=7,framealpha=0.6)
    ax.set_title(title+'   [sqrt-log gauge: rho=1-sqrt(u/34), u=-log10(|z|/Rtop); blue=state, red=transverse]',fontsize=7.5)
    fig.tight_layout(); fig.savefig(path,dpi=160); plt.close(fig)
def layout(N):
    if N==4: return None  # 3D tetrahedron
    if N%2==1:
        ang=2*np.pi*np.arange(N)/N
        return np.c_[np.cos(ang),np.sin(ang)]
    ang=2*np.pi*np.arange(N-1)/(N-1)
    P=np.zeros((N,2)); P[:N-1]=np.c_[np.cos(ang),np.sin(ang)]
    return P
def centered_gram(z,N):
    """シリーズ確立済み定義：d_ij² = z_ij²、B = −1/2 J D² J（J = I − 11ᵀ/N）"""
    E=edges(N); D2=np.zeros((N,N),complex)
    for val,(i,j) in zip(z*z,E): D2[i,j]=D2[j,i]=val
    J=np.eye(N)-np.ones((N,N))/N
    return -0.5*J@D2@J
def takagi(B):
    """Autonne–Takagi 分解 B = Σ σ_k v_k v_kᵀ（σ_k ≥ 0 降順、v_k 正規直交）。
    実対称 2N×2N 埋め込み M=[[ReB,ImB],[ImB,−ReB]] の正固有対 (σ,(x;y)) → v=x+iy が B conj(v)=σv を満たす。
    位相自由度は ±1 のみ（先頭最大成分の実部 ≥0 に固定）→ 座標は本質的に一意＝事実の幾何。"""
    Nn=B.shape[0]; P=B.real; Q=B.imag
    M=np.block([[P,Q],[Q,-P]])
    w,V=np.linalg.eigh(M)
    idx=np.argsort(w)[::-1][:Nn]
    sig=[]; vecs=[]
    for j in idx:
        if w[j]<=1e-14: continue
        v=V[:Nn,j]+1j*V[Nn:,j]; v=v/np.linalg.norm(v)
        k=int(np.argmax(np.abs(v)))
        if v[k].real<0 or (v[k].real==0 and v[k].imag<0): v=-v
        sig.append(float(w[j])); vecs.append(v)
    return np.array(sig),(np.array(vecs).T if vecs else np.zeros((Nn,0)))
def simplex_dim(scales):
    """単体の次元読出し：rank（非零 Takagi 軸数）と実効次元 d_eff=(Σσ)²/Σσ²（σ=Takagi 値=r_k²）"""
    sg=scales**2; sg=sg[sg>1e-12*max(1e-300,sg.max())]
    return (len(sg), float((sg.sum())**2/(sg*sg).sum())) if len(sg) else (0,0.0)
def takagi_scales(z,N):
    return np.sqrt(np.maximum(np.linalg.svd(centered_gram(z,N),compute_uv=False),0.0))
def takagi_coords(z,N,kmax=3):
    sig,Vk=takagi(centered_gram(z,N))
    k=min(kmax,Vk.shape[1]); X=Vk[:,:k]*np.sqrt(sig[:k])[None,:]
    return sig,X   # X[i,k] = 頂点 i の第 k 正準軸の複素座標、(x_i−x_j)·(x_i−x_j)=z_ij²
def geometry_fig(N,v,Z,r0,path,title,axlim=None):
    """事実の幾何：複素シンプレックスの上位 3 正準軸平面に頂点を実座標で描画＋軸スケールスペクトル。
    辺は状態の位相で彩色（細線）。axlim は N ごとに全局面共通（絶対スケール）。"""
    E=edges(N); ph=(np.angle(gauge(v,Z)))%(2*np.pi)
    sig,X=takagi_coords(Z,N,kmax=3)
    pv=v.real/np.linalg.norm(v.real); qv=v.imag-(v.imag@pv)*pv; qv/=np.linalg.norm(qv)
    Zp=Z-pv*(pv@Z)-qv*(qv@Z)
    rf=takagi_scales(Z,N); rp=takagi_scales(Zp,N)
    kshow=X.shape[1]
    fig,axs=plt.subplots(1,kshow+1,figsize=(4.3*(kshow+1),4.6))
    if axlim is None: axlim=1.1*max(1e-12,float(np.abs(X).max()))
    for k in range(kshow):
        ax=axs[k]
        segs=[[(X[i,k].real,X[i,k].imag),(X[j,k].real,X[j,k].imag)] for (i,j) in E]
        ax.add_collection(LineCollection(segs,colors=[CMAP(x/(2*np.pi)) for x in ph],linewidths=0.9,alpha=0.55))
        ax.scatter(X[:,k].real,X[:,k].imag,color='k',s=28,zorder=3)
        for i in range(N): ax.annotate(str(i),(X[i,k].real,X[i,k].imag),textcoords='offset points',xytext=(5,5),fontsize=8)
        ax.set_xlim(-axlim,axlim); ax.set_ylim(-axlim,axlim); ax.set_aspect('equal')
        ax.axhline(0,color='k',lw=0.3,alpha=0.3); ax.axvline(0,color='k',lw=0.3,alpha=0.3); ax.grid(True,alpha=0.15)
        ax.set_title(f'canonical axis {k+1} (complex plane)\nr_{k+1}=sqrt(sigma)={np.sqrt(sig[k]):.3e}',fontsize=8)
    ax=axs[kshow]; kk=np.arange(1,len(rf)+1)
    ax.semilogy(kk,np.maximum(rf,1e-40),'o-',color='tab:blue',label='r_k of FULL state Z')
    ax.semilogy(np.arange(1,len(rp)+1),np.maximum(rp,1e-40),'^--',color='tab:red',label='r_k of TRANSVERSE Z_perp')
    ax.set_xlabel('axis k'); ax.set_ylabel('r_k'); ax.set_ylim(1e-36,10); ax.grid(True,which='both',alpha=.25); ax.legend(fontsize=7)
    ax.set_title('axis-scale spectrum\n(blue/red convention = complex-plane fig)',fontsize=8)
    rk_,de_=simplex_dim(rf); rkp_,dep_=simplex_dim(rp)
    fig.suptitle(title+f'   [Takagi coords, (x_i-x_j)·(x_i-x_j)=z_ij²]   dim: full rank={rk_} d_eff={de_:.2f} | perp rank={rkp_} d_eff={dep_:.2f} (max N-1={N-1})',fontsize=9)
    fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(path,dpi=160); plt.close(fig)
def plane_basis(N,v):
    A=adjacency(N); K=A*np.imag(np.conj(v)[:,None]*v[None,:])
    w,U=np.linalg.eigh(1j*K)
    pos=np.where(w>1e-12)[0]
    planes=[]
    for j in pos:  # 平面 = (+σ, −σ) の対。−σ 側は複素共役固有ベクトル
        planes.append((float(w[j]),U[:,j]))
    kern=np.where(np.abs(w)<=1e-12)[0]
    return planes,U[:,kern]
st_rows=[]; an_rows=[]
for N in [4,5,6,7,8]:
    dd=os.path.join(ROOT,'data',f'N{N}'); fd=os.path.join(ROOT,'figures',f'N{N}'); os.makedirs(fd,exist_ok=True)
    pz=np.load(os.path.join(dd,'parent_v.npz')); v=pz['v']; r0=float(pz['r'])
    S=np.load(os.path.join(dd,'states_treatment.npz'))['Z']
    ts=np.genfromtxt(os.path.join(dd,'treatment_linear124_amplitude_aware_timeseries.csv'),delimiter=',',names=True)
    f=ts['H_perp']/ts['H_total']; steps=len(f)-1
    stages,kind=stage_taus(f,steps)
    # 従来図
    fig,ax=plt.subplots(figsize=(10,6))
    ax.semilogy(np.arange(len(f)),np.maximum(f,1e-40),lw=0.9)
    for lab,tau in stages: ax.axvline(tau,color='r',ls=':',alpha=0.6); ax.annotate(lab.split('_')[0],(tau,1e-38),color='r',fontsize=8,rotation=90)
    ax.set_xlabel('step'); ax.set_ylabel('H_perp/H_total'); ax.set_title(f'N={N} handmade ({str(pz["design"])}) : H_perp fraction'); ax.grid(True,which='both',alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(fd,f'N{N}_Hperp_frac.png'),dpi=170); plt.close(fig)
    for col,name in [('PR_over_M','PR_over_M'),('amp_std','amp_std'),('abs_ZT_Z','abs_ZT_Z')]:
        fig,ax=plt.subplots(figsize=(9,4.5))
        (ax.semilogy if col=='abs_ZT_Z' else ax.plot)(np.arange(len(f)),np.maximum(ts[col],1e-30) if col=='abs_ZT_Z' else ts[col],lw=0.9)
        ax.set_xlabel('step'); ax.set_ylabel(name); ax.set_title(f'N={N} handmade: {name}'); ax.grid(True,alpha=.25)
        fig.tight_layout(); fig.savefig(os.path.join(fd,f'N{N}_{name}.png'),dpi=160); plt.close(fig)
    # 平面占有
    planes,kern=plane_basis(N,v); sub=np.arange(0,len(S),20)
    occ=np.zeros((len(sub),len(planes)))
    for a,t in enumerate(sub):
        Z=S[t]
        for b,(sg,u) in enumerate(planes):
            occ[a,b]=abs(np.vdot(u,Z))**2+abs(np.vdot(np.conj(u),Z))**2
    Hs=(np.abs(S[sub])**2).sum(axis=1)
    with open(os.path.join(dd,'plane_occupation.csv'),'w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['step']+[f'plane{b}_sigma{planes[b][0]:.6f}' for b in range(len(planes))])
        for a,t in enumerate(sub): w.writerow([int(t)]+list(occ[a]/Hs[a]))
    fig,ax=plt.subplots(figsize=(10,6))
    for b,(sg,u) in enumerate(planes):
        ax.semilogy(sub,np.maximum(occ[:,b]/Hs,1e-40),lw=0.9,label=f'σ/r²={sg/(r0*r0):.3f}')
    ax.set_xlabel('step'); ax.set_ylabel('plane occupation / H'); ax.set_title(f'N={N} handmade: plane occupation (basis: iK(v), step0)')
    ax.legend(fontsize=7,ncol=2); ax.grid(True,which='both',alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(fd,f'N{N}_plane_occupation.png'),dpi=170); plt.close(fig)
    # Takagi 軸スケール時系列（シリーズ定義の幾何、40 step 刻み）
    sub2=np.arange(0,len(S),40); RF=np.array([takagi_scales(S[t],N) for t in sub2]); 
    pv=v.real/np.linalg.norm(v.real); qv=v.imag-(v.imag@pv)*pv; qv/=np.linalg.norm(qv)
    RP=np.array([takagi_scales(S[t]-pv*(pv@S[t])-qv*(qv@S[t]),N) for t in sub2])
    fig,axs=plt.subplots(1,2,figsize=(12,5),sharey=True)
    for kk in range(RF.shape[1]):
        axs[0].semilogy(sub2,np.maximum(RF[:,kk],1e-40),lw=0.9)
        axs[1].semilogy(sub2,np.maximum(RP[:,kk],1e-40),lw=0.9)
    axs[0].set_title(f'N={N}: full-simplex canonical axes r_k'); axs[1].set_title('perp-component axes r_k')
    for a in axs: a.set_xlabel('step'); a.grid(True,which='both',alpha=.25)
    axs[0].set_ylabel('r_k'); fig.tight_layout(); fig.savefig(os.path.join(fd,f'N{N}_takagi_axes_timeseries.png'),dpi=170); plt.close(fig)
    with open(os.path.join(dd,'takagi_axes.csv'),'w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['step']+[f'r{k+1}_full' for k in range(RF.shape[1])]+[f'r{k+1}_perp' for k in range(RP.shape[1])])
        for a2,t in enumerate(sub2): w.writerow([int(t)]+list(RF[a2])+list(RP[a2]))
    # 単体対称性の実測：設計の対称 g（奇数 C=Z_N 巡回／偶数 A=色回転（極固定）＋位相 2π/(N−1) のねじ型）
    gperm=([(i+1)%N for i in range(N)] if N%2==1 else [(i+1)%(N-1) for i in range(N-1)]+[N-1])
    def sym_dev_t(t):
        E2=edges(N); Dm=np.zeros((N,N),complex); z=S[t]
        for val,(i,j) in zip(z*z,E2): Dm[i,j]=Dm[j,i]=val
        Dg=Dm[np.ix_(gperm,gperm)]; ip=np.vdot(Dm,Dg)
        phf=np.exp(1j*np.angle(ip)) if abs(ip)>0 else 1.0
        return float(np.linalg.norm(Dg-phf*Dm)/np.linalg.norm(Dm)), float(np.angle(ip))
    SYD=np.array([sym_dev_t(t) for t in sub2])
    fig,ax=plt.subplots(figsize=(10,5))
    ax.semilogy(sub2,np.maximum(SYD[:,0],1e-18),lw=1.0,color='tab:purple')
    phpred=(2*np.pi/(N-1)) if N%2==0 else 0.0
    ax.set_xlabel('step'); ax.set_ylabel('min_phi ||D2(g)-e^{i phi} D2|| / ||D2||')
    ax.set_title(f'N={N}: simplex symmetry deviation (screw: perm + phase; predicted phi={phpred:.3f})')
    ax.grid(True,which='both',alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(fd,f'N{N}_symmetry_deviation.png'),dpi=170); plt.close(fig)
    with open(os.path.join(dd,'symmetry_deviation.csv'),'w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['step','screw_deviation','optimal_phase']); 
        for a2,t in enumerate(sub2): w.writerow([int(t),SYD[a2,0],SYD[a2,1]])
    # 球面・ヌル錐診断：頂点ごとのエルミート半径 R_i(t) と ヌル錐逸脱 |x_i·x_i|/|B|max(t)
    def sphere_null(t):
        E2=edges(N); D2=np.zeros((N,N),complex)
        z=S[t]
        for val,(i,j) in zip(z*z,E2): D2[i,j]=D2[j,i]=val
        J2=np.eye(N)-np.ones((N,N))/N; B2=-0.5*J2@D2@J2
        sg2,V2=takagi(B2); X2=(V2*np.sqrt(sg2)[None,:]) if V2.shape[1] else np.zeros((N,0))
        R=np.linalg.norm(X2,axis=1); nd=np.abs(np.diag(B2))/max(np.abs(B2).max(),1e-300)
        return R,nd
    SN=np.array([np.concatenate(sphere_null(t)) for t in sub2])
    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for i in range(N):
        axs[0].plot(sub2,SN[:,i],lw=0.9,label=f'v{i}')
        axs[1].semilogy(sub2,np.maximum(SN[:,N+i],1e-18),lw=0.9)
    axs[0].set_ylabel('vertex Hermitian radius R_i'); axs[0].set_title(f'N={N}: sphere (all R_i equal = spherical)')
    axs[1].set_ylabel('|x_i·x_i| / |B|max'); axs[1].set_title('null-cone deviation per vertex (0 = local closure)')
    for a in axs: a.set_xlabel('step'); a.grid(True,which='both',alpha=.25)
    axs[0].legend(fontsize=6,ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(fd,f'N{N}_sphere_nullcone.png'),dpi=170); plt.close(fig)
    with open(os.path.join(dd,'sphere_nullcone.csv'),'w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['step']+[f'R_v{i}' for i in range(N)]+[f'nulldev_v{i}' for i in range(N)])
        for a2,t in enumerate(sub2): w.writerow([int(t)]+list(SN[a2]))
    dims=np.array([[*simplex_dim(RF[a2]),*simplex_dim(RP[a2])] for a2 in range(len(sub2))])
    fig,ax=plt.subplots(figsize=(10,5))
    ax.plot(sub2,dims[:,1],label='d_eff full simplex',color='tab:blue')
    ax.plot(sub2,dims[:,3],label='d_eff perp component',color='tab:red',ls='--')
    ax.axhline(N-1,color='gray',ls=':',lw=0.8); ax.annotate(f'N-1={N-1}',(sub2[-1]*0.02,N-1+0.05),fontsize=8,color='gray')
    ax.set_xlabel('step'); ax.set_ylabel('effective complex dimension d_eff'); ax.set_ylim(0,N-0.5)
    ax.set_title(f'N={N}: simplex effective dimension (d_eff=(Σσ)²/Σσ², σ=Takagi values of B)')
    ax.legend(); ax.grid(True,alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(fd,f'N{N}_simplex_dimension.png'),dpi=170); plt.close(fig)
    with open(os.path.join(dd,'simplex_dimension.csv'),'w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['step','rank_full','d_eff_full','rank_perp','d_eff_perp'])
        for a2,t in enumerate(sub2): w.writerow([int(t),int(dims[a2,0]),dims[a2,1],int(dims[a2,2]),dims[a2,3]])
    # 5 局面図
    axlim=1.1*max(float(np.abs(takagi_coords(S[tau],N,3)[1]).max()) for _,tau in stages)
    for lab,tau in stages:
        Z=S[tau]
        complex_plane_fig(v,Z,r0,os.path.join(fd,f'N{N}_{lab}_tau{tau}_complexplane.png'),f'N={N} {lab} τ={tau}  (f={f[tau]:.2e})')
        geometry_fig(N,v,Z,r0,os.path.join(fd,f'N{N}_{lab}_tau{tau}_geometry.png'),f'N={N} {lab} τ={tau}',axlim=axlim)
        st_rows.append([N,str(pz['design']),kind,lab,tau,float(f[tau])])
    gf=json.load(open(os.path.join(dd,'summary.json')))['treatment']
    an_rows.append([N,str(pz['design']),kind,float(np.nanmax(f)),int(np.nanargmax(f)),float(f[-1]),
                    (gf['growth_fit'] or {}).get('slope_ln_Hperp_per_step'),gf['onset_Hperp_fraction_gt_0.05'],gf['Htotal_max_abs_drift']])
    print(f"N={N}: kind={kind} stages={[(l,t) for l,t in stages]} max f={np.nanmax(f):.2e} final={f[-1]:.2e}")
os.makedirs(os.path.join(ROOT,'results'),exist_ok=True)
with open(os.path.join(ROOT,'results','stage_taus.csv'),'w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['N','design','kind','stage','tau','Hperp_frac']); w.writerows(st_rows)
with open(os.path.join(ROOT,'results','analysis_summary.csv'),'w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['N','design','kind','max_frac','max_frac_step','final_frac','growth_rate_per_step','onset_frac_gt5pct','Htotal_drift']); w.writerows(an_rows)
print("PASS3 OK")

# ---- 合成インフレーション図（全 N 重ね描き、従来の比較図と同型）----
def composite():
    fig,ax=plt.subplots(figsize=(11,6.5))
    colors={4:'tab:gray',5:'tab:green',6:'tab:blue',7:'tab:olive',8:'tab:red'}
    for N in [4,5,6,7,8]:
        ts=np.genfromtxt(os.path.join(ROOT,'data',f'N{N}','treatment_linear124_amplitude_aware_timeseries.csv'),delimiter=',',names=True)
        f=ts['H_perp']/ts['H_total']; design='A' if N%2==0 else 'C'
        ax.semilogy(np.arange(len(f)),np.maximum(f,1e-40),lw=1.1,color=colors[N],label=f'N={N} (design {design})')
    ax.set_xlabel('step'); ax.set_ylabel('H_perp/H_total'); ax.set_ylim(1e-36,3)
    ax.set_title('Handmade parents N=4..8 (seedless, 2pi/124, amplitude-aware K): even/odd stripe of inflation')
    ax.legend(); ax.grid(True,which='both',alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT,'figures','composite_Hperp_frac_N4toN8.png'),dpi=180); plt.close(fig)
composite(); print('composite figure written')
