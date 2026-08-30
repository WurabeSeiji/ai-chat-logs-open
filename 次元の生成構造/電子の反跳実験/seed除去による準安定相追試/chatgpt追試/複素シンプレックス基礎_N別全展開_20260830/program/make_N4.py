#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""複素シンプレックス基礎：N=4 の全計算展開と図化（make_N3.py と同じ要領）。
設計：3 つの完全マッチング（対辺ペア）に位相 0/60/120°。スケール規約 r²=1/15 を全 N で共通に使う。"""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=os.path.join(ROOT,'figures')
N=4; E=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]; M=6
r=1.0/np.sqrt(15.0)
# 完全マッチング（対辺ペア）: {(0,3),(1,2)}=色0, {(0,2),(1,3)}=色1, {(0,1),(2,3)}=色2
color={(0,3):0,(1,2):0,(0,2):1,(1,3):1,(0,1):2,(2,3):2}
theta=np.array([np.radians(60.0*color[e]) for e in E])
z=r*np.exp(1j*theta)
print(f"r = 1/√15 = {r:.9f}   ‖v‖ = √(6r²) = √(6/15) = {np.sqrt(6)*r:.9f}")
print("色（マッチング）と位相:")
for e,c in sorted(color.items(),key=lambda x:(x[1],x[0])): print(f"  組{e}: 色{c} → θ={60*c}°")
for k,e in enumerate(E):
    print(f"組{e} θ={np.degrees(theta[k]):5.1f}°: a={z[k].real:+.9f} b={z[k].imag:+.9f} |z|={abs(z[k]):.9f}")
d2=z*z
print("d²=z²:")
for k,e in enumerate(E): print(f"  {e}: {d2[k].real:+.9f}{d2[k].imag:+.9f}i (角度 {np.degrees(np.angle(d2[k]))%360:.0f}°)")
print(f"大域閉塞 Σz² = {d2.sum().real:+.2e}{d2.sum().imag:+.2e}i  （各角度 0/120/240° が 2 本ずつ → 2r²(1+ω+ω²)=0）")
print("局所閉塞（頂点ごと）:")
for i in range(N):
    ssum=sum(d2[k] for k,(a,b) in enumerate(E) if i in (a,b))
    print(f"  頂点{i}: {ssum.real:+.2e}{ssum.imag:+.2e}i  {'成立' if abs(ssum)<1e-15 else '破れ'}  （各頂点は 3 色を 1 本ずつ持つ → r²(1+ω+ω²)=0）")
print("和則（各組の隣接 4 本）:")
for k,(a,b) in enumerate(E):
    adj=[j for j,(c,d) in enumerate(E) if j!=k and ({a,b}&{c,d})]
    s2=sum(np.sin(theta[j]-theta[k])**2 for j in adj); s22=sum(np.sin(2*(theta[j]-theta[k])) for j in adj)
    print(f"  組{E[k]}: 隣接={[E[j] for j in adj]}  Σsin²φ={s2:.6f} (=N−1=3)  Σsin2φ={s22:+.1e}")
A=np.zeros((M,M))
for a2 in range(M):
    for b2 in range(M):
        if a2!=b2 and set(E[a2])&set(E[b2]): A[a2,b2]=1.0
K=A*np.imag(np.conj(z)[:,None]*z[None,:]); hv=1j*(K@z); mu=(np.vdot(z,hv)/np.vdot(z,z)).real
print(f"自己無撞着: μ={mu:+.9f}  μ/r²={mu/r**2:+.6f} (=−(N−1)=−3)  残差={np.linalg.norm(hv-mu*z)/np.linalg.norm(z):.1e}")
D=np.zeros((N,N),complex)
for val,(i,j) in zip(d2,E): D[i,j]=D[j,i]=val
J=np.eye(N)-np.ones((N,N))/N; Bm=-0.5*J@D@J
print("B = −½JD²J =")
for row in Bm: print("   [",", ".join(f"{x.real:+.6f}{x.imag:+.6f}i" for x in row),"]")
sv=np.linalg.svd(Bm,compute_uv=False)
print(f"Takagi 値 σ = {np.round(sv,9).tolist()}  軸スケール √σ = {np.round(np.sqrt(sv),9).tolist()}")
print(f"rank = {(sv>1e-12).sum()} = N−1 = 3（N=3 と違い最大次元）")
Mre=np.block([[Bm.real,Bm.imag],[Bm.imag,-Bm.real]]); w,V=np.linalg.eigh(Mre)
idx=np.argsort(w)[::-1][:N]; X=[]
for j in idx:
    if w[j]<=1e-14: continue
    vv=V[:N,j]+1j*V[N:,j]; vv/=np.linalg.norm(vv); X.append(np.sqrt(w[j])*vv)
X=np.array(X).T
err=max(abs(((X[i]-X[j])@(X[i]-X[j]))-D[i,j]) for i,j in E)
print(f"埋め込み（3 複素軸）: 距離再現誤差 max = {err:.1e}")
print("頂点のエルミート半径:", np.round(np.linalg.norm(X,axis=1),9).tolist())
nulldev=[abs(X[i]@X[i]) for i in range(N)]
print("頂点の x·x（光円錐なら 0）:", [f"{x:.1e}" for x in nulldev])
cols=['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown']
# 図1: z 複素平面（マッチング色で）
fig,ax=plt.subplots(figsize=(5.6,5.6))
th=np.linspace(0,2*np.pi,200); ax.plot(r*np.cos(th),r*np.sin(th),'--',color='gray',lw=0.8)
mc={0:'tab:red',1:'tab:blue',2:'tab:green'}
seen=set()
for k,e in enumerate(E):
    c=mc[color[e]]
    ax.annotate('',xy=(z[k].real,z[k].imag),xytext=(0,0),arrowprops=dict(arrowstyle='->',color=c,lw=2,alpha=0.8))
    if color[e] not in seen:
        ax.annotate(f'matching {color[e]} ({60*color[e]}°) ×2',(z[k].real*1.28,z[k].imag*1.28),fontsize=8,color=c,ha='center'); seen.add(color[e])
ax.set_xlim(-.36,.36); ax.set_ylim(-.36,.36); ax.set_aspect('equal'); ax.grid(alpha=.2)
ax.axhline(0,color='k',lw=.4); ax.axvline(0,color='k',lw=.4)
ax.set_xlabel('a = Re z'); ax.set_ylabel('b = Im z'); ax.set_title('N=4: 6 waves = 3 matchings x 2, phases 0/60/120 deg',fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG,'N4_fig1_z_complex_plane.png'),dpi=170); plt.close(fig)
# 図2: d² 閉塞（6 本の首尾接続六角形）
fig,ax=plt.subplots(figsize=(5.6,5.6))
order=np.argsort(np.angle(d2))
pos=0+0j
for k in order:
    c=mc[color[E[k]]]
    ax.annotate('',xy=((pos+d2[k]).real,(pos+d2[k]).imag),xytext=(pos.real,pos.imag),arrowprops=dict(arrowstyle='->',color=c,lw=2))
    pos+=d2[k]
ax.scatter([0],[0],color='k',s=40,zorder=3); ax.annotate('start = end (Σz²=0)',(0.004,0.004),fontsize=9)
ax.set_aspect('equal'); ax.grid(alpha=.2); lim=0.16
ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
ax.set_xlabel('Re d²'); ax.set_ylabel('Im d²'); ax.set_title('N=4 closure: six d² vectors close (2 per angle)',fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG,'N4_fig2_d2_closure.png'),dpi=170); plt.close(fig)
# 図3: 幾何（3 複素軸パネル＋正四面体の骨格は模式で対辺色分け）
fig,axs=plt.subplots(1,3,figsize=(13,4.6))
for kk in range(3):
    ax=axs[kk]; c1=X[:,kk]
    for k,(i,j) in enumerate(E):
        ax.plot([c1[i].real,c1[j].real],[c1[i].imag,c1[j].imag],color=mc[color[E[k]]],lw=1.6,alpha=.75)
    ax.scatter(c1.real,c1.imag,color='k',s=40,zorder=3)
    for i in range(N): ax.annotate(f'v{i}',(c1[i].real,c1[i].imag),textcoords='offset points',xytext=(6,6),fontsize=9)
    ax.set_aspect('equal'); ax.grid(alpha=.2); ax.set_title(f'canonical axis {kk+1} (r={np.sqrt(sv[kk]):.4f})',fontsize=9)
    lim=0.30; ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
fig.suptitle('N=4 factual geometry: vertices in 3 equal complex axes (rank 3, round)',fontsize=10)
fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(os.path.join(FIG,'N4_fig3_geometry_axes.png'),dpi=170); plt.close(fig)
print("figures written")
