#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""複素シンプレックス基礎：N=3 の全計算展開と図化。
波の値 z_e（大きさ r 共通、位相 0/60/120°）→ |z| → d²=z² → 閉塞 → 和則 → B → Takagi 軸 → 埋め込み座標。
スケールは自由（比だけが物理）なので、数値例として r²=1/15 を採る（‖v‖²=3r²=1/5）。"""
import os, json
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=os.path.join(ROOT,'figures')
N=3; E=[(0,1),(0,2),(1,2)]; M=3
r=1.0/np.sqrt(15.0)
NORM=np.sqrt(3.0)*r
theta=np.radians([120.0,0.0,60.0])   # 辺 (0,1)=120°, (0,2)=0°, (1,2)=60°（実親と同じ割り当て、全体回転は除く）
z=r*np.exp(1j*theta)
print(f"r = 1/√15 = {r:.9f}   ‖v‖ = √(3r²) = √(1/5) = {NORM:.9f}")
for k,(e,th) in enumerate(zip(E,np.degrees(theta))):
    print(f"辺{e} θ={th:5.1f}°: a={z[k].real:+.9f} b={z[k].imag:+.9f} |z|=√(a²+b²)=√({z[k].real**2:.9f}+{z[k].imag**2:.9f})=√{abs(z[k])**2:.9f}={abs(z[k]):.9f}")
d2=z*z
print("d²=z²:", [f"{e}: {d2[k].real:+.9f}{d2[k].imag:+.9f}i (|d²|={abs(d2[k]):.9f}, 角度={np.degrees(np.angle(d2[k]))%360:.1f}°)" for k,e in enumerate(E)])
print(f"閉塞 Σz² = {d2.sum().real:+.2e}{d2.sum().imag:+.2e}i")
print("局所閉塞（頂点ごと Σ_j z_ij²）:")
for i in range(N):
    s=sum(d2[k] for k,(a,b) in enumerate(E) if i in (a,b))
    print(f"  頂点{i}: {s.real:+.9f}{s.imag:+.9f}i  |・|/r²={abs(s)/r**2:.4f}")
# 和則
print("和則（各辺の隣接 2 本）:")
for k,(a,b) in enumerate(E):
    adj=[j for j,(c,d) in enumerate(E) if j!=k and ({a,b}&{c,d})]
    s2=sum(np.sin(theta[j]-theta[k])**2 for j in adj); s22=sum(np.sin(2*(theta[j]-theta[k])) for j in adj)
    print(f"  辺{E[k]}: Σsin²φ={s2:.6f}  Σsin2φ={s22:+.2e}   （N−1=2 でなく 3/2 ！）")
# 生成子 K と μ
A=np.zeros((M,M))
for a2 in range(M):
    for b2 in range(M):
        if a2!=b2 and set(E[a2])&set(E[b2]): A[a2,b2]=1.0
K=A*np.imag(np.conj(z)[:,None]*z[None,:]); hv=1j*(K@z); mu=(np.vdot(z,hv)/np.vdot(z,z)).real
res=np.linalg.norm(hv-mu*z)/np.linalg.norm(z)
print(f"自己無撞着: μ={mu:+.9f}  μ/r²={mu/r**2:+.6f} (=−3/2)  残差={res:.1e}")
# B と Takagi
D=np.zeros((N,N),complex)
for val,(i,j) in zip(d2,E): D[i,j]=D[j,i]=val
J=np.eye(N)-np.ones((N,N))/N; Bm=-0.5*J@D@J
print("B = −½JD²J ="); 
for row in Bm: print("   [",", ".join(f"{x.real:+.6f}{x.imag:+.6f}i" for x in row),"]")
s=np.linalg.svd(Bm,compute_uv=False)
print(f"Takagi 値 σ = {np.round(s,9).tolist()}  → 軸スケール r_k=√σ = {np.round(np.sqrt(s),9).tolist()}")
print(f"rank（>1e-12）= {(s>1e-12).sum()} ＝ 1 ！（N−1=2 でない：N=3 の例外）")
# 座標（Autonne–Takagi）
Mre=np.block([[Bm.real,Bm.imag],[Bm.imag,-Bm.real]]); w,V=np.linalg.eigh(Mre)
idx=np.argsort(w)[::-1][:N]; X=[]
for j in idx:
    if w[j]<=1e-14: continue
    vv=V[:N,j]+1j*V[N:,j]; vv/=np.linalg.norm(vv); X.append(np.sqrt(w[j])*vv)
X=np.array(X).T
print("埋め込み座標（第 1 正準軸の複素平面、頂点 0,1,2）:", np.round(X[:,0],6).tolist())
err=max(abs(((X[i]-X[j])@(X[i]-X[j]))-D[i,j]) for i,j in E)
print(f"距離再現誤差 max = {err:.1e}")
# ---- 図 1: z の複素平面
fig,ax=plt.subplots(figsize=(5.6,5.6))
th=np.linspace(0,2*np.pi,200); ax.plot(r*np.cos(th),r*np.sin(th),'--',color='gray',lw=0.8)
cols=['tab:red','tab:blue','tab:green']
for k,e in enumerate(E):
    ax.annotate('',xy=(z[k].real,z[k].imag),xytext=(0,0),arrowprops=dict(arrowstyle='->',color=cols[k],lw=2))
    ax.annotate(f'z{e} θ={np.degrees(theta[k]):.0f}°',(z[k].real*1.15,z[k].imag*1.15),fontsize=9,color=cols[k],ha='center')
ax.set_xlim(-.36,.36); ax.set_ylim(-.36,.36); ax.set_aspect('equal'); ax.grid(alpha=.2)
ax.axhline(0,color='k',lw=.4); ax.axvline(0,color='k',lw=.4)
ax.set_xlabel('a = Re z'); ax.set_ylabel('b = Im z'); ax.set_title(f'N=3: wave values z_e (|z|={r:.4f}, phases 0/60/120 deg)')
fig.tight_layout(); fig.savefig(os.path.join(FIG,'N3_fig1_z_complex_plane.png'),dpi=170); plt.close(fig)
# ---- 図 2: d² の閉塞（ベクトルを首尾接続すると閉じた三角形）
fig,ax=plt.subplots(figsize=(5.6,5.6))
p=0+0j
for k,e in enumerate(E):
    ax.annotate('',xy=((p+d2[k]).real,(p+d2[k]).imag),xytext=(p.real,p.imag),arrowprops=dict(arrowstyle='->',color=cols[k],lw=2))
    mid=p+d2[k]/2; ax.annotate(f'd²{e}',(mid.real*1.0+0.004,mid.imag+0.004),fontsize=9,color=cols[k])
    p+=d2[k]
ax.scatter([0],[0],color='k',s=40,zorder=3); ax.annotate('start = end (Σz²=0)',(0.004,-0.012),fontsize=9)
ax.set_aspect('equal'); ax.grid(alpha=.2); lim=0.09
ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
ax.set_xlabel('Re d²'); ax.set_ylabel('Im d²'); ax.set_title('N=3 closure: the three d² vectors sum to zero',fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG,'N3_fig2_d2_closure_triangle.png'),dpi=170); plt.close(fig)
# ---- 図 3: 埋め込み（第 1 正準軸の複素平面に 3 頂点、辺の色=位相）
fig,ax=plt.subplots(figsize=(5.6,5.6))
c1=X[:,0]
for k,(i,j) in enumerate(E):
    ax.plot([c1[i].real,c1[j].real],[c1[i].imag,c1[j].imag],color=cols[k],lw=2,alpha=.8)
ax.scatter(c1.real,c1.imag,color='k',s=45,zorder=3)
for i in range(N): ax.annotate(f'v{i}',(c1[i].real,c1[i].imag),textcoords='offset points',xytext=(7,7),fontsize=10)
ax.set_aspect('equal'); ax.grid(alpha=.2); ax.axhline(0,color='k',lw=.4); ax.axvline(0,color='k',lw=.4)
ax.set_xlabel('Re x (canonical axis 1)'); ax.set_ylabel('Im x'); ax.set_title('N=3 geometry: equilateral triangle in ONE complex axis (rank 1)',fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG,'N3_fig3_geometry_rank1.png'),dpi=170); plt.close(fig)
print("figures written")
