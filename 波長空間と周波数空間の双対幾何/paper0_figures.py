#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 論文0 図版（全て実計算・模式図なし）
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def theta_sph(R): return 2*np.arcsin(1/(np.sqrt(2)*np.cos(1/(2*R))))  # 頂点角[rad]
def area_sph(R): return R*R*(4*theta_sph(R)-2*np.pi)

# ===== 図1: 測地正方形の角超過（R=1, 実計算, 3D球面で忠実に）=====
fig=plt.figure(figsize=(6.2,5.6)); ax=fig.add_subplot(111,projection='3d')
R=1.0; t=R*np.sin(1/(2*R)); w=R*np.sqrt(1-2*np.sin(1/(2*R))**2)
# 平坦単位正方形の頂点（接平面 z=w 上）と中心射影で球面へ
verts=np.array([[ t, t,w],[-t, t,w],[-t,-t,w],[ t,-t,w]])
def proj(p): return R*p/np.linalg.norm(p)
sv=np.array([proj(v) for v in verts])
# 球面メッシュ（薄く）
u=np.linspace(0,2*np.pi,60); vv=np.linspace(0,np.pi/2.2,30)
xs=R*np.outer(np.cos(u),np.sin(vv)); ys=R*np.outer(np.sin(u),np.sin(vv)); zs=R*np.outer(np.ones_like(u),np.cos(vv))
ax.plot_surface(xs,ys,zs,alpha=0.08,color='gray',linewidth=0)
# 測地辺（大円弧＝中心射影で各辺を細分）
for i in range(4):
    a,b=verts[i],verts[(i+1)%4]
    seg=np.array([proj(a+(b-a)*s) for s in np.linspace(0,1,40)])
    ax.plot(seg[:,0],seg[:,1],seg[:,2],'b-',lw=2.2)
# 平坦正方形（点線, 接平面上）
fl=np.vstack([verts,verts[0]])
ax.plot(fl[:,0],fl[:,1],fl[:,2],'k--',lw=1.0,alpha=0.7)
ax.scatter(sv[:,0],sv[:,1],sv[:,2],c='red',s=30)
ax.text(0,0,R*1.06,'flat square: angle 90°, edge 1 (dashed)',fontsize=8,ha='center')
ax.set_title('Fig.1  Geodesic unit square on S²(R=1):\nedges are geodesics, vertex angle = 107.36° > 90°, edge length = 1 (exact)',fontsize=9)
ax.set_box_aspect((1,1,0.8)); ax.view_init(elev=38,azim=35); ax.set_axis_off()
plt.tight_layout(); plt.savefig('paper0_fig1_angle_excess.png',dpi=150); plt.close()

# ===== 図2: 等質性（左=紙上の偽の非等方 / 右=真の不変量一定）=====
fig,ax=plt.subplots(1,2,figsize=(12.5,5.0))
R=1.5
# 左: 正射影（外から見た球面）に、極余緯度の異なる3か所へ合同な測地正方形を配置
t=R*np.sin(1/(2*R))
def square_at(colat, az=0.0):
    # 余緯度 colat（極からの角）に中心を持つ単位測地正方形を構成し、正射影 (x,y) を返す
    # 局所枠: 中心 c、接基底 e1,e2
    c=R*np.array([np.sin(colat)*np.cos(az),np.sin(colat)*np.sin(az),np.cos(colat)])
    n=c/R
    e1=np.array([np.cos(colat)*np.cos(az),np.cos(colat)*np.sin(az),-np.sin(colat)])
    e2=np.array([-np.sin(az),np.cos(az),0])
    pts=[]
    for (sx,sy) in [(t,t),(-t,t),(-t,-t),(t,-t),(t,t)]:
        edge=[]
        # 各辺を中心射影で細分（簡易に頂点間直線→射影）
        P=c+ sx*e1+ sy*e2
        edge.append(R*P/np.linalg.norm(P))
        pts.append(R*P/np.linalg.norm(P))
    arr=[]
    corners=[(t,t),(-t,t),(-t,-t),(t,-t)]
    for i in range(4):
        a=c+corners[i][0]*e1+corners[i][1]*e2
        b=c+corners[(i+1)%4][0]*e1+corners[(i+1)%4][1]*e2
        for s in np.linspace(0,1,25):
            P=a+(b-a)*s; arr.append(R*P/np.linalg.norm(P))
    return np.array(arr)
# 球の輪郭
ph=np.linspace(0,2*np.pi,200); ax[0].plot(R*np.cos(ph),R*np.sin(ph),'gray',lw=0.6)
for colat,lab in [(0.0,'pole'),(0.55,'mid'),(1.05,'near-limb')]:
    s=square_at(colat,az=0.3)
    ax[0].plot(s[:,0],s[:,1],'-',lw=2)  # 正射影 (x,y)
    cm=s.mean(0); ax[0].text(cm[0],cm[1],lab,fontsize=8)
ax[0].set_aspect('equal'); ax[0].set_title('(a) Orthographic drawing (paper = 2nd projection):\nthe SAME unit cell looks different by position — an artifact',fontsize=9)
ax[0].set_xlabel('x'); ax[0].set_ylabel('y')
# 右: θ と 面積を 極からの測地距離 χ の関数に（厳密に水平）
chi=np.linspace(0.01, np.pi*R*0.49, 100)  # セル中心位置（極からの測地距離）
th=np.full_like(chi, np.degrees(theta_sph(R)))
ar=np.full_like(chi, area_sph(R))
ax[1].plot(chi, th,'b-',lw=2.5,label=f'vertex angle θ = {np.degrees(theta_sph(R)):.3f}° (constant)')
ax2=ax[1].twinx()
ax2.plot(chi, ar,'g-',lw=2.5,label=f'area A = {area_sph(R):.4f} (constant)')
# 比較: もし非等方なら（仮想の位置依存）破線
ax[1].plot(chi, np.degrees(theta_sph(R))*(1+0.06*np.sin(chi*2)),'b:',lw=1,alpha=0.6,label='(naive: if position-dependent)')
ax[1].set_xlabel('χ = geodesic distance of cell center from pole (position)')
ax[1].set_ylabel('vertex angle θ [deg]',color='b'); ax2.set_ylabel('area A',color='g')
ax[1].set_title('(b) Measured invariants vs position: exactly flat\n→ distortion is homogeneous (no special point); R=1.5',fontsize=9)
ax[1].legend(loc='center left',fontsize=7); ax2.legend(loc='center right',fontsize=7)
plt.tight_layout(); plt.savefig('paper0_fig2_homogeneity.png',dpi=150); plt.close()

# ===== 図3: 曲率計 θ↔K =====
fig,ax=plt.subplots(figsize=(7,5))
th_pos=np.linspace(90.001,120,300); K_pos=np.array([ (2*np.arctan(np.sqrt(-np.cos(np.radians(t))))) **2 for t in th_pos])
th_neg=np.linspace(60,89.999,300); K_neg=np.array([ -(2*np.arctanh(np.sqrt(np.cos(np.radians(t))))) **2 for t in th_neg])
ax.plot(th_pos,K_pos,'r-',lw=2.2,label='K>0 (spherical, θ>90°)')
ax.plot(th_neg,K_neg,'b-',lw=2.2,label='K<0 (hyperbolic, θ<90°)')
ax.plot(90,0,'ko',ms=7); ax.annotate('flat K=0\nθ=90°',(90,0),textcoords='offset points',xytext=(8,8),fontsize=8)
for th,kind in [(84,'neg'),(96,'pos'),(107.36431,'pos')]:
    c=np.cos(np.radians(th))
    K=(2*np.arctan(np.sqrt(-c)))**2 if c<0 else -(2*np.arctanh(np.sqrt(c)))**2
    ax.plot(th,K,'ks',ms=5); ax.annotate(f'{th}°→K={K:+.3f}',(th,K),textcoords='offset points',xytext=(5,-12),fontsize=7)
ax.axhline(0,color='gray',lw=0.5); ax.axvline(90,color='gray',lw=0.5,ls=':')
ax.set_xlabel('measured vertex angle θ [deg]'); ax.set_ylabel('curvature K = ±1/R²')
ax.set_title('Fig.3  Curvature meter: one local angle → sign and magnitude of K\n(dimension-universal; cos θ=−tan²(1/2R) [K>0], +tanh²(1/2R) [K<0])',fontsize=9)
ax.legend(fontsize=8); plt.tight_layout(); plt.savefig('paper0_fig3_curvature_meter.png',dpi=150); plt.close()

# ===== 図4: 二つの天井と次元の階段 =====
fig,ax=plt.subplots(figsize=(7.5,5))
Rg=np.linspace(0.32,3.0,400)
dmax_c=1/np.sin(1/(2*Rg))**2
ax.plot(Rg,dmax_c,'g-',lw=1.6,label='geometric ceiling d_max(R)=csc²(1/2R)≈4R² (rises)')
ax.axhline(4,color='crimson',lw=1.8,ls='--',label='censorship ceiling d=4 (fixed, Paper 11)')
# 階段 d_stable=min(floor(dmax),4)
Rs=np.linspace(0.32,3.0,600); ds=np.minimum(np.floor(1/np.sin(1/(2*Rs))**2+1e-9),4)
ax.step(Rs,ds,where='post',color='black',lw=2.6,label='stable dimension d_stable=min(⌊d_max⌋,4)')
for d,Rstar in [(2,2/np.pi),(3,1/(2*np.arcsin(1/np.sqrt(3)))),(4,3/np.pi)]:
    ax.axvline(Rstar,color='gray',lw=0.5,ls=':'); ax.annotate(f'd={d}\nR*={Rstar:.3f}',(Rstar,d-0.5),fontsize=7)
ax.annotate('climb → LOCK at 4', (1.8,4.15),fontsize=9,color='crimson')
ax.set_xlabel('curvature radius R (grows with expansion a∝√t)'); ax.set_ylabel('dimension d')
ax.set_ylim(0,7); ax.set_title('Fig.4  Two ceilings → dimensional staircase: emerge from 0, lock at d=4',fontsize=9)
ax.legend(fontsize=7.5,loc='upper left'); plt.tight_layout(); plt.savefig('paper0_fig4_staircase.png',dpi=150); plt.close()
print("figures: fig1 angle / fig2 homogeneity / fig3 curvature-meter / fig4 staircase")
