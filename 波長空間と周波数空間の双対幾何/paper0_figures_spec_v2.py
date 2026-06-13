#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 論文0 図版 A〜D（仕様書 paper0_figures_spec_for_claudecode.md を厳密実装、R=3固定、全て実計算）
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=3.0

# ===== 図A: 1+1次元 円周→線分（等間隔保存、両パネル同一スケール）=====
# 横並びでなく縦並びにし、x軸スケールを厳密一致させる（弧長=線分間隔が目で一致するように）
fig=plt.figure(figsize=(11,7.2))
ax_c=fig.add_subplot(2,1,1)   # 円（上）
ax_l=fig.add_subplot(2,1,2)   # 線分（下）, x軸は測地距離で円と同尺
ph=np.linspace(0,2*np.pi,400); ax_c.plot(R*np.cos(ph),R*np.sin(ph),'k-',lw=1.2)
ks=np.arange(12); pts=np.array([[R*np.cos(np.deg2rad(30*k)),R*np.sin(np.deg2rad(30*k))] for k in ks])
ax_c.scatter(pts[:,0],pts[:,1],c='black',s=40,zorder=3)
for k in ks:  # 各点に弧長目盛り（測地距離 = R*角）を添える
    pass
ax_c.annotate('arc length π/2 between neighbors',(R*np.cos(np.deg2rad(15)),R*np.sin(np.deg2rad(15))),textcoords='offset points',xytext=(8,4),fontsize=8)
ax_c.plot(R*np.cos(np.pi),R*np.sin(np.pi),'rx',ms=11,mew=2); ax_c.annotate('cut here (180°)',(-R,0),textcoords='offset points',xytext=(-8,12),fontsize=8,color='r')
ax_c.set_aspect('equal'); ax_c.set_title('(top) circle S¹(R=3): 12 points at 30° (geodesic arc π/2 between neighbors)',fontsize=9)
ax_c.set_xlabel('x'); ax_c.set_ylabel('y')
# 線分: 測地距離を横軸に。点間隔 π/2、全長 2πR=6π。x軸範囲を [−R, 6π−R] 等にして弧長尺を円と揃える
xseg=np.array([np.pi/2*k for k in range(12)])
ax_l.plot([xseg[0],xseg[-1]+np.pi/2],[0,0],'k-',lw=1.2)
ax_l.scatter(xseg,np.zeros_like(xseg),c='black',s=40,zorder=3)
# 各区間に弧長 π/2 の寸法線
for k in range(11):
    ax_l.annotate('',(xseg[k+1],0.35),(xseg[k],0.35),arrowprops=dict(arrowstyle='<->',lw=0.6,color='tab:blue'))
ax_l.annotate('each gap = π/2 ≈ 1.571 (identical to the arc length above)',(xseg[3],0.45),fontsize=8,color='tab:blue')
ax_l.annotate('total length = circumference 2πR = 6π ≈ 18.85',(xseg[5],-0.5),fontsize=8)
ax_l.set_ylim(-1.0,1.0); ax_l.set_yticks([])
ax_l.set_xlabel('geodesic distance along the unrolled line (same length scale as the circle radius axis)')
ax_l.set_title('(bottom) unrolled to a straight line: every gap is exactly π/2 — spacing is unchanged, NOT shrunk',fontsize=9)
ax_l.set_aspect('equal')   # ★ 円と同じ等積スケール → 間隔が目で厳密一致
fig.suptitle('Fig.A  The 1+1D map (R=3): geodesic-length-preserving — spacing is invariant (reference, zero distortion)',fontsize=10)
plt.tight_layout(); plt.savefig('paper0_figA_1d_reference.png',dpi=150); plt.close()

# ===== 図B: 2+1次元 展開図（平面図+側面図, かまぼこ膨らみ）=====
fig,ax=plt.subplots(2,1,figsize=(11,6.0),gridspec_kw={'height_ratios':[1.4,1]})
wband=np.pi/2  # 帯の測地幅 π/2
# 平面図: 長さ方向12点×幅方向（膨らみは載せない=正射影でwを潰す）
xseg=np.array([np.pi/2*k for k in range(12)])
for yv in np.linspace(-wband/2,wband/2,4):
    ax[0].plot([0,2*np.pi*R],[yv,yv],'gray',lw=0.6)
for xv in xseg:
    ax[0].plot([xv,xv],[-wband/2,wband/2],'gray',lw=0.6)
ax[0].scatter(np.repeat(xseg,2),np.tile([-wband/2,wband/2],12),c='black',s=14)
ax[0].set_ylim(-wband,wband); ax[0].set_aspect('equal'); ax[0].set_xlabel('length direction (geodesic)'); ax[0].set_ylabel('width π/2')
ax[0].set_title('(top) plan view: width π/2 × length 2πR=6π — nearly uniform grid (bulge projected out)',fontsize=9)
# 側面図: かまぼこ profile（幅方向断面の w 高さ）。長さ軸と共有
chi=np.linspace(-wband/2,wband/2,100)
h_prof=R*(np.cos(0)-np.cos(chi/R)); h_prof=R*(1-np.cos(np.abs(chi)/R))  # 中心で0…ではなく縁で持ち上がり? 仕様: 中心最大
# 仕様: 帯中心(幅0)でw最大、縁(±π/4)で基準。中心からの弧頂持ち上がり = R(1-cos(χ/R))、χは縁からの距離
h_prof = R*(1-np.cos((wband/2 - np.abs(chi))/R))  # 縁で0, 中心で最大 h
hmax=R*(1-np.cos((wband/2)/R))
# 側面図は「幅方向断面」を長さ方向の代表位置に描く。ここでは幅方向プロファイルを1枚示す。
ax[1].fill_between(chi, 0, h_prof, color='tab:blue', alpha=0.25)
ax[1].plot(chi,h_prof,'b-',lw=2)
ax[1].annotate(f'bulge height h = R(1−cos(π/12)) = {hmax:.5f}',(0,hmax),textcoords='offset points',xytext=(10,2),fontsize=8)
ax[1].set_xlabel('width direction (geodesic, span π/2)'); ax[1].set_ylabel('curvature direction w')
ax[1].set_title('(side) bulge profile: kamaboko height h≈0.102 (small vs width π/2). SAME for all d≥2 (§2 dim-independence)',fontsize=9)
fig.suptitle('Fig.B  Unfolded 2+1D map (R=3): plan grid stays uniform, distortion isolated into the side bulge',fontsize=10)
plt.tight_layout(); plt.savefig('paper0_figB_band_unfold.png',dpi=150); plt.close()

# ===== 図C: 測地正方形の角超過（R=3 主 ｜ R=1 視認用、2パネル）=====
fig,axc=plt.subplots(1,2,figsize=(11,5.6))
def draw_geo_square(ax,Rv,title):
    t=Rv*np.sin(1/(2*Rv)); w=Rv*np.sqrt(1-2*np.sin(1/(2*Rv))**2)
    th=np.degrees(2*np.arcsin(1/(np.sqrt(2)*np.cos(1/(2*Rv)))))
    Ar=Rv*Rv*(4*(2*np.arcsin(1/(np.sqrt(2)*np.cos(1/(2*Rv)))))-2*np.pi)
    fl=np.array([[0.5,0.5],[-0.5,0.5],[-0.5,-0.5],[0.5,-0.5],[0.5,0.5]])
    ax.plot(fl[:,0],fl[:,1],'--',color='gray',lw=1.3,label='flat unit square (90 deg)')
    verts=np.array([[ t, t,w],[-t, t,w],[-t,-t,w],[ t,-t,w]])
    proj=lambda p: Rv*p/np.linalg.norm(p)
    sv=np.array([proj(v) for v in verts]); scale=0.5/sv[0,0]
    first=True
    for i in range(4):
        aa,bb=verts[i],verts[(i+1)%4]
        seg=np.array([proj(aa+(bb-aa)*ss) for ss in np.linspace(0,1,60)])[:,:2]*scale
        ax.plot(seg[:,0],seg[:,1],'b-',lw=2.6,label=('geodesic square (%.2f deg)'%th) if first else None); first=False
    ax.annotate('angle %.2f deg (=90+%.2f)\narea %.4f (+%.1f%%)\nedge length 1 (geodesic)'%(th,th-90,Ar,(Ar-1)*100),(0,-0.66),fontsize=8,ha='center')
    ax.set_aspect('equal'); ax.set_xlim(-0.95,0.95); ax.set_ylim(-0.95,0.95)
    ax.set_title(title,fontsize=9); ax.legend(fontsize=7.5,loc='upper right')
    return th,Ar
th3,A3=draw_geo_square(axc[0],3.0,'(a) R=3 (mild, the series regime): angle 91.62 deg,\nbarely visible - distortion ~1/R^2 is small')
th1,A1=draw_geo_square(axc[1],1.0,'(b) R=1 (strong, for visibility): angle 107.36 deg,\narea +21% - same effect, magnified')
fig.suptitle('Fig.C  Geodesic unit square: edges are geodesics, yet interior angle > 90 deg (grows as R shrinks)',fontsize=10)
plt.tight_layout(); plt.savefig('paper0_figC_angle_excess.png',dpi=150); plt.close()

# ===== 図D: c_d 棒グラフ =====
fig,ax=plt.subplots(figsize=(6,4))
ds=[1,2,3,4,5]; cd=[d*(d-1)/12 for d in ds]
bars=ax.bar(ds,cd,color=['gray','tab:blue','tab:blue','tab:red','tab:blue'])
for d,c in zip(ds,cd): ax.text(d,c+0.03,f'{c:.4f}' if c>0 else '0',ha='center',fontsize=8)
ax.set_xlabel('dimension d'); ax.set_ylabel('curvature-excess coefficient c_d')
ax.set_title('Fig.D  c_d = d(d−1)/12: zero at d=1 (1D logic wave curvature-exact),\npositive only for d≥2 (distortion in coupled geometry)',fontsize=9)
ax.set_xticks(ds); plt.tight_layout(); plt.savefig('paper0_figD_cd.png',dpi=150); plt.close()

print("=== 検算照合（作図に使った厳密値 vs 仕様書）===")
print(f"図A 点間 π/2 = {np.pi/2:.7f}, 線分長 2πR = {2*np.pi*R:.7f}（修正: 円を1点で切って開く=全周。仕様書の3πは半周で12点が入らない不整合）")
print(f"図B 膨らみ高 h = {hmax:.7f} (仕様 0.1022225)  一致={abs(hmax-0.1022225)<1e-6}")
print(f"図C R=3: θ = {th3:.5f}° (仕様 91.62171), 面積 = {A3:.7f} (仕様 1.0189503) / R=1: θ = {th1:.5f}°, 面積 = {A1:.5f}")
print(f"図D c_d = {cd} (仕様 [0,1/6,1/2,1,5/3] for d=1..5)")
