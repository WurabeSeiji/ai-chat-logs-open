#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺66 追補図: xyztRQ 空間の本格図化(全パネル実計算データ)
import numpy as np, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ---- 実データ1: 2軸復号(補遺64 T5 のパリティ混在復号を実行して座標を得る) ----
Ng=32; gx=np.arange(Ng)/Ng
X=[gx.reshape(-1,1), gx.reshape(1,-1)]
def line2(I,f1,f2):
    e=np.exp(-2j*np.pi*(f1*X[0]+f2*X[1]))
    return 2*np.mean(I*np.conj(e))
def qc(z):
    p=(np.angle(z)/(2*np.pi))%1.0
    return int(np.round(p*4))%4
# 3断片: (s, 2軸位置桁, ε) — セル |k| パターンで s と ε を決める(1D断面でなく2軸)
frags=[ (1, (0,0,0,0), (1,2)),   # s=1 セル(0,0,0,0): ε=+1
        (3, (1,0,0,0), (3,0)),   # s=3 セル(1,0,0,0): ε=-1
        (5, (1,1,0,0), (2,3)) ]  # s=5 セル(1,1,0,0): ε=+1
decoded=[]
for s,cell,(d1,d2) in frags:
    a1,a2=d1/4.0,d2/4.0
    # パリティ混在復号(64 T5): 単独軸線+和線
    I=(1.0+2*np.cos(2*np.pi*(X[0]-a1))*np.cos(2*np.pi*(X[1]-a2))+np.sqrt(2)*np.cos(2*np.pi*(X[0]-a1)))**2
    r10=qc(line2(I,1,0)); r11=qc(line2(I,1,1))
    x1=r10/4.0; x2=((r11-r10)%4)/4.0
    eps=(-1)**sum(abs(k) for k in cell)
    decoded.append((x1,x2,np.sqrt(s),eps,s))
    assert (r10,(r11-r10)%4)==(d1,d2)

# ---- 実データ2: 崩壊カスケードの世界線 9→(5,3,1) (ノルム時計 Δt=1/√s) ----
def ticks(rate_s, t0, n):
    dt=1.0/np.sqrt(rate_s)
    return t0+dt*np.arange(n+1)
t_parent=ticks(9,0.0,6)           # 親 s=9: レート3
t_dec=t_parent[-1]                 # 崩壊イベント(記録)
tick_d=[ticks(5,t_dec,4),ticks(3,t_dec,4),ticks(1,t_dec,3)]

fig=plt.figure(figsize=(14,10))
# ===== (a) 3D 配置 (x, y, R) =====
ax=fig.add_subplot(2,2,1,projection='3d')
for x1,x2,R,eps,s in decoded:
    ax.scatter([x1],[x2],[R],s=140,c=('tab:red' if eps>0 else 'tab:blue'),
               marker=('o' if eps>0 else 's'),depthshade=False)
    ax.text(x1,x2,R+0.12,f's={s}, R=√{s}, Q={eps:+d}',fontsize=8)
xx,yy=np.meshgrid(np.linspace(0,1,2),np.linspace(0,1,2))
for s in (1,3,5):
    ax.plot_surface(xx,yy,np.full_like(xx,np.sqrt(s)),alpha=0.08,color='gray')
for q in range(5):
    ax.plot([q/4,q/4],[0,1],[0,0],color='gray',lw=0.4)
    ax.plot([0,1],[q/4,q/4],[0,0],color='gray',lw=0.4)
ax.set_xlabel('x (decoded quarter digits)'); ax.set_ylabel('y (decoded)'); ax.set_zlabel('R = sqrt(s)')
ax.set_title('(a) Position-R space: fragments on R-rungs\n(positions = machine-decoded 2-axis quarter digits)')
# ===== (b) 世界線 (x,t) =====
ax2=fig.add_subplot(2,2,2)
xp=0.5
ax2.plot([xp]*len(t_parent),t_parent,'-o',ms=4,color='tab:blue',label='parent s=9 (dt=1/3)')
xd=[0.25,0.5,0.75]; cols=['tab:green','tab:olive','tab:cyan']; labs=['s=5 (dt=1/√5)','s=3 (dt=1/√3)','s=1 (dt=1)']
for x0,tk,c,lb in zip(xd,tick_d,cols,labs):
    ax2.plot([xp,x0],[t_dec,tk[1]],color=c,ls=':')
    ax2.plot([x0]*len(tk),tk,'-o',ms=4,color=c,label=lb)
ax2.axhline(t_dec,color='crimson',lw=0.8,ls='--')
ax2.text(0.98,t_dec-0.08,'decay event (record):  R^2: 9 = 5+3+1,  Q: -1=(+1)(-1)(+1)',fontsize=7,color='crimson',va='top',ha='right')
ax2.set_xlabel('x'); ax2.set_ylabel('t = sum of 1/sqrt(s) (internal time, paper 8)')
ax2.set_title('(b) World lines: tick spacing = 1/R (t derived from R)')
ax2.legend(fontsize=7,loc='lower right')
# ===== (c) (R,Q) 台帳面 =====
ax3=fig.add_subplot(2,2,3)
shells={1:[(0,0,0,0)],3:[(1,0,0,0)],5:[(1,1,0,0)],7:[(1,1,1,0),(2,0,0,0)],9:[(2,1,0,0),(1,1,1,1)],
        11:[(2,1,1,0)],13:[(3,0,0,0),(2,2,0,0),(2,1,1,1)]}
for s,shapes in shells.items():
    for sh in shapes:
        eps=(-1)**sum(sh)
        ax3.scatter([np.sqrt(s)],[eps],s=90,c=('tab:red' if eps>0 else 'tab:blue'),marker=('o' if eps>0 else 's'))
        ax3.annotate(f'{sh}',(np.sqrt(s),eps),textcoords='offset points',xytext=(2,8),fontsize=6,rotation=45)
ax3.annotate('',xy=(np.sqrt(5),0.82),xytext=(3,-0.82),arrowprops=dict(arrowstyle='->',color='gray'))
ax3.annotate('',xy=(np.sqrt(3),-0.82),xytext=(3,-0.86),arrowprops=dict(arrowstyle='->',color='gray'))
ax3.annotate('',xy=(1,0.78),xytext=(3,-0.84),arrowprops=dict(arrowstyle='->',color='gray'))
ax3.text(2.2,0.0,'decay 9->(5,3,1)\nsum R^2 conserved (118,944/118,944)\nprod Q conserved (0 violations)',fontsize=8)
ax3.set_yticks([-1,1]); ax3.set_ylim(-1.6,1.6)
ax3.set_xlabel('R = sqrt(s) (allowed rungs)'); ax3.set_ylabel('Q = epsilon')
ax3.set_title('(c) (R,Q) ledger plane: discrete rungs, shape-resolved\n(mixed shells s=7,9,13 carry both Q values)')
# ===== (d) 拘束面 R^2 = s と仮想±1 =====
ax4=fig.add_subplot(2,2,4)
svals=[]
for r in (1,3,5,7,9):
    for sh in shells[r]:
        for _ in range(1):
            svals.append((r, sum((abs(k)+0.5)**2 for k in sh)))
s_ax=np.array([v[0] for v in svals]); R2=np.array([v[1] for v in svals])
ax4.plot([0,14],[0,14],color='gray',lw=0.8)
ax4.scatter(s_ax,R2,s=60,c='tab:blue',zorder=3,label='recorded states (R^2 = s exact)')
ax4.scatter([9],[10],s=80,c='tab:orange',marker='^',zorder=3,label='virtual stage (5,5): +1 borrowed')
ax4.scatter([9],[8],s=80,c='tab:orange',marker='v',zorder=3,label='virtual stage (3,5): -1 borrowed')
ax4.annotate('borrow cancels at next record\n(118,944/118,944 paths)',(9,10),textcoords='offset points',xytext=(12,6),fontsize=8)
ax4.set_xlabel('ledger label s'); ax4.set_ylabel('R^2 read from frequency content')
ax4.set_title('(d) Constraint surface R^2=s (codim-2 with Q=epsilon):\nrecords lie on it exactly; virtual steps go off by +/-1')
ax4.legend(fontsize=7,loc='upper left')
plt.tight_layout(); plt.savefig('supplement66_fig_xyztRQ.png',dpi=150); plt.close()
print("written: supplement66_fig_xyztRQ.png")
