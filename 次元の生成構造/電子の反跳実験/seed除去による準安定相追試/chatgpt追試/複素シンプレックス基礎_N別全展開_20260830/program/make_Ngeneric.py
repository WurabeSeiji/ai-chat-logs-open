#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""複素シンプレックス基礎：N=5..16 の全計算展開・図 3 枚・md 断片の生成。
設計＝偶数 N: 円法の完全マッチング分解（色 c に θ=c·180°/(N−1)）／奇数 N: 距離クラス分解（距離 d に θ=(d−1)·180°/n, n=(N−1)/2）。
スケール規約 r²=1/15（全 N 共通、§0.2）。"""
import os, sys
import numpy as np
from fractions import Fraction
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG=os.path.join(ROOT,'figures'); SEC=os.path.join(ROOT,'sections')
os.makedirs(SEC,exist_ok=True)
N=int(sys.argv[1])
E=[(i,j) for i in range(N) for j in range(i+1,N)]; M=len(E)
r=1.0/np.sqrt(15.0)
if N%2==0:
    n=N-1; col={}
    for rr in range(n):
        col[tuple(sorted((rr,N-1)))]=rr
        for k in range(1,N//2): col[tuple(sorted(((rr-k)%n,(rr+k)%n)))]=rr
    ncls=n; step=180.0/(N-1); kind='even'
else:
    nn=(N-1)//2; col={}
    for d in range(1,nn+1):
        for i in range(N): col[tuple(sorted((i,(i+d)%N)))]=d-1
    ncls=nn; step=180.0/nn; kind='odd'
theta=np.array([np.radians(step*col[e]) for e in E]); z=r*np.exp(1j*theta); d2=z*z
A=np.zeros((M,M))
for a in range(M):
    for b in range(M):
        if a!=b and set(E[a])&set(E[b]): A[a,b]=1.0
K=A*np.imag(np.conj(z)[:,None]*z[None,:]); hv=1j*(K@z); mu=(np.vdot(z,hv)/np.vdot(z,z)).real
res=float(np.linalg.norm(hv-mu*z)/np.linalg.norm(z))
s2=[]; s22=[]
for k2,(a,b) in enumerate(E):
    adj=[j for j,(c,d) in enumerate(E) if j!=k2 and ({a,b}&{c,d})]
    s2.append(sum(np.sin(theta[j]-theta[k2])**2 for j in adj)); s22.append(abs(sum(np.sin(2*(theta[j]-theta[k2])) for j in adj)))
glob=abs(d2.sum()); loc=max(abs(sum(d2[k] for k,(a,b) in enumerate(E) if i in (a,b))) for i in range(N))
D=np.zeros((N,N),complex)
for val,(i,j) in zip(d2,E): D[i,j]=D[j,i]=val
J=np.eye(N)-np.ones((N,N))/N; Bm=-0.5*J@D@J
sv=np.linalg.svd(Bm,compute_uv=False); rank=int((sv>1e-12).sum()); axes=np.sqrt(np.maximum(sv,0))
Mre=np.block([[Bm.real,Bm.imag],[Bm.imag,-Bm.real]]); w,V=np.linalg.eigh(Mre)
idx=np.argsort(w)[::-1][:N]; X=[]
for j in idx:
    if w[j]<=1e-14: continue
    vv=V[:N,j]+1j*V[N:,j]; vv/=np.linalg.norm(vv); X.append(np.sqrt(w[j])*vv)
X=np.array(X).T
err=max(abs(((X[i]-X[j])@(X[i]-X[j]))-D[i,j]) for i,j in E)
R=np.linalg.norm(X,axis=1); nulldev=max(abs(X[i]@X[i]) for i in range(N))
ax_pos=np.sort(np.sqrt(sv[sv>1e-12]))[::-1]
uax=[]; 
for x in ax_pos:
    if not uax or abs(uax[-1][0]-x)>1e-8*ax_pos[0]: uax.append([x,1])
    else: uax[-1][1]+=1
round_flag=(len(uax)==1)
uax=[u[0] for u in uax]
d2real=bool(np.abs(d2.imag).max()<1e-15)
# ---- 図
mcmap=plt.get_cmap('tab10')
def ccolor(c): return mcmap(c%10)
fig,ax=plt.subplots(figsize=(5.6,5.6))
th=np.linspace(0,2*np.pi,200); ax.plot(r*np.cos(th),r*np.sin(th),'--',color='gray',lw=0.8)
for c in range(ncls):
    zz=r*np.exp(1j*np.radians(step*c)); cnt=sum(1 for e in E if col[e]==c)
    ax.annotate('',xy=(zz.real,zz.imag),xytext=(0,0),arrowprops=dict(arrowstyle='->',color=ccolor(c),lw=2))
    ax.annotate(f'{step*c:.0f}° ×{cnt}',(zz.real*1.22,zz.imag*1.22),fontsize=8,color=ccolor(c),ha='center')
ax.set_xlim(-.36,.36); ax.set_ylim(-.36,.36); ax.set_aspect('equal'); ax.grid(alpha=.2)
ax.axhline(0,color='k',lw=.4); ax.axvline(0,color='k',lw=.4)
ax.set_xlabel('a = Re z'); ax.set_ylabel('b = Im z'); ax.set_title(f'N={N}: {M} waves in {ncls} classes, step {step:.1f} deg',fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG,f'N{N}_fig1_z_complex_plane.png'),dpi=170); plt.close(fig)
fig,ax=plt.subplots(figsize=(5.6,5.6))
pos=0+0j
for k in np.argsort(np.angle(d2)):
    ax.annotate('',xy=((pos+d2[k]).real,(pos+d2[k]).imag),xytext=(pos.real,pos.imag),arrowprops=dict(arrowstyle='->',color=ccolor(col[E[k]]),lw=1.6,alpha=.85))
    pos+=d2[k]
ax.scatter([0],[0],color='k',s=40,zorder=3); ax.annotate('start = end',(0.004,0.004),fontsize=9)
ax.set_aspect('equal'); ax.grid(alpha=.2); ax.autoscale()
ax.set_xlabel('Re d²'); ax.set_ylabel('Im d²'); ax.set_title(f'N={N} closure: {M} d² vectors close',fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG,f'N{N}_fig2_d2_closure.png'),dpi=170); plt.close(fig)
kshow=min(3,X.shape[1])
fig,axs=plt.subplots(1,kshow+1,figsize=(4.2*(kshow+1),4.4))
lim=1.15*float(np.abs(X[:,:kshow]).max())
for kk in range(kshow):
    ax=axs[kk]; c1=X[:,kk]
    for k,(i,j) in enumerate(E):
        ax.plot([c1[i].real,c1[j].real],[c1[i].imag,c1[j].imag],color=ccolor(col[E[k]]),lw=1.0,alpha=.55)
    ax.scatter(c1.real,c1.imag,color='k',s=30,zorder=3)
    for i in range(N): ax.annotate(str(i),(c1[i].real,c1[i].imag),textcoords='offset points',xytext=(5,5),fontsize=8)
    ax.set_aspect('equal'); ax.grid(alpha=.2); ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
    ax.set_title(f'axis {kk+1}: r={axes[kk]:.4f}',fontsize=9)
ax=axs[kshow]; kk2=np.arange(1,len(axes)+1)
ax.plot(kk2,axes,'o-'); ax.set_xlabel('axis k'); ax.set_ylabel('scale √σ'); ax.set_title(f'axis spectrum (rank {rank})',fontsize=9); ax.grid(alpha=.25); ax.set_ylim(0,max(axes)*1.2)
fig.suptitle(f'N={N} factual geometry (Takagi coordinates)',fontsize=10)
fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(FIG,f'N{N}_fig3_geometry_axes.png'),dpi=170); plt.close(fig)
# ---- md 断片
sec=N-1  # N=5 → ## 4
edge_note=f'円周の隣どうしの {N} 本' ; diag=M-N
if kind=='even':
    design=(f'**構成（偶数 N の設計＝完全マッチング分解）**。{N} 人で全員に相手がいるペアの組ませ方（完全マッチング、1 組 {N//2} 本）が、'
            f'ちょうど {N-1} セットで全 {M} 組を重複なく覆う（総当たり戦の {N-1} 節）。セット（色）c = 0..{N-2} の波に位相 θ = c×{step:.6g}°'
            f'（= c×180°/{N-1}）を与える。各頂点はどの色の波もちょうど 1 本ずつ持つ。')
    locreason=(f'各頂点は {N-1} 色を 1 本ずつ持つ。d² の角度は色 c で c×{2*step:.6g}° = c×360°/{N-1}、つまり **1 の {N-1} 乗根がちょうど 1 個ずつ**。'
               f'冪根の総和は零なので局所和 = r²×({N-1} 乗根の和) = 0。')
    globreason=f'全体では各角度が {N//2} 本ずつ → {N//2}×r²×（{N-1} 乗根の和）= 0。'
else:
    nn=(N-1)//2
    design=(f'**構成（奇数 N の設計＝距離クラス分解）**。頂点を円周に並べ、組 (i,j) の「距離」を d = |i−j| を N で折り返した値（1〜{nn}）とする。'
            f'距離 d のクラスは各頂点にちょうど 2 本ずつ属する（{N} 本、全 {nn} クラスで {M} 組を覆う）。クラス d の波に位相 θ = (d−1)×{step:.6g}°（= (d−1)×180°/{nn}）を与える。')
    locreason=(f'各頂点は {nn} クラスを 2 本ずつ持つ。d² の角度は (d−1)×{2*step:.6g}° = (d−1)×360°/{nn}、つまり **1 の {nn} 乗根が 2 個ずつ**。'
               f'冪根の総和は零なので局所和 = 2r²×({nn} 乗根の和) = 0。'+('' if nn>1 else '（※ 乗根が 1 種類しかない場合は和が零にならない——それが N=3。）'))
    globreason=f'全体では各角度が {N} 本ずつ → {N}×r²×（{nn} 乗根の和）= 0。'
cls_rows=[]
for c in range(ncls):
    ph=step*c; a_=r*np.cos(np.radians(ph)); b_=r*np.sin(np.radians(ph))
    mem=[str(e) for e in E if col[e]==c]; 
    cls_rows.append(f"| {c} | {ph:.6g}° | {a_:+.9f} | {b_:+.9f} | {len(mem)} | {' '.join(mem)} |")
d2_rows=[]
for c in range(ncls):
    ph2=(2*step*c)%360; val=d2[[k for k,e in enumerate(E) if col[e]==c][0]]
    d2_rows.append(f"| {c} | {ph2:.6g}° | {val.real:+.9f}{val.imag:+.9f}i |")
axes_txt='、'.join(f'{x:.9f}' for x in uax)
mult='（全軸等長＝丸い）' if round_flag else f'（等長でない。相異なる値 {len(uax)} 種：長い軸あり）'
Rs=np.sort(R)[::-1]; uR=[]
for x in Rs:
    if not uR or abs(uR[-1]-x)>1e-8*Rs[0]: uR.append(x)
Rtxt='、'.join(f'{x:.9f}' for x in uR)
mu_frac=Fraction(N-1,15)
special=''
if d2real: special='\n**特記（この N だけ）**：位相が 0°/90° のため d² が**全て実数**（+r² と −r²）。この配置は実の擬ユークリッド空間に埋め込める（符号数 (2,2)）。詳しくは本節末尾の補足。\n'
if N in (9,15): special+=f'\n**特記**：N={N} では一部の距離クラスが 1 つの輪でなく複数の輪に分かれる（距離が N と公約数を持つため）。設計の成立に必要なのは「各頂点に各クラスが 2 本ずつ」だけなので、結論は変わらない（残差 {res:.0e} で検証済み）。\n'
frag=f'''## {sec}. N={N} の全展開

### {sec}.0 N={N} の複素シンプレックスの定義

頂点 {N} 個 {{0,…,{N-1}}}。2 点の組は M = {N}×{N-1}/2 = {M} 個。各組に波 z が 1 個（実数 {2*M} 個が状態の全て）。2 乗距離は d² = z²。**N={N} の複素シンプレックスとは、{N} 頂点とこの {M} 個の複素 2 乗距離の組**。

### {sec}.1 数え上げ（N={N}）

| 量 | 値 | 計算 |
|---|---|---|
| 頂点 | {N} | — |
| 波の総数 M | {M} | {N}×{N-1}/2 |
| 辺 | {N} | {edge_note} |
| 対角線 | {diag} | M − N = {M} − {N} |
| 面（三角形） | {N*(N-1)*(N-2)//6} | C({N},3) |
| 位相の種類 | {ncls}（θ = k×{step:.6g}°） | §{sec}.2 の構成 |
| d² の角度の種類 | {ncls}（k×{2*step:.6g}°） | 位相の 2 倍 |
| rank | {rank}{' = N−1' if rank==N-1 else ''} | §{sec}.4 の計算結果 |

### {sec}.2 波の値 (a, b) と構成、自己無撞着の検証

{design}
スケールは規約 r² = 1/15（§0.2。‖v‖² = {M}r² = {Fraction(M,15)}）。同じクラスの波は同じ (a, b) を持つので、クラスごとに列挙する：

| クラス | θ | a = r·cosθ | b = r·sinθ | 本数 | 所属する組 |
|---|---|---|---|---|---|
{chr(10).join(cls_rows)}

**自己無撞着の検証【定理・実測】**：和則（§0.1）を全 {M} 組で数値検査——実部条件 Σsin²φ の値は全組で {np.mean(s2):.6f}（= N−1 = {N-1}、ばらつき {max(abs(x-(N-1)) for x in s2):.0e}）、虚部条件 |Σsin2φ| ≤ {max(s22):.0e}。方程式の残差 {res:.1e}。よって μ/r² = −(N−1) = −{N-1}、**μ = −{mu_frac} = {mu:+.9f}**。
※この構成の一意性は証明していない【未証明】。
{special}
### {sec}.3 2 乗距離（クラスごと）

| クラス | d² の角度 | d² の値 |
|---|---|---|
{chr(10).join(d2_rows)}

### {sec}.4 閉塞の判定・B 行列と形

- **大域閉塞【定理・実測】**：{globreason}実測 |Σz²| = {glob:.0e}。
- **局所閉塞【定理・実測】**：{locreason}実測 max|局所和| = {loc:.0e}。よって**全頂点が光円錐上**（埋め込みで max|x·x| = {nulldev:.0e}）。
- **B = −½JD²J の Takagi 軸【実測・機械精度】**：軸スケールの相異なる値は {axes_txt} {mult}。rank = {rank}。頂点のエルミート半径：{Rtxt}。埋め込みの距離再現誤差 max = {err:.0e}。

### {sec}.5 N={N} のまとめ（種別）

- 【選択＋実測】{'完全マッチング分解' if kind=='even' else '距離クラス分解'}に位相を等間隔で配った等モジュラー状態は自己無撞着（残差 {res:.0e}）。
- 【定理・実測】大域閉塞・局所閉塞とも成立。全頂点が光円錐上。
- 【定理】μ/r² = −(N−1) = −{N-1}。
- 【実測（機械精度）】rank = {rank}{'（最大次元）' if rank==N-1 else ''}、軸は{'全て等長（丸い）' if round_flag else '等長でない（1 本だけ長い軸を持つ）'}。
- 【規約の帰結】大きさに制限はない（スケール対称性）。

### 図（N={N}）

![N={N} 波の複素平面](figures/N{N}_fig1_z_complex_plane.png)
![N={N} d² の閉塞](figures/N{N}_fig2_d2_closure.png)
![N={N} 幾何（Takagi 座標）](figures/N{N}_fig3_geometry_axes.png)
'''
open(os.path.join(SEC,f'N{N}.md'),'w').write(frag)
print(f"N={N}: 残差={res:.1e} μ={mu:+.6f} rank={rank} 丸い={round_flag} d²実={d2real} 断片と図を出力")
