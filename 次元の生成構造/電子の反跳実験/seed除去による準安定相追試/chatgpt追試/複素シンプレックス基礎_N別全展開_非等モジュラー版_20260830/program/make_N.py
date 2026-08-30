#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""複素シンプレックス基礎（非等モジュラー版）：N=3..16 の全計算展開・図 3 枚・md 断片の生成。
等モジュラー版 make_N3.py / make_N4.py / make_Ngeneric.py と同じ量を同じ順に計算し、同じ様式の図と節を出す。
違いは状態だけ：state_provider.state(N)（クラス重み付き族 or 多様体上の代表点）。等モジュラー版の値も同時に計算して並記する。"""
import os, sys
import numpy as np
from fractions import Fraction
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,os.path.dirname(__file__))
from state_provider import *
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG=os.path.join(ROOT,'figures'); SEC=os.path.join(ROOT,'sections'); os.makedirs(SEC,exist_ok=True); os.makedirs(FIG,exist_ok=True)
N=int(sys.argv[1])
E=edges(N); M=len(E); A=adjacency(N)
z,kind,col,q,step=state(N); z0=equimodular(N)
cls=np.array([col[e] for e in E]); rbar=np.sqrt(RBAR2)
theta=np.angle(z); amp=np.abs(z); d2=z*z
def deg(x):
    d=np.degrees(x)%360.0
    return 0.0 if abs(d-360.0)<1e-6 else d
def analyze(v):
    K=K_of(N,v,A); hv=1j*(K@v); mu=(np.vdot(v,hv)/np.vdot(v,v)).real; res=float(np.linalg.norm(hv-mu*v)/np.linalg.norm(v))
    th=np.angle(v); a2=np.abs(v)**2; s2=[]; s22=[]
    for k2,(a,b) in enumerate(E):
        adj=[j for j,(c,d) in enumerate(E) if j!=k2 and ({a,b}&{c,d})]
        s2.append(sum(a2[j]*np.sin(th[j]-th[k2])**2 for j in adj)); s22.append(abs(sum(a2[j]*np.sin(2*(th[j]-th[k2])) for j in adj)))
    dd=v*v; glob=abs(dd.sum()); locs=[abs(sum(dd[k] for k,(a,b) in enumerate(E) if i in (a,b))) for i in range(N)]
    D=np.zeros((N,N),complex)
    for val,(i,j) in zip(dd,E): D[i,j]=D[j,i]=val
    J=np.eye(N)-np.ones((N,N))/N; Bm=-0.5*J@D@J
    sv=np.linalg.svd(Bm,compute_uv=False); rank=int((sv>1e-12).sum()); axes=np.sqrt(np.maximum(sv,0))
    Mre=np.block([[Bm.real,Bm.imag],[Bm.imag,-Bm.real]]); w,V=np.linalg.eigh(Mre)
    idx=np.argsort(w)[::-1][:N]; X=[]
    for j in idx:
        if w[j]<=1e-14: continue
        vv=V[:N,j]+1j*V[N:,j]; vv/=np.linalg.norm(vv); X.append(np.sqrt(w[j])*vv)
    X=np.array(X).T
    err=max(abs(((X[i]-X[j])@(X[i]-X[j]))-D[i,j]) for i,j in E)
    R=np.linalg.norm(X,axis=1); nulls=[abs(X[i]@X[i]) for i in range(N)]
    ax_pos=np.sort(np.sqrt(sv[sv>1e-12]))[::-1]; uax=[]
    for x in ax_pos:
        if not uax or abs(uax[-1][0]-x)>1e-8*ax_pos[0]: uax.append([x,1])
        else: uax[-1][1]+=1
    return dict(mu=float(mu),res=res,s2=s2,s22=s22,glob=glob,locs=locs,Bm=Bm,sv=sv,rank=rank,axes=axes,X=X,err=err,R=R,nulls=nulls,uax=[u[0] for u in uax],round=(len(uax)==1),d2real=bool(np.abs(dd.imag).max()<1e-15))
Z=analyze(z); Z0=analyze(z0)
def fmt(x,n=9): return f'{x:+.{n}f}'
def uniq(vals,tol=1e-8):
    out=[]
    for x in sorted(vals,reverse=True):
        if not out or abs(out[-1]-x)>tol*max(abs(out[0]),1e-300): out.append(x)
    return out
n_phase=len(uniq(list((np.degrees(theta)%360).round(9)),tol=1e-7)); n_amp=len(uniq(list(amp))); n_d2ang=len(uniq(list((np.degrees(np.angle(d2))%360).round(9)),tol=1e-7))
# ---------------- 図（等モジュラー版と同じ様式） ----------------
mcmap=plt.get_cmap('tab10')
def ccolor(c): return mcmap(c%10)
fig,ax=plt.subplots(figsize=(5.6,5.6))
th=np.linspace(0,2*np.pi,200); ax.plot(rbar*np.cos(th),rbar*np.sin(th),'--',color='gray',lw=0.8)
if kind=='class':
    for c in range(q):
        k=[kk for kk,e in enumerate(E) if col[e]==c][0]; cnt=int((cls==c).sum())
        ax.annotate('',xy=(z[k].real,z[k].imag),xytext=(0,0),arrowprops=dict(arrowstyle='->',color=ccolor(c),lw=2))
        ax.annotate(f'{step*c:.0f}° ×{cnt} |z|={amp[k]:.3f}',(z[k].real*1.25,z[k].imag*1.25),fontsize=7,color=ccolor(c),ha='center')
    ttl=f'N={N}: {M} waves in {q} classes, step {step:.1f} deg, class-weighted |z|'
else:
    for k,e in enumerate(E):
        ax.annotate('',xy=(z[k].real,z[k].imag),xytext=(0,0),arrowprops=dict(arrowstyle='->',color=ccolor(col[e]),lw=1.6,alpha=.85))
        if N<=4: ax.annotate(f'z{e} θ={deg(theta[k]):.1f}°',(z[k].real*1.15,z[k].imag*1.15),fontsize=8,color=ccolor(col[e]),ha='center')
    ttl=f'N={N}: {M} waves, generic self-consistent state (|z| and phases all differ)'
ax.set_xlim(-.4,.4); ax.set_ylim(-.4,.4); ax.set_aspect('equal'); ax.grid(alpha=.2)
ax.axhline(0,color='k',lw=.4); ax.axvline(0,color='k',lw=.4)
ax.set_xlabel('a = Re z'); ax.set_ylabel('b = Im z'); ax.set_title(ttl,fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(FIG,f'N{N}_fig1_z_complex_plane.png'),dpi=170); plt.close(fig)
fig,ax=plt.subplots(figsize=(5.6,5.6))
pos=0+0j; path=[0j]
for k in np.argsort(np.angle(d2)):
    ax.annotate('',xy=((pos+d2[k]).real,(pos+d2[k]).imag),xytext=(pos.real,pos.imag),arrowprops=dict(arrowstyle='->',color=ccolor(col[E[k]]),lw=1.6,alpha=.85))
    pos+=d2[k]; path.append(pos)
Pp=np.array(path); lim=1.12*max(abs(Pp.real).max(),abs(Pp.imag).max(),0.05)
ax.scatter([0],[0],color='k',s=40,zorder=3); ax.annotate('start = end',(0.02*lim,0.02*lim),fontsize=9)
ax.set_aspect('equal'); ax.grid(alpha=.2); ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
ax.set_xlabel('Re d²'); ax.set_ylabel('Im d²'); ax.set_title(f'N={N} closure: {M} d² vectors of unequal length close',fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG,f'N{N}_fig2_d2_closure.png'),dpi=170); plt.close(fig)
X=Z['X']; axes=Z['axes']; rank=Z['rank']
kshow=min(3,X.shape[1])
fig,axs=plt.subplots(1,kshow+1,figsize=(4.2*(kshow+1),4.4))
if kshow+1==1: axs=[axs]
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
ax.plot(kk2,axes,'o-',label='non-equimodular'); ax.plot(kk2,Z0['axes'],'s--',color='gray',label='equimodular version'); ax.legend(fontsize=7)
ax.set_xlabel('axis k'); ax.set_ylabel('scale √σ'); ax.set_title(f'axis spectrum (rank {rank})',fontsize=9); ax.grid(alpha=.25); ax.set_ylim(0,max(axes.max(),Z0['axes'].max())*1.2)
fig.suptitle(f'N={N} factual geometry (Takagi coordinates), non-equimodular state',fontsize=10)
fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(FIG,f'N{N}_fig3_geometry_axes.png'),dpi=170); plt.close(fig)
# ---------------- md 断片 ----------------
sec=N-1
diag=M-N; faces=N*(N-1)*(N-2)//6
mu_txt=f'{Fraction(Z["mu"]).limit_denominator(1000)}'
axes_txt='、'.join(f'{x:.9f}' for x in Z['uax']); axes0_txt='、'.join(f'{x:.9f}' for x in Z0['uax'])
mult='（全軸等長＝丸い）' if Z['round'] else f'（等長でない。相異なる値 {len(Z["uax"])} 種）'
mult0='丸い' if Z0['round'] else f'相異なる値 {len(Z0["uax"])} 種'
Rtxt='、'.join(f'{x:.9f}' for x in uniq(list(Z['R'])))
locmax=max(Z['locs']); loc_ok=locmax<1e-12
# 表
if kind=='class':
    rows=[]
    for c in range(q):
        k=[kk for kk,e in enumerate(E) if col[e]==c][0]; mem=[f'({e[0]},{e[1]})' for e in E if col[e]==c]
        rows.append(f"| {c} | {step*c:.6g}° | {amp[k]**2:.9f} | {amp[k]:.9f} | {fmt(z[k].real)} | {fmt(z[k].imag)} | {len(mem)} | {' '.join(mem)} |")
    wave_table='| クラス | θ | a_c = \\|z\\|² | \\|z\\| | a | b | 本数 | 所属する組 |\n|---|---|---|---|---|---|---|---|\n'+'\n'.join(rows)
    d2rows=[]
    for c in range(q):
        k=[kk for kk,e in enumerate(E) if col[e]==c][0]
        d2rows.append(f"| {c} | {(2*step*c)%360:.6g}° | {amp[k]**2:.9f} | {fmt(d2[k].real)}{fmt(d2[k].imag)}i |")
    d2_table='| クラス | d² の角度 | \\|d²\\| = a_c | d² の値 |\n|---|---|---|---|\n'+'\n'.join(d2rows)
else:
    rows=[]
    for k,e in enumerate(E):
        rows.append(f"| {e} | {col[e]}（{step*col[e]:.0f}°） | {deg(theta[k]):.3f}° | {amp[k]:.9f} | {fmt(z[k].real)} | {fmt(z[k].imag)} |")
    wave_table='| 組 | 等モジュラー版のクラス（θ） | θ | \\|z\\| | a | b |\n|---|---|---|---|---|---|\n'+'\n'.join(rows)
    d2rows=[]
    for k,e in enumerate(E):
        d2rows.append(f"| {e} | {deg(2*theta[k]):.3f}° | {amp[k]**2:.9f} | {fmt(d2[k].real)}{fmt(d2[k].imag)}i |")
    d2_table='| 組 | d² の角度 | \\|d²\\| = \\|z\\|² | d² の値 |\n|---|---|---|---|\n'+'\n'.join(d2rows)
# 構成の説明
if kind=='class':
    ac=[amp[[kk for kk,e in enumerate(E) if col[e]==c][0]]**2 for c in range(q)]
    omega_sum=abs(sum(ac[c]*np.exp(2j*np.pi*c/q) for c in range(q)))
    if N%2==0:
        design=(f'**構成（偶数 N の設計＝完全マッチング分解、位相は等モジュラー版と同一）**。{N} 人で全員に相手がいるペアの組ませ方（完全マッチング、1 組 {N//2} 本）が、'
                f'ちょうど {N-1} セットで全 {M} 組を重複なく覆う。セット（色）c = 0..{N-2} の波に位相 θ = c×{step:.6g}°（= c×180°/{N-1}）を与える。各頂点はどの色の波もちょうど 1 本ずつ持つ。')
        adjn=2; locreason=(f'各頂点は {N-1} 色を 1 本ずつ持つ。d² は色 c で a_c·ω^c（ω = e^{{i·360°/{N-1}}}）なので局所和 = Σ_c a_c ω^c = 0（§{sec}.2 の条件そのもの）。')
        globreason=f'全体では各色が {N//2} 本ずつ → {N//2}×Σ_c a_c ω^c = 0。'
    else:
        nn=q
        design=(f'**構成（奇数 N の設計＝距離クラス分解、位相は等モジュラー版と同一）**。頂点を円周に並べ、組 (i,j) の「距離」を d = |i−j| を N で折り返した値（1〜{nn}）とする。'
                f'距離 d のクラスは各頂点にちょうど 2 本ずつ属する（{N} 本、全 {nn} クラスで {M} 組を覆う）。クラス d の波に位相 θ = (d−1)×{step:.6g}°（= (d−1)×180°/{nn}）を与える。')
        adjn=4; locreason=(f'各頂点は {nn} クラスを 2 本ずつ持つ。d² はクラス c で a_c·ω^c（ω = e^{{i·360°/{nn}}}）なので局所和 = 2Σ_c a_c ω^c = 0（§{sec}.2 の条件そのもの）。')
        globreason=f'全体では各クラスが {N} 本ずつ → {N}×Σ_c a_c ω^c = 0。'
    amp_design=(f'\n**振幅（等モジュラー版との唯一の違い）**。クラス c の振幅二乗を a_c = r̄²(1 + 0.6·cos(4πc/{q}))、r̄² = 1/15 とする（§0.2 の選択 (A)）。'
                f'値は {"、".join(f"a_{c} = {ac[c]:.9f}" for c in range(q))}（平均は r̄² = 1/15 = 0.066666667、‖v‖² = {M}r̄² = {Fraction(M,15)} で等モジュラー版と同じ）。'
                f'\n\n**この振幅が自己無撞着を保つ理由（導出）**。組 e（クラス c）から見て他の各クラス c′ に隣接組が {adjn} 本、同じクラスの隣接組は位相差 0 で寄与しない。一般和則（§0.1）に代入し、実部条件と虚部条件を 1 本の複素条件にまとめると'
                f' ω^{{−c}}·Σ_{{c′}} a_{{c′}} ω^{{c′}} = Σ_{{c′}} a_{{c′}} + μ/{adjn//2}（全ての c で同じ実数）。q = {q} ≥ 3 でこれが全 c について実になるのは **Σ_c a_c ω^c = 0** のときだけで、そのとき μ = −{adjn//2}·Σ_c a_c = −(N−1)·r̄²。'
                f'a_c = r̄²(1+0.6cos(4πc/q)) は三角和の直交性で Σ_c a_c ω^c = 0 を厳密に満たす（実測 |Σ a_c ω^c| = {omega_sum:.0e}）。つまり**必要なのは「振幅二乗を重みとした 1 の {q} 乗根の重心が原点」であって、振幅が等しいことではない**。')
    uniq_note='※この族は q−3 = '+str(q-3)+' 個の実パラメータを持つ連続族で、等モジュラー版はその 1 点（ε=0）である。'
else:
    if N==3:
        design=('**構成（多様体上の代表点）**。N=3 では 3 本の波が互いに隣接する。等モジュラー版の状態（位相 0/60/120°、大きさ共通）から、自己無撞着写像 F(v) = i·K(v)v − μ(v)v のヤコビアン核（スケール・全体回転を除いて 2 次元）の固定方向へ 0.3‖v‖ 動かし、Gauss–Newton で F = 0 に戻した状態（§0.2 の選択 (B)）。'
                '大きさも位相差も 3 本で全て異なる。')
    elif N==4:
        design=('**構成（多様体上の代表点）**。等モジュラー版の状態（3 ペアに 0/60/120°）から、自己無撞着写像のヤコビアン核（スケール・全体回転を除いて 5 次元）の固定方向へ 0.3‖v‖ 動かし、Gauss–Newton で F = 0 に戻した状態（§0.2 の選択 (B)）。'
                'N=4 では局所閉塞を保ったまま非等モジュラーにする道がない（{S_i=0, W_i 均等} の解集合はスケール分の 1 次元しかない）ので、この状態は**局所閉塞を破る**（§'+str(sec)+'.4）。')
    else:
        design=(f'**構成（多様体上の代表点）**。クラス重み付き族（§0.2 (A)）は q = {q} ≤ 3 のため点に退化する（Σ_c a_c ω^c = 0 が a_c 全部等しいを強制）。そこで等モジュラー版の状態から、ヤコビアン核の固定方向へ 0.3‖v‖ 動かし、'
                f'**局所閉塞 S_i = 0 と頂点重み W_i = Σ_k|z_ik|² の均等**を拘束にした Newton で戻した状態（§0.2 の選択 (B)。§0.1 の定理「局所閉塞＋頂点重み均等 ⇒ 自己無撞着、μ = −W」により自動的に自己無撞着）。')
    amp_design=''
    uniq_note='※この状態は等モジュラー点を通る自己無撞着解多様体（次元は §0.2 参照）上の 1 点で、一意性はない。'
special=''
if Z['d2real']: special='\n**特記（この N だけ）**：d² が**全て実数**。\n'
if Z0['d2real'] and not Z['d2real']: special=f'\n**特記**：等モジュラー版では位相が 0°/90° で d² が全て実数（実の擬ユークリッド空間 R^{{2,2}} に置けた）が、この状態では位相が {n_phase} 種類に散り、d² は複素になる（max|Im d²| = {np.abs(d2.imag).max():.3f}）。**実埋め込みは等モジュラー点だけの性質**である。\n'
if N in (9,15): special+=f'\n**特記**：N={N} では一部の距離クラスが複数の輪に分かれるが、必要なのは「各頂点に各クラスが 2 本ずつ」だけなので結論は変わらない（残差 {Z["res"]:.0e}）。\n'
if loc_ok: loc_line=f'- **局所閉塞【{"定理・" if kind=="class" else "選択＋"}実測】**：{locreason if kind=="class" else "拘束に入れた分枝を選んだので成立。"}実測 max|局所和| = {locmax:.0e}。よって**全頂点が複素ヌル錐上**（埋め込みで max|x·x| = {max(Z["nulls"]):.0e}）。'
else: loc_line=f'- **局所閉塞【実測】：破れる**。頂点ごとの |局所和|/r̄² = {"、".join(f"{x*15:.3f}" for x in Z["locs"])}（等モジュラー版は全頂点 0）。頂点は複素ヌル錐上に**ない**（|x·x|/r̄² = {"、".join(f"{x*15:.3f}" for x in Z["nulls"])}）。'
s2=Z['s2']; s22=Z['s22']
frag=f'''## {sec}. N={N} の全展開

### {sec}.0 N={N} の複素シンプレックスの定義

頂点 {N} 個 {{0,…,{N-1}}}。2 点の組は M = {N}×{N-1}/2 = {M} 個。各組に波 z が 1 個（実数 {2*M} 個が状態の全て）。2 乗距離は d² = z²。**N={N} の複素シンプレックスとは、{N} 頂点とこの {M} 個の複素 2 乗距離の組**。

### {sec}.1 数え上げ（N={N}）

| 量 | 値 | 計算 |
|---|---|---|
| 頂点 | {N} | — |
| 波の総数 M | {M} | {N}×{N-1}/2 |
| 辺 | {N} | 円周の隣どうしの {N} 本 |
| 対角線 | {diag} | M − N = {M} − {N} |
| 面（三角形） | {faces} | C({N},3) |
| 位相の種類 | {n_phase}{f'（θ = k×{step:.6g}°）' if kind=='class' else '（全組で異なる）' if n_phase==M else ''} | §{sec}.2 の構成（等モジュラー版：{q}） |
| 大きさ \\|z\\| の種類 | {n_amp} | §{sec}.2（等モジュラー版：1） |
| d² の角度の種類 | {n_d2ang} | 位相の 2 倍 |
| rank | {Z['rank']}{' = N−1' if Z['rank']==N-1 else ''} | §{sec}.4 の計算結果（等モジュラー版：{Z0['rank']}） |

### {sec}.2 波の値 (a, b) と構成、自己無撞着の検証

{design}{amp_design}
スケールは規約 〈|z|²〉 = r̄² = 1/15（§0.2。‖v‖² = {M}r̄² = {Fraction(M,15)}、等モジュラー版と同じ）。{'同じクラスの波は同じ (a, b) を持つので、クラスごとに列挙する' if kind=='class' else '全組を列挙する'}：

{wave_table}

**自己無撞着の検証【定理・実測】**：一般和則（§0.1）を全 {M} 組で数値検査——実部条件 Σ_f |z_f|² sin²φ の値は全組で {np.mean(s2):.9f}（ばらつき {max(abs(x-np.mean(s2)) for x in s2):.0e}。= −μ）、虚部条件 |Σ_f |z_f|² sin2φ| ≤ {max(s22):.0e}。方程式の残差 {Z['res']:.1e}。よって **μ = {Z['mu']:+.9f}**{f' = −(N−1)·r̄² = −{Fraction(N-1,15)}' if abs(Z['mu']+(N-1)/15)<1e-9 else ''}——**等モジュラー版と同じ値**（{Z0['mu']:+.9f}）。
{uniq_note}
{special}
### {sec}.3 2 乗距離{'（クラスごと）' if kind=='class' else '（組ごと）'}

{d2_table}

{'**角度は等モジュラー版と同じ、大きさ |d²| = a_c がクラスごとに違う。**' if kind=='class' else '**角度も大きさも組ごとに違う。**'}

### {sec}.4 閉塞の判定・B 行列と形

- **大域閉塞【定理・実測】**：{globreason if kind=='class' else '自己無撞着（μ≠0）から従う（§0.2 の定理）。'}実測 |Σz²| = {Z['glob']:.0e}。
{loc_line}
- **B = −½JD²J の Takagi 軸【実測・機械精度】**：軸スケールの相異なる値は {axes_txt} {mult}。**等モジュラー版**：{axes0_txt}（{mult0}）。rank = {Z['rank']}。頂点のエルミート半径：{Rtxt}。埋め込みの距離再現誤差 max = {Z['err']:.0e}。

### {sec}.5 N={N} のまとめ（種別）

- 【選択＋実測】{'クラス重み付き（振幅二乗 a_c = r̄²(1+0.6cos(4πc/q))）の状態は自己無撞着' if kind=='class' else '多様体上の代表点は自己無撞着'}（残差 {Z['res']:.0e}）。等モジュラーは**使っていない**。
- 【定理・実測】大域閉塞は成立。局所閉塞は{'成立（全頂点が複素ヌル錐上）' if loc_ok else '**破れる**'}。
- 【定理】μ = {Z['mu']:+.9f}{f' = −(N−1)·〈|z|²〉' if abs(Z['mu']+(N-1)/15)<1e-9 else ''}：等モジュラー版と同一（μ は平均振幅二乗だけで決まる）。
- 【実測（機械精度）】rank = {Z['rank']}{'（最大次元）' if Z['rank']==N-1 else ''}、軸は{'全て等長（丸い）' if Z['round'] else '等長でない'}（等モジュラー版：{mult0}）。**等モジュラー版との差は形（Takagi 軸）だけ**{'' if loc_ok else '、および局所閉塞の破れ'}。
- 【規約の帰結】大きさに制限はない（スケール対称性）。

### 図（N={N}）

![N={N} 波の複素平面](figures/N{N}_fig1_z_complex_plane.png)

![N={N} d² の閉塞](figures/N{N}_fig2_d2_closure.png)

![N={N} 幾何（Takagi 座標）](figures/N{N}_fig3_geometry_axes.png)
'''
open(os.path.join(SEC,f'N{N}.md'),'w').write(frag)
# 一覧表用
import json
json.dump(dict(N=N,M=M,q=q,kind=kind,n_phase=n_phase,n_amp=n_amp,n_d2ang=n_d2ang,rank=Z['rank'],round=Z['round'],round0=Z0['round'],loc_ok=bool(loc_ok),mu=Z['mu'],res=Z['res'],axes=[float(x) for x in Z['uax']],axes0=[float(x) for x in Z0['uax']],d2real=Z['d2real'],d2real0=Z0['d2real'],amp_min=float(amp.min()),amp_max=float(amp.max())),open(os.path.join(SEC,f'N{N}.json'),'w'))
print(f"N={N} {kind}: 残差={Z['res']:.1e} μ={Z['mu']:+.9f} 位相{n_phase}種 大きさ{n_amp}種 rank={Z['rank']} 丸い={Z['round']}(等モジュラー {Z0['round']}) 局所閉塞={'成立' if loc_ok else '破れ'} d²実={Z['d2real']}")
