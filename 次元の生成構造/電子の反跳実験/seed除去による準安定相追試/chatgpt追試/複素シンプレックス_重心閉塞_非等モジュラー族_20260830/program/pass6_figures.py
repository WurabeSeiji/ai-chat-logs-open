#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス6：図化（fig8: 乱数均衡親の走行、fig9: 全親の ρ−1 帯）（読出し専用）。`--static` で走行に依存しない図のみ。
fig1: N=6 の d² ベクトル連鎖（等モジュラー vs ε=0.6）——長さが違っても重心が原点＝閉じる
fig2: N=3..16 の q と族次元 q−3（丸い N を強調）
fig3: 族の 1 パラメータ路に沿う Takagi 軸（N=6,8,9,11）
fig4: ヤコビアン核次元（非自明）と N²−4N−1、クラス一様族次元
fig5: 走行 H⊥/H の時間発展（N6/N8/N9、ε 別）＋ |z| の std/mean
fig6: 予測 λ_f vs 実測 λ
fig7: 最終状態の |z_e| 分布（初期との比較）"""
import os, csv, json, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0,os.path.dirname(__file__))
from common import *
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=os.path.join(ROOT,'figures')
STATIC='--static' in sys.argv
cm=plt.get_cmap('tab10'); R2=1.0/15.0
# fig1
fig,axs=plt.subplots(1,3,figsize=(14,4.6))
for ax,(eps,title) in zip(axs[:2],[(0.0,'equimodular: 15 d² vectors, 3 per class'),(0.6,'class-weighted (k=2, ε=0.6): lengths differ')]):
    th,cls,q=phases(6); a=R2*cos_mode(q,2,eps); v=state(6,a); d2=v*v
    order=np.argsort(cls,kind='stable'); pos=0j; path=[0j]
    for e in order:
        ax.annotate('',xy=(pos.real+d2[e].real,pos.imag+d2[e].imag),xytext=(pos.real,pos.imag),arrowprops=dict(arrowstyle='->',color=cm(cls[e]),lw=1.6)); pos+=d2[e]; path.append(pos)
    P=np.array(path); m=1.15*max(abs(P.real).max(),abs(P.imag).max()); ax.set_xlim(-m,m); ax.set_ylim(-m,m)
    ax.scatter([0],[0],color='k',s=40,zorder=3); ax.annotate(f'start = end, |Σd²| = {abs(d2.sum()):.0e}',(0.003,0.003),fontsize=8)
    ax.set_aspect('equal'); ax.set_title(f'N=6 {title}',fontsize=9); ax.set_xlabel('Re d²'); ax.set_ylabel('Im d²'); ax.grid(alpha=.3)
ax=axs[2]; th,cls,q=phases(6); a=R2*cos_mode(q,2,0.6); v=state(6,a); d2=v*v; E=edges(6)
for c in range(q):
    idx=[e for e in range(len(E)) if cls[e]==c]; ax.scatter(d2[idx].real,d2[idx].imag,color=cm(c),s=40,label=f'class {c}: a_c={a[c]:.4f}')
for c in range(q): ax.plot([0,R2*np.cos(2*np.pi*c/q)],[0,R2*np.sin(2*np.pi*c/q)],color=cm(c),lw=.8,ls='--')
ax.scatter([0],[0],color='k',s=30); ax.set_aspect('equal'); ax.legend(fontsize=7); ax.set_title('d²=a_c ω^c: weighted roots of unity, Σ a_c ω^c = 0',fontsize=9); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig1_centroid_closure_N6.png'),dpi=150); plt.close()
# fig2
T=list(csv.DictReader(open(os.path.join(ROOT,'results','family_table.csv'))))
Ns=[int(r['N']) for r in T]; qs=[int(r['q']) for r in T]; dims=[int(r['dim_measured']) for r in T]; rnd=[r['eq_round']=='True' for r in T]
fig,ax=plt.subplots(figsize=(9,4)); x=np.array(Ns)
ax.bar(x-0.2,qs,0.4,label='q (phase classes)',color='#9ecae1'); ax.bar(x+0.2,dims,0.4,label='family dim (q−3, scale removed)',color='#fd8d3c')
for xi,r_,d in zip(x,rnd,dims):
    if r_: ax.annotate('round',(xi,max(d,0)+0.4),ha='center',fontsize=8,color='green')
ax.set_xticks(Ns); ax.set_xlabel('N'); ax.set_ylabel('count'); ax.legend(); ax.set_title('phase classes q and dimension of the class-weighted self-consistent family; "round" = all Takagi axes equal (N=4,5,7)',fontsize=9); ax.grid(axis='y',alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig2_family_dimension.png'),dpi=150); plt.close()
# fig3
fig,axs=plt.subplots(1,4,figsize=(16,3.8))
for ax,N in zip(axs,[6,8,9,11]):
    p=np.genfromtxt(os.path.join(ROOT,'data',f'family_path_N{N}.csv'),delimiter=',',names=True)
    for i in range(N-1): ax.plot(p['eps'],p[f'axis{i+1}'],lw=1.4)
    ax.axvline(0,color='k',ls=':',lw=.8); ax.set_title(f'N={N} (q={phases(N)[2]}): Takagi axes along a_c=1+ε cos(4πc/q)',fontsize=9); ax.set_xlabel('ε'); ax.set_ylabel('axis scale √σ_k'); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig3_takagi_along_family.png'),dpi=150); plt.close()
# fig4
Jt=list(csv.DictReader(open(os.path.join(ROOT,'results','jacobian_nullity.csv'))))
Ns=np.array([int(r['N']) for r in Jt]); nn=np.array([int(r['nontrivial']) for r in Jt]); cf=np.array([int(r['class_family_dim']) for r in Jt])
fig,ax=plt.subplots(figsize=(8,4)); ax.plot(Ns,nn,'o-',label='Jacobian nullity − 2 (measured)'); ax.plot(Ns,Ns**2-4*Ns-1,'--',color='gray',label='N²−4N−1'); ax.plot(Ns,cf,'s-',color='#fd8d3c',label='class-weighted family dim (q−3)')
ax.set_xlabel('N'); ax.set_ylabel('dimension'); ax.set_xticks(Ns); ax.legend(); ax.grid(alpha=.3); ax.set_title('kernel of DF at the equimodular point vs. the class-weighted family',fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig4_jacobian_nullity.png'),dpi=150); plt.close()
print("static figures done")
if STATIC: raise SystemExit
pred={r['tag']:r for r in csv.DictReader(open(os.path.join(ROOT,'results','parents_predictions.csv')))}
def load(t):
    f=os.path.join(ROOT,'data',t,'treatment_linear124_amplitude_aware_timeseries.csv'); f=f if os.path.exists(f) else f+'.gz'   # gzip 済みでも読める
    return np.genfromtxt(f,delimiter=',',names=True) if os.path.exists(f) else None
# fig5
for N in [6,8,9]:
    tags=[t for t in pred if t.startswith(f'N{N}_')]
    fig,axs=plt.subplots(1,2,figsize=(13,4.3))
    for i,t in enumerate(tags):
        a=load(t)
        if a is None: continue
        lab=f"ε={pred[t]['eps']} k={pred[t]['k']} (pred λ={float(pred[t]['pred_lambda_f']):.5f})"
        axs[0].semilogy(a['step'],np.maximum(a['H_perp']/a['H_total'],1e-40),lw=1.2,color=cm(i),label=lab)
        axs[1].plot(a['step'],a['amp_std']/np.sqrt(a['H_total']/len(a['step'])*0+a['H_total'])*np.sqrt(N*(N-1)/2),lw=1.2,color=cm(i),label=lab)
    axs[0].set_xlabel('step'); axs[0].set_ylabel('H⊥/H'); axs[0].set_title(f'N={N}: deviation from parent (log)',fontsize=9); axs[0].legend(fontsize=7); axs[0].grid(alpha=.3)
    axs[1].set_xlabel('step'); axs[1].set_ylabel('std|z| / rms|z|'); axs[1].set_title(f'N={N}: amplitude spread (0 = equimodular)',fontsize=9); axs[1].legend(fontsize=7); axs[1].grid(alpha=.3)
    plt.tight_layout(); plt.savefig(os.path.join(FIG,f'fig5_dynamics_N{N}.png'),dpi=150); plt.close()
# fig6, fig7
D=list(csv.DictReader(open(os.path.join(ROOT,'results','dynamics_summary.csv'))))
fig,ax=plt.subplots(figsize=(5,5))
for r in D:
    if r['lambda_measured']: ax.scatter(float(r['pred_lambda_f']),float(r['lambda_measured']),s=50); ax.annotate(r['tag'],(float(r['pred_lambda_f']),float(r['lambda_measured'])),fontsize=7,xytext=(4,4),textcoords='offset points')
lim=[0,0.021]; ax.plot(lim,lim,'--',color='gray'); ax.set_xlim(lim); ax.set_ylim(lim); ax.set_xlabel('predicted λ_f (co-rotating monodromy, before run)'); ax.set_ylabel('measured λ (ln H⊥ slope)'); ax.grid(alpha=.3); ax.set_title('growth rate: prediction vs measurement',fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig6_lambda_pred_vs_meas.png'),dpi=150); plt.close()
fig,axs=plt.subplots(1,3,figsize=(15,4))
for ax,N in zip(axs,[6,8,9]):
    tags=[t for t in pred if t.startswith(f'N{N}_')]
    for i,t in enumerate(tags):
        f=os.path.join(ROOT,'data',t,'states_treatment.npz')
        if not os.path.exists(f): continue
        Z=np.load(f)['Z']; z0=np.abs(Z[0]); zT=np.abs(Z[-1]); idx=np.argsort(z0)
        ax.plot(np.sort(zT)[::-1],'o-',ms=3,color=cm(i),label=f"ε={pred[t]['eps']} k={pred[t]['k']} final"); ax.plot(np.sort(z0)[::-1],':',color=cm(i),lw=1)
    ax.set_xlabel('edge (sorted by |z|)'); ax.set_ylabel('|z_e|'); ax.set_title(f'N={N}: |z_e| at step 40000 (solid) vs step 0 (dotted)',fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig7_final_amplitudes.png'),dpi=150); plt.close()
# fig8: 乱数均衡親（対称性なし）の走行
rp=os.path.join(ROOT,'results','balanced_random_parents.csv')
if os.path.exists(rp):
    RB=list(csv.DictReader(open(rp)))
    fig,axs=plt.subplots(1,4,figsize=(17,3.8))
    for ax,N in zip(axs,[5,6,7,8]):
        for i,r in enumerate([r for r in RB if int(r['N'])==N]):
            a=load(r['tag'])
            if a is None: continue
            ax.semilogy(a['step'],np.maximum(a['H_perp']/a['H_total'],1e-40),lw=1.1,color=cm(i),label=f"s{r['seed']} pred λ={float(r['lambda_f']):.4f}")
        # 対称親（等モジュラー）の参照
        ref={5:None,6:'N6_eps0.00_k2',7:None,8:'N8_eps0.00_k2'}[N]
        if ref: a=load(ref); ax.semilogy(a['step'],np.maximum(a['H_perp']/a['H_total'],1e-40),'k--',lw=1,label='symmetric equimodular parent')
        ax.set_title(f'N={N}: random balanced parents (S_i=0, W_i=W), no symmetry',fontsize=9); ax.set_xlabel('step'); ax.set_ylabel('H⊥/H'); ax.legend(fontsize=6); ax.grid(alpha=.3)
    plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig8_random_balanced_dynamics.png'),dpi=150); plt.close()
    # fig9: ρ−1 の帯
    MC=list(csv.DictReader(open(os.path.join(ROOT,'results','monodromy_calibration.csv'))))
    fig,ax=plt.subplots(figsize=(8,4))
    ax.scatter([int(r['N']) for r in MC],[float(r['rho_minus_1']) for r in MC],marker='s',s=60,color='k',label='symmetric equimodular parent (1-factor / distance classes)')
    ax.scatter([int(r['N'])+0.15 for r in pred.values() if r['eps'] not in ('0.0','0.00')],[float(r['pred_rho'])-1 for r in pred.values() if r['eps'] not in ('0.0','0.00')],marker='^',s=40,color='#fd8d3c',label='class-weighted family member')
    ax.scatter([int(r['N'])-0.15 for r in RB],[float(r['rho_minus_1']) for r in RB],marker='o',s=30,color='#3182bd',label='random balanced parent (no symmetry)')
    ax.axhspan(1e-9,1e-3,color='green',alpha=.08); ax.text(4.1,3e-4,'neutral band (floor in 40000 steps)',fontsize=8,color='green')
    ax.set_yscale('log'); ax.set_ylim(1e-9,1e-1); ax.set_xlabel('N'); ax.set_ylabel('ρ − 1 (co-rotating monodromy)'); ax.legend(fontsize=7,loc='lower right'); ax.grid(alpha=.3); ax.set_title('instability of self-consistent parents: symmetric vs generic',fontsize=9)
    plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig9_rho_bands.png'),dpi=150); plt.close()
print("PASS6 OK")
