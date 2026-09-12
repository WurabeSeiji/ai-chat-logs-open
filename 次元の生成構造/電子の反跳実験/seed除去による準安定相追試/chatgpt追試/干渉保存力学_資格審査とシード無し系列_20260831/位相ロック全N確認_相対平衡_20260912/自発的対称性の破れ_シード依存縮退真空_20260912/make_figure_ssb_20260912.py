#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SSB本体の実データ図（模式図でなく実測）。
図左: N=6 同一床＋異シード5本→各シードのロック複素平面(実測z_f)を重ね描き＝同型(等半径)だが位相配置は相異(縮退真空)。
図中: 各シードのΔ(1step回転)が全て-2π/Nに一致(実測)。
図右: 全N6-40バッチ(ssb_batch_summary.csv)のN依存＝真空差(中央)vs N。
入力: full_N3_N40_sweep N=6床, 厳密one_step; ssb_batch_summary.csv。"""
import math, os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
GEN=os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','make_parent型初期値_自己無撞着構造_20260907'))
def adj(N):
    ea,eb=np.triu_indices(N,1);M=len(ea);A=np.zeros((M,M))
    for e in range(M):
        s=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]);s[e]=False;A[e,s]=1.0
    return A
def one_step(z,A,den):
    u=np.exp(1j*np.angle(z));H=A*(np.conj(u)[:,None]*u[None,:]);np.fill_diagonal(H,0)
    H=1j*np.imag(H);w,V=np.linalg.eigh(H);ph=np.exp(-1j*(2*math.pi/den)*w);return V@(ph*(V.conj().T@z))
N=6;A=adj(N);M=N*(N-1)//2;den=N;STEPS=800;ideal=-2*math.pi/N
Z0=np.asarray(np.load(GEN+f'/full_N3_N40_sweep/states/hm_N{N}_den_{N}_states_500.npz')['Z'][0],np.complex128)
finals=[];deltas=[];ratios=[]
for seed in range(5):
    rng=np.random.default_rng(1000*N+seed)
    z=Z0+1e-8*(rng.standard_normal(M)+1j*rng.standard_normal(M));z=z/np.linalg.norm(z)*np.linalg.norm(Z0)
    for _ in range(STEPS): z=one_step(z,A,den)
    # 剛体回転除去(wave0基準)で真空の"形"を比較
    w=z*np.conj(z[0]/abs(z[0]));finals.append(w)
    z1=one_step(z,A,den);deltas.append(float(np.angle(np.vdot(z,z1))));r=np.abs(z);ratios.append(r.max()/r.min())
fig,ax=plt.subplots(1,3,figsize=(16,5.2))
cols=plt.cm.tab10(np.arange(5))
for i,w in enumerate(finals):
    ax[0].plot(w.real,w.imag,'o',ms=5,color=cols[i],alpha=.8,label=f'seed{i}')
th=np.linspace(0,2*np.pi,200);r0=np.abs(finals[0]).mean();ax[0].plot(r0*np.cos(th),r0*np.sin(th),'k:',lw=.8)
ax[0].set_aspect('equal');ax[0].set_xlabel('Re w (rigid-rot removed)');ax[0].set_ylabel('Im w')
ax[0].set_title('(a) N=6 locked vacua: equal radius, DIFFERENT phase config per seed');ax[0].legend(fontsize=7)
ax[1].plot(range(5),np.array(deltas),'o',ms=9,label='Δ per seed (measured)')
ax[1].axhline(ideal,ls='--',c='k',label=f'-2π/N={ideal:.5f}')
ax[1].set_xlabel('seed');ax[1].set_ylabel('1-step rotation Δ [rad]');ax[1].set_ylim(ideal-0.01,ideal+0.01)
ax[1].set_title(f'(b) all seeds lock to same Δ=-2π/N (amp ratio {min(ratios):.4f}-{max(ratios):.4f})');ax[1].legend(fontsize=8)
rows=list(csv.DictReader(open('全N6_40_SSBバッチ_20260912/ssb_batch_summary.csv')))
Nn=[int(r['N']) for r in rows];vac=[float(r['vacuum_spread_med_deg']) for r in rows]
ax[2].plot(Nn,vac,'o-');ax[2].set_xlabel('N');ax[2].set_ylabel('vacuum phase-config spread (median deg)')
ax[2].set_title('(c) all N=6-40: SSB 35/35, seed-dependent vacuum spread vs N')
fig.suptitle('Spontaneous symmetry breaking (measured, exact one_step): same floor + tiny seeds -> degenerate vacua')
fig.tight_layout();fig.savefig('fig_ssb_degenerate_vacua.png',dpi=120);print('wrote fig_ssb_degenerate_vacua.png')
