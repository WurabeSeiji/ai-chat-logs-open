#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""単一波EOM＋90°必然性の実測図（読出しのみ・物理無変更）。
図1: sin(2Δφ)の零点=k·90°（第2高調波=振幅駆動の零点が90°格子の理由）＋床/ロックの個別項分布(実測)。
図2: EOM検証(実測 (exp(εK)z-z)/ε の振幅/位相成分 vs 導出式) の散布一致(N=6)。
入力: full_N3_N40_sweep の den=N npz（N=6床/ロック）。"""
import math, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.linalg import expm
GEN=os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','make_parent型初期値_自己無撞着構造_20260907'))
def adj(N):
    ea,eb=np.triu_indices(N,1);M=len(ea);A=np.zeros((M,M))
    for e in range(M):
        s=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]);s[e]=False;A[e,s]=1.0
    return A
def K_of(z,A):
    ph=np.angle(z);d=ph[None,:]-ph[:,None];return A*np.sin(d)
N=6;A=adj(N);M=N*(N-1)//2
S=np.asarray(np.load(GEN+f'/full_N3_N40_sweep/states/hm_N{N}_den_{N}_states_500.npz')['Z'],np.complex128)
Z0=S[0];Zf=S[-1]
# 図1
fig,ax=plt.subplots(1,2,figsize=(13,5))
psi=np.linspace(0,360,2000)
ax[0].plot(psi,np.sin(2*np.radians(psi)),lw=1.5)
for k in range(5): ax[0].axvline(90*k,ls=':',c='r',alpha=.6)
ax[0].axhline(0,c='gray',lw=.5)
ax[0].set_xlabel('phase difference Δφ [deg]');ax[0].set_ylabel('amplitude drive sin(2Δφ)')
ax[0].set_title('(1) amplitude drive = 2nd harmonic; zeros exactly at k·90°')
ax[0].set_xticks([0,90,180,270,360])
# 床とロックの個別項|sin(2Δφ)|分布（実測、隣接辺のみ）
def terms(z):
    ph=np.angle(z);d=ph[None,:]-ph[:,None];return np.abs((A*np.sin(2*d))[A>0])
ax[1].hist(terms(Z0),bins=30,alpha=.6,label='floor step0 (90° lattice): all ~0 (termwise)')
ax[1].hist(terms(Zf),bins=30,alpha=.6,label='locked (tilted vacuum): O(1), sum cancels')
ax[1].set_xlabel('|sin(2Δφ)| per adjacent edge pair (measured)');ax[1].set_ylabel('count')
ax[1].set_title('(1b) floor=termwise-zero (unique 90°) vs locked=cancellation');ax[1].legend(fontsize=8)
fig.tight_layout();fig.savefig('fig_90deg_necessity.png',dpi=120);print('wrote fig_90deg_necessity.png')
# 図2 EOM検証散布
eps=1e-6
def amp_eom(z):
    r=np.abs(z);ph=np.angle(z);d=ph[None,:]-ph[:,None];return 0.5*np.sum(A*r[None,:]*np.sin(2*d),1)
def ph_eom(z):
    r=np.abs(z);ph=np.angle(z);d=ph[None,:]-ph[:,None];return np.sum(A*r[None,:]*np.sin(d)**2,1)
LA=[];RA=[];LP=[];RP=[]
for t in np.linspace(1,len(S)-2,10).astype(int):
    z=S[t];zk=expm(eps*K_of(z,A))@z;r=np.abs(z)
    LA+=list((np.abs(zk)-r)/eps);RA+=list(amp_eom(z))
    LP+=list(r*np.angle(zk*np.conj(z))/eps);RP+=list(ph_eom(z))
fig,ax=plt.subplots(1,2,figsize=(12,5.5))
for a,(L,R,tt) in zip(ax,[(LA,RA,'amplitude  dr/dτ = ½Σ r_f sin(2Δφ)'),(LP,RP,'phase  r dφ/dτ = Σ r_f sin²Δφ')]):
    L=np.array(L);R=np.array(R);a.plot(R,L,'o',ms=3,alpha=.5)
    lo,hi=min(R.min(),L.min()),max(R.max(),L.max());a.plot([lo,hi],[lo,hi],'k--',lw=1)
    a.set_xlabel('derived RHS');a.set_ylabel('measured LHS (exp(εK)z)');a.set_title(tt)
fig.suptitle(f'(2) single-wave EOM verification N={N} (measured vs derived, y=x)')
fig.tight_layout();fig.savefig('fig_single_wave_eom_verify.png',dpi=120);print('wrote fig_single_wave_eom_verify.png')
