#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""インフレ曲線＝線形不安定 の実測図（読出しのみ）。N依存も示す。
図: (左)|μ_max| と onset の N依存、(中)解析ランプ傾き2log10|μ| vs 実測 vs N、(右)N=6のH⊥/H実測ランプに解析傾き重畳。"""
import math, os
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
rv=lambda z:np.concatenate([z.real,z.imag]);cp=lambda v,M:v[:M]+1j*v[M:]
Ns=list(range(6,17)); amax=[];slope_pred=[];slope_meas=[];onset=[]
hp6=None
for N in Ns:
    A=adj(N);M=N*(N-1)//2
    S=np.asarray(np.load(GEN+f'/full_N3_N40_sweep/states/hm_N{N}_den_{N}_states_500.npz')['Z'],np.complex128)
    Z0=S[0];n=2*M;base=rv(one_step(Z0,A,1.0*N));J=np.zeros((n,n));eps=1e-8
    for k in range(n):
        v=rv(Z0).copy();v[k]+=eps;J[:,k]=(rv(one_step(cp(v,M),A,N))-base)/eps
    a=float(np.max(np.abs(np.linalg.eigvals(J))));amax.append(a);slope_pred.append(2*math.log10(a))
    p=Z0.real/np.linalg.norm(Z0.real);q=Z0.imag-np.dot(Z0.imag,p)*p;q/=np.linalg.norm(q)
    hpv=np.array([(lambda zp:np.vdot(zp,zp).real/np.vdot(z,z).real)(z-p*np.dot(p,z)-q*np.dot(q,z)) for z in S])
    lo=int(np.argmax(hpv>1e-18));hi=int(np.argmax(hpv>1e-8));slope_meas.append(float(np.polyfit(np.arange(lo,hi),np.log10(hpv[lo:hi]),1)[0]))
    onset.append(int(np.argmax(hpv>1e-3)))
    if N==6: hp6=hpv
fig,ax=plt.subplots(1,3,figsize=(16,5))
ax[0].plot(Ns,amax,'o-',label='|μ_max| (floor Jacobian)');ax[0].set_xlabel('N');ax[0].set_ylabel('|μ_max|');ax[0].axhline(1,ls=':',c='r')
ax2=ax[0].twinx();ax2.plot(Ns,onset,'s--',c='tab:orange',label='onset (measured)');ax2.set_ylabel('onset step')
ax[0].set_title('(a) instability rate |μ_max| ↓ and onset ↑ with N');ax[0].legend(loc='upper right',fontsize=8);ax2.legend(loc='center right',fontsize=8)
ax[1].plot(Ns,slope_pred,'o-',label='analytic 2log10|μ_max| (floor eigenvalue)')
ax[1].plot(Ns,slope_meas,'s-',label='measured ramp slope')
ax[1].set_xlabel('N');ax[1].set_ylabel('log10(H⊥/H) slope /step');ax[1].set_title('(b) analytic (floor) vs measured; gap=non-normal correction');ax[1].legend(fontsize=8)
st=np.arange(len(hp6));ax[2].semilogy(st,np.maximum(hp6,1e-33),label='N=6 measured H⊥/H')
lo=int(np.argmax(hp6>1e-25));x=np.arange(lo,lo+60);ax[2].semilogy(x,hp6[lo]*10**(slope_pred[0]*(x-lo)),'r--',label='analytic slope 2log10|μ_max|')
ax[2].set_xlabel('step');ax[2].set_ylabel('H⊥/H');ax[2].set_title('(c) N=6: exponential ramp, analytic slope overlaid');ax[2].legend(fontsize=8);ax[2].set_ylim(1e-33,3)
fig.tight_layout();fig.savefig('fig_inflation_linear_instability.png',dpi=120);print('wrote fig_inflation_linear_instability.png')
