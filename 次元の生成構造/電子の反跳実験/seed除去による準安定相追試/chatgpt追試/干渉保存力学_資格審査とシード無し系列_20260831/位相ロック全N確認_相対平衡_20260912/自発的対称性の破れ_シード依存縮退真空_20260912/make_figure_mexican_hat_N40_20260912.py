#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=40 メキシカンハット型SSB（実測データから再構成・模式図でない）。N=6版と同一手順。
床=full_N3_N40_sweep のZ0、秩序変数r(τ)=√(H⊥/H)は N40_2000step検証走行の2000step実測、
λ=床ヤコビアン最大固有値のln、r_vac=飽和実測、縮退真空5本は厳密one_stepでシード走行。"""
import math, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
HERE=os.path.dirname(os.path.abspath(__file__))
GEN=os.path.abspath(os.path.join(HERE,'..','..','make_parent型初期値_自己無撞着構造_20260907'))
PL=os.path.abspath(os.path.join(HERE,'..'))
def adj(N):
    ea,eb=np.triu_indices(N,1);M=len(ea);A=np.zeros((M,M))
    for e in range(M):
        s=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]);s[e]=False;A[e,s]=1.0
    return A
def one_step(z,A,den):
    u=np.exp(1j*np.angle(z));H=A*(np.conj(u)[:,None]*u[None,:]);np.fill_diagonal(H,0)
    H=1j*np.imag(H);w,V=np.linalg.eigh(H);ph=np.exp(-1j*(2*math.pi/den)*w);return V@(ph*(V.conj().T@z))
rv=lambda z:np.concatenate([z.real,z.imag]);cp=lambda v,M:v[:M]+1j*v[M:]
N=40;A=adj(N);M=N*(N-1)//2;den=N
# 秩序変数は 2000step 実測（本シリーズ走行）
S=np.asarray(np.load(PL+'/N40_2000step検証走行_20260912/results/hm_N40_den_40_states_2000.npz')['Z'],np.complex128)
Z0=S[0];p=Z0.real/np.linalg.norm(Z0.real);q=Z0.imag-np.dot(Z0.imag,p)*p;q/=np.linalg.norm(q)
Pperp=np.eye(M)-np.outer(p,p)-np.outer(q,q)
def perpfrac(z):
    a=np.dot(p,z);b=np.dot(q,z);zp=z-p*a-q*b;return np.vdot(zp,zp).real/np.vdot(z,z).real
hp=np.array([perpfrac(z) for z in S]);r_traj=np.sqrt(np.maximum(hp,0));r_vac=float(np.sqrt(np.median(hp[-100:])))
# λ 床ヤコビアン（M=780→2M=1560, ~数分）
base=rv(one_step(Z0,A,den));n=2*M;J=np.zeros((n,n));eps=1e-8
for k in range(n):
    v=rv(Z0).copy();v[k]+=eps;J[:,k]=(rv(one_step(cp(v,M),A,den))-base)/eps
lam=math.log(float(np.max(np.abs(np.linalg.eigvals(J)))))
def V(r): return -0.5*lam*r**2 + (lam/(4*r_vac**2))*r**4
# 縮退真空5本（厳密one_step, N=40はlock~1300step→1600step走行）
perps=[]
for seed in range(5):
    rng=np.random.default_rng(1000*N+seed)
    z=Z0+1e-8*(rng.standard_normal(M)+1j*rng.standard_normal(M));z=z/np.linalg.norm(z)*np.linalg.norm(Z0)
    for _ in range(1600): z=one_step(z,A,den)
    zp=Pperp@z;perps.append(np.concatenate([zp.real,zp.imag]))
Pm=np.array(perps);Pc=Pm-Pm.mean(0);U,s,Vt=np.linalg.svd(Pc,full_matrices=False);coords=Pc@Vt[:2].T;ang=np.arctan2(coords[:,1],coords[:,0])
fig=plt.figure(figsize=(14,6))
axa=fig.add_subplot(1,2,1)
rr=np.linspace(0,1.15*r_vac,300);axa.plot(rr,V(rr),'b-',lw=2,label='reconstructed V(r) (Landau, measured λ,r_vac)')
axa.plot([0],[V(0)],'r^',ms=12,label='floor (symmetric, unstable max)')
axa.plot([r_vac],[V(r_vac)],'go',ms=11,label='vacuum (broken, min)');axa.plot([-r_vac],[V(r_vac)],'go',ms=11)
tt=np.arange(len(r_traj));sc=axa.scatter(r_traj,V(r_traj),c=tt,cmap='viridis',s=8,zorder=5)
axa.set_xlabel('order parameter  r = √(H⊥/H)');axa.set_ylabel('effective potential V(r)')
axa.set_title(f'(a) N=40: rolls from unstable floor (r=0) to vacuum r_vac={r_vac:.3f}\nλ=ln|μ_max|={lam:.3f} (measured)')
axa.legend(fontsize=8);plt.colorbar(sc,ax=axa,label='step τ')
axb=fig.add_subplot(1,2,2,projection='3d')
R=np.linspace(0,1.15*r_vac,60);TH=np.linspace(0,2*np.pi,80);Rg,Tg=np.meshgrid(R,TH);X=Rg*np.cos(Tg);Y=Rg*np.sin(Tg);Z=V(Rg)
axb.plot_surface(X,Y,Z,cmap='coolwarm',alpha=.6,linewidth=0,antialiased=True)
axb.scatter([0],[0],[V(0)],c='r',s=60,label='floor (unstable apex)')
for a in ang: axb.scatter([r_vac*np.cos(a)],[r_vac*np.sin(a)],[V(r_vac)],c='k',s=45)
axb.scatter([r_vac*np.cos(ang[0])],[r_vac*np.sin(ang[0])],[V(r_vac)],c='k',s=45,label='seed-selected vacua (measured, brim)')
axb.set_xlabel('perp dir 1');axb.set_ylabel('perp dir 2');axb.set_zlabel('V')
axb.set_title('(b) Mexican-hat: symmetric unstable apex → degenerate vacua ring\n(seeds land at different brim angles = spontaneous choice)')
axb.legend(fontsize=8);axb.view_init(elev=28,azim=35)
fig.suptitle('N=40 Mexican-hat SSB reconstructed from measured order-parameter dynamics (not schematic)')
fig.tight_layout();fig.savefig('fig_mexican_hat_ssb_N40.png',dpi=120);print('wrote fig_mexican_hat_ssb_N40.png lam=%.3f r_vac=%.3f 真空角=%s'%(lam,r_vac,np.round(np.degrees(ang),1)))
