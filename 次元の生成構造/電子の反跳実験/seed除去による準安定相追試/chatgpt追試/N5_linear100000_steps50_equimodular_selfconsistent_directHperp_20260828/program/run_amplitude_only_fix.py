#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N=5 linear100000 seedless + make_parent amplitude-normalization removed experiment.
Baseline is copied from the previously verified linear124 experiment.
Based on the amplitude-aware linear124 version. Also removes the remaining explicit initial amplitude normalization.
The interaction-generator correction remains:
    K_ij = sin(theta_j-theta_i)
        -> Im(conj(Z_i) * Z_j)
for edge pairs sharing a vertex.
The exact exponential rotation exp((2*pi/124) K) is unchanged.
Starting from the prior seedless/all-three-fixes experiment, the only new change is removal of v/=norm(v) inside make_parent. Other make_parent logic is unchanged.
"""
import os, math, csv, json, sys, hashlib
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import original_engine as eng

ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA=os.path.join(ROOT,'data'); FIG=os.path.join(ROOT,'figures')
os.makedirs(DATA,exist_ok=True); os.makedirs(FIG,exist_ok=True)
N=5; STEPS=50; SEED=0; DELTA=0.0; L=100000


def adjacency_mask(s):
    m=s.m
    A=np.zeros((m,m),dtype=float)
    for i in range(m):
        share=(s.ea==s.ea[i])|(s.ea==s.eb[i])|(s.eb==s.ea[i])|(s.eb==s.eb[i])
        share[i]=False
        A[i,share]=1.0
    return A


def dense_K_phase_only(s,z):
    """Exact baseline generator from the verified linear124 experiment."""
    m=s.m; K=np.zeros((m,m)); th=np.angle(z)
    for i in range(m):
        share=(s.ea==s.ea[i])|(s.ea==s.eb[i])|(s.eb==s.ea[i])|(s.eb==s.eb[i])
        share[i]=False
        K[i,share]=np.sin(th[share]-th[i])
    return K


def dense_K_amplitude_aware(s,z):
    """ONE change: retain complex amplitudes in interaction coefficient.

    For connected edge-wave pair (i,j):
        K_ij = Im(conj(z_i) z_j)
             = |z_i||z_j| sin(theta_j-theta_i)
    This is real and antisymmetric by construction.
    """
    A=adjacency_mask(s)
    pair=np.imag(np.conj(z)[:,None] * z[None,:])
    K=A*pair
    np.fill_diagonal(K,0.0)
    return K


def exp_linear_step(s,z,mode):
    K = dense_K_phase_only(s,z) if mode=='baseline_phase_only' else dense_K_amplitude_aware(s,z)
    # diagnostic gate: generator must be real antisymmetric
    asym=np.linalg.norm(K+K.T)
    if asym > 1e-10:
        raise RuntimeError(f'K antisymmetry failure: {asym}')
    H=1j*K
    w,V=np.linalg.eigh(H)
    return V @ (np.exp(-1j*(2*math.pi/L)*w)*(V.conj().T@z)), K


def metrics(s,Z,p,q,t,K=None):
    htot=float(np.vdot(Z,Z).real)
    Zperp=Z-p*(p@Z)-q*(q@Z)              # 変更1: 親平面 (p,q) への直交成分を直接作る（差し引き htot-hpar は丸め床 1e-15 で埋もれる）
    hperp=float(np.vdot(Zperp,Zperp).real)
    hpar=htot-hperp
    pr=eng.participation_ratio(Z)
    ztz=abs(complex(Z@Z))
    # 判断2(A6-b): σ₁ は「この step を実際に回した生成子 K」の最大固有値（位相のみ読出し器 set_theta は使わない）
    amp=np.abs(Z)
    if K is None:
        knorm=float('nan'); kasym=float('nan'); sig1=float('nan')
    else:
        knorm=float(np.linalg.norm(K,2)); kasym=float(np.linalg.norm(K+K.T))
        sig1=float(np.max(np.abs(np.linalg.eigvalsh(1j*K))))
    return [t,hpar,hperp,htot,pr,pr/s.m,ztz,sig1,float(amp.min()),float(amp.max()),float(amp.std()),knorm,kasym]


def run(mode,v,g,p,q):
    s=eng.LowRankSystem(N)
    Z=v.copy()  # SEEDLESS: remove external perturbation DELTA*g; make_parent unchanged
    rows=[]; states=np.empty((STEPS+1,s.m),dtype=np.complex128)
    Kprev=None
    dphi=[]  # 判断2(A6-b): 実測位相進み arg<Z_t, Z_{t+1}>（ANGLE·σ₁ との検算用）
    for t in range(STEPS+1):
        states[t]=Z
        rows.append(metrics(s,Z,p,q,t,Kprev))
        if t==STEPS: break
        Zn,Kprev=exp_linear_step(s,Z,mode)
        dphi.append(float(np.angle(np.vdot(Z,Zn)))); Z=Zn
    a=np.asarray(rows,float); dphi=np.array(dphi+[np.nan])
    # 列追加：measured_phase_advance（step t→t+1）、angle_times_sigma1（同 step の K の σ₁ × 2π/L；K は直前 step のものなので 1 step ずらして対応）
    sig_next=np.r_[a[1:,7],np.nan]
    return np.c_[a,dphi,(2*math.pi/L)*sig_next],states

# EXACT SAME parent and seed construction as prior linear124 comparison
s0=eng.LowRankSystem(N)
rng=np.random.default_rng(40260721+1000*N+SEED)
v,res,sig=eng.make_parent(s0,rng)
g=None  # 判断5(S1): 外部 seed g は計算しない（zero_closure_kernel_seed の定義はエンジンに残す）
p=v.real/np.linalg.norm(v.real)
q=v.imag-(v.imag@p)*p; q/=np.linalg.norm(q)

# 変更2: baseline（位相のみ K）枝は実行しない
treat,states_t=run('treatment_amplitude_aware',v,g,p,q)

headers=['step','H_parallel','H_perp','H_total','PR','PR_over_M','abs_ZT_Z','sigma1_of_step_generator','amp_min','amp_max','amp_std','K_spectral_norm','K_antisym_error','measured_phase_advance','angle_times_sigma1_next']
for name,a in [('treatment_linear124_amplitude_aware',treat)]:
    with open(os.path.join(DATA,name+'_timeseries.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(headers); w.writerows(a)
np.savez_compressed(os.path.join(DATA,'states_treatment.npz'),Z=states_t)

# Reproduction gate against previously saved linear124 baseline
ref_path=os.path.join(ROOT,'data','previous_treatment_with_normalization.csv')
repro_max_abs=float('nan')

# summaries
def first_cross(a,col,threshold):
    ix=np.where(a[:,col]>threshold)[0]
    return int(ix[0]) if len(ix) else None

def frac_onset(a,thr=0.05):
    f=a[:,2]/a[:,3]; ix=np.where(f>thr)[0]
    return int(ix[0]) if len(ix) else None

def growth_fit(a,lo=1e-10,hi=1e-3):
    y=a[:,2]; mask=(y>lo)&(y<hi)&np.isfinite(y)
    if mask.sum()<3: return None
    x=a[mask,0]; ly=np.log(y[mask])
    slope,inter=np.polyfit(x,ly,1)
    pred=slope*x+inter
    ssr=((ly-pred)**2).sum(); sst=((ly-ly.mean())**2).sum()
    return {'slope_ln_Hperp_per_step':float(slope),'intercept':float(inter),'R2':float(1-ssr/sst),'n':int(mask.sum()),'step_min':int(x.min()),'step_max':int(x.max())}

def summarize(a):
    frac=a[:,2]/a[:,3]
    return {
      'onset_Hperp_fraction_gt_0.05': frac_onset(a,0.05),
      'max_Hperp': float(np.nanmax(a[:,2])),
      'max_Hperp_step': int(a[np.nanargmax(a[:,2]),0]),
      'max_Hperp_fraction': float(np.nanmax(frac)),
      'max_Hperp_fraction_step': int(a[np.nanargmax(frac),0]),
      'final_Hparallel':float(a[-1,1]),'final_Hperp':float(a[-1,2]),'final_Htotal':float(a[-1,3]),
      'final_PR_over_M':float(a[-1,5]),'final_abs_ZT_Z':float(a[-1,6]),
      'final_amp_min':float(a[-1,8]),'final_amp_max':float(a[-1,9]),'final_amp_std':float(a[-1,10]),
      'Htotal_max_abs_drift':float(np.nanmax(np.abs(a[:,3]-a[0,3]))),
      'growth_fit':growth_fit(a)
    }

summary={
 'experiment':'N5 linear124 seedless + make_parent v normalization removed; no external DELTA*g; no initialization Z normalization; linear exponential rotation; amplitude-aware K',
 'N':N,'M':10,'L':L,'steps':STEPS,'seed':SEED,'delta':DELTA,
 'parent_residual':float(res),'parent_sigma':[float(x) for x in sig],
 'fixes':['make_parent v/=norm(v) removed; all other make_parent logic unchanged','removed external perturbation seed DELTA*g: Z0=v','removed explicit Z/=norm(Z) at initialization','linear exponential rotation retained','amplitude-aware K retained: Im(conj(Z_i)*Z_j)','DECISION A5: zero_closure_generic normalization removed (unused)','DECISION A6(b): sigma1 read from the generator actually used in each branch; measured phase advance recorded','DECISION S1: zero_closure_kernel_seed not executed (definition kept)','DECISION R3(ii): validate_against_dense rewritten for amplitude-aware raw exp rotation (tool only)'],
  'treatment':summarize(treat),
}
with open(os.path.join(DATA,'summary.json'),'w') as f: json.dump(summary,f,indent=2,ensure_ascii=False)

# key-step comparison
key_steps=sorted(set([0,1,2,3,5,10,15,20,25,30,40,50]))
with open(os.path.join(DATA,'key_steps.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['step','treat_Hperp','treat_Hperp_frac','treat_PR_M','treat_amp_std'])
    for t in key_steps:
        a=treat[t]
        w.writerow([t,a[2],a[2]/a[3],a[5],a[10]])

# plots (match original H_perp log style and additional diagnostics)
import matplotlib.pyplot as plt
plt.figure(figsize=(10,6))
plt.semilogy(treat[:,0],np.maximum(treat[:,2],1e-40),label='treatment: amplitude-aware K (H_perp = |Z - P Z|^2, direct)')
plt.xlabel('step'); plt.ylabel('H_perp'); plt.title('N=5 linear100000: H_perp (single-change control)'); plt.legend(); plt.grid(True,which='both',alpha=.25); plt.tight_layout(); plt.savefig(os.path.join(FIG,'N5_Hperp_baseline_vs_amplitude_aware.png'),dpi=180); plt.close()

plt.figure(figsize=(10,6))
plt.plot(treat[:,0],treat[:,5],label='amplitude-aware')
plt.xlabel('step'); plt.ylabel('PR / M'); plt.title('N=5 linear100000: participation ratio'); plt.legend(); plt.grid(True,alpha=.25); plt.tight_layout(); plt.savefig(os.path.join(FIG,'N5_PR_baseline_vs_amplitude_aware.png'),dpi=180); plt.close()

plt.figure(figsize=(10,6))
plt.plot(treat[:,0],treat[:,10],label='amplitude-aware')
plt.xlabel('step'); plt.ylabel('std(|Z_m|)'); plt.title('N=5 linear100000: amplitude dispersion'); plt.legend(); plt.grid(True,alpha=.25); plt.tight_layout(); plt.savefig(os.path.join(FIG,'N5_amplitude_std_compare.png'),dpi=180); plt.close()

plt.figure(figsize=(10,6))
plt.semilogy(treat[:,0],np.maximum(treat[:,6],1e-30),label='amplitude-aware')
plt.xlabel('step'); plt.ylabel('|Z^T Z|'); plt.title('N=5 linear100000: squared-closure residual'); plt.legend(); plt.grid(True,which='both',alpha=.25); plt.tight_layout(); plt.savefig(os.path.join(FIG,'N5_closure_residual_compare.png'),dpi=180); plt.close()

print(json.dumps(summary,indent=2,ensure_ascii=False))
