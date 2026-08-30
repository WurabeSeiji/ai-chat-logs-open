#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""力学走行（データ収集のみ、図なし）。手作り親パッケージ（../手作り親_設計点走行_N4toN8_20260829/program/pass2_run.py）の
byte 同一コピーに 3 行の変更のみ：(1) 引数を親タグ（例 N6_eps0.60_k2, random_N5_s0。N はタグ中の N\d+ から取得）に (2) 親の読込先 data/<tag>/ (3) summary の experiment 文字列。
力学は N8 テンプレートと同一：振幅込み K・厳密線形回転 exp((2π/124)K)・種なし Z=v.copy()・正規化なし・直接読出し H⊥・STEPS=40000。"""
import os, math, csv, json, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import original_engine as eng

ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import re; TAG=sys.argv[1]; N=int(re.search(r"N(\d+)",TAG).group(1)); STEPS=40000; SEED=0; DELTA=0.0; L=124
DATA=os.path.join(ROOT,'data',TAG)

def adjacency_mask(s):
    m=s.m
    A=np.zeros((m,m),dtype=float)
    for i in range(m):
        share=(s.ea==s.ea[i])|(s.ea==s.eb[i])|(s.eb==s.ea[i])|(s.eb==s.eb[i])
        share[i]=False
        A[i,share]=1.0
    return A

def dense_K_amplitude_aware(s,z):
    A=adjacency_mask(s)
    pair=np.imag(np.conj(z)[:,None] * z[None,:])
    K=A*pair
    np.fill_diagonal(K,0.0)
    return K

def exp_linear_step(s,z):
    K=dense_K_amplitude_aware(s,z)
    asym=np.linalg.norm(K+K.T)
    if asym > 1e-10:
        raise RuntimeError(f'K antisymmetry failure: {asym}')
    H=1j*K
    w,V=np.linalg.eigh(H)
    return V @ (np.exp(-1j*(2*math.pi/L)*w)*(V.conj().T@z)), K

def metrics(s,Z,p,q,t,K=None):
    htot=float(np.vdot(Z,Z).real)
    Zperp=Z-p*(p@Z)-q*(q@Z)
    hperp=float(np.vdot(Zperp,Zperp).real)
    hpar=htot-hperp
    pr=eng.participation_ratio(Z)
    ztz=abs(complex(Z@Z))
    amp=np.abs(Z)
    if K is None:
        knorm=float('nan'); kasym=float('nan'); sig1=float('nan')
    else:
        knorm=float(np.linalg.norm(K,2)); kasym=float(np.linalg.norm(K+K.T))
        sig1=float(np.max(np.abs(np.linalg.eigvalsh(1j*K))))
    return [t,hpar,hperp,htot,pr,pr/s.m,ztz,sig1,float(amp.min()),float(amp.max()),float(amp.std()),knorm,kasym]

def run(v,p,q):
    s=eng.LowRankSystem(N)
    Z=v.copy()  # SEEDLESS
    rows=[]; states=np.empty((STEPS+1,s.m),dtype=np.complex128)
    Kprev=None; dphi=[]
    for t in range(STEPS+1):
        states[t]=Z
        rows.append(metrics(s,Z,p,q,t,Kprev))
        if t==STEPS: break
        Zn,Kprev=exp_linear_step(s,Z)
        dphi.append(float(np.angle(np.vdot(Z,Zn)))); Z=Zn
    a=np.asarray(rows,float); dphi=np.array(dphi+[np.nan])
    sig_next=np.r_[a[1:,7],np.nan]
    return np.c_[a,dphi,(2*math.pi/L)*sig_next],states

pz=np.load(os.path.join(DATA,'parent_v.npz'))
v=pz['v']; res=float(pz['residual']); sig=pz['sigma']
p=v.real/np.linalg.norm(v.real)
q=v.imag-(v.imag@p)*p; q/=np.linalg.norm(q)
treat,states_t=run(v,p,q)

headers=['step','H_parallel','H_perp','H_total','PR','PR_over_M','abs_ZT_Z','sigma1_of_step_generator','amp_min','amp_max','amp_std','K_spectral_norm','K_antisym_error','measured_phase_advance','angle_times_sigma1_next']
with open(os.path.join(DATA,'treatment_linear124_amplitude_aware_timeseries.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(headers); w.writerows(treat)
np.savez_compressed(os.path.join(DATA,'states_treatment.npz'),Z=states_t)

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
    return {'onset_Hperp_fraction_gt_0.05':frac_onset(a,0.05),
      'max_Hperp':float(np.nanmax(a[:,2])),'max_Hperp_step':int(a[np.nanargmax(a[:,2]),0]),
      'max_Hperp_fraction':float(np.nanmax(frac)),'max_Hperp_fraction_step':int(a[np.nanargmax(frac),0]),
      'final_Hparallel':float(a[-1,1]),'final_Hperp':float(a[-1,2]),'final_Htotal':float(a[-1,3]),
      'final_PR_over_M':float(a[-1,5]),'final_abs_ZT_Z':float(a[-1,6]),
      'final_amp_min':float(a[-1,8]),'final_amp_max':float(a[-1,9]),'final_amp_std':float(a[-1,10]),
      'Htotal_max_abs_drift':float(np.nanmax(np.abs(a[:,3]-a[0,3]))),
      'growth_fit':growth_fit(a)}
summary={'experiment':f'{TAG} class-weighted parent (design {str(pz["design"])}) linear124 seedless; dynamics identical to N8 template / handmade-parent pass2_run.py','N':N,'M':N*(N-1)//2,'L':L,'steps':STEPS,'seed':SEED,'delta':DELTA,
 'parent_residual':res,'parent_sigma':[float(x) for x in sig],'parent_design':str(pz['design']),'parent_r':float(pz['r']),
 'treatment':summarize(treat)}
with open(os.path.join(DATA,'summary.json'),'w') as f: json.dump(summary,f,indent=2,ensure_ascii=False)
key_steps=sorted(set([0,25,50,75,100,125,150,200,300,500,750,1000,1500,2000,3000,4000,5000,7500,10000,15000,20000,25000,30000,35000,40000]))
with open(os.path.join(DATA,'key_steps.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['step','treat_Hperp','treat_Hperp_frac','treat_PR_M','treat_amp_std'])
    for t in key_steps:
        a=treat[t]; w.writerow([t,a[2],a[2]/a[3],a[5],a[10]])
print(f"RUN {TAG} done: final H_perp/H={treat[-1,2]/treat[-1,3]:.3e} max frac={summary['treatment']['max_Hperp_fraction']:.3e}")
