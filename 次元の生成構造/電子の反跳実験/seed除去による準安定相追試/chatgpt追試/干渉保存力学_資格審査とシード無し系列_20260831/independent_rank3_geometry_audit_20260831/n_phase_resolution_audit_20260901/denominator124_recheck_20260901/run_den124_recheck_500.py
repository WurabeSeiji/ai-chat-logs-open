#!/usr/bin/env python3
import os, math, csv, time
import numpy as np

IN_DIR='/mnt/data/hm_rank3_inputs'
OUT_DIR='/mnt/data/denominator124_recheck_20260901'
os.makedirs(OUT_DIR, exist_ok=True)
STEPS=500
DEN=124

def edges(N): return np.triu_indices(N,k=1)
def adjacency(N):
    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),float)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e])
        share[e]=False; A[e,share]=1.0
    return A

def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H

def step(z,A,d):
    H=H_of(z,A); w,V=np.linalg.eigh(H)
    return V@(np.exp(-1j*d*w)*(V.conj().T@z))

def plane(v):
    p=v.real/np.linalg.norm(v.real)
    q=v.imag-(v.imag@p)*p; q/=np.linalg.norm(q)
    return p,q

def metrics(z,p,q):
    h=float(np.vdot(z,z).real)
    zp=z-p*(p@z)-q*(q@z)
    hp=float(np.vdot(zp,zp).real)
    return h,hp/h,float(abs(z@z)/h)

rows=[]; sums=[]
delta=2*math.pi/DEN
for N in range(3,17):
    Z0=np.load(os.path.join(IN_DIR,f'hm_N{N}_states_treatment.npz'))['Z'][0].astype(np.complex128)
    A=adjacency(N); p,q=plane(Z0); z=Z0.copy()
    states=np.empty((STEPS+1,len(z)),np.complex128); hpv=[]; rr=[]; t0=time.time()
    for t in range(STEPS+1):
        states[t]=z; h,hp,clo=metrics(z,p,q); hpv.append(hp)
        r=dict(N=N,denominator=DEN,delta=delta,step=t,H_total=h,Hperp_frac=hp,global_closure=clo)
        rr.append(r); rows.append(r)
        if t<STEPS: z=step(z,A,delta)
    hpv=np.asarray(hpv); ix=np.where(hpv>0.05)[0]
    sums.append(dict(N=N,denominator=DEN,delta=delta,runtime_sec=time.time()-t0,
                     onset_Hperp_005=(int(ix[0]) if len(ix) else ''),
                     max_Hperp_frac=float(hpv.max()),final_Hperp_frac=float(hpv[-1]),
                     H_drift_abs=float(abs(rr[-1]['H_total']-rr[0]['H_total'])),
                     initial_closure=rr[0]['global_closure'],final_closure=rr[-1]['global_closure']))
    np.savez_compressed(os.path.join(OUT_DIR,f'hm_N{N}_den_124_states_500.npz'),Z=states,N=N,denominator=DEN,delta=delta)

with open(os.path.join(OUT_DIR,'timeseries_den124_500.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
with open(os.path.join(OUT_DIR,'summary_den124_500.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(sums[0].keys())); w.writeheader(); w.writerows(sums)
