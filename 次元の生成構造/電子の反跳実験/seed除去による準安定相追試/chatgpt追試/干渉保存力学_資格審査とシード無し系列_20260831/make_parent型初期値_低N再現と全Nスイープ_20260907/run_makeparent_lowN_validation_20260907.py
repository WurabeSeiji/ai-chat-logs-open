#!/usr/bin/env python3
import os, math, csv, json
import numpy as np

BASE=os.path.dirname(os.path.abspath(__file__))
PARENT_DIR=os.path.join(BASE,'parents_actual_N3_N7')
OUT=os.path.join(BASE,'lowN_validation','results')
os.makedirs(OUT,exist_ok=True)
STEPS=500
N_LIST=[3,4,5,6,7]
OFFSETS=(-2,-1,0,1,2)

def edges(N):
    a,b=np.triu_indices(N,k=1); return a.astype(np.int64),b.astype(np.int64)

def adjacency(N):
    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    return A

def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)

def one_step(z,A,den):
    H=H_of(np.exp(1j*np.angle(z)),A); H=(1j*np.imag(H)).astype(np.complex128,copy=False)
    w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)

def plane(v):
    p=v.real.astype(np.float64,copy=True); p/=np.linalg.norm(p)
    q=v.imag.astype(np.float64,copy=True); q-=np.dot(q,p)*p; q/=np.linalg.norm(q)
    return p,q

def metrics(z,p,q):
    h=np.vdot(z,z).real; zp=z-p*np.dot(p,z)-q*np.dot(q,z); hp=np.vdot(zp,zp).real
    return float(hp/h),float(h),float(abs(z@z)/h)

rows=[]; summaries=[]
for N in N_LIST:
    z0=np.array(np.load(os.path.join(PARENT_DIR,f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],dtype=np.complex128,copy=True)
    A=adjacency(N); p,q=plane(z0)
    pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0] + [(124,'124')]
    for den,label in pairs:
        z=z0.copy(); vals=np.empty(STEPS+1); states=np.empty((STEPS+1,z.size),np.complex128); closures=np.empty(STEPS+1); htot=np.empty(STEPS+1)
        for t in range(STEPS+1):
            states[t]=z; vals[t],htot[t],closures[t]=metrics(z,p,q)
            if t<STEPS: z=one_step(z,A,den)
        np.savez_compressed(os.path.join(OUT,f'hm_N{N}_den_{den}_states_500.npz'),Z=states,N=np.int64(N),denominator=np.int64(den),steps=np.int64(STEPS))
        rows.extend((N,label,den,t,vals[t],htot[t],closures[t]) for t in range(STEPS+1))
        ix=np.flatnonzero(vals>0.05)
        summaries.append((N,label,den,int(ix[0]) if ix.size else -1,float(vals[0]),float(vals[1]),float(vals[10]),float(vals[100]),float(vals[-1]),float(vals.max()),float(closures[0]),float(closures[-1])))
with open(os.path.join(OUT,'timeseries.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','series','denominator','step','Hperp_frac','H_total','global_closure']); w.writerows(rows)
with open(os.path.join(OUT,'summary.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','series','denominator','onset_gt_0.05','initial','step1','step10','step100','final','max','closure0','closure500']); w.writerows(summaries)
print('LOW-N VALIDATION DONE')
