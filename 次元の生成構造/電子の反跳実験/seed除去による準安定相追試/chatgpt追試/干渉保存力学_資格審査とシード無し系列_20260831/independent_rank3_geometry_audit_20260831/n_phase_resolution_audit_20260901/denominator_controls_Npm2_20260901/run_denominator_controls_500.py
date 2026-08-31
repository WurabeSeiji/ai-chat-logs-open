#!/usr/bin/env python3
import os, math, json, csv, time
import numpy as np

IN_DIR='/mnt/data/hm_rank3_inputs'
OUT_DIR='/mnt/data/denominator_controls_20260901'
STEPS=500
OFFSETS=[-2,-1,0,1,2]
os.makedirs(OUT_DIR,exist_ok=True)
os.makedirs(os.path.join(OUT_DIR,'figures'),exist_ok=True)

def edges(N):
    return np.triu_indices(N,k=1)

def adjacency(N):
    ea,eb=edges(N); M=len(ea)
    A=np.zeros((M,M),float)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e])
        share[e]=False; A[e,share]=1.0
    return A

def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H

def step(z,A,delta):
    H=H_of(z,A); w,V=np.linalg.eigh(H)
    return V @ (np.exp(-1j*delta*w)*(V.conj().T@z))

def plane_basis(v):
    p=v.real/np.linalg.norm(v.real)
    q=v.imag-(v.imag@p)*p; q/=np.linalg.norm(q)
    return p,q

def metrics(z,p,q):
    h=float(np.vdot(z,z).real)
    zp=z-p*(p@z)-q*(q@z)
    hp=float(np.vdot(zp,zp).real)
    return h,hp/h,float(abs(z@z)/h)

all_rows=[]; summaries=[]
for off in OFFSETS:
    scheme=f'N{off:+d}' if off else 'N'
    print('SCHEME',scheme,flush=True)
    for N in range(3,17):
        den=N+off
        if den <= 0: raise RuntimeError((N,off,den))
        delta=2*math.pi/den
        path=os.path.join(IN_DIR,f'hm_N{N}_states_treatment.npz')
        Z0=np.load(path)['Z'][0].astype(np.complex128)
        A=adjacency(N); p,q=plane_basis(Z0)
        z=Z0.copy(); states=np.empty((STEPS+1,len(z)),np.complex128)
        rows=[]; t0=time.time()
        for t in range(STEPS+1):
            states[t]=z
            h,hp,clo=metrics(z,p,q)
            r=dict(offset=off,scheme=scheme,N=N,denominator=den,delta=delta,step=t,H_total=h,Hperp_frac=hp,global_closure=clo)
            rows.append(r); all_rows.append(r)
            if t<STEPS: z=step(z,A,delta)
        rt=time.time()-t0
        hpv=np.array([r['Hperp_frac'] for r in rows]); ix=np.where(hpv>0.05)[0]
        summaries.append(dict(offset=off,scheme=scheme,N=N,denominator=den,delta=delta,runtime_sec=rt,
                              onset_Hperp_005=(int(ix[0]) if len(ix) else ''),max_Hperp_frac=float(hpv.max()),
                              final_Hperp_frac=float(hpv[-1]),H_drift_abs=float(abs(rows[-1]['H_total']-rows[0]['H_total'])),
                              initial_closure=rows[0]['global_closure'],final_closure=rows[-1]['global_closure']))
        np.savez_compressed(os.path.join(OUT_DIR,f'hm_N{N}_den_{den}_offset_{off:+d}_states_500.npz'),Z=states,N=N,denominator=den,delta=delta,offset=off)
        print(' ',N,'den',den,'rt',round(rt,2),'onset',summaries[-1]['onset_Hperp_005'],'final',hpv[-1],flush=True)

with open(os.path.join(OUT_DIR,'timeseries_metrics_denominator_controls_500.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(all_rows[0].keys())); w.writeheader(); w.writerows(all_rows)
with open(os.path.join(OUT_DIR,'summary_denominator_controls_500.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summaries[0].keys())); w.writeheader(); w.writerows(summaries)
with open(os.path.join(OUT_DIR,'manifest.json'),'w') as f:
    json.dump({'steps':STEPS,'N_range':[3,16],'offsets':OFFSETS,'delta':'2*pi/(N+offset)',
               'input':'original hm initial state Z[0] from hm_rank3_inputs; complex128 unchanged',
               'interaction':'H_ef=A_ef*conj(z_e)*z_f; z_next=exp(-i delta H(z))z',
               'seed':'none beyond original initial-condition floating-point content','normalization':False},f,indent=2)
print('DONE')
