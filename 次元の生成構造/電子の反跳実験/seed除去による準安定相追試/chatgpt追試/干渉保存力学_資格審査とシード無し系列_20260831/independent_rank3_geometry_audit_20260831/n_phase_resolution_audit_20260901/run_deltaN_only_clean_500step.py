#!/usr/bin/env python3
import os, math, json, csv, time
import numpy as np

IN_DIR='/mnt/data/hm_rank3_inputs'
OUT_DIR='/mnt/data/deltaN_only_clean_20260901'
STEPS=500


def edges(N):
    a,b=np.triu_indices(N,k=1)
    return a,b

def adjacency(N):
    ea,eb=edges(N); M=len(ea)
    A=np.zeros((M,M),float)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e])
        share[e]=False
        A[e,share]=1.0
    return A

def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:])
    np.fill_diagonal(H,0.0)
    return H

def step(z,A,delta):
    H=H_of(z,A)
    w,V=np.linalg.eigh(H)
    return V @ (np.exp(-1j*delta*w)*(V.conj().T@z))

def plane_basis(v):
    p=v.real.copy(); npn=np.linalg.norm(p)
    if npn<1e-14:
        p=v.imag.copy(); npn=np.linalg.norm(p)
    p/=npn
    q=v.imag-(v.imag@p)*p
    nq=np.linalg.norm(q)
    if nq<1e-14:
        # deterministic orthogonal fallback
        q=np.zeros_like(p); q[np.argmin(np.abs(p))]=1.0; q-= (q@p)*p; nq=np.linalg.norm(q)
    q/=nq
    return p,q

def phase_grid_resid(z,N):
    th=np.angle(z)
    u=np.mean(np.exp(1j*N*th))
    if abs(u)<1e-300:
        return 1.0, float('nan')
    phi=np.angle(u)/N
    resid=float(np.mean(1.0-np.cos(N*(th-phi))))
    return resid,float(phi)

def metrics(z,p,q,N,A):
    h=float(np.vdot(z,z).real)
    zp=z-p*(p@z)-q*(q@z)
    hp=float(np.vdot(zp,zp).real)
    a2=np.abs(z)**2
    pr=float(h*h/np.sum(a2*a2))
    clos=float(abs(z@z)/h)
    loc=[]
    ea,eb=edges(N); zz=z*z
    S=np.zeros(N,complex)
    for e,(i,j) in enumerate(zip(ea,eb)):
        S[i]+=zz[e]; S[j]+=zz[e]
    local=float(np.max(np.abs(S))/h)
    grid,phi=phase_grid_resid(z,N)
    H=H_of(z,A)
    ev=np.linalg.eigvalsh(H)
    return dict(H_total=h,Hperp_frac=hp/h,global_closure=clos,local_closure=local,
                PR_over_M=pr/len(z),amp_min=float(np.abs(z).min()),amp_max=float(np.abs(z).max()),
                amp_std=float(np.abs(z).std()),phase_ZN_resid=grid,phase_coset_offset=phi,
                H_spectral_radius=float(np.max(np.abs(ev))))

def run_one(N,delta,label):
    path=os.path.join(IN_DIR,f'hm_N{N}_states_treatment.npz')
    Z0=np.load(path)['Z'][0].astype(np.complex128)
    A=adjacency(N); p,q=plane_basis(Z0)
    states=np.empty((STEPS+1,len(Z0)),complex)
    rows=[]; z=Z0.copy(); t0=time.time()
    for t in range(STEPS+1):
        states[t]=z
        m=metrics(z,p,q,N,A); m.update(N=N,step=t,scheme=label,delta=delta)
        rows.append(m)
        if t<STEPS: z=step(z,A,delta)
    runtime=time.time()-t0
    np.savez_compressed(os.path.join(OUT_DIR,f'hm_N{N}_{label}_states_500.npz'),Z=states,delta=delta,N=N)
    return rows,runtime

def summarize(rows,runtime):
    r0=rows[0]; rf=rows[-1]
    hp=np.array([r['Hperp_frac'] for r in rows])
    ix=np.where(hp>0.05)[0]
    onset=int(ix[0]) if len(ix) else None
    return dict(N=rf['N'],scheme=rf['scheme'],delta=rf['delta'],runtime_sec=runtime,
                onset_Hperp_005=onset,max_Hperp_frac=float(hp.max()),final_Hperp_frac=rf['Hperp_frac'],
                H_drift_abs=abs(rf['H_total']-r0['H_total']),
                initial_global_closure=r0['global_closure'],final_global_closure=rf['global_closure'],
                initial_local_closure=r0['local_closure'],final_local_closure=rf['local_closure'],
                initial_phase_ZN_resid=r0['phase_ZN_resid'],final_phase_ZN_resid=rf['phase_ZN_resid'],
                initial_PR_over_M=r0['PR_over_M'],final_PR_over_M=rf['PR_over_M'],
                final_amp_ratio=(rf['amp_max']/max(rf['amp_min'],1e-300)))

all_rows=[]; summaries=[]
for N in range(3,17):
    print('N',N,flush=True)
    label,delta='deltaN',2*math.pi/N
    rows,rt=run_one(N,delta,label)
    all_rows.extend(rows); summaries.append(summarize(rows,rt))
    print(' ',label,'runtime',rt,'final hp',rows[-1]['Hperp_frac'],flush=True)

# CSV
fields=list(all_rows[0].keys())
with open(os.path.join(OUT_DIR,'timeseries_metrics.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(all_rows)
fields2=list(summaries[0].keys())
with open(os.path.join(OUT_DIR,'summary.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields2); w.writeheader(); w.writerows(summaries)
with open(os.path.join(OUT_DIR,'manifest.json'),'w') as f:
    json.dump(dict(steps=STEPS,input_dir=IN_DIR,schemes={'deltaN':'2*pi/N'},
                   interaction='H_ef=A_ef*conj(z_e)*z_f, z_next=exp(-i delta H(z)) z',
                   no_phase_projection=True,no_normalization=True),f,indent=2)
print('DONE')
