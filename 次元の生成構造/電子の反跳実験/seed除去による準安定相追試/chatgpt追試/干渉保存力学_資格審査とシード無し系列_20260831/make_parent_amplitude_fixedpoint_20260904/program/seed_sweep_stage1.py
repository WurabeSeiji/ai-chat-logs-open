#!/usr/bin/env python3
import csv, math, sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'source_snapshot'))
import original_engine as eng
BASE=40260721; ITERS=400; BETA=0.5
NS=range(3,11); DELTAS=range(20)

def run(N,delta):
    seed=BASE+1000*N+delta
    s=eng.LowRankSystem(N); rng=np.random.default_rng(seed)
    theta=rng.uniform(0,2*np.pi,s.m)
    first_stable=None; prev_norm=None; prev_sig2=None; max_iderr=0.0
    v=y=None; sig=None; chi=None
    trace=[]
    for it in range(ITERS):
        s.set_theta(theta); ev,EV=np.linalg.eig(s.J@s.G); idx=int(np.argmin(ev.imag))
        y=EV[:,idx].astype(complex); v=s.w(y); sig=float(-ev[idx].imag)
        chi=float(np.real(1j*np.vdot(y,s.J@y))); norm2=float(np.vdot(v,v).real)
        max_iderr=max(max_iderr,abs(norm2-sig*chi))
        if prev_norm is not None and abs(norm2-prev_norm)<1e-12 and abs(sig*sig-prev_sig2)<1e-12 and first_stable is None:
            first_stable=it+1
        else:
            if prev_norm is not None and (abs(norm2-prev_norm)>=1e-12 or abs(sig*sig-prev_sig2)>=1e-12): first_stable=None
        prev_norm=norm2; prev_sig2=sig*sig
        tn=np.angle(v); theta=np.angle((1-BETA)*np.exp(1j*theta)+BETA*np.exp(1j*tn))
        if it>=ITERS-10: trace.append((norm2,sig*sig,chi))
    # corrected self-consistency using theta=angle(v)
    s.set_theta(np.angle(v)); kv=s.kmatvec(v); iv=1j*kv; mu=np.vdot(v,iv)/np.vdot(v,v); cres=float(np.linalg.norm(iv-mu*v)/np.linalg.norm(v))
    return dict(N=N,delta=delta,seed=seed,M=s.m,norm2=norm2,r2=norm2/s.m,Nr2=N*norm2/s.m,
                sigma2=sig*sig,chi=chi,identity_error=max_iderr,corrected_residual=cres,
                last10_norm2_span=max(x[0] for x in trace)-min(x[0] for x in trace),
                last10_sigma2_span=max(x[1] for x in trace)-min(x[1] for x in trace))

rows=[]
for N in NS:
    for d in DELTAS:
        r=run(N,d); rows.append(r)
    vals=np.array([r['norm2'] for r in rows if r['N']==N])
    print(f'N={N} norm2 range {vals.min():.15g} .. {vals.max():.15g} unique~ {len(np.unique(np.round(vals,12)))}')

out=HERE/'results'/'seed_sweep_N3_N10_20seeds.csv'
with open(out,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
# cluster report by norm2 rounded 1e-10
with open(HERE/'results'/'seed_branch_summary_N3_N10.csv','w',newline='') as f:
    fields=['N','branch_id','count','norm2_mean','norm2_min','norm2_max','Nr2_mean','sigma2_mean','corrected_residual_max']
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
    for N in NS:
        rr=[r for r in rows if r['N']==N]; rr.sort(key=lambda x:x['norm2']); clusters=[]
        for r in rr:
            if not clusters or abs(r['norm2']-np.mean([x['norm2'] for x in clusters[-1]]))>1e-9: clusters.append([r])
            else: clusters[-1].append(r)
        for bi,c in enumerate(clusters,1):
            w.writerow(dict(N=N,branch_id=bi,count=len(c),norm2_mean=np.mean([x['norm2'] for x in c]),norm2_min=min(x['norm2'] for x in c),norm2_max=max(x['norm2'] for x in c),Nr2_mean=np.mean([x['Nr2'] for x in c]),sigma2_mean=np.mean([x['sigma2'] for x in c]),corrected_residual_max=max(x['corrected_residual'] for x in c)))
