#!/usr/bin/env python3
import csv, sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'source_snapshot'))
import original_engine as eng
BASE=40260721; ITERS=2000; BETA=.5
NS=range(3,11); DELTAS=range(20)

def run(N,d):
 seed=BASE+1000*N+d; s=eng.LowRankSystem(N); rng=np.random.default_rng(seed); th=rng.uniform(0,2*np.pi,s.m)
 v=y=None; sig=chi=None
 for it in range(ITERS):
  s.set_theta(th); ev,EV=np.linalg.eig(s.J@s.G); idx=int(np.argmin(ev.imag)); y=EV[:,idx].astype(complex); v=s.w(y); sig=float(-ev[idx].imag); chi=float(np.real(1j*np.vdot(y,s.J@y)))
  tn=np.angle(v); th=np.angle((1-BETA)*np.exp(1j*th)+BETA*np.exp(1j*tn))
 s.set_theta(np.angle(v)); kv=s.kmatvec(v); iv=1j*kv; mu=np.vdot(v,iv)/np.vdot(v,v); cres=float(np.linalg.norm(iv-mu*v)/np.linalg.norm(v)); n2=float(np.vdot(v,v).real)
 return dict(N=N,delta=d,seed=seed,iterations=ITERS,norm2=n2,r2=n2/s.m,Nr2=N*n2/s.m,sigma2=sig*sig,chi=chi,corrected_residual=cres,converged=(cres<1e-9))
rows=[]
for N in NS:
 for d in DELTAS: rows.append(run(N,d))
 rr=[x for x in rows if x['N']==N]; conv=[x for x in rr if x['converged']]
 print(f'N={N} converged {len(conv)}/20 range '+(f"{min(x['norm2'] for x in conv):.15g}..{max(x['norm2'] for x in conv):.15g}" if conv else 'none'))
with open(HERE/'results'/'seed_sweep_N3_N10_20seeds_2000iter.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
with open(HERE/'results'/'converged_branch_summary_N3_N10_2000iter.csv','w',newline='') as f:
 fields=['N','branch_id','count','norm2_mean','Nr2_mean','sigma2_mean','residual_max','seed_deltas']; w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
 for N in NS:
  rr=sorted([x for x in rows if x['N']==N and x['converged']],key=lambda x:x['norm2']); clusters=[]
  for r in rr:
   if not clusters or abs(r['norm2']-np.mean([x['norm2'] for x in clusters[-1]]))>1e-9: clusters.append([r])
   else: clusters[-1].append(r)
  for bi,c in enumerate(clusters,1):
   w.writerow(dict(N=N,branch_id=bi,count=len(c),norm2_mean=np.mean([x['norm2'] for x in c]),Nr2_mean=np.mean([x['Nr2'] for x in c]),sigma2_mean=np.mean([x['sigma2'] for x in c]),residual_max=max(x['corrected_residual'] for x in c),seed_deltas=';'.join(str(x['delta']) for x in c)))
