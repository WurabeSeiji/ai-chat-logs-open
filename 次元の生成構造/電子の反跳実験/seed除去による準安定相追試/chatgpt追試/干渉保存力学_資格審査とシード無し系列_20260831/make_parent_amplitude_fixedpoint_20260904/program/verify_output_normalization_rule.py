#!/usr/bin/env python3
import csv, sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'source_snapshot'))
import original_engine as eng
BASE=40260721; ITERS=2000; BETA=.5
rows=[]
for N in range(3,11):
 for d in range(20):
  s=eng.LowRankSystem(N); rng=np.random.default_rng(BASE+1000*N+d); th=rng.uniform(0,2*np.pi,s.m)
  for it in range(ITERS):
   s.set_theta(th); ev,EV=np.linalg.eig(s.J@s.G); idx=int(np.argmin(ev.imag)); y=EV[:,idx].astype(complex); v=s.w(y); sig_gen=float(-ev[idx].imag); phi=np.angle(v); th=np.angle(.5*np.exp(1j*th)+.5*np.exp(1j*phi))
  s.set_theta(np.angle(v)); kv=s.kmatvec(v); iv=1j*kv; mu=np.vdot(v,iv)/np.vdot(v,v); cres=float(np.linalg.norm(iv-mu*v)/np.linalg.norm(v))
  ev2=np.linalg.eigvals(s.J@s.G); sig_self=float(-ev2[np.argmin(ev2.imag)].imag)
  yself=(1j/sig_self)*(s.J@s.wt(v)); yn2=float(np.vdot(yself,yself).real)
  rows.append(dict(N=N,delta=d,converged=(cres<1e-9),corrected_residual=cres,sigma_gen=sig_gen,sigma_self=sig_self,sigma_ratio=sig_self/sig_gen,yself_norm2=yn2,norm2=float(np.vdot(v,v).real)))
conv=[r for r in rows if r['converged']]
print('converged',len(conv),'of',len(rows))
print('sigma_ratio max dev from 2',max(abs(r['sigma_ratio']-2) for r in conv))
print('yself_norm2 max dev from 1/4',max(abs(r['yself_norm2']-.25) for r in conv))
for N in range(3,11):
 c=[r for r in conv if r['N']==N]
 print(N,len(c),max(abs(r['sigma_ratio']-2) for r in c) if c else None,max(abs(r['yself_norm2']-.25) for r in c) if c else None)
with open(HERE/'results'/'output_normalization_rule_N3_N10_20seeds_2000iter.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
