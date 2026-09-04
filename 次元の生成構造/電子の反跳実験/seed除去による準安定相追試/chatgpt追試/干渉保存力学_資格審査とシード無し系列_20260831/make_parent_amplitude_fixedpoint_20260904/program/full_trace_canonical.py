#!/usr/bin/env python3
import csv, sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'source_snapshot'))
import original_engine as eng
BASE=40260721; ITERS=400; BETA=.5
rows=[]; summary=[]
for N in range(3,17):
    s=eng.LowRankSystem(N); rng=np.random.default_rng(BASE+1000*N); theta=rng.uniform(0,2*np.pi,s.m)
    vals=[]
    for it in range(ITERS):
        s.set_theta(theta); ev,EV=np.linalg.eig(s.J@s.G); idx=int(np.argmin(ev.imag)); y=EV[:,idx].astype(complex); v=s.w(y)
        sig=float(-ev[idx].imag); chi=float(np.real(1j*np.vdot(y,s.J@y))); norm2=float(np.vdot(v,v).real)
        tn=np.angle(v); d=np.angle(np.exp(1j*(tn-theta)))
        rows.append(dict(N=N,iteration=it+1,sigma=sig,sigma2=sig*sig,chi=chi,norm2=norm2,r2=norm2/s.m,Nr2=N*norm2/s.m,identity_error=abs(norm2-sig*chi),phase_output_vs_input_rms=float(np.sqrt(np.mean(d*d)))))
        vals.append((norm2,sig*sig))
        theta=np.angle((1-BETA)*np.exp(1j*theta)+BETA*np.exp(1j*tn))
    # earliest t such that all subsequent norm2 and sigma2 stay within 1e-10 of final
    fn,fs=vals[-1]; onset=None
    for k in range(ITERS):
        if max(abs(a-fn) for a,b in vals[k:])<1e-10 and max(abs(b-fs) for a,b in vals[k:])<1e-10:
            onset=k+1; break
    summary.append(dict(N=N,stable_from_iteration=onset,final_norm2=fn,final_sigma2=fs,max_tail_norm2_dev=max(abs(a-fn) for a,b in vals[(onset-1 if onset else 0):]),max_tail_sigma2_dev=max(abs(b-fs) for a,b in vals[(onset-1 if onset else 0):])))
    print(N,onset,fn,fs)
with open(HERE/'results'/'canonical_full_trace_N3_N16.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
with open(HERE/'results'/'canonical_invariant_onset_N3_N16.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=summary[0].keys()); w.writeheader(); w.writerows(summary)
