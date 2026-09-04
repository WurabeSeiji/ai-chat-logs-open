#!/usr/bin/env python3
import csv, sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'source_snapshot'))
import original_engine as eng
BASE=40260721; ITERS=2000
base=list(csv.DictReader(open(HERE/'results'/'output_normalization_rule_N3_N10_20seeds_2000iter.csv')))
sel=[(int(r['N']),int(r['delta'])) for r in base if r['converged']=='True']
def step(s,t):
    s.set_theta(t); ev,V=np.linalg.eig(s.J@s.G); j=int(np.argmin(ev.imag)); v=s.w(V[:,j].astype(complex)); p=np.angle(v)
    tn=np.angle(.5*np.exp(1j*t)+.5*np.exp(1j*p)); return tn,p,float(-ev[j].imag)
def spread(z): return float(np.max(np.abs(z-z[0])))
def wrap(x): return np.angle(np.exp(1j*x))
rows=[]
for N,dseed in sel:
    s=eng.LowRankSystem(N); t=np.random.default_rng(BASE+1000*N+dseed).uniform(0,2*np.pi,s.m)
    for _ in range(ITERS): t,_,_=step(s,t)
    t0=t.copy(); t1,p0,sg0=step(s,t0); t2,p1,sg1=step(s,t1)
    # projective means: equality modulo edgewise pi and one global phase, represented by exp(2i difference) common across edges
    th2cycle=spread(np.exp(2j*(t2-t0)))
    phfixed=spread(np.exp(2j*(p1-p0)))
    mid0=float(np.max(np.abs(np.exp(2j*t1)-np.exp(1j*(t0+p0)))))
    mid1=float(np.max(np.abs(np.exp(2j*t2)-np.exp(1j*(t1+p1)))))
    d=wrap(t1-t0); delta=wrap(p0-t0)
    step6=spread(np.exp(6j*d))
    cube=spread(np.exp(3j*delta))
    delta_mid=float(np.max(np.abs(np.exp(1j*delta)-np.exp(2j*d))))
    # direct derived identity: use projective quantities Rtheta=e^{2i(t2-t0)}, Rphi=e^{2i(p1-p0)}.
    # From midpoint: exp(6i(t1-t0)) differs from Rtheta/Rphi only by a global factor, so spreads should agree.
    derived=np.exp(6j*(t1-t0)) * np.exp(2j*(p1-p0)) / np.exp(2j*(t2-t0))
    derived_global_spread=spread(derived)
    rows.append(dict(N=N,delta_seed=dseed,sigma0=sg0,sigma1=sg1,
        theta_projective_two_cycle_spread=th2cycle,
        phi_projective_fixed_spread=phfixed,
        midpoint0_error=mid0,midpoint1_error=mid1,
        sixth_step_spread=step6,cube_delta_spread=cube,
        delta_equals_2step_error=delta_mid,
        derived_global_identity_spread=derived_global_spread))
out=HERE/'results'/'projective_cycle_proof_converged128.csv'
with open(out,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('cases',len(rows))
for k in ['theta_projective_two_cycle_spread','phi_projective_fixed_spread','midpoint0_error','midpoint1_error','sixth_step_spread','cube_delta_spread','delta_equals_2step_error','derived_global_identity_spread']:
    v=[r[k] for r in rows]; print(k,'max',max(v),'mean',sum(v)/len(v))
