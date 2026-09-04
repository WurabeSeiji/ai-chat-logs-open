#!/usr/bin/env python3
import csv,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE/'source_snapshot'))
import original_engine as eng
base_rows=list(csv.DictReader(open(HERE/'results'/'output_normalization_rule_N3_N10_20seeds_2000iter.csv')))
sel=[(int(r['N']),int(r['delta'])) for r in base_rows if r['converged']=='True']
BASE=40260721; ITERS=2000
rows=[]
for N,d in sel:
 s=eng.LowRankSystem(N); rng=np.random.default_rng(BASE+1000*N+d); th=rng.uniform(0,2*np.pi,s.m)
 for it in range(ITERS):
  s.set_theta(th);ev,EV=np.linalg.eig(s.J@s.G);idx=int(np.argmin(ev.imag));v=s.w(EV[:,idx].astype(complex));phi=np.angle(v);th=np.angle(.5*np.exp(1j*th)+.5*np.exp(1j*phi))
 # recompute current theta tuple
 s.set_theta(th);ev,EV=np.linalg.eig(s.J@s.G);idx=int(np.argmin(ev.imag));v=s.w(EV[:,idx].astype(complex));phi=np.angle(v);sg=float(-ev[idx].imag)
 A=eng._adjacency(s); pairs=np.argwhere(np.triu(A,1)>0); M=s.m; eye=np.eye(M)
 Kt=np.column_stack([s.kmatvec(eye[:,j]) for j in range(M)])
 s.set_theta(phi); Kp=np.column_stack([s.kmatvec(eye[:,j]) for j in range(M)])
 kp=Kp@v; mu=np.vdot(v,1j*kp)/np.vdot(v,v); selfres=float(np.linalg.norm(1j*kp-mu*v)/np.linalg.norm(v))
 materr=float(np.linalg.norm(Kp+2*Kt)/max(np.linalg.norm(Kp),1e-30))
 wrap=lambda x: np.angle(np.exp(1j*x))
 dth=np.array([wrap(th[f]-th[e]) for e,f in pairs]); dph=np.array([wrap(phi[f]-phi[e]) for e,f in pairs])
 qth=float(np.min(np.abs(np.sin(dth)[:,None]-np.array([-0.5,0,0.5])[None,:]),axis=1).max())
 qph=float(np.min(np.abs(np.sin(dph)[:,None]-np.array([-1,0,1])[None,:]),axis=1).max())
 delta=wrap(phi-th); u3=np.exp(3j*delta); cubesp=float(np.max(np.abs(u3-u3[0])))
 thnext=np.angle(.5*np.exp(1j*th)+.5*np.exp(1j*phi)); midpoint=float(np.max(np.abs(np.exp(2j*thnext)-np.exp(1j*(th+phi)))))
 rows.append(dict(N=N,delta_seed=d,self_residual=selfres,sigma_gen=sg,Kphi_plus_2Ktheta_rel=materr,sin_theta_quant_error=qth,sin_phi_quant_error=qph,cube_delta_spread=cubesp,midpoint_identity_error=midpoint))
out=HERE/'results'/'phase_quantization_converged128.csv'
with open(out,'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print('cases',len(rows))
for k in ['self_residual','Kphi_plus_2Ktheta_rel','sin_theta_quant_error','sin_phi_quant_error','cube_delta_spread','midpoint_identity_error']:
 print(k,'max',max(r[k] for r in rows),'mean',sum(r[k] for r in rows)/len(rows))
