#!/usr/bin/env python3
"""All-run post-collection diagnostic: Q2-anchored real 2-frame, no fitted basis.
This is an algebraic coordinate readout, NOT a proof of two independent waves.
Run after analyze_survey.py. No integration or mutation of inputs.
"""
from pathlib import Path
import json, math, os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
from scipy import sparse
ROOT=Path(os.environ.get('SURVEY_OUTPUT_DIR',str(Path(__file__).resolve().parent)))
ROOT.mkdir(parents=True,exist_ok=True)
R=json.loads((ROOT/'all_results.json').read_text())
for r in R:
 out=ROOT/'per_run'/f'{r["run_id"]}_frame.json'
 if out.exists(): continue
 Z=np.load(r['raw_path'],allow_pickle=False)['Z']; H0=r['H0']; Q0=complex(r['Q20_re'],r['Q20_im']); chi=abs(Q0)/H0
 if r['analytic_fixed'] or 1-chi<1e-12:
  out.write_text(json.dumps({'run_id':r['run_id'],'frame_valid':False,'reason':'real-state fixed-point sector'})); continue
 alpha=np.angle(Q0)/2; Zg=Z*np.exp(-1j*alpha)
 aa=math.sqrt((H0+abs(Q0))/2);bb=math.sqrt((H0-abs(Q0))/2)
 p=Zg.real/aa;q=Zg.imag/bb
 cpp=np.sum(p[:-1]*p[1:],axis=1);cpq=np.sum(p[:-1]*q[1:],axis=1)
 cqp=np.sum(q[:-1]*p[1:],axis=1);cqq=np.sum(q[:-1]*q[1:],axis=1)
 phi=np.arctan2(cpq-cqp,cpp+cqq)
 gramerr=max(float(np.max(abs(np.sum(p*p,axis=1)-1))),float(np.max(abs(np.sum(q*q,axis=1)-1))),float(np.max(abs(np.sum(p*q,axis=1)))))
 plane_leak=np.sqrt(np.maximum(0,1-(cpp*cpp+cpq*cpq+cqp*cqp+cqq*cqq)/2))
 L=r['L'];M=r['M'];ea,eb=np.triu_indices(L,1)
 Inc=sparse.csr_matrix((np.ones(2*M),(np.r_[np.arange(M),np.arange(M)],np.r_[ea,eb])),shape=(M,L))
 def A(Y):return (Inc@(Inc.T@Y.T)).T-2*Y
 u=np.exp(1j*np.angle(Zg))
 KZ=(u.conj()*A(u*Zg)-u*A(u.conj()*Zg))/(2j)
 Kp=KZ.real/aa;Kq=KZ.imag/bb
 omega=np.sum(p*Kq,axis=1)
 offspeed=np.sqrt(np.maximum(0,(np.sum(Kp*Kp+Kq*Kq,axis=1)-2*omega**2)/2))
 expected=r['den']
 dt=2*np.pi/r['den']; phase_error=np.angle(np.exp(1j*(phi-dt*omega[:-1])))
 result={'run_id':r['run_id'],'frame_valid':True,'frame_gram_error_max':gramerr,'windows':{}}
 for name,lo,hi in [('W1',0,1024),('W2',1024,2048),('W3',2048,3072),('W4',3072,4096)]:
  w=omega[lo:hi];pp=phi[lo:hi]; err=phase_error[lo:hi]
  result['windows'][name]={'omega_generator_mean':float(w.mean()),'omega_generator_std':float(w.std()),
   'generator_cycles_per_step_mean':float(w.mean()/r['den']),
   'principal_phase_mean':float(pp.mean()),'principal_phase_std':float(pp.std()),
   'sampled_cycles_per_step':float(pp.mean()/(2*np.pi)),
   'sampled_principal_period':float(2*np.pi/abs(pp.mean())) if abs(pp.mean())>1e-12 else None,
   'off_plane_step_rms':float(np.sqrt(np.mean(plane_leak[lo:hi]**2))),
   'off_plane_generator_rms':float(np.sqrt(np.mean(offspeed[lo:hi]**2))),
   'phase_prediction_error_rms':float(np.sqrt(np.mean(err**2))),
   'polar_orientation_flip_fraction':float(np.mean((cpp[lo:hi]*cqq[lo:hi]-cpq[lo:hi]*cqp[lo:hi])<0))}
 np.savez_compressed(ROOT/'timeseries'/f'{r["run_id"]}_frame.npz',phi_principal=phi,omega_generator=omega,plane_leak=plane_leak,offspeed=offspeed,phase_prediction_error=phase_error)
 out.write_text(json.dumps(result,indent=2));print('FRAME',r['run_id'],flush=True)
allframes=[json.loads((ROOT/'per_run'/f'{r["run_id"]}_frame.json').read_text()) for r in R]
(ROOT/'all_frames.json').write_text(json.dumps(allframes,indent=2))
print('ALL_99_FRAMES_COMPLETE',flush=True)
