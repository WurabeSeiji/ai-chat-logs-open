#!/usr/bin/env python3
"""Instantaneous generator spectrum through an exact low-rank factorization.
K=C S^T-S C^T, C=diag(cos theta) B_inc,S=diag(sin theta) B_inc.
Thin QR gives all nonzero eigenvalues in <=2L dimensions. Zero eigenspace
is aggregated, not arbitrarily split. Reads saved states at fixed 512-step grid.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import json,numpy as np
from scipy.linalg import qr,eigh
P=Path('/mnt/data/survey_analysis_20260915');R=json.loads((P/'all_results.json').read_text())
O=P/'instantaneous';O.mkdir(exist_ok=True)
for r in R:
 op=O/(r['run_id']+'.json')
 if op.exists() and 'dominant_gap_relative' in op.read_text():continue
 Z=np.load(r['raw_path'])['Z'];L=r['L'];M=r['M'];ea,eb=np.triu_indices(L,1)
 B=np.zeros((M,L));B[np.arange(M),ea]=1;B[np.arange(M),eb]=1
 J=np.zeros((2*L,2*L),complex);J[:L,L:]=1j*np.eye(L);J[L:,:L]=-1j*np.eye(L)
 out=[]; prev=None
 for t in range(0,4097,512):
  z=Z[t];theta=np.angle(z);W=np.concatenate((np.cos(theta)[:,None]*B,np.sin(theta)[:,None]*B),axis=1)
  Q,T=qr(W,mode='economic',check_finite=False);h=T@J@T.T;h=(h+h.conj().T)/2
  e,V=eigh(h,check_finite=False);coeff=V.conj().T@(Q.T@z);power=abs(coeff)**2/np.vdot(z,z).real
  valid=abs(e)>1e-10*max(1,float(max(abs(e))));null=float(max(0,1-power[valid].sum()))
  if valid.any():
   ix=np.where(valid)[0][np.argmax(power[valid])];ip=int(np.argmin(abs(e+e[ix])))
   u=Q@V[:,ix];pairweight=float(power[ix]+power[ip]); sub=Q@V[:,[ix,ip]]
   alignment=float(np.linalg.norm(prev.conj().T@sub,'fro')**2/2) if prev is not None else None
   pred=(1+np.sqrt(max(0,1-r['chi']**2)))/2
   out.append({'t':t,'dominant_weight':float(power[ix]),'opposite_weight':float(power[ip]),'pair_weight':pairweight,'null_weight':null,
               'dominant_eigenvalue_H':float(e[ix]),'cycles_per_step':float(-e[ix]/r['den']),
               'principal_cycles_per_step':float((-e[ix]/r['den']+.5)%1-.5),
               'plane_overlap_previous':alignment,'single_ellipse_dominant_weight':float(pred),
               'nonzero_rank':int(valid.sum()),'dominant_gap_relative':float(np.min(np.delete(abs(e-e[ix]),ix))/max(1,abs(e[ix])))})
   prev=sub
  else:out.append({'t':t,'dominant_weight':0.,'opposite_weight':0.,'pair_weight':0.,'null_weight':1.,'nonzero_rank':0})
 op.write_text(json.dumps({'run_id':r['run_id'],'samples':out},indent=2));print(r['run_id'],flush=True)
ALL=[json.loads((O/(r['run_id']+'.json')).read_text()) for r in R]
(P/'all_instantaneous.json').write_text(json.dumps(ALL,indent=2));print('ALL99_INSTANTANEOUS_COMPLETE')
