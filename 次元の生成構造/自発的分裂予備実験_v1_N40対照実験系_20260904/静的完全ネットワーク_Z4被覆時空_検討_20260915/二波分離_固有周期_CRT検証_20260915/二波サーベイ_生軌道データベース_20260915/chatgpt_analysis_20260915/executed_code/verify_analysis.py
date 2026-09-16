#!/usr/bin/env python3
"""Independent checks on the derived metrics and local spectral explanation.
Reads existing trajectories only. No new dynamics run.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import numpy as np,json
from scipy.linalg import eigh,svdvals
P=Path('/mnt/data/survey_analysis_20260915')
R=json.loads((P/'all_results.json').read_text()); RR={r['run_id']:r for r in R}
checks=[]
for rid in ['L8_ma2_mb4_den8','L12_ma2_mb3_den12','L32_ma2_mb4_den32']:
 r=RR[rid]; Z=np.load(r['raw_path'])['Z']; X=Z[3072:4096].T
 s=svdvals(X);evals=s*s/1024
 ev=np.array(r['windows']['W4']['eigenvalues'])
 ec=svdvals(X-X.mean(axis=1)[:,None])**2/1024
 ec0=np.array(r['windows']['W4']['centered_eigenvalues'])
 checks.append({'run_id':rid,'svd_vs_covariance_max_scaled':float(np.max(abs(evals-ev))/r['H0']),
                'centered_svd_vs_covariance_max_scaled':float(np.max(abs(ec-ec0))/r['H0'])})
spec=[]
for L in [12,24,32]:
 for den in [L,40]:
  rid=f'L{L}_ma2_mb4_den{den}';r=RR[rid];Z=np.load(r['raw_path'])['Z']; t=3584;z=Z[t]
  ea,eb=np.triu_indices(L,1);Inc=np.zeros((len(ea),L));Inc[np.arange(len(ea)),ea]=1;Inc[np.arange(len(ea)),eb]=1
  A=Inc@Inc.T-2*np.eye(len(ea));u=np.exp(1j*np.angle(z));K=A*np.imag(u.conj()[:,None]*u[None,:]);H=1j*K
  e,V=eigh(H,check_finite=False);weight=abs(V.conj().T@z)**2/np.vdot(z,z).real
  idx=np.argsort(weight)[::-1][:4]
  local=V@(np.exp(-2j*np.pi*e/den)*(V.conj().T@z))
  spec.append({'run_id':rid,'sample_t':t,'local_step_residual':float(np.linalg.norm(local-Z[t+1])/np.linalg.norm(z)),
       'spectral_max_abs':float(max(abs(e))),
       'top_components':[{'eigenvalue_H':float(e[i]),'weight':float(weight[i]),'signed_cycles_per_step':float(-e[i]/den),'principal_cycles':float((-e[i]/den+.5)%1-.5)} for i in idx]})
(P/'independent_analysis_checks.json').write_text(json.dumps({'svd_checks':checks,'local_spectra':spec},indent=2))
print(json.dumps({'svd_checks':checks,'local_spectra':spec},indent=2))
