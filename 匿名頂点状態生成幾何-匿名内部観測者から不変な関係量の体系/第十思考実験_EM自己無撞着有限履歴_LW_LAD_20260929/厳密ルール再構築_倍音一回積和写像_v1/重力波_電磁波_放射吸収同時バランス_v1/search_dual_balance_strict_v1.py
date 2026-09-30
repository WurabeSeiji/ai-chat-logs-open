#!/usr/bin/env python3
import json, math, argparse
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import strict_harmonic_rn_tensor_v1 as rn


def init_mixed(layout,r,q,x0,sigma,amp,eta):
    K=layout.K
    s=np.zeros((layout.nvar,K),complex)
    g=amp*np.exp(-0.5*((r-x0)/sigma)**2)
    for i in range(layout.nx):
        s[layout.field_var(0,i),:]=g[i]
        s[layout.field_var(1,i),:]=eta*g[i]
    qp=[1.0,q,q*q,q*q*q,q*q*q*q]
    for p,val in enumerate(qp): s[layout.q_var(p),:]=val
    return rn.dft_encode(s).reshape(-1)


def run_precompiled(q,x0,eta,M,pi,pj,layout,r,dt=0.12,steps=700,sigma=1.5,amp=1e-4):
    z=init_mixed(layout,r,q,x0,sigma,amp,eta)
    flux=np.zeros((2,2),float)
    max_qerr=0.0; max_nonreal=0.0
    for n in range(steps):
        old=z
        z=rn.transition(old,M,pi,pj)
        if not np.isfinite(z).all(): raise FloatingPointError('nonfinite state')
        qread,qp,qerr,fields,ft,S=rn.readout(z,old,layout,r,dt)
        max_qerr=max(max_qerr,qerr)
        if n in (0,steps-1): max_nonreal=max(max_nonreal,float(np.max(np.abs(rn.decode_state(z,layout).imag))))
        for ch in (0,1):
            flux[ch,0]+=max(0.0,-S[ch,1])*dt
            flux[ch,1]+=max(0.0,S[ch,-2])*dt
    deltas=[]
    for ch in (0,1):
        den=flux[ch].sum()
        deltas.append(float((flux[ch,1]-flux[ch,0])/den) if den>0 else float('nan'))
    return {'q':q,'x0':x0,'eta':eta,
            'flux':{'ch0_inner':float(flux[0,0]),'ch0_outer':float(flux[0,1]),
                    'ch1_inner':float(flux[1,0]),'ch1_outer':float(flux[1,1])},
            'delta':[deltas[0],deltas[1]],
            'objective_norm':float(math.hypot(*deltas)),
            'audit':{'q_power_max_error':max_qerr,'max_decoded_imaginary':max_nonreal}}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--q',type=float,required=True); ap.add_argument('--out',required=True)
    ap.add_argument('--nx',type=int,default=61); ap.add_argument('--K',type=int,default=4); ap.add_argument('--steps',type=int,default=700)
    a=ap.parse_args()
    r=np.linspace(2.4,62.4,a.nx); layout=rn.Layout(nx=a.nx,K=a.K); dt=0.12
    M,pi,pj,nrules=rn.compile_sparse_quadratic(layout,r,dt)
    cache={}; history=[]
    def ev(x0,eta):
        key=(round(float(x0),10),round(float(eta),10))
        if key not in cache:
            rr=run_precompiled(a.q,float(x0),float(eta),M,pi,pj,layout,r,dt=dt,steps=a.steps)
            cache[key]=rr; history.append(rr)
        return cache[key]
    # coarse deterministic scan: enough to seed nonlinear solve, no randomness.
    best=None
    for x0 in np.linspace(2.65,3.35,8):
        for eta in np.linspace(-2.0,2.0,9):
            rr=ev(x0,eta)
            # reject effectively absent channel to avoid meaningless delta
            f=rr['flux']; e0=f['ch0_inner']+f['ch0_outer']; e1=f['ch1_inner']+f['ch1_outer']
            if e0<1e-14 or e1<1e-14: continue
            if best is None or rr['objective_norm']<best['objective_norm']: best=rr
    if best is None: raise RuntimeError('no nonzero two-channel seed')
    def fun(x):
        rr=ev(x[0],x[1]); return np.array(rr['delta'])
    sol=least_squares(fun,[best['x0'],best['eta']],bounds=([2.4,-5.0],[4.0,5.0]),xtol=2e-8,ftol=2e-8,gtol=2e-8,max_nfev=35,verbose=0)
    final=ev(sol.x[0],sol.x[1])
    result={'model':'strict harmonic-only RN coupled simultaneous two-channel flux-balance search',
            'persistent_state':'complex harmonic coefficients only',
            'transition':'one fixed sparse quadratic interaction array; unchanged during search and across q cases',
            'search_variables':'initial-state x0 and channel amplitude ratio eta only',
            'balance_definition':'delta_ch=(F_outer-F_inner)/(F_outer+F_inner); seek delta0=delta1=0',
            'parameters':{'q_initial_state':a.q,'nx':a.nx,'K':a.K,'dt':dt,'steps':a.steps,'sigma':1.5,'amp_ch0':1e-4},
            'interaction':{'shape':list(M.shape),'nnz':int(M.nnz),'monomial_pairs':int(len(pi)),'sample_rules':int(nrules)},
            'coarse_best':best,'optimizer':{'success':bool(sol.success),'status':int(sol.status),'message':sol.message,'nfev':int(sol.nfev),'x':[float(v) for v in sol.x],'cost':float(sol.cost)},
            'final':final,'all_evaluations':history}
    Path(a.out).write_text(json.dumps(result,indent=2,ensure_ascii=False))
    print(json.dumps({'q':a.q,'coarse_best':best,'optimizer':result['optimizer'],'final':final,'evaluations':len(history)},indent=2,ensure_ascii=False))

if __name__=='__main__': main()
