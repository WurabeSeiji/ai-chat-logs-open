#!/usr/bin/env python3
import json, math, argparse
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import strict_harmonic_rn_tensor_v1 as rn

def init_two_pulses(layout,r,q,x0,x1,eta,sigma=1.5,amp=1e-4):
    K=layout.K; s=np.zeros((layout.nvar,K),complex)
    g0=amp*np.exp(-0.5*((r-x0)/sigma)**2)
    g1=(eta*amp)*np.exp(-0.5*((r-x1)/sigma)**2)
    for i in range(layout.nx):
        s[layout.field_var(0,i),:]=g0[i]
        s[layout.field_var(1,i),:]=g1[i]
    qp=[1.0,q,q*q,q*q*q,q*q*q*q]
    for p,val in enumerate(qp): s[layout.q_var(p),:]=val
    return rn.dft_encode(s).reshape(-1)

def run_precompiled(q,x0,x1,eta,M,pi,pj,layout,r,dt=0.12,steps=700):
    z=init_two_pulses(layout,r,q,x0,x1,eta)
    flux=np.zeros((2,2),float); max_qerr=0.; max_nonreal=0.
    for n in range(steps):
        old=z; z=rn.transition(old,M,pi,pj)
        qread,qp,qerr,fields,ft,S=rn.readout(z,old,layout,r,dt)
        max_qerr=max(max_qerr,qerr)
        if n in (0,steps-1): max_nonreal=max(max_nonreal,float(np.max(np.abs(rn.decode_state(z,layout).imag))))
        for ch in (0,1):
            flux[ch,0]+=max(0.0,-S[ch,1])*dt
            flux[ch,1]+=max(0.0,S[ch,-2])*dt
    delta=[]
    for ch in (0,1):
        den=flux[ch].sum(); delta.append(float((flux[ch,1]-flux[ch,0])/den) if den else float('nan'))
    return {'q':q,'x0':x0,'x1':x1,'eta':eta,'flux':{'ch0_inner':float(flux[0,0]),'ch0_outer':float(flux[0,1]),'ch1_inner':float(flux[1,0]),'ch1_outer':float(flux[1,1])},'delta':delta,'objective_norm':float(math.hypot(*delta)),'audit':{'q_power_max_error':max_qerr,'max_decoded_imaginary':max_nonreal}}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--q',type=float,required=True); ap.add_argument('--out',required=True); ap.add_argument('--eta',type=float,default=-1.0)
    a=ap.parse_args(); nx=61; K=4; dt=.12; steps=700
    r=np.linspace(2.4,62.4,nx); layout=rn.Layout(nx=nx,K=K); M,pi,pj,nrules=rn.compile_sparse_quadratic(layout,r,dt)
    hist=[]; cache={}
    def ev(x0,x1):
        key=(round(float(x0),10),round(float(x1),10))
        if key not in cache:
            rr=run_precompiled(a.q,float(x0),float(x1),a.eta,M,pi,pj,layout,r,dt,steps); cache[key]=rr; hist.append(rr)
        return cache[key]
    best=None
    grid=np.linspace(2.6,3.5,7)
    for x0 in grid:
      for x1 in grid:
        rr=ev(x0,x1)
        if best is None or rr['objective_norm']<best['objective_norm']: best=rr
    def fun(x): return np.array(ev(x[0],x[1])['delta'])
    sol=least_squares(fun,[best['x0'],best['x1']],bounds=([2.4,2.4],[4.2,4.2]),xtol=1e-9,ftol=1e-9,gtol=1e-9,max_nfev=45)
    final=ev(sol.x[0],sol.x[1])
    out={'model':'strict harmonic-only RN coupled simultaneous two-channel balance with two-pulse initial state','persistent_state':'complex harmonic coefficients only','transition':'same fixed sparse quadratic interaction array','search_variables':'initial x0_ch0 and x0_ch1 only; eta fixed as initial-state choice','balance_definition':'delta_ch=(F_outer-F_inner)/(F_outer+F_inner)','parameters':{'q_initial_state':a.q,'eta_initial_state':a.eta,'nx':nx,'K':K,'dt':dt,'steps':steps,'sigma':1.5,'amp_ch0':1e-4},'interaction':{'shape':list(M.shape),'nnz':int(M.nnz),'monomial_pairs':int(len(pi)),'sample_rules':int(nrules)},'coarse_best':best,'optimizer':{'success':bool(sol.success),'message':sol.message,'nfev':int(sol.nfev),'x':[float(v) for v in sol.x],'cost':float(sol.cost)},'final':final,'all_evaluations':hist}
    Path(a.out).write_text(json.dumps(out,indent=2,ensure_ascii=False)); print(json.dumps({'q':a.q,'optimizer':out['optimizer'],'final':final,'evaluations':len(hist)},indent=2,ensure_ascii=False))
if __name__=='__main__': main()
