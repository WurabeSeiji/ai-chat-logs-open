#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""E3-A0 fixed-point v1: minimal reciprocal perturbation of two published 2-value Kepler sectors.

Key point tested:
  raw hierarchical coordinates exchange the local area readout at O(eps), while
  the exact normal modes remain independent 2-value laws.  In the nondegenerate
  limit, the slow/outer normal-mode law parameter shifts only at O(eps^2).

This is a diagnostic candidate, not a derived 3-state law.
"""
import csv, json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def om(u,v): return float(u[0]*v[1]-u[1]*v[0])
def rot(t): return np.array([[math.cos(t),-math.sin(t)],[math.sin(t),math.cos(t)]],float)

def module(n,e,u0):
    th=2*math.pi/n
    T=np.diag([math.sqrt(1-e), math.sqrt(1+e)])
    S=T@rot(th)@np.linalg.inv(T)
    x0=T@np.array([math.cos(u0),math.sin(u0)])
    xm1=np.linalg.inv(S)@x0
    return 2*math.cos(th),xm1,x0

def conic_fit(z):
    r=np.abs(z); A=np.column_stack([np.ones(len(z)),z.real,z.imag])
    c,*_=np.linalg.lstsq(A,r,rcond=None)
    pred=A@c
    resid=float(np.sqrt(np.mean((r-pred)**2))/max(float(np.mean(r)),1e-30))
    ev=np.array([-c[1],-c[2]])
    return resid,float(np.linalg.norm(ev)),float(math.atan2(ev[1],ev[0]))

def simulate(n_in,n_out,eps,u_in,u_out,steps=16000,e_in=.35,e_out=.25):
    ti,im1,i0=module(n_in,e_in,u_in); to,om1,o0=module(n_out,e_out,u_out)
    xi=np.empty((steps+1,2)); xo=np.empty((steps+1,2)); xi[0]=i0; xo[0]=o0
    xi[1]=ti*i0+eps*o0-im1; xo[1]=to*o0+eps*i0-om1
    for k in range(1,steps):
        xi[k+1]=ti*xi[k]+eps*xo[k]-xi[k-1]
        xo[k+1]=to*xo[k]+eps*xi[k]-xo[k-1]
    qi=np.array([om(xi[k],xi[k+1]) for k in range(steps)])
    qo=np.array([om(xo[k],xo[k+1]) for k in range(steps)])
    qt=qi+qo

    # Diagonalize the hierarchy-space coupling. Each eigen-coordinate is again a 2D two-value law.
    H=np.array([[ti,eps],[eps,to]])
    lam,V=np.linalg.eigh(H) # columns eigenvectors, ascending lambda
    # mode[j,k,:] = V[0,j] Xin + V[1,j] Xout
    modes=[]
    for j in range(2): modes.append(V[0,j]*xi+V[1,j]*xo)
    # outer-dominant eigenmode
    jout=int(np.argmax(np.abs(V[1,:])))
    jin=1-jout
    mout=modes[jout]; minm=modes[jin]
    qmo=np.array([om(mout[k],mout[k+1]) for k in range(steps)])
    qmi=np.array([om(minm[k],minm[k+1]) for k in range(steps)])

    # Use full sequences for exact conic diagnostic. For raw external, fit over first 8 unperturbed radial periods
    # so the metric measures local departure rather than long beat accumulation.
    win=max(64,int(8*n_out/2))
    zraw=(xo[:win,0]+1j*xo[:win,1])**2
    zmout=(mout[:,0]+1j*mout[:,1])**2
    raw_resid,raw_e,_=conic_fit(zraw)
    mode_resid,mode_e,_=conic_fit(zmout)

    scale=max(abs(qt[0]),1e-30)
    mix_angle=.5*math.atan2(2*eps,to-ti) if abs(to-ti)+abs(eps)>0 else 0.0
    return {
      'n_in':n_in,'n_out':n_out,'eps':eps,'tau_in':ti,'tau_out':to,
      'delta_tau':to-ti,'lambda_in_eff':float(lam[jin]),'lambda_out_eff':float(lam[jout]),
      'tau_out_shift':float(lam[jout]-to),'tau_in_shift':float(lam[jin]-ti),
      'mix_angle_abs':abs(float(mix_angle)),
      'q_total_drift':float(np.max(np.abs(qt-qt[0]))/scale),
      'q_in_mod':float(np.std(qi)/max(abs(np.mean(qi)),1e-30)),
      'q_out_mod':float(np.std(qo)/max(abs(np.mean(qo)),1e-30)),
      'q_mode_out_drift':float(np.max(np.abs(qmo-qmo[0]))/max(abs(qmo[0]),1e-30)),
      'q_mode_in_drift':float(np.max(np.abs(qmi-qmi[0]))/max(abs(qmi[0]),1e-30)),
      'raw_outer_conic_resid':raw_resid,'raw_outer_ecc':raw_e,
      'mode_outer_conic_resid':mode_resid,'mode_outer_ecc':mode_e,
      'mode_stable':bool(np.max(np.abs(lam))<=2+1e-12),
    }

def pfit(rows,key,epsmax=None):
    pts=[]
    for r in rows:
        if r['eps']<=0: continue
        if epsmax is not None and r['eps']>epsmax: continue
        y=abs(r[key])
        if y>1e-14 and np.isfinite(y): pts.append((r['eps'],y))
    if len(pts)<3:return None
    s,b=np.polyfit(np.log([x for x,y in pts]),np.log([y for x,y in pts]),1)
    return {'slope':float(s),'prefactor':float(math.exp(b)),'n':len(pts)}

def suite(name,nin,nout,epslist,phases,steps):
    raw=[]
    for e in epslist:
      for j,p in enumerate(phases):
        r=simulate(nin,nout,e,.37+p,1.11,steps=steps);r['suite']=name;r['phase']=j;raw.append(r)
    agg=[]
    for e in epslist:
      ss=[r for r in raw if r['eps']==e]
      def med(k):return float(np.median([r[k] for r in ss]))
      def mx(k):return float(np.max([r[k] for r in ss]))
      agg.append({'suite':name,'eps':e,'phases':len(ss),'stable':all(r['mode_stable'] for r in ss),
        'q_total_drift_max':mx('q_total_drift'),'q_in_mod_med':med('q_in_mod'),'q_out_mod_med':med('q_out_mod'),
        'raw_conic_resid_med':med('raw_outer_conic_resid'),'mode_conic_resid_max':mx('mode_outer_conic_resid'),
        'mode_q_drift_max':max(mx('q_mode_out_drift'),mx('q_mode_in_drift')),
        'mix_angle':med('mix_angle_abs'),'tau_out_shift':med('tau_out_shift')})
    small=max(epslist)*.3
    fits={k:pfit(agg,k,small) for k in ['q_in_mod_med','q_out_mod_med','raw_conic_resid_med','mix_angle','tau_out_shift']}
    return raw,agg,fits

def main():
    phases=[0,.43,.91,1.37,1.93,2.51,3.07,3.61]
    eps_h=[0,1e-5,3e-5,1e-4,3e-4,1e-3,2e-3,4e-3,6e-3,8e-3]
    eps_d=[0,1e-5,3e-5,1e-4,3e-4,1e-3,3e-3,1e-2,2e-2]
    rh,ah,fh=suite('hierarchical_nonres_31_127',31,127,eps_h,phases,steps=20000)
    rd,ad,fd=suite('degenerate_31_31',31,31,eps_d,phases,steps=12000)
    out={'experiment':'E3-A0-fixed-point-v1','status':'diagnostic candidate, not derived 3-state law',
         'recurrence':'Xi[k+1]+Xi[k-1]=tau_i Xi[k]+eps Xo[k]; Xo[k+1]+Xo[k-1]=tau_o Xo[k]+eps Xi[k]',
         'hierarchical_fits':fh,'degenerate_fits':fd,'aggregates':ah+ad,'raw':rh+rd}
    jp=HERE/'E3A0_hierarchical_two_kepler_perturbation_full_v1.json'; cp=HERE/'E3A0_hierarchical_two_kepler_perturbation_aggregates_v1.csv'
    json.dump(out,open(jp,'w'),ensure_ascii=False,indent=2)
    cols=list((ah+ad)[0].keys()); f=open(cp,'w',newline='');w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(ah+ad);f.close()
    print(json.dumps({'hierarchical_fits':fh,'degenerate_fits':fd,'hierarchical':ah,'degenerate':ad,'json':str(jp),'csv':str(cp)},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
