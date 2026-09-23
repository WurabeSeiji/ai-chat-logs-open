#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 4 numerical experiment B: generate a gravity orbit only by X_{k+1}=S_k X_k.

The S_k family is area preserving (det S_k=1).  The angular gauge uses the Nth-order
binomial recurrence for Schwarzschild dphi/dchi.  Position readout is the same quadratic
map as Thought Experiment 2: (x,y)=(a^2-b^2,2ab).
"""
import argparse,csv,json,math
from pathlib import Path
import numpy as np

A=30.0; B=27.0

def coeffs(N):
    c=[1.0]
    for m in range(N): c.append(c[-1]*(2*m+1)/(m+1))
    return np.array(c)

def gN(chi,p,e,c):
    u=(3.0+e*np.cos(chi))/p
    # Horner
    y=np.zeros_like(np.asarray(chi,dtype=float))+c[-1]
    for a in c[-2::-1]: y=a+u*y
    return y

def R(th):
    c,s=math.cos(th),math.sin(th); return np.array([[c,-s],[s,c]],float)

def Phi(X):
    a,b=X[...,0],X[...,1]
    return np.column_stack((a*a-b*b,2*a*b))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--order',type=int,default=12); ap.add_argument('--steps-per-radial-orbit',type=int,default=4000); ap.add_argument('--radial-orbits',type=float,default=3.0); args=ap.parse_args()
    out=Path(__file__).resolve().parent/'results'; out.mkdir(exist_ok=True)
    e=math.sqrt(1-(B/A)**2); p=A*(1-e*e); rp=A*(1-e)
    nsteps=int(round(args.steps_per_radial_orbit*args.radial_orbits)); dchi=2*math.pi/args.steps_per_radial_orbit
    chi=np.arange(nsteps+1,dtype=float)*dchi
    c=coeffs(args.order)
    # midpoint quadrature for angle increment, generated from recurrence coefficients
    mid=(chi[:-1]+chi[1:])*0.5
    dphi=gN(mid,p,e,c)*dchi
    phi=np.zeros(nsteps+1); phi[1:]=np.cumsum(dphi)
    r=p/(1+e*np.cos(chi))
    theta=0.5*phi
    X=np.zeros((nsteps+1,2)); X[0]=[math.sqrt(r[0]),0.0]
    det_err=[]; map_err=[]
    for k in range(nsteps):
        rho=math.sqrt(r[k+1]/r[k])
        S=R(theta[k+1]) @ np.diag([rho,1/rho]) @ R(-theta[k])
        X[k+1]=S@X[k]
        det_err.append(abs(np.linalg.det(S)-1.0))
        target=math.sqrt(r[k+1])*np.array([math.cos(theta[k+1]),math.sin(theta[k+1])])
        map_err.append(float(np.linalg.norm(X[k+1]-target)))
    Y=Phi(X); rr=np.hypot(Y[:,0],Y[:,1]); ph=np.unwrap(np.arctan2(Y[:,1],Y[:,0]))
    path=out/f'sn_gravity_C1_N{args.order}.csv'
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.writer(f); w.writerow(['k','chi','x','y','r','phi','X0','X1'])
        for k in range(nsteps+1): w.writerow([k,f'{chi[k]:.15g}',f'{Y[k,0]:.15g}',f'{Y[k,1]:.15g}',f'{rr[k]:.15g}',f'{ph[k]:.15g}',f'{X[k,0]:.15g}',f'{X[k,1]:.15g}'])
    # per-radial-orbit angle and precession
    idx=np.arange(0,nsteps+1,args.steps_per_radial_orbit,dtype=int)
    idx=idx[idx<=nsteps]; peri_phi=ph[idx]
    prec=float(np.mean(np.diff(peri_phi)-2*math.pi)) if len(peri_phi)>1 else None
    meta={'model':'X_{k+1}=S_k X_k; det(S_k)=1; quadratic readout Phi', 'a':A,'b':B,'e':e,'p':p,
          'order_N':args.order,'coefficients':c.tolist(),'steps_per_radial_orbit':args.steps_per_radial_orbit,
          'radial_orbits':args.radial_orbits,'precession_per_orbit_rad':prec,
          'max_det_S_minus_1':float(max(det_err,default=0.0)),'max_state_map_error':float(max(map_err,default=0.0)),
          'output_file':path.name}
    with open(out/f'sn_gravity_C1_N{args.order}_meta.json','w',encoding='utf-8') as f: json.dump(meta,f,ensure_ascii=False,indent=2)
    print(json.dumps(meta,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
