#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 4 numerical experiment A: reproduce the Paper 3 C1 strong-field PN orbit.

Uses the same equations/integrator as generate_paper3_reference_orbit.py, but writes only
one clean relative-orbit time series for direct comparison with the S(n) model.
"""
import csv, json, math
from pathlib import Path
import numpy as np
import generate_paper3_reference_orbit as p3

A=30.0; B=27.0; Q=1.0; N_ORBITS=3.0

def write_csv(path, cols):
    names=list(cols)
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.writer(f); w.writerow(names)
        for row in zip(*(cols[k] for k in names)): w.writerow([f"{float(v):.15g}" for v in row])

def main():
    out=Path(__file__).resolve().parent/'results'; out.mkdir(exist_ok=True)
    ma=Q/(1+Q); mb=1/(1+Q); nu=ma*mb
    e=math.sqrt(1-(B/A)**2); p=A*(1-e*e); T0=2*math.pi*A**1.5
    init,_,info=p3.initial_state(A,B,nu,ma,mb)
    ts, st, dv, stop=p3.integrate(A,B,nu,ma,mb,N_ORBITS*T0,rr=True,init=init)
    # sample densely on uniform coordinate-time grid for comparison/plotting
    tq=np.linspace(0.0,ts[-1],12001)
    y=p3.interpolate(ts,st,dv,tq)
    x,yy,vx,vy,ta,tb=(y[:,i] for i in range(6))
    r=np.hypot(x,yy); phi=np.unwrap(np.arctan2(yy,x))
    write_csv(out/'reference_pn_C1.csv',{'t':tq,'x':x,'y':yy,'r':r,'phi':phi})
    tp,pp=p3.periapsis_passages(tq,x,yy)
    meta={
      'source':'Paper 3 generator copy: Newtonian + 1PN + 2PN conservative + 2.5PN radiation reaction',
      'a':A,'b':B,'q_mass_ratio':Q,'nu':nu,'e0':e,'p0':p,'orbits_requested':N_ORBITS,
      'stop_reason':stop,'initial_condition':info,
      'periapsis_times':tp.tolist(),'periapsis_phi':pp.tolist(),
      'precession_per_orbit_rad': float((pp[-1]-pp[0])/(len(pp)-1)-2*math.pi) if len(pp)>1 else None,
      'reference_file':'reference_pn_C1.csv'
    }
    with open(out/'reference_pn_C1_meta.json','w',encoding='utf-8') as f: json.dump(meta,f,ensure_ascii=False,indent=2)
    print(json.dumps(meta,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
