#!/usr/bin/env python3
import math, json, sys
import numpy as np
sys.path.insert(0, '/mnt/data')
import experiment03_core as e3

# Scale gauges chosen only to keep exponentials in floating-point range.
LAM = {'chi':0.1, 'p':0.01, 'e':1.0, 'phi':0.1, 't':1.0e-5}
KEYS = ('chi','p','e','phi','t')


def enc(vals):
    return np.array([math.exp(LAM[k]*vals[i]) for i,k in enumerate(KEYS)], dtype=float)


def dec(Z):
    return tuple(math.log(float(Z[i]))/LAM[k] for i,k in enumerate(KEYS))


def multiplicative_step_oracle(Z,dchi,nu,c):
    # Audit/proof-of-equivalence shell only: log decode is used to evaluate
    # the existing effective PN/RK factor generator. Persistent state update
    # itself is purely multiplicative. A strict physical implementation must
    # next replace this decoder-dependent factor generator.
    chi,p,e,phi,t = dec(Z)
    p1,e1,dt = e3.rk4_elements(chi,p,e,dchi,nu)
    cm=chi+0.5*dchi
    pm=0.5*(p+p1); em=0.5*(e+e1)
    dphi=e3.gN_scalar(cm,pm,em,c)*dchi
    factors=np.array([
        math.exp(LAM['chi']*dchi),
        math.exp(LAM['p']*(p1-p)),
        math.exp(LAM['e']*(e1-e)),
        math.exp(LAM['phi']*dphi),
        math.exp(LAM['t']*dt),
    ],dtype=float)
    return Z*factors


def additive_step(state,dchi,nu,c):
    chi,p,e,phi,t=state
    p1,e1,dt=e3.rk4_elements(chi,p,e,dchi,nu)
    cm=chi+0.5*dchi
    pm=0.5*(p+p1); em=0.5*(e+e1)
    dphi=e3.gN_scalar(cm,pm,em,c)*dchi
    return (chi+dchi,p1,e1,phi+dphi,t+dt)


def run_case(p0,e0,q=1.0,order=12,steps=4000,radial_orbits=3.0):
    nsteps=int(round(steps*radial_orbits)); dchi=2*math.pi/steps
    c=e3.coeffs(order); ma=q/(1+q); mb=1/(1+q); nu=ma*mb
    a=(0.0,p0,e0,0.0,0.0); Z=enc(a)
    worst={k:0.0 for k in KEYS}
    worst_r=worst_x=worst_y=0.0
    for _ in range(nsteps):
        a=additive_step(a,dchi,nu,c)
        Z=multiplicative_step_oracle(Z,dchi,nu,c)
        b=dec(Z)
        for i,k in enumerate(KEYS): worst[k]=max(worst[k],abs(a[i]-b[i]))
        chia,pa,ea,phia,_=a; chib,pb,eb,phib,_=b
        ra=pa/(1+ea*math.cos(chia)); rb=pb/(1+eb*math.cos(chib))
        xa=ra*math.cos(phia); ya=ra*math.sin(phia)
        xb=rb*math.cos(phib); yb=rb*math.sin(phib)
        worst_r=max(worst_r,abs(ra-rb)); worst_x=max(worst_x,abs(xa-xb)); worst_y=max(worst_y,abs(ya-yb))
    return {'p0':p0,'e0':e0,'q':q,'max_abs_state':worst,'max_abs_r':worst_r,'max_abs_x':worst_x,'max_abs_y':worst_y}


def main():
    cases=[(p,e,1.0) for p in (24.3,60.0,120.0) for e in (0.20,0.45,0.55)]
    cases += [(60.0,0.45,q) for q in (2.0,4.0,10.0)]
    rows=[run_case(*x) for x in cases]
    summary={'lambda':LAM,'cases':rows,'worst':{}}
    for k in KEYS:
        summary['worst'][k]=max(r['max_abs_state'][k] for r in rows)
    for k in ('r','x','y'):
        summary['worst'][k]=max(r['max_abs_'+k] for r in rows)
    print(json.dumps(summary,indent=2))
    with open('/mnt/data/five_state_multiplicative_audit.json','w') as f: json.dump(summary,f,indent=2)

if __name__=='__main__': main()
