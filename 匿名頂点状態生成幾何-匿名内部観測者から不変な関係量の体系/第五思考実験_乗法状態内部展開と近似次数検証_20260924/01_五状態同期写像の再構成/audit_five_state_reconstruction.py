#!/usr/bin/env python3
import math, json
import numpy as np
import sys
sys.path.insert(0, '/mnt/data')
import experiment03_core as e3


def original_run(p0,e0,qmass=1.0,order=12,steps=4000,radial_orbits=3.0):
    nsteps=int(round(steps*radial_orbits))
    dchi=2.0*math.pi/steps
    c=e3.coeffs(order)
    ma=qmass/(1+qmass); mb=1/(1+qmass); nu=ma*mb
    chi=np.arange(nsteps+1,dtype=float)*dchi
    pp=np.empty(nsteps+1); ee=np.empty(nsteps+1); phi=np.empty(nsteps+1); tosc=np.empty(nsteps+1)
    pp[0],ee[0],phi[0],tosc[0]=p0,e0,0.0,0.0
    for k in range(nsteps):
        p1,e1,dt=e3.rk4_elements(chi[k],pp[k],ee[k],dchi,nu)
        pp[k+1],ee[k+1]=p1,e1
        tosc[k+1]=tosc[k]+dt
        cm=chi[k]+0.5*dchi
        pm=0.5*(pp[k]+pp[k+1]); em=0.5*(ee[k]+ee[k+1])
        phi[k+1]=phi[k]+e3.gN_scalar(cm,pm,em,c)*dchi
    r=pp/(1.0+ee*np.cos(chi)); theta=0.5*phi
    X=np.zeros((nsteps+1,2)); X[0]=[math.sqrt(r[0]),0.0]
    for k in range(nsteps):
        rho=math.sqrt(r[k+1]/r[k])
        S=e3.rot(theta[k+1]) @ np.diag([rho,1.0/rho]) @ e3.rot(-theta[k])
        X[k+1]=S@X[k]
    Y=e3.Phi(X)
    return {'chi':chi,'p':pp,'e':ee,'phi':phi,'t':tosc,'r':r,'X':X,'Y':Y}


def step5(state,dchi,nu,c):
    # ONLY persistent input: chi,p,e,phi,t. All other values are same-step temporaries.
    chi,p,e,phi,t=state
    p1,e1,dt=e3.rk4_elements(chi,p,e,dchi,nu)
    cm=chi+0.5*dchi
    pm=0.5*(p+p1); em=0.5*(e+e1)
    dphi=e3.gN_scalar(cm,pm,em,c)*dchi
    return (chi+dchi, p1, e1, phi+dphi, t+dt)


def derive_from_state(state):
    chi,p,e,phi,t=state
    r=p/(1.0+e*math.cos(chi))
    theta=0.5*phi
    X=np.array([math.sqrt(r)*math.cos(theta), math.sqrt(r)*math.sin(theta)],dtype=float)
    Y=np.array([X[0]*X[0]-X[1]*X[1], 2.0*X[0]*X[1]],dtype=float)
    return r,X,Y


def five_state_run(p0,e0,qmass=1.0,order=12,steps=4000,radial_orbits=3.0):
    nsteps=int(round(steps*radial_orbits)); dchi=2.0*math.pi/steps; c=e3.coeffs(order)
    ma=qmass/(1+qmass); mb=1/(1+qmass); nu=ma*mb
    state=(0.0,p0,e0,0.0,0.0)
    # History below is write-only logging; it is never read by step5.
    H=np.empty((nsteps+1,5)); R=np.empty(nsteps+1); X=np.empty((nsteps+1,2)); Y=np.empty((nsteps+1,2))
    H[0]=state; R[0],X[0],Y[0]=derive_from_state(state)
    for k in range(nsteps):
        state=step5(state,dchi,nu,c)
        H[k+1]=state
        R[k+1],X[k+1],Y[k+1]=derive_from_state(state)
    return {'chi':H[:,0],'p':H[:,1],'e':H[:,2],'phi':H[:,3],'t':H[:,4],'r':R,'X':X,'Y':Y}


def compare_case(p,e,q,steps=4000,radial_orbits=3.0):
    a=original_run(p,e,q,steps=steps,radial_orbits=radial_orbits)
    b=five_state_run(p,e,q,steps=steps,radial_orbits=radial_orbits)
    out={'p0':p,'e0':e,'q':q}
    for key in ['chi','p','e','phi','t','r','X','Y']:
        diff=np.asarray(a[key])-np.asarray(b[key])
        out[f'max_abs_{key}']=float(np.max(np.abs(diff)))
        out[f'rms_{key}']=float(np.sqrt(np.mean(diff*diff)))
    # compare observables reconstructed directly from 5-state vs original X-recursion Phi(X)
    out['max_abs_xy']=out['max_abs_Y']
    out['max_abs_r']=out['max_abs_r']
    return out

if __name__=='__main__':
    cases=[]
    for p in (24.3,60.0,120.0):
        for e in (0.20,0.45,0.55):
            cases.append((p,e,1.0))
    for q in (2.0,4.0,10.0):
        cases.append((60.0,0.45,q))
    results=[compare_case(*c) for c in cases]
    summary={
      'cases':results,
      'worst':{k:max(r[k] for r in results) for k in results[0] if k.startswith('max_abs_')}
    }
    print(json.dumps(summary,indent=2))
    with open('/mnt/data/five_state_audit_results.json','w') as f: json.dump(summary,f,indent=2)
