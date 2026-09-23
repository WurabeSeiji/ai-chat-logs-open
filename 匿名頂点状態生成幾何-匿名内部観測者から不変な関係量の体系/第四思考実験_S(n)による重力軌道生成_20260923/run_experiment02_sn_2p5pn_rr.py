#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 4 experiment 02: S(n) gravity orbit with local 2.5PN radiation reaction.

Conservative geometry:
    r = p/(1+e cos chi)
    dphi/dchi = sum_{m=0}^N c_m ((3+e cos chi)/p)^m
    c_0=1, c_{m+1}=((2m+1)/(m+1)) c_m

Dissipative update:
    Use the same local 2.5PN radiation-reaction acceleration term as Paper 3,
    evaluated on the Newtonian osculating ellipse.  Because a_RR is 2.5PN,
    Newtonian orbital relations are sufficient for this leading dissipative update;
    conservative corrections to that evaluation enter at higher PN order.

The osculating elements p,e are integrated versus chi by RK4.  The internal state
is then advanced with the same det(S_k)=1 map used in experiment 01.
"""
import argparse, csv, json, math
from pathlib import Path
import numpy as np

A=30.0; B=27.0; Q=1.0


def coeffs(N):
    c=[1.0]
    for m in range(N):
        c.append(c[-1]*(2*m+1)/(m+1))
    return np.array(c,dtype=float)


def gN_scalar(chi,p,e,c):
    u=(3.0+e*math.cos(chi))/p
    y=float(c[-1])
    for a in c[-2::-1]:
        y=float(a)+u*y
    return y


def rot(th):
    c,s=math.cos(th),math.sin(th)
    return np.array([[c,-s],[s,c]],float)


def Phi(X):
    a,b=X[...,0],X[...,1]
    return np.column_stack((a*a-b*b,2*a*b))


def rr_rates_chi(chi,p,e,nu):
    """Return dp/dchi, de/dchi, dt/dchi from local 2.5PN RR term.

    Newtonian osculating relations (G=c=M=1):
      r=p/(1+e cos chi), h=sqrt(p),
      rdot=e sin chi/sqrt(p), vt=(1+e cos chi)/sqrt(p),
      dt/dchi=r^2/h.
    Paper-3 RR acceleration:
      a_RR = C (an*n + bv*v), C=(8/5) nu / r^3,
      an=rdot(18 v^2 + 2/(3r) - 25 rdot^2),
      bv=-(6 v^2 - 2/r - 15 rdot^2).
    """
    ce=math.cos(chi); se=math.sin(chi)
    den=1.0+e*ce
    r=p/den
    h=math.sqrt(p)
    rdot=e*se/h
    vt=den/h
    v2=rdot*rdot+vt*vt
    invr=1.0/r
    C=1.6*nu*invr**3
    an=rdot*(18.0*v2+(2.0/3.0)*invr-25.0*rdot*rdot)
    bv=-(6.0*v2-2.0*invr-15.0*rdot*rdot)
    Edot=C*(an*rdot+bv*v2)
    pdot=2.0*C*bv*p
    E=(e*e-1.0)/(2.0*p)
    edot=(p*Edot+E*pdot)/e
    dt_dchi=r*r/h
    return pdot*dt_dchi, edot*dt_dchi, dt_dchi, {
        'r':r,'rdot':rdot,'vt':vt,'v2':v2,'Edot':Edot,'pdot':pdot,'edot':edot,
        'C':C,'an':an,'bv':bv
    }


def rk4_elements(chi,p,e,hchi,nu):
    def f(x,pp,ee):
        dp,de,dt,_=rr_rates_chi(x,pp,ee,nu)
        return np.array([dp,de,dt],float)
    y=np.array([p,e,0.0],float)
    k1=f(chi,y[0],y[1])
    y2=y+0.5*hchi*k1
    k2=f(chi+0.5*hchi,y2[0],y2[1])
    y3=y+0.5*hchi*k2
    k3=f(chi+0.5*hchi,y3[0],y3[1])
    y4=y+hchi*k3
    k4=f(chi+hchi,y4[0],y4[1])
    dy=(hchi/6.0)*(k1+2*k2+2*k3+k4)
    return p+dy[0], e+dy[1], dy[2]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--order',type=int,default=12)
    ap.add_argument('--steps-per-radial-orbit',type=int,default=4000)
    ap.add_argument('--radial-orbits',type=float,default=3.0)
    args=ap.parse_args()
    out=Path(__file__).resolve().parent/'results'; out.mkdir(exist_ok=True)

    ma=Q/(1+Q); mb=1/(1+Q); nu=ma*mb
    e0=math.sqrt(1-(B/A)**2); p0=A*(1-e0*e0)
    nsteps=int(round(args.steps_per_radial_orbit*args.radial_orbits))
    dchi=2*math.pi/args.steps_per_radial_orbit
    c=coeffs(args.order)

    chi=np.arange(nsteps+1,dtype=float)*dchi
    p=np.empty(nsteps+1); e=np.empty(nsteps+1); phi=np.empty(nsteps+1); t=np.empty(nsteps+1)
    p[0]=p0; e[0]=e0; phi[0]=0.0; t[0]=0.0

    # Evolve p,e under local 2.5PN RR; evolve phi with conservative S(n) angular gauge.
    for k in range(nsteps):
        pk1,ek1,dt=rk4_elements(chi[k],p[k],e[k],dchi,nu)
        if not (pk1>0.0 and 0.0<ek1<1.0):
            raise RuntimeError(f'Invalid osculating elements at step {k}: p={pk1}, e={ek1}')
        p[k+1]=pk1; e[k+1]=ek1; t[k+1]=t[k]+dt
        # midpoint using averaged evolving elements
        cm=chi[k]+0.5*dchi; pm=0.5*(p[k]+p[k+1]); em=0.5*(e[k]+e[k+1])
        phi[k+1]=phi[k]+gN_scalar(cm,pm,em,c)*dchi

    r=p/(1.0+e*np.cos(chi))
    theta=0.5*phi
    X=np.zeros((nsteps+1,2)); X[0]=[math.sqrt(r[0]),0.0]
    det_err=[]; map_err=[]
    for k in range(nsteps):
        rho=math.sqrt(r[k+1]/r[k])
        S=rot(theta[k+1]) @ np.diag([rho,1.0/rho]) @ rot(-theta[k])
        X[k+1]=S@X[k]
        det_err.append(abs(np.linalg.det(S)-1.0))
        target=math.sqrt(r[k+1])*np.array([math.cos(theta[k+1]),math.sin(theta[k+1])])
        map_err.append(float(np.linalg.norm(X[k+1]-target)))

    Y=Phi(X); rr=np.hypot(Y[:,0],Y[:,1]); ph=np.unwrap(np.arctan2(Y[:,1],Y[:,0]))
    stem=f'experiment02_sn_2p5pn_rr_C1_N{args.order}'
    path=out/f'{stem}.csv'
    with open(path,'w',encoding='utf-8',newline='') as f:
        w=csv.writer(f); w.writerow(['k','chi','t_osc','p','e','x','y','r','phi','X0','X1'])
        for k in range(nsteps+1):
            w.writerow([k,f'{chi[k]:.15g}',f'{t[k]:.15g}',f'{p[k]:.15g}',f'{e[k]:.15g}',
                        f'{Y[k,0]:.15g}',f'{Y[k,1]:.15g}',f'{rr[k]:.15g}',f'{ph[k]:.15g}',
                        f'{X[k,0]:.15g}',f'{X[k,1]:.15g}'])

    idx=np.arange(0,nsteps+1,args.steps_per_radial_orbit,dtype=int); idx=idx[idx<=nsteps]
    peri_phi=ph[idx]
    prec=float(np.mean(np.diff(peri_phi)-2*math.pi)) if len(peri_phi)>1 else None
    meta={
        'experiment':'02',
        'model':'S(n) conservative Schwarzschild angular gauge + local 2.5PN RR osculating-element recurrence',
        'a0':A,'b0':B,'q_mass_ratio':Q,'nu':nu,'e0':e0,'p0':p0,
        'order_N':args.order,'coefficients':c.tolist(),
        'steps_per_radial_orbit':args.steps_per_radial_orbit,'radial_orbits':args.radial_orbits,
        'rr_source':'same 2.5PN local radiation-reaction acceleration term as Paper 3 (Kidder 1995 Eq. 2.2 implementation)',
        'rr_element_update':'Newtonian osculating p=h^2, e^2=1+2Ep; RK4 in chi',
        'p_final':float(p[-1]),'e_final':float(e[-1]),'r_peri_initial':float(r[0]),
        'nominal_r_at_final_peri_phase':float(r[-1]),
        'precession_per_radial_cycle_rad':prec,
        'max_det_S_minus_1':float(max(det_err,default=0.0)),
        'max_state_map_error':float(max(map_err,default=0.0)),
        'output_file':path.name,
    }
    with open(out/f'{stem}_meta.json','w',encoding='utf-8') as f:
        json.dump(meta,f,ensure_ascii=False,indent=2)
    print(json.dumps(meta,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
