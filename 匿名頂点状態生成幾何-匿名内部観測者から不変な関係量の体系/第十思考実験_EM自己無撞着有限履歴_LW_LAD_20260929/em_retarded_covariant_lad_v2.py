#!/usr/bin/env python3
"""Two-charge retarded Maxwell + covariant Lorentz-Dirac solver scaffold.

Flat spacetime, metric (-,+,+,+), units c=1 and 4*pi*eps0=1.
Mutual fields: exact Lienard-Wiechert retarded fields for point charges.
Self-force: covariant Lorentz-Dirac equation
    m a^mu = f_ext^mu + (2/3) q^2 ( d a^mu/dtau - a^2 u^mu )
which is evolved as a first-order system in (x,u,a):
    du^mu/dtau = a^mu
    da^mu/dtau = (a^mu - f_ext^mu/m)/tau0 + a^2 u^mu
    tau0 = 2 q^2/(3m)

This is the exact classical point-charge LAD equation once the renormalized point-charge
model is accepted; it can exhibit runaway/preacceleration pathologies. No LL reduction is
used in the primary dynamics.

A full Fourier encoding of the finite history window is computed as an exact reversible
representation layer; no dipole radiation-power formula is used.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np

C=1.0; KE=1.0
ETA=np.diag([-1.0,1.0,1.0,1.0])

def mdot(a,b): return float(a @ ETA @ b)

def normalize_u(u):
    # derive physical gamma from spatial proper velocity, preserving future direction
    us=u[1:]
    u0=math.sqrt(1.0+float(np.dot(us,us)))
    return np.r_[u0,us]

def project_a(u,a):
    # enforce u.a=0 for u.u=-1: a_perp = a + (u.a) u
    return a + mdot(u,a)*u

def v_from_u(u): return u[1:]/u[0]

def gamma_from_v(v): return 1.0/math.sqrt(1.0-float(np.dot(v,v)))

def interp(times,vals,tq):
    if tq<=times[0]: return vals[0].copy()
    if tq>=times[-1]: return vals[-1].copy()
    j=int(np.searchsorted(times,tq)); t0,t1=times[j-1],times[j]
    w=(tq-t0)/(t1-t0); return (1-w)*vals[j-1]+w*vals[j]

def retarded_time(times,r_src,r_obs,t_obs,max_iter=80,tol=1e-13):
    lo,hi=times[0],min(t_obs,times[-1]); tr=hi
    for _ in range(max_iter):
        rs=interp(times,r_src,tr); new=t_obs-np.linalg.norm(r_obs-rs)
        new=min(max(new,lo),hi)
        if abs(new-tr)<tol: return new
        tr=new
    return tr

def lw_fields(q,r_obs,t_obs,times,r_src,v_src,a3_src):
    tr=retarded_time(times,r_src,r_obs,t_obs)
    rs=interp(times,r_src,tr); vs=interp(times,v_src,tr); acc=interp(times,a3_src,tr)
    Rv=r_obs-rs; R=float(np.linalg.norm(Rv))
    if R<1e-11: raise FloatingPointError('mutual LW singularity')
    n=Rv/R; beta=vs; b2=float(np.dot(beta,beta)); gam=1/math.sqrt(max(1e-30,1-b2))
    kap=1-float(np.dot(n,beta));
    if abs(kap)<1e-12: raise FloatingPointError('LW kappa singularity')
    bdot=acc
    Enear=KE*q*(n-beta)/(gam*gam*kap**3*R**2)
    Erad=KE*q*np.cross(n,np.cross(n-beta,bdot))/(kap**3*R)
    E=Enear+Erad; B=np.cross(n,E)
    return E,B,tr,Enear,Erad

def force3_to_four(v,F3):
    gam=gamma_from_v(v)
    return np.r_[gam*float(np.dot(F3,v)), gam*F3]

def four_accel_to_three(u,a4):
    # a4 spatial = gamma^2 a3 + gamma^4(v.a3)v ; invert via known identities
    gam=u[0]; v=v_from_u(u)
    # a3 = (a4_sp - v*a4_0)/gamma^2
    return (a4[1:]-v*a4[0])/(gam*gam)

def fourier_encode(samples): return np.fft.fft(samples,axis=0)/samples.shape[0]
def fourier_decode(c): return np.fft.ifft(c*c.shape[0],axis=0).real

def deriv(state, i, p, times, r_hist, v_hist, a3_hist):
    # state=(r,u,a4), evaluate mutual retarded external force and exact LAD ODE
    r,u,a4=state
    j=1-i
    E,B,tr,Enear,Erad=lw_fields(p[j]['q'],r,times[-1],times,r_hist[j],v_hist[j],a3_hist[j])
    v=v_from_u(u); F3=p[i]['q']*(E+np.cross(v,B)); f4=force3_to_four(v,F3)
    tau0=(2.0/3.0)*KE*p[i]['q']**2/p[i]['m']
    if tau0==0: adot4=np.zeros(4)
    else:
        a2=mdot(a4,a4)
        adot4=(a4-f4/p[i]['m'])/tau0 + a2*u
    gamma=u[0]
    # convert proper-time derivatives to coordinate-time derivatives
    drdt=v
    dudt=a4/gamma
    dadt=adot4/gamma
    diag=dict(retarded_time=tr,E_near_norm=float(np.linalg.norm(Enear)),E_rad_norm=float(np.linalg.norm(Erad)),F3=F3.tolist(),a2=a2 if tau0 else 0.0)
    return (drdt,dudt,dadt),diag

def rk4_one(state,i,p,dt,times,rh,vh,ah):
    def add(s,k,fac): return (s[0]+fac*k[0],s[1]+fac*k[1],s[2]+fac*k[2])
    k1,d1=deriv(state,i,p,times,rh,vh,ah)
    k2,_=deriv(add(state,k1,dt/2),i,p,times,rh,vh,ah)
    k3,_=deriv(add(state,k2,dt/2),i,p,times,rh,vh,ah)
    k4,_=deriv(add(state,k3,dt),i,p,times,rh,vh,ah)
    out=tuple(s + dt*(a+2*b+2*c+d)/6 for s,a,b,c,d in zip(state,k1,k2,k3,k4))
    r,u,a4=out; u=normalize_u(u); a4=project_a(u,a4)
    return (r,u,a4),d1

def simulate(args):
    p=[dict(q=args.q1,m=args.m1),dict(q=args.q2,m=args.m2)]
    n=args.history; dt=args.dt; times=np.linspace(-(n-1)*dt,0,n)
    r=np.zeros((2,n,3)); u=np.zeros((2,n,4)); a4=np.zeros((2,n,4)); a3=np.zeros((2,n,3)); v=np.zeros((2,n,3))
    r[0,:,0]=-args.separation/2; r[1,:,0]=args.separation/2
    for i,vy in enumerate([args.v1,args.v2]):
        vv=np.array([0.,vy,0.]); g=gamma_from_v(vv); u[i,:,:]=np.r_[g,g*vv]; v[i,:,:]=vv
    # initialize acceleration history from mutual Lorentz force only (physical branch seed)
    for i in range(2):
        j=1-i; dist=np.linalg.norm(r[i,-1]-r[j,-1]); E=KE*p[j]['q']*(r[i,-1]-r[j,-1])/dist**3
        F3=p[i]['q']*E; f4=force3_to_four(v[i,-1],F3); a4[i,:,:]=project_a(u[i,-1],f4/p[i]['m']); a3[i,:,:]=four_accel_to_three(u[i,-1],a4[i,-1])
    diagnostics=[]
    for step in range(args.steps):
        states=[]; per=[]
        for i in range(2):
            st=(r[i,-1].copy(),u[i,-1].copy(),a4[i,-1].copy())
            ns,di=rk4_one(st,i,p,dt,times,r,v,a3); states.append(ns); per.append(di)
        tnew=times[-1]+dt; times=np.r_[times[1:],tnew]
        for i,(ri,ui,ai) in enumerate(states):
            r[i]=np.vstack([r[i,1:],ri]); u[i]=np.vstack([u[i,1:],ui]); a4[i]=np.vstack([a4[i,1:],ai])
            vi=v_from_u(ui); a3i=four_accel_to_three(ui,ai)
            v[i]=np.vstack([v[i,1:],vi]); a3[i]=np.vstack([a3[i,1:],a3i])
        # Fourier/harmonic representation of complete finite worldline state window
        flat=np.concatenate([r.transpose(1,0,2).reshape(n,6),u.transpose(1,0,2).reshape(n,8),a4.transpose(1,0,2).reshape(n,8)],axis=1)
        cc=fourier_encode(flat); back=fourier_decode(cc); ferr=float(np.max(np.abs(back-flat)))
        diagnostics.append(dict(step=step+1,time=float(tnew),separation=float(np.linalg.norm(r[1,-1]-r[0,-1])),harmonic_reconstruction_maxerr=ferr,particles=per,
                                u_norm=[mdot(u[i,-1],u[i,-1]) for i in range(2)],ua_orth=[mdot(u[i,-1],a4[i,-1]) for i in range(2)]))
        if not np.isfinite(flat).all(): raise FloatingPointError('nonfinite state')
    return dict(parameters=vars(args),final=dict(r=r[:,-1].tolist(),v=v[:,-1].tolist(),u=u[:,-1].tolist(),a4=a4[:,-1].tolist()),diagnostics=diagnostics)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--q1',type=float,default=0.002); ap.add_argument('--q2',type=float,default=-0.002)
    ap.add_argument('--m1',type=float,default=1.0); ap.add_argument('--m2',type=float,default=1.0)
    ap.add_argument('--separation',type=float,default=10.0); ap.add_argument('--v1',type=float,default=0.002); ap.add_argument('--v2',type=float,default=-0.002)
    ap.add_argument('--history',type=int,default=512); ap.add_argument('--dt',type=float,default=1e-5); ap.add_argument('--steps',type=int,default=1000)
    ap.add_argument('--out',default='em_covariant_lad_result.json')
    a=ap.parse_args(); res=simulate(a); Path(a.out).write_text(json.dumps(res,ensure_ascii=False,indent=2));
    d=res['diagnostics']; print(json.dumps(dict(final=res['final'],max_harmonic_err=max(x['harmonic_reconstruction_maxerr'] for x in d),last=d[-1]),ensure_ascii=False,indent=2))
if __name__=='__main__': main()
