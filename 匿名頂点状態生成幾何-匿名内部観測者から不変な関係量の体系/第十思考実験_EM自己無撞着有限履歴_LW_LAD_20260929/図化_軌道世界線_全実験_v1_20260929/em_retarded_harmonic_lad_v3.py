#!/usr/bin/env python3
"""Harmonic-primary two-charge retarded Maxwell + covariant Lorentz-Dirac solver.

Flat spacetime, metric (-,+,+,+), units c=1 and 4*pi*eps0=1.

Persistent finite-history state is NOT a time-sample history array.  It is the
full set of DFT/Fourier coefficients of the finite worldline window for each
particle and each component of (r,u,a4).  History samples are decoded only as a
readout operation when the retarded Lienard-Wiechert solver needs them.

The sliding-window write is performed directly in harmonic coefficient space:

    C'_m = exp(+i*2*pi*m/N) * [ C_m + (y - x_oldest)/N ]

where x_oldest = sum_m C_m, and y is the new sample.  Thus the persistent state
is harmonic coefficients plus the scalar current time; no external persistent
history array is used after initialization.

Mutual fields: exact retarded Lienard-Wiechert point-charge fields.
Self-force: covariant Lorentz-Dirac equation
    m a^mu = f_ext^mu + (2/3) q^2 ( d a^mu/dtau - a^2 u^mu )
No dipole radiation-power approximation is used.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np

C=1.0; KE=1.0
ETA=np.diag([-1.0,1.0,1.0,1.0])

def mdot(a,b): return float(a @ ETA @ b)

def normalize_u(u):
    us=u[1:]
    return np.r_[math.sqrt(1.0+float(np.dot(us,us))),us]

def project_a(u,a):
    return a + mdot(u,a)*u

def v_from_u(u): return u[1:]/u[0]

def gamma_from_v(v): return 1.0/math.sqrt(1.0-float(np.dot(v,v)))

def four_accel_to_three(u,a4):
    gam=u[0]; v=v_from_u(u)
    return (a4[1:]-v*a4[0])/(gam*gam)

def force3_to_four(v,F3):
    gam=gamma_from_v(v)
    return np.r_[gam*float(np.dot(F3,v)),gam*F3]

# ---------- harmonic history state ----------

def encode_samples(samples):
    """samples shape (N,D) -> complex coefficients shape (N,D)."""
    return np.fft.fft(samples,axis=0)/samples.shape[0]

def decode_samples(coeff):
    """Exact grid-point readout of the full-mode harmonic state."""
    return np.fft.ifft(coeff*coeff.shape[0],axis=0).real

def oldest_from_coeff(coeff):
    """x_0 = sum_m C_m for DFT convention used here."""
    return np.sum(coeff,axis=0).real

def harmonic_shift_append(coeff,new_sample):
    """Shift finite history left by one and append new_sample, in coefficient space."""
    n=coeff.shape[0]
    old=oldest_from_coeff(coeff)
    phase=np.exp(2j*np.pi*np.arange(n)/n)[:,None]
    return phase*(coeff + (new_sample-old)[None,:]/n)

def read_history(coeff):
    """Readout only; returned array is ephemeral, not persistent state."""
    return decode_samples(coeff)

# ---------- retarded Maxwell ----------

def interp(times,vals,tq):
    if tq<=times[0]: return vals[0].copy()
    if tq>=times[-1]: return vals[-1].copy()
    j=int(np.searchsorted(times,tq)); t0,t1=times[j-1],times[j]
    w=(tq-t0)/(t1-t0)
    return (1-w)*vals[j-1]+w*vals[j]

def retarded_time(times,r_src,r_obs,t_obs,max_iter=80,tol=1e-13):
    lo,hi=times[0],min(t_obs,times[-1]); tr=hi
    for _ in range(max_iter):
        rs=interp(times,r_src,tr)
        new=t_obs-np.linalg.norm(r_obs-rs)
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
    kap=1-float(np.dot(n,beta))
    if abs(kap)<1e-12: raise FloatingPointError('LW kappa singularity')
    Enear=KE*q*(n-beta)/(gam*gam*kap**3*R**2)
    Erad=KE*q*np.cross(n,np.cross(n-beta,acc))/(kap**3*R)
    E=Enear+Erad; B=np.cross(n,E)
    return E,B,tr,Enear,Erad

# ---------- covariant LAD ----------

def derive_readout(coeffs,t_end,dt):
    """Decode harmonic persistent state to ephemeral arrays needed by LW interpolation."""
    n=coeffs.shape[1]
    times=t_end-(n-1)*dt + np.arange(n)*dt
    flat=np.stack([read_history(coeffs[i]) for i in range(2)],axis=0)  # (2,N,11)
    r=flat[:,:,:3]; u=flat[:,:,3:7]; a4=flat[:,:,7:11]
    gam=u[:,:,0]
    v=u[:,:,1:4]/gam[:,:,None]
    a3=(a4[:,:,1:4]-v*a4[:,:,0,None])/(gam[:,:,None]**2)
    return times,r,u,a4,v,a3,flat

def deriv(state,i,p,t_obs,times,r_hist,v_hist,a3_hist):
    r,u,a4=state; j=1-i
    E,B,tr,Enear,Erad=lw_fields(p[j]['q'],r,t_obs,times,r_hist[j],v_hist[j],a3_hist[j])
    v=v_from_u(u); F3=p[i]['q']*(E+np.cross(v,B)); f4=force3_to_four(v,F3)
    tau0=(2.0/3.0)*KE*p[i]['q']**2/p[i]['m']
    if tau0==0:
        adot4=np.zeros(4); a2=0.0
    else:
        a2=mdot(a4,a4)
        adot4=(a4-f4/p[i]['m'])/tau0 + a2*u
    gam=u[0]
    drdt=v; dudt=a4/gam; dadt=adot4/gam
    diag=dict(retarded_time=float(tr),E_near_norm=float(np.linalg.norm(Enear)),
              E_rad_norm=float(np.linalg.norm(Erad)),F3=F3.tolist(),a2=float(a2))
    return (drdt,dudt,dadt),diag

def rk4_one(state,i,p,dt,t_obs,times,rh,vh,ah):
    def add(s,k,fac): return tuple(sx+fac*kx for sx,kx in zip(s,k))
    k1,d1=deriv(state,i,p,t_obs,times,rh,vh,ah)
    k2,_=deriv(add(state,k1,dt/2),i,p,t_obs,times,rh,vh,ah)
    k3,_=deriv(add(state,k2,dt/2),i,p,t_obs,times,rh,vh,ah)
    k4,_=deriv(add(state,k3,dt),i,p,t_obs,times,rh,vh,ah)
    out=tuple(s + dt*(a+2*b+2*c+d)/6 for s,a,b,c,d in zip(state,k1,k2,k3,k4))
    r,u,a4=out; u=normalize_u(u); a4=project_a(u,a4)
    return (r,u,a4),d1

def initial_coeffs(args,p):
    n=args.history
    flat=np.zeros((2,n,11),float)
    flat[0,:,0]=-args.separation/2; flat[1,:,0]=args.separation/2
    for i,vy in enumerate([args.v1,args.v2]):
        vv=np.array([0.,vy,0.]); g=gamma_from_v(vv)
        flat[i,:,3:7]=np.r_[g,g*vv]
    # seed acceleration from mutual instantaneous Lorentz force; constant over prehistory
    for i in range(2):
        j=1-i
        ri=flat[i,-1,:3]; rj=flat[j,-1,:3]
        dist=np.linalg.norm(ri-rj); E=KE*p[j]['q']*(ri-rj)/dist**3
        vi=v_from_u(flat[i,-1,3:7]); F3=p[i]['q']*E; f4=force3_to_four(vi,F3)
        aseed=project_a(flat[i,-1,3:7],f4/p[i]['m'])
        flat[i,:,7:11]=aseed
    coeffs=np.stack([encode_samples(flat[i]) for i in range(2)],axis=0)
    return coeffs,flat

def simulate(args):
    p=[dict(q=args.q1,m=args.m1),dict(q=args.q2,m=args.m2)]
    n=args.history; dt=args.dt
    coeffs,initial_flat=initial_coeffs(args,p)
    t_end=0.0
    diagnostics=[]
    max_write_read_err=0.0
    for step in range(args.steps):
        times,rh,uh,a4h,vh,a3h,flat=derive_readout(coeffs,t_end,dt)
        # current sample is newest grid value read from harmonic state
        states=[]; per=[]
        for i in range(2):
            st=(rh[i,-1].copy(),uh[i,-1].copy(),a4h[i,-1].copy())
            ns,di=rk4_one(st,i,p,dt,t_end,times,rh,vh,a3h)
            states.append(ns); per.append(di)
        t_new=t_end+dt
        # write directly into harmonic state, no persistent history array
        for i,(ri,ui,ai) in enumerate(states):
            sample=np.r_[ri,ui,ai]
            coeffs[i]=harmonic_shift_append(coeffs[i],sample)
        t_end=t_new
        # validation only: read back grid and compare newest to just-written sample
        rb=np.stack([read_history(coeffs[i]) for i in range(2)],axis=0)
        wrerr=max(float(np.max(np.abs(rb[i,-1]-np.r_[states[i][0],states[i][1],states[i][2]]))) for i in range(2))
        max_write_read_err=max(max_write_read_err,wrerr)
        diagnostics.append(dict(step=step+1,time=float(t_end),
            separation=float(np.linalg.norm(rb[1,-1,:3]-rb[0,-1,:3])),
            harmonic_write_read_maxerr=wrerr,particles=per,
            u_norm=[mdot(rb[i,-1,3:7],rb[i,-1,3:7]) for i in range(2)],
            ua_orth=[mdot(rb[i,-1,3:7],rb[i,-1,7:11]) for i in range(2)]))
        if not np.isfinite(coeffs).all(): raise FloatingPointError('nonfinite harmonic state')
    final_rb=np.stack([read_history(coeffs[i]) for i in range(2)],axis=0)
    return dict(parameters=vars(args),persistent_state='full complex harmonic coefficients only',
                harmonic_state_shape=list(coeffs.shape),
                final=dict(r=final_rb[:,-1,:3].tolist(),
                           u=final_rb[:,-1,3:7].tolist(),
                           a4=final_rb[:,-1,7:11].tolist(),
                           v=[v_from_u(final_rb[i,-1,3:7]).tolist() for i in range(2)]),
                max_harmonic_write_read_err=max_write_read_err,
                diagnostics=diagnostics)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--q1',type=float,default=0.002); ap.add_argument('--q2',type=float,default=-0.002)
    ap.add_argument('--m1',type=float,default=1.0); ap.add_argument('--m2',type=float,default=1.0)
    ap.add_argument('--separation',type=float,default=10.0); ap.add_argument('--v1',type=float,default=0.002); ap.add_argument('--v2',type=float,default=-0.002)
    ap.add_argument('--history',type=int,default=512); ap.add_argument('--dt',type=float,default=1e-5); ap.add_argument('--steps',type=int,default=1000)
    ap.add_argument('--out',default='em_harmonic_lad_v3.json')
    a=ap.parse_args(); res=simulate(a)
    Path(a.out).write_text(json.dumps(res,ensure_ascii=False,indent=2))
    d=res['diagnostics']
    print(json.dumps(dict(persistent_state=res['persistent_state'],harmonic_state_shape=res['harmonic_state_shape'],
        final=res['final'],max_harmonic_write_read_err=res['max_harmonic_write_read_err'],last=d[-1]),ensure_ascii=False,indent=2))
if __name__=='__main__': main()
