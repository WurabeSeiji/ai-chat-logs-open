#!/usr/bin/env python3
import math
import numpy as np

H = 2.0*math.pi/4000.0
ORDER=12

def coeffs(N):
    c=[1.0]
    for m in range(N):
        c.append(c[-1]*(2*m+1)/(m+1))
    return np.asarray(c,float)
C=coeffs(ORDER)

def gN(x,p,e):
    u=(3.0+e*math.cos(x))/p
    y=float(C[-1])
    for a in C[-2::-1]:
        y=float(a)+u*y
    return y

def rr(x,p,e,nu):
    ce,se=math.cos(x),math.sin(x)
    den=1.0+e*ce
    r=p/den
    h=math.sqrt(p)
    rdot=e*se/h
    vt=den/h
    v2=rdot*rdot+vt*vt
    invr=1.0/r
    Crr=1.6*nu*invr**3
    an=rdot*(18.0*v2+(2.0/3.0)*invr-25.0*rdot*rdot)
    bv=-(6.0*v2-2.0*invr-15.0*rdot*rdot)
    Edot=Crr*(an*rdot+bv*v2)
    pdot=2.0*Crr*bv*p
    E=(e*e-1.0)/(2.0*p)
    edot=(p*Edot+E*pdot)/e
    dt_dchi=r*r/h
    return np.array([pdot*dt_dchi,edot*dt_dchi,dt_dchi],float)

def ref(s,nu):
    x,p,e,phi,t=map(float,s)
    y=np.array([p,e,0.0])
    k1=rr(x,y[0],y[1],nu)
    y2=y+0.5*H*k1
    k2=rr(x+0.5*H,y2[0],y2[1],nu)
    y3=y+0.5*H*k2
    k3=rr(x+0.5*H,y3[0],y3[1],nu)
    y4=y+H*k3
    k4=rr(x+H,y4[0],y4[1],nu)
    dy=(H/6.0)*(k1+2*k2+2*k3+k4)
    p1,e1,dt=p+dy[0],e+dy[1],dy[2]
    cm=x+0.5*H
    pm=0.5*(p+p1); em=0.5*(e+e1)
    phi1=phi+gN(cm,pm,em)*H
    return np.array([x+H,p1,e1,phi1,t+dt])

def inline(s,nu):
    x,p,e,phi,t=map(float,s)
    # Entire RK4 update is a function of old x,p,e only.
    k1=rr(x,p,e,nu)
    k2=rr(x+0.5*H,p+0.5*H*k1[0],e+0.5*H*k1[1],nu)
    k3=rr(x+0.5*H,p+0.5*H*k2[0],e+0.5*H*k2[1],nu)
    k4=rr(x+H,p+H*k3[0],e+H*k3[1],nu)
    d=(H/6.0)*(k1+2*k2+2*k3+k4)
    dphi=H*gN(x+0.5*H,p+0.5*d[0],e+0.5*d[1])
    return np.array([x+H,p+d[0],e+d[1],phi+dphi,t+d[2]])

def main():
    rng=np.random.default_rng(20260924)
    max_abs=np.zeros(5)
    max_rel=np.zeros(5)
    worst=None
    for _ in range(100000):
        q=float(np.exp(rng.uniform(math.log(1.0),math.log(10.0))))
        nu=q/(1+q)**2
        s=np.array([
            rng.uniform(0,6*math.pi),
            rng.uniform(20,150),
            rng.uniform(0.1,0.65),
            rng.uniform(-10,30),
            rng.uniform(0,5000),
        ])
        a=ref(s,nu); b=inline(s,nu)
        d=np.abs(a-b)
        r=d/np.maximum(1.0,np.abs(a))
        if d.max()>max_abs.max(): worst=(s.copy(),nu,a.copy(),b.copy(),d.copy())
        max_abs=np.maximum(max_abs,d); max_rel=np.maximum(max_rel,r)
    print('max_abs',max_abs)
    print('max_rel',max_rel)
    print('worst',worst)

    # Pure algebraic exponential identity test, no log decode involved in equality.
    lam=np.array([1e-3,1e-3,1e-2,1e-3,1e-6])
    max_exp=np.zeros(5)
    for _ in range(100000):
        q=float(np.exp(rng.uniform(math.log(1.0),math.log(10.0))))
        nu=q/(1+q)**2
        s=np.array([rng.uniform(0,6*math.pi),rng.uniform(20,150),rng.uniform(.1,.65),rng.uniform(-10,30),rng.uniform(0,5000)])
        nxt=inline(s,nu)
        delta=nxt-s
        lhs=np.exp(lam*nxt)
        rhs=np.exp(lam*s)*np.exp(lam*delta)
        max_exp=np.maximum(max_exp,np.abs(lhs-rhs))
    print('exp_identity_max_abs',max_exp)

if __name__=='__main__': main()
