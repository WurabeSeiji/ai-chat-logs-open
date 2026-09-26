#!/usr/bin/env python3
import math, json
import numpy as np

P_VALUES=(24.3,60.0,120.0)
E_VALUES=(0.20,0.45,0.55)
Q_VALUES=(1.0,2.0,4.0,10.0)
ORDER=12
STEPS=4000
ORBITS=3.0
H=2.0*math.pi/STEPS


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


def increments_current_only(s,nu):
    # IMPORTANT: every RHS below is a function only of the current five-state s.
    x,p,e,phi,t = map(float,s)
    k1=rr(x,p,e,nu)
    k2=rr(x+0.5*H,
          p+0.5*H*k1[0],
          e+0.5*H*k1[1],nu)
    k3=rr(x+0.5*H,
          p+0.5*H*k2[0],
          e+0.5*H*k2[1],nu)
    k4=rr(x+H,
          p+H*k3[0],
          e+H*k3[1],nu)
    d=(H/6.0)*(k1+2.0*k2+2.0*k3+k4)
    # Original code used p_{n+1},e_{n+1} in the midpoint.
    # Substitution gives midpoint directly from current state via d(s):
    pm=p+0.5*d[0]
    em=e+0.5*d[1]
    dphi=gN(x+0.5*H,pm,em)*H
    return np.array([H,d[0],d[1],dphi,d[2]],float)


def step_inline(s,nu):
    # single synchronous array update: no component of s_next is used on the RHS
    return s + increments_current_only(s,nu)


def step_reference_sequential(s,nu):
    # literal source-code semantics, but from the same current state x rather than a pre-built chi array
    x,p,e,phi,t=map(float,s)
    k1=rr(x,p,e,nu)
    y2=np.array([p,e,0.0])+0.5*H*k1
    k2=rr(x+0.5*H,y2[0],y2[1],nu)
    y3=np.array([p,e,0.0])+0.5*H*k2
    k3=rr(x+0.5*H,y3[0],y3[1],nu)
    y4=np.array([p,e,0.0])+H*k3
    k4=rr(x+H,y4[0],y4[1],nu)
    d=(H/6.0)*(k1+2*k2+2*k3+k4)
    p1=p+d[0]; e1=e+d[1]; t1=t+d[2]
    pm=0.5*(p+p1); em=0.5*(e+e1)
    phi1=phi+gN(x+0.5*H,pm,em)*H
    x1=x+H
    return np.array([x1,p1,e1,phi1,t1],float)


def exp_step(s,nu,lambdas=None):
    # The simultaneous additive state update represented as ONE elementwise product.
    # This proves the accumulator update can be exponentiated. It does NOT claim delta(s)
    # is already closed in encoded variables without a logarithmic decode.
    if lambdas is None:
        lambdas=np.array([1e-3,1e-3,1e-2,1e-3,1e-6])
    d=increments_current_only(s,nu)
    z=np.exp(lambdas*s)
    z1=z*np.exp(lambdas*d)
    s1=np.log(z1)/lambdas
    return s1


def jacobian(fun,s,nu):
    s=np.asarray(s,float)
    J=np.zeros((5,5))
    for j in range(5):
        eps=1e-7*max(1.0,abs(s[j]))
        sp=s.copy(); sm=s.copy(); sp[j]+=eps; sm[j]-=eps
        J[:,j]=(fun(sp,nu)-fun(sm,nu))/(2*eps)
    return J


def run_case(p0,e0,q):
    nu=q/(1+q)**2
    s0=np.array([0.0,p0,e0,0.0,0.0],float)
    a=s0.copy(); b=s0.copy(); m=s0.copy()
    nsteps=int(round(STEPS*ORBITS))
    max_ref=np.zeros(5)
    max_mul=np.zeros(5)
    for _ in range(nsteps):
        a1=step_reference_sequential(a,nu)
        b1=step_inline(b,nu)
        m1=exp_step(m,nu)
        max_ref=np.maximum(max_ref,np.abs(a1-b1))
        max_mul=np.maximum(max_mul,np.abs(b1-m1))
        a,b,m=a1,b1,m1
    return {
        'p0':p0,'e0':e0,'q':q,'nu':nu,
        'max_abs_sequential_vs_inline':max_ref.tolist(),
        'max_abs_inline_vs_exp_product_roundtrip':max_mul.tolist(),
        'final_inline':b.tolist(),
    }

cases=[]
for p in P_VALUES:
    for e in E_VALUES:
        cases.append(run_case(p,e,1.0))
for q in (2.0,4.0,10.0):
    cases.append(run_case(60.0,0.45,q))

# Test whether one fixed linear map x' = Mx+b (equiv. fixed monomial map in exp coordinates)
# can represent the 5-state map. Varying Jacobian rules this out for these coordinates.
probe_states=[
    np.array([0.1,24.3,0.2,0.3,10.0]),
    np.array([1.2,60.0,0.45,2.0,100.0]),
    np.array([2.4,120.0,0.55,5.0,300.0]),
]
Js=[jacobian(step_inline,s,0.25) for s in probe_states]
jac_diff=max(float(np.max(np.abs(Js[i]-Js[j]))) for i in range(len(Js)) for j in range(i))

result={
 'state_order':['chi','p','e','phi','t'],
 'n_cases':len(cases),
 'global_max_sequential_vs_inline':np.max(np.array([c['max_abs_sequential_vs_inline'] for c in cases]),axis=0).tolist(),
 'global_max_inline_vs_exp_product_roundtrip':np.max(np.array([c['max_abs_inline_vs_exp_product_roundtrip'] for c in cases]),axis=0).tolist(),
 'fixed_linear_map_test':{
    'criterion':'A fixed monomial map in exponential coordinates implies a constant Jacobian in the additive coordinates (up to homogeneous bias).',
    'max_abs_jacobian_difference_between_probes':jac_diff,
    'jacobians':[J.tolist() for J in Js],
 },
 'cases':cases,
}
with open('/mnt/data/five_state_full_inline_audit.json','w') as f: json.dump(result,f,indent=2)
print(json.dumps({k:v for k,v in result.items() if k!='cases' and k!='fixed_linear_map_test'},indent=2))
print('jacobian_diff',jac_diff)
