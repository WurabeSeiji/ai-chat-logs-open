#!/usr/bin/env python3
import math, json
import numpy as np

SEED=20260924
SAMPLES=100000

def coeffs(N):
    c=[1.0]
    for m in range(N): c.append(c[-1]*(2*m+1)/(m+1))
    return np.array(c,float)

def rr_original(U,P,E,nu):
    inv=1.0/U
    c=(U+inv)/2
    s=(U-inv)/(2j)
    den=1+E*c
    r=P/den
    hp=P**0.5
    rdot=E*s/hp
    vt=den/hp
    v2=rdot*rdot+vt*vt
    invr=1/r
    C=1.6*nu*invr**3
    an=rdot*(18*v2+(2/3)*invr-25*rdot*rdot)
    bv=-(6*v2-2*invr-15*rdot*rdot)
    Edot=C*(an*rdot+bv*v2)
    pdot=2*C*bv*P
    En=(E*E-1)/(2*P)
    edot=(P*Edot+En*pdot)/E
    dt=r*r/hp
    return np.array([pdot*dt, edot*dt, dt],complex)

def rr_laurent(U,P,E,nu):
    Ui=1/U
    fp=-(nu/(P**1.5))*(
        6*E**3*(U**3+Ui**3)
        +20*E**2*(U**2+Ui**2)
        +(18/5*E**3+112/5*E)*(U+Ui)
        +56/5*E**2+64/5
    )
    fe=-(nu/(P**2.5))*(
        5/4*E**4*(U**5+Ui**5)
        +55/6*E**3*(U**4+Ui**4)
        +(117/20*E**4+316/15*E**2)*(U**3+Ui**3)
        +(22*E**3+56/3*E)*(U**2+Ui**2)
        +(5/2*E**4+404/15*E**2+32/5)*(U+Ui)
        +121/15*E**3+304/15*E
    )
    ft=P**1.5/(1+E*(U+Ui)/2)**2
    return np.array([fp,fe,ft],complex)

def g_recurrence(U,P,E,N=12):
    u=(3+E*(U+1/U)/2)/P
    c=coeffs(N)
    y=complex(c[-1])
    for a in c[-2::-1]: y=a+u*y
    return y

def g_explicit12(U,P,E):
    u=(3+E*(U+1/U)/2)/P
    c=[1,1,3/2,5/2,35/8,63/8,231/16,429/16,6435/128,12155/128,46189/256,88179/256,676039/1024]
    return sum(c[m]*u**m for m in range(13))

def main():
    rng=np.random.default_rng(SEED)
    max_abs=np.zeros(3); max_rel=np.zeros(3); max_g=0.0; max_c=0.0
    worst=None
    for _ in range(SAMPLES):
        chi=rng.uniform(-8*math.pi,8*math.pi)
        U=np.exp(1j*chi)
        P=rng.uniform(20,160)
        E=rng.uniform(0.05,0.70)
        q=float(np.exp(rng.uniform(math.log(1),math.log(10))))
        nu=q/(1+q)**2
        a=rr_original(U,P,E,nu); b=rr_laurent(U,P,E,nu)
        d=np.abs(a-b); rel=d/np.maximum(1,np.abs(a))
        if d.max()>max_abs.max(): worst={'chi':float(chi),'P':float(P),'E':float(E),'q':float(q),'original':[[float(z.real),float(z.imag)] for z in a],'laurent':[[float(z.real),float(z.imag)] for z in b],'abs_diff':d.tolist()}
        max_abs=np.maximum(max_abs,d); max_rel=np.maximum(max_rel,rel)
        max_g=max(max_g,abs(g_recurrence(U,P,E)-g_explicit12(U,P,E)))
    cr=coeffs(12)
    ce=np.array([1,1,3/2,5/2,35/8,63/8,231/16,429/16,6435/128,12155/128,46189/256,88179/256,676039/1024],float)
    max_c=float(np.max(np.abs(cr-ce)))
    out={
      'seed':SEED,'samples':SAMPLES,
      'ranges':{'chi':'[-8pi,8pi]','P':'[20,160]','E':'[0.05,0.70]','q_mass_ratio':'log-uniform [1,10]'},
      'max_abs_rr_components':max_abs.tolist(),
      'max_rel_rr_components':max_rel.tolist(),
      'max_abs_g12_recurrence_vs_explicit':float(max_g),
      'max_abs_coefficient_difference':max_c,
      'g12_coefficients':ce.tolist(),
      'worst_rr_sample':worst,
      'interpretation':'Floating-point equivalence test of the explicit Laurent formulas against the original trigonometric RR formulas. This is an algebraic verification, not an independent physical validation.'
    }
    with open('/mnt/data/expanded_rr_g12_verification.json','w') as f: json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
if __name__=='__main__': main()
