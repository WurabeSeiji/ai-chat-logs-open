import math, ast, inspect
import numpy as np

HSTEP = 2.0*math.pi/4000.0
ORDER=12
TAU0=1.0e4
PH=complex(math.cos(HSTEP), math.sin(HSTEP))
PH2=complex(math.cos(HSTEP/2), math.sin(HSTEP/2))

def coeffs(N):
    c=[1.0]
    for m in range(N): c.append(c[-1]*(2*m+1)/(m+1))
    return np.array(c,float)
C=coeffs(ORDER)

def trigU(U):
    inv=1.0/U
    return (0.5*(U+inv)).real, ((U-inv)/(2j)).real

def rr(U,P,E,nu):
    ce,se=trigU(U)
    den=1.0+E*ce
    r=P/den
    hp=math.sqrt(P)
    rdot=E*se/hp
    vt=den/hp
    v2=rdot*rdot+vt*vt
    invr=1.0/r
    Crr=1.6*nu*invr**3
    an=rdot*(18*v2+(2/3)*invr-25*rdot*rdot)
    bv=-(6*v2-2*invr-15*rdot*rdot)
    Edot=Crr*(an*rdot+bv*v2)
    pdot=2*Crr*bv*P
    En=(E*E-1)/(2*P)
    edot=(P*Edot+En*pdot)/E
    dt_dchi=r*r/hp
    return np.array([pdot*dt_dchi,edot*dt_dchi,dt_dchi],float)

def gN(U,P,E):
    ce,_=trigU(U)
    u=(3+E*ce)/P
    y=float(C[-1])
    for a in C[-2::-1]: y=float(a)+u*y
    return y

def step_mul(Z,nu):
    # persistent state remains [U, P, E, Hphi, Q]
    U=complex(Z[0]); P=float(Z[1].real); E=float(Z[2].real); Hphi=complex(Z[3]); Q=float(Z[4].real)
    # all stage expressions are functions of old U,P,E only
    k1=rr(U,P,E,nu)
    k2=rr(U*PH2,P+0.5*HSTEP*k1[0],E+0.5*HSTEP*k1[1],nu)
    k3=rr(U*PH2,P+0.5*HSTEP*k2[0],E+0.5*HSTEP*k2[1],nu)
    k4=rr(U*PH,P+HSTEP*k3[0],E+HSTEP*k3[1],nu)
    d=(HSTEP/6)*(k1+2*k2+2*k3+k4)
    dphi=HSTEP*gN(U*PH2,P+0.5*d[0],E+0.5*d[1])
    fac=np.array([
        PH,
        1+d[0]/P,
        1+d[1]/E,
        complex(math.cos(dphi), math.sin(dphi)),
        math.exp(d[2]/TAU0)
    ],complex)
    return np.asarray(Z,complex)*fac

def step_ref(x,nu):
    chi,p,e,phi,t=map(float,x)
    U=complex(math.cos(chi),math.sin(chi))
    k1=rr(U,p,e,nu)
    k2=rr(U*PH2,p+0.5*HSTEP*k1[0],e+0.5*HSTEP*k1[1],nu)
    k3=rr(U*PH2,p+0.5*HSTEP*k2[0],e+0.5*HSTEP*k2[1],nu)
    k4=rr(U*PH,p+HSTEP*k3[0],e+HSTEP*k3[1],nu)
    d=(HSTEP/6)*(k1+2*k2+2*k3+k4)
    dphi=HSTEP*gN(U*PH2,p+0.5*d[0],e+0.5*d[1])
    return np.array([chi+HSTEP,p+d[0],e+d[1],phi+dphi,t+d[2]],float)

def init(x):
    chi,p,e,phi,t=map(float,x)
    return np.array([complex(math.cos(chi),math.sin(chi)),p,e,complex(math.cos(phi),math.sin(phi)),math.exp(t/TAU0)],complex)

def readout(Z):
    U,P,E,Hphi,Q=Z
    return np.array([math.atan2(U.imag,U.real), P.real, E.real, math.atan2(Hphi.imag,Hphi.real), TAU0*math.log(Q.real)],float)

def static_audit():
    src=inspect.getsource(step_mul)
    tree=ast.parse(src)
    bad=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.Call):
            nm=node.func.id if isinstance(node.func,ast.Name) else (node.func.attr if isinstance(node.func,ast.Attribute) else '')
            if nm in {'log','log1p','angle','unwrap','atan2','arctan2'}: bad.append((nm,node.lineno))
    return bad

def run_case(p0,e0,q,steps=12000):
    nu=q/(1+q)**2
    x=np.array([0.,p0,e0,0.,0.])
    z=init(x)
    max_phase_chi=max_phase_phi=max_p=max_e=max_t=max_r=max_xy=0.0
    for _ in range(steps):
        x=step_ref(x,nu)
        z=step_mul(z,nu)
        U,P,E,Hp,Q=z
        # compare phase directly (no branch issue)
        uref=complex(math.cos(x[0]),math.sin(x[0]))
        href=complex(math.cos(x[3]),math.sin(x[3]))
        max_phase_chi=max(max_phase_chi,abs(U-uref))
        max_phase_phi=max(max_phase_phi,abs(Hp-href))
        max_p=max(max_p,abs(P.real-x[1])); max_e=max(max_e,abs(E.real-x[2]))
        tr=TAU0*math.log(Q.real)
        max_t=max(max_t,abs(tr-x[4]))
        # physical orbit readout from multiplicative state without feeding back
        ce=0.5*(U+1/U); r=P.real/(1+E.real*ce.real)
        xs=r*Hp.real; ys=r*Hp.imag
        rr=x[1]/(1+x[2]*math.cos(x[0])); xr=rr*math.cos(x[3]); yr=rr*math.sin(x[3])
        max_r=max(max_r,abs(r-rr)); max_xy=max(max_xy,abs(xs-xr),abs(ys-yr))
    return np.array([max_phase_chi,max_p,max_e,max_phase_phi,max_t,max_r,max_xy])

print('forbidden_calls_in_update',static_audit())
cases=[]
for p in (24.3,60.,120.):
  for e in (.20,.45,.55): cases.append((p,e,1.0))
for q in (2.,4.,10.): cases.append((60.,.45,q))
g=np.zeros(7)
for c in cases:
    r=run_case(*c); g=np.maximum(g,r); print(c, r)
print('GLOBAL',g)

if __name__ == '__main__':
    import csv
    data=np.genfromtxt('/mnt/data/experiment02_sn_2p5pn_rr_C1_N12.csv',delimiter=',',names=True)
    x=np.array([0.,24.3,0.4358898943540673,0.,0.])
    z=init(x)
    errs={k:0.0 for k in ['p','e','t','r','x','y','phi_phase']}
    for i,row in enumerate(data):
        U,P,E,Hp,Q=z
        ce=(0.5*(U+1/U)).real
        r=P.real/(1+E.real*ce)
        xx=r*Hp.real; yy=r*Hp.imag
        href=complex(math.cos(float(row['phi'])),math.sin(float(row['phi'])))
        errs['p']=max(errs['p'],abs(P.real-row['p']))
        errs['e']=max(errs['e'],abs(E.real-row['e']))
        errs['t']=max(errs['t'],abs(TAU0*math.log(Q.real)-row['t_osc']))
        errs['r']=max(errs['r'],abs(r-row['r']))
        errs['x']=max(errs['x'],abs(xx-row['x']))
        errs['y']=max(errs['y'],abs(yy-row['y']))
        errs['phi_phase']=max(errs['phi_phase'],abs(Hp-href))
        if i+1<len(data): z=step_mul(z,0.25)
    print('CSV_C1',errs)
