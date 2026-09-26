#!/usr/bin/env python3
import math, json, csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

OUT=Path('/mnt/data')
NU=0.25
NS=(12,16,20,24)


def coeffs(N):
    c=[1.0]
    for m in range(N):
        c.append(c[-1]*(2*m+1)/(m+1))
    return np.array(c,float)

COEFF={N:coeffs(N) for N in NS}


def gN_numeric(U,P,E,N):
    C=0.5*(U+1.0/U)
    u=(3.0+E*C)/P
    c=COEFF[N]
    y=complex(c[-1])
    for a in c[-2::-1]:
        y=complex(a)+u*y
    return y


def rr_numeric(U,P,E,nu=NU):
    Cc=0.5*(U+1.0/U)
    Ss=(U-1.0/U)/(2j)
    den=1.0+E*Cc
    r=P/den
    sqP=P**0.5
    rdot=E*Ss/sqP
    vt=den/sqP
    v2=rdot*rdot+vt*vt
    invr=den/P
    Crr=1.6*nu*invr**3
    an=rdot*(18.0*v2+(2.0/3.0)*invr-25.0*rdot*rdot)
    bv=-(6.0*v2-2.0*invr-15.0*rdot*rdot)
    Edot=Crr*(an*rdot+bv*v2)
    pdot=2.0*Crr*bv*P
    Eorb=(E*E-1.0)/(2.0*P)
    edot=(P*Edot+Eorb*pdot)/E
    dt=r*r/sqP
    return pdot*dt, edot*dt, dt


def rhs_complex(x,z,N,tau0):
    U,P,E,H,Q=z
    dp,de,dt=rr_numeric(U,P,E)
    g=gN_numeric(U,P,E,N)
    return np.array([1j*U, dp, de, 1j*g*H, (dt/tau0)*Q],complex)

class TPS:
    __slots__=('c','n')
    def __init__(self,c,n=None):
        a=np.asarray(c,dtype=complex)
        if n is None: n=len(a)-1
        self.n=n
        if len(a)<n+1:
            b=np.zeros(n+1,complex); b[:len(a)]=a; a=b
        else:
            a=a[:n+1].copy()
        self.c=a
    @staticmethod
    def const(v,n):
        a=np.zeros(n+1,complex); a[0]=v; return TPS(a,n)
    def _coerce(self,o):
        return o if isinstance(o,TPS) else TPS.const(o,self.n)
    def __add__(self,o):
        o=self._coerce(o); return TPS(self.c+o.c,self.n)
    __radd__=__add__
    def __neg__(self): return TPS(-self.c,self.n)
    def __sub__(self,o): return self+(-self._coerce(o))
    def __rsub__(self,o): return self._coerce(o)-self
    def __mul__(self,o):
        o=self._coerce(o); n=self.n; a=np.zeros(n+1,complex)
        for k in range(n+1):
            a[k]=sum(self.c[j]*o.c[k-j] for j in range(k+1))
        return TPS(a,n)
    __rmul__=__mul__
    def inv(self):
        n=self.n; a=self.c; b=np.zeros(n+1,complex)
        b[0]=1.0/a[0]
        for k in range(1,n+1):
            b[k]=-(sum(a[j]*b[k-j] for j in range(1,k+1)))/a[0]
        return TPS(b,n)
    def __truediv__(self,o): return self*self._coerce(o).inv()
    def __rtruediv__(self,o): return self._coerce(o)*self.inv()
    def pow_int(self,m):
        if m==0: return TPS.const(1,self.n)
        if m<0: return self.pow_int(-m).inv()
        out=TPS.const(1,self.n); base=self; k=m
        while k:
            if k&1: out=out*base
            base=base*base; k//=2
        return out
    def pow_real(self,alpha):
        # a0^alpha * (1 + delta)^alpha, delta[0]=0; finite exact truncation
        n=self.n; a0=self.c[0]
        delta=self/a0-1.0
        out=TPS.const(1,n)
        term=TPS.const(1,n)
        binom=1.0
        for m in range(1,n+1):
            term=term*delta
            binom*= (alpha-(m-1))/m
            out=out+binom*term
        return (a0**alpha)*out


def rr_tps(U,P,E):
    Cc=(U+1.0/U)*0.5
    Ss=(U-1.0/U)/(2j)
    den=1.0+E*Cc
    r=P/den
    sqP=P.pow_real(0.5)
    rdot=E*Ss/sqP
    vt=den/sqP
    v2=rdot*rdot+vt*vt
    invr=den/P
    Crr=1.6*NU*invr.pow_int(3)
    an=rdot*(18.0*v2+(2.0/3.0)*invr-25.0*rdot*rdot)
    bv=-(6.0*v2-2.0*invr-15.0*rdot*rdot)
    Edot=Crr*(an*rdot+bv*v2)
    pdot=2.0*Crr*bv*P
    Eorb=(E*E-1.0)/(2.0*P)
    edot=(P*Edot+Eorb*pdot)/E
    dt=r*r/sqP
    return pdot*dt, edot*dt, dt


def gN_tps(U,P,E,N):
    Cc=(U+1.0/U)*0.5
    u=(3.0+E*Cc)/P
    c=COEFF[N]
    y=TPS.const(c[-1],U.n)
    for a in c[-2::-1]:
        y=a+u*y
    return y


def taylor_step(z,h,order,N,tau0):
    # Direct Taylor coefficient recurrence on the SAME five state variables.
    # No RK stages. No extra physical state. No log/arg in dynamics.
    a=np.zeros((5,order+1),complex)
    a[:,0]=z
    for k in range(order):
        series=[TPS(a[j],order) for j in range(5)]
        U,P,E,H,Q=series
        dp,de,dt=rr_tps(U,P,E)
        g=gN_tps(U,P,E,N)
        F=[1j*U, dp, de, 1j*g*H, (dt/tau0)*Q]
        for j in range(5):
            a[j,k+1]=F[j].c[k]/(k+1)
    powers=np.array([h**k for k in range(order+1)],complex)
    return a@powers


def integrate_taylor(z0,orbits,steps_per_orbit,order,N,tau0):
    h=2*math.pi/steps_per_orbit
    z=z0.copy()
    for _ in range(int(round(orbits*steps_per_orbit))):
        z=taylor_step(z,h,order,N,tau0)
    return z


def ref_complex(z0,orbits,N,tau0):
    sol=solve_ivp(lambda x,z: rhs_complex(x,z,N,tau0), (0,2*math.pi*orbits), z0,
                  method='DOP853', rtol=2e-13, atol=2e-15)
    if not sol.success: raise RuntimeError(sol.message)
    return sol.y[:,-1]


def state_error(z,ref,p0,tau0):
    U,P,E,H,Q=z; Ur,Pr,Er,Hr,Qr=ref
    # final-observable errors; log only here, never dynamics
    phi_err=abs(np.angle(H/Hr))
    t= tau0*math.log(abs(Q))
    tr=tau0*math.log(abs(Qr))
    vals={
      'U_abs':abs(U-Ur),
      'P_rel':abs(P-Pr)/p0,
      'E_abs':abs(E-Er),
      'phi_abs':phi_err,
      't_rel_orbit':abs(t-tr)/tau0,
    }
    vals['combined']=max(vals.values())
    return {k:float(v) for k,v in vals.items()}

# real RR for angular study

def rr_real(chi,p,e):
    ce,se=math.cos(chi),math.sin(chi)
    den=1+e*ce; r=p/den; hp=math.sqrt(p)
    rdot=e*se/hp; vt=den/hp; v2=rdot*rdot+vt*vt; invr=den/p
    C=1.6*NU*invr**3
    an=rdot*(18*v2+(2/3)*invr-25*rdot*rdot)
    bv=-(6*v2-2*invr-15*rdot*rdot)
    Edot=C*(an*rdot+bv*v2); pdot=2*C*bv*p
    Eorb=(e*e-1)/(2*p); edot=(p*Edot+Eorb*pdot)/e
    dt=r*r/hp
    return pdot*dt,edot*dt,dt

def gN_real(chi,p,e,N):
    u=(3+e*math.cos(chi))/p
    c=COEFF[N]; y=c[-1]
    for a in c[-2::-1]: y=a+u*y
    return y

def gexact_real(chi,p,e):
    u=(3+e*math.cos(chi))/p
    return (1-2*u)**-0.5


def angular_study_case(p0,e0,orbits=3.0):
    # y=[p,e,t,phi_exact,phi12,phi16,phi20,phi24]
    def rhs(x,y):
        p,e=y[0],y[1]
        dp,de,dt=rr_real(x,p,e)
        return [dp,de,dt,gexact_real(x,p,e)] + [gN_real(x,p,e,N) for N in NS]
    y0=[p0,e0,0,0]+[0]*len(NS)
    sol=solve_ivp(rhs,(0,2*math.pi*orbits),y0,method='DOP853',rtol=2e-13,atol=2e-15,dense_output=True)
    yf=sol.y[:,-1]
    phi_exact=yf[3]
    out={'p0':p0,'e0':e0,'orbits':orbits,'final_p':float(yf[0]),'final_e':float(yf[1]),'phi_exact':float(phi_exact),'N':{}}
    xs=np.linspace(0,2*math.pi*orbits,12001)
    ys=sol.sol(xs); ps=ys[0]; es=ys[1]
    umax=float(np.max((3+es*np.cos(xs))/ps))
    out['u_max']=umax
    for j,N in enumerate(NS):
        phiN=yf[4+j]
        # instantaneous max absolute / relative factor error sampled on same p,e path
        c=COEFF[N]
        us=(3+es*np.cos(xs))/ps
        # Horner vectorized
        gn=np.full_like(us,c[-1])
        for aa in c[-2::-1]: gn=aa+us*gn
        ge=(1-2*us)**-0.5
        absmax=float(np.max(np.abs(gn-ge)))
        relmax=float(np.max(np.abs((gn-ge)/ge)))
        cnext=coeffs(N+1)[-1]
        first_omitted=float(cnext*(umax**(N+1)))
        out['N'][str(N)]={
            'phi_final':float(phiN),
            'phi_abs_error':float(abs(phiN-phi_exact)),
            'g_abs_max_error':absmax,
            'g_rel_max_error':relmax,
            'first_omitted_at_umax':first_omitted,
        }
    return out


def main():
    # angular truncation: full 3x3 grid, 3 radial orbits
    angular=[]
    for p0 in (24.3,60.0,120.0):
        for e0 in (0.20,0.45,0.55):
            angular.append(angular_study_case(p0,e0,3.0))

    # integration order study: representative strong/mid/weak cases, N=24 fixed
    integ=[]
    cases=[(24.3,0.45),(60.0,0.45),(120.0,0.45)]
    step_counts=(8,16,32,64)
    orders=(4,6,8)
    for p0,e0 in cases:
        a0=p0/(1-e0*e0); tau0=2*math.pi*a0**1.5
        z0=np.array([1+0j,p0+0j,e0+0j,1+0j,1+0j],complex)
        ref=ref_complex(z0,1.0,24,tau0)
        cres={'p0':p0,'e0':e0,'N':24,'tau0':tau0,'orders':{}}
        for order in orders:
            rows=[]
            for s in step_counts:
                z=integrate_taylor(z0,1.0,s,order,24,tau0)
                er=state_error(z,ref,p0,tau0)
                rows.append({'steps_per_orbit':s,'h':2*math.pi/s,**er})
            # empirical rates from combined error
            for i in range(1,len(rows)):
                e1=rows[i-1]['combined']; e2=rows[i]['combined']
                rows[i]['empirical_order']=float(math.log(e1/e2,2)) if e2>0 and e1>0 else None
            cres['orders'][str(order)]=rows
        integ.append(cres)

    result={
      'definition':{
        'physical_state':['U=e^{i chi}','P=p','E=e','H=e^{i phi}','Q=e^{t/tau0}'],
        'state_dimension':5,
        'interaction_law_count':1,
        'angular_orders':list(NS),
        'integration_orders':[4,6,8],
        'note':'Taylor order changes coefficient depth only; physical state dimension and vector-field law are unchanged.'
      },
      'angular_truncation':angular,
      'integration_order':integ,
    }
    with open(OUT/'order_without_more_S_results.json','w') as f: json.dump(result,f,indent=2)

    # flat CSVs
    with open(OUT/'angular_order_convergence.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['p0','e0','orbits','u_max','N','phi_abs_error','g_abs_max_error','g_rel_max_error','first_omitted_at_umax'])
        for c in angular:
            for N,v in c['N'].items():
                w.writerow([c['p0'],c['e0'],c['orbits'],c['u_max'],N,v['phi_abs_error'],v['g_abs_max_error'],v['g_rel_max_error'],v['first_omitted_at_umax']])
    with open(OUT/'integration_order_convergence.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['p0','e0','N','order','steps_per_orbit','h','combined','U_abs','P_rel','E_abs','phi_abs','t_rel_orbit','empirical_order'])
        for c in integ:
            for o,rows in c['orders'].items():
                for r in rows:
                    w.writerow([c['p0'],c['e0'],c['N'],o,r['steps_per_orbit'],r['h'],r['combined'],r['U_abs'],r['P_rel'],r['E_abs'],r['phi_abs'],r['t_rel_orbit'],r.get('empirical_order')])
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
