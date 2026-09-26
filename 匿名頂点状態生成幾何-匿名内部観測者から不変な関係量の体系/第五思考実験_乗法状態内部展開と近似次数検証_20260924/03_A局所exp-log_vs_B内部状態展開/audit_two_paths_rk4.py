#!/usr/bin/env python3
import math, json
import numpy as np

H = 2.0*math.pi/4000.0
ORDER=12
TAU0=1.0e4
PH=complex(math.cos(H),math.sin(H))
PH2=complex(math.cos(H/2),math.sin(H/2))


def coeffs(N):
    c=[1.0]
    for m in range(N): c.append(c[-1]*(2*m+1)/(m+1))
    return np.array(c,float)
C=coeffs(ORDER)

def trigU(U):
    inv=1/U
    return float((0.5*(U+inv)).real), float(((U-inv)/(2j)).real)

def rr(U,P,E,nu):
    ce,se=trigU(U)
    den=1.0+E*ce
    r=P/den
    hp=math.sqrt(P)
    rdot=E*se/hp
    vt=den/hp
    v2=rdot*rdot+vt*vt
    invr=1/r
    Cr=1.6*nu*invr**3
    an=rdot*(18*v2+(2/3)*invr-25*rdot*rdot)
    bv=-(6*v2-2*invr-15*rdot*rdot)
    Edot=Cr*(an*rdot+bv*v2)
    pdot=2*Cr*bv*P
    En=(E*E-1)/(2*P)
    edot=(P*Edot+En*pdot)/E
    dt=r*r/hp
    return np.array([pdot*dt,edot*dt,dt],float)

def gN(U,P,E):
    ce,_=trigU(U)
    u=(3+E*ce)/P
    y=float(C[-1])
    for a in C[-2::-1]: y=float(a)+u*y
    return y

def ref_step(x,nu):
    U,P,E,Hphi,Q=x
    k1=rr(U,P,E,nu)
    k2=rr(U*PH2,P+0.5*H*k1[0],E+0.5*H*k1[1],nu)
    k3=rr(U*PH2,P+0.5*H*k2[0],E+0.5*H*k2[1],nu)
    k4=rr(U*PH,P+H*k3[0],E+H*k3[1],nu)
    d=(H/6)*(k1+2*k2+2*k3+k4)
    dphi=H*gN(U*PH2,P+0.5*d[0],E+0.5*d[1])
    return np.array([U*PH,P+d[0],E+d[1],Hphi*complex(math.cos(dphi),math.sin(dphi)),Q*math.exp(d[2]/TAU0)],complex)

# ---------- A: local exp/multiply/log only around additive RK4 combinations ----------
def ladd(a,b,lam=1e-2):
    # Represents a+b as log(exp(lam*a)*exp(lam*b))/lam.
    return math.log(math.exp(lam*a)*math.exp(lam*b))/lam

def lsum(vals,lam=1e-2):
    # left fold using local exp/multiply/log
    s=0.0
    for v in vals:
        s=ladd(s,float(v),lam)
    return s

def step_A(x,nu):
    U,P,E,Hphi,Q=x
    P=float(P.real); E=float(E.real); Q=float(Q.real)
    k1=rr(U,P,E,nu)
    P2=ladd(P,0.5*H*k1[0]); E2=ladd(E,0.5*H*k1[1])
    k2=rr(U*PH2,P2,E2,nu)
    P3=ladd(P,0.5*H*k2[0]); E3=ladd(E,0.5*H*k2[1])
    k3=rr(U*PH2,P3,E3,nu)
    P4=ladd(P,H*k3[0]); E4=ladd(E,H*k3[1])
    k4=rr(U*PH,P4,E4,nu)
    d=np.array([
        (H/6)*lsum([k1[0],2*k2[0],2*k3[0],k4[0]]),
        (H/6)*lsum([k1[1],2*k2[1],2*k3[1],k4[1]]),
        (H/6)*lsum([k1[2],2*k2[2],2*k3[2],k4[2]]),
    ])
    Pm=ladd(P,0.5*d[0]); Em=ladd(E,0.5*d[1])
    dphi=H*gN(U*PH2,Pm,Em)
    Pn=ladd(P,d[0]); En=ladd(E,d[1])
    # Hphi and Q are multiplicative states already.
    return np.array([U*PH,Pn,En,Hphi*complex(math.cos(dphi),math.sin(dphi)),Q*math.exp(d[2]/TAU0)],complex)

# ---------- B: explicit synchronous microstate, 11-phase internal cycle ----------
# Microstate fields kept in dict for readability; *all* are explicit states.
# phase q is one-hot length 11, so control is also explicit state.

def init_B(U,P,E,Hphi,Q):
    return {
      'U':complex(U),'P':float(complex(P).real),'E':float(complex(E).real),'H':complex(Hphi),'Q':float(complex(Q).real),
      'k1':np.zeros(3),'k2':np.zeros(3),'k3':np.zeros(3),'k4':np.zeros(3),
      'Us':complex(U),'Ps':float(complex(P).real),'Es':float(complex(E).real),
      'd':np.zeros(3),'Um':complex(U),'Pm':float(complex(P).real),'Em':float(complex(E).real),'dphi':0.0,
      'q':np.array([1.0]+[0.0]*10)
    }

def blend(q,cands):
    # All cands were computed from OLD microstate. q selects without hidden branch state.
    x=cands[0]*q[0]
    for j in range(1,len(cands)): x=x+cands[j]*q[j]
    return x

def micro_B(s,nu):
    U,P,E,Hh,Q=s['U'],s['P'],s['E'],s['H'],s['Q']
    k1,k2,k3,k4=s['k1'],s['k2'],s['k3'],s['k4']
    Us,Ps,Es=s['Us'],s['Ps'],s['Es']
    d=s['d']; Um,Pm,Em=s['Um'],s['Pm'],s['Em']; dphi=s['dphi']; q=s['q']

    # Candidate values for each internal phase; every candidate reads OLD s only.
    # 0 compute k1
    # 1 form stage2
    # 2 compute k2
    # 3 form stage3
    # 4 compute k3
    # 5 form stage4
    # 6 compute k4
    # 7 combine d
    # 8 form midpoint
    # 9 compute dphi
    # 10 commit macro state and reset work registers
    k1c=[rr(U,P,E,nu),k1,k1,k1,k1,k1,k1,k1,k1,k1,np.zeros(3)]
    k2c=[k2,k2,rr(Us,Ps,Es,nu),k2,k2,k2,k2,k2,k2,k2,np.zeros(3)]
    k3c=[k3,k3,k3,k3,rr(Us,Ps,Es,nu),k3,k3,k3,k3,k3,np.zeros(3)]
    k4c=[k4,k4,k4,k4,k4,k4,rr(Us,Ps,Es,nu),k4,k4,k4,np.zeros(3)]

    Usc=[Us,U*PH2,Us,U*PH2,Us,U*PH,Us,Us,Us,Us,U*PH]
    Psc=[Ps,P+0.5*H*k1[0],Ps,P+0.5*H*k2[0],Ps,P+H*k3[0],Ps,Ps,Ps,Ps,P]
    Esc=[Es,E+0.5*H*k1[1],Es,E+0.5*H*k2[1],Es,E+H*k3[1],Es,Es,Es,Es,E]

    dc=[d,d,d,d,d,d,d,(H/6)*(k1+2*k2+2*k3+k4),d,d,np.zeros(3)]
    Umc=[Um,Um,Um,Um,Um,Um,Um,Um,U*PH2,Um,U*PH]
    Pmc=[Pm,Pm,Pm,Pm,Pm,Pm,Pm,Pm,P+0.5*d[0],Pm,P+d[0]]
    Emc=[Em,Em,Em,Em,Em,Em,Em,Em,E+0.5*d[1],Em,E+d[1]]
    dpc=[dphi,dphi,dphi,dphi,dphi,dphi,dphi,dphi,dphi,H*gN(Um,Pm,Em),0.0]

    Uc=[U]*10+[U*PH]
    Pc=[P]*10+[P+d[0]]
    Ec=[E]*10+[E+d[1]]
    Hc=[Hh]*10+[Hh*complex(math.cos(dphi),math.sin(dphi))]
    Qc=[Q]*10+[Q*math.exp(d[2]/TAU0)]

    out={}
    out['U']=blend(q,Uc); out['P']=float(blend(q,Pc)); out['E']=float(blend(q,Ec)); out['H']=blend(q,Hc); out['Q']=float(blend(q,Qc))
    out['k1']=blend(q,k1c); out['k2']=blend(q,k2c); out['k3']=blend(q,k3c); out['k4']=blend(q,k4c)
    out['Us']=blend(q,Usc); out['Ps']=float(blend(q,Psc)); out['Es']=float(blend(q,Esc))
    out['d']=blend(q,dc); out['Um']=blend(q,Umc); out['Pm']=float(blend(q,Pmc)); out['Em']=float(blend(q,Emc)); out['dphi']=float(blend(q,dpc))
    out['q']=np.roll(q,1)
    return out

def macro_B(s,nu):
    # 11 explicit internal synchronous transitions, no hidden work variables across them.
    for _ in range(11): s=micro_B(s,nu)
    return s

def coreB(s):
    return np.array([s['U'],s['P'],s['E'],s['H'],s['Q']],complex)


def compare_case(p0,e0,qmass,nsteps=12000):
    nu=qmass/(1+qmass)**2
    x0=np.array([1+0j,p0,e0,1+0j,1+0j],complex)
    r=x0.copy(); a=x0.copy(); b=init_B(*x0)
    ma=np.zeros(5); mb=np.zeros(5)
    for _ in range(nsteps):
        r=ref_step(r,nu)
        a=step_A(a,nu)
        b=macro_B(b,nu)
        bb=coreB(b)
        ma=np.maximum(ma,np.abs(a-r))
        mb=np.maximum(mb,np.abs(bb-r))
    return ma,mb,r,a,coreB(b)

def csv_check():
    path='/mnt/data/experiment02_sn_2p5pn_rr_C1_N12.csv'
    data=np.genfromtxt(path,delimiter=',',names=True)
    e0=float(data['e'][0]); p0=float(data['p'][0]); nu=0.25
    a=np.array([1+0j,p0,e0,1+0j,1+0j],complex)
    b=init_B(*a)
    ea={'A_p':0,'A_e':0,'A_t':0,'A_xy':0,'B_p':0,'B_e':0,'B_t':0,'B_xy':0}
    for i,row in enumerate(data):
        for label,z in [('A',a),('B',coreB(b))]:
            U,P,E,Hh,Q=z
            ce=trigU(U)[0]; rr0=P.real/(1+E.real*ce)
            xx=rr0*Hh.real; yy=rr0*Hh.imag
            ea[label+'_p']=max(ea[label+'_p'],abs(P.real-row['p']))
            ea[label+'_e']=max(ea[label+'_e'],abs(E.real-row['e']))
            tt=TAU0*math.log(Q.real)
            ea[label+'_t']=max(ea[label+'_t'],abs(tt-row['t_osc']))
            ea[label+'_xy']=max(ea[label+'_xy'],abs(xx-row['x']),abs(yy-row['y']))
        if i+1<len(data):
            a=step_A(a,nu); b=macro_B(b,nu)
    return ea

def main():
    cases=[]
    for p in (24.3,60.,120.):
        for e in (.20,.45,.55): cases.append((p,e,1.0))
    for q in (2.,4.,10.): cases.append((60.,.45,q))
    ga=np.zeros(5); gb=np.zeros(5); rows=[]
    for c in cases:
        ma,mb,*_=compare_case(*c)
        ga=np.maximum(ga,ma); gb=np.maximum(gb,mb)
        rows.append({'case':c,'A':ma.tolist(),'B':mb.tolist()})
        print(c,'A',ma,'B',mb)
    print('GLOBAL_A',ga)
    print('GLOBAL_B',gb)
    cc=csv_check(); print('CSV',cc)
    out={'global_A':ga.tolist(),'global_B':gb.tolist(),'cases':rows,'csv':cc,'method_B_microphases':11,'method_B_explicit_phase_state':True}
    with open('/mnt/data/two_path_rk4_audit.json','w') as f: json.dump(out,f,indent=2)
if __name__=='__main__': main()
