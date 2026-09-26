import math, csv, json, hashlib
from pathlib import Path
import numpy as np
from numba import njit

HERE=Path(__file__).resolve().parent
REF=np.genfromtxt(HERE/'reference_raw.csv',delimiter=',',names=True)
META=json.loads((HERE/'reference_metadata.json').read_text())
HSTEP=2*math.pi/4000.0; PH=complex(math.cos(HSTEP),math.sin(HSTEP)); TAU0=1e4
N=0.25; C=1.09; D=0.36; R0=50.0; RF=20.0
A=(4/3)*N*D*C; B=(64/5)*N*C*C

@njit
def rates(P):
    rootC=math.sqrt(C); rootP=math.sqrt(P)
    dP=-(4/3)*N*D*rootC/rootP -(64/5)*N*C*rootC/(P*rootP)
    dt=P*rootP/rootC
    return dP,dt

@njit
def run(maxsteps=600000,sample_every=100):
    # fixed allocation upper bound; trim at return
    ns=maxsteps//sample_every+3
    steps=np.empty(ns,np.int64); Ps=np.empty(ns); ts=np.empty(ns); hrs=np.empty(ns); his=np.empty(ns); urs=np.empty(ns); uis=np.empty(ns)
    P=R0; Q=1.0; Hh=1+0j; U=1+0j; step=0; j=0
    while P>RF and step<=maxsteps:
        if step%sample_every==0:
            steps[j]=step; Ps[j]=P; ts[j]=TAU0*math.log(Q); hrs[j]=Hh.real; his[j]=Hh.imag; urs[j]=U.real; uis[j]=U.imag; j+=1
        k1p,k1t=rates(P)
        k2p,k2t=rates(P+0.5*HSTEP*k1p)
        k3p,k3t=rates(P+0.5*HSTEP*k2p)
        k4p,k4t=rates(P+HSTEP*k3p)
        dP=(HSTEP/6)*(k1p+2*k2p+2*k3p+k4p)
        dt=(HSTEP/6)*(k1t+2*k2t+2*k3t+k4t)
        P=P+dP; Q=Q*math.exp(dt/TAU0); Hh=Hh*PH; U=U*PH; step+=1
    return steps[:j],Ps[:j],ts[:j],hrs[:j],his[:j],urs[:j],uis[:j],step,P,Q,Hh,U

def Ft(x):
    return x**3/(3*A)-B*x*x/(2*A*A)+B*B*x/(A**3)-B**3/(A**4)*math.log(A*x+B)
def Fp(x):
    return math.sqrt(C)*(2*x**1.5/(3*A)-2*B*math.sqrt(x)/(A*A)+2*B**1.5/(A**2.5)*math.atan(math.sqrt(A*x/B)))
def tref(P): return Ft(R0)-Ft(P)
def pref(P): return Fp(R0)-Fp(P)

steps,Ps,ts,hrs,his,urs,uis,ncross,Plast,Qlast,Hlast,Ulast=run()
# build by interpolation one-way
rt=np.asarray(REF['t'],float); rr=np.asarray(REF['r'],float); rp=np.asarray(REF['phi'],float)
out=HERE/'charged_eight_state_benchmark_sampled.csv'
head=['step','r_state','t_state','x_state','y_state','U_re','U_im','H_re','H_im','E_state','N_state','C_state','D_state','N_drift','C_drift','D_drift','t_closed_at_r','phi_closed_at_r','abs_t_closed_residual','abs_H_closed_residual','r_saved_at_t','x_saved_at_t','y_saved_at_t','abs_r_saved_residual','xy_saved_residual']
maxvals=[0,0,0,0]
with out.open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f); w.writerow(head)
    for i,st in enumerate(steps):
        P=float(Ps[i]); t=float(ts[i]); hr=float(hrs[i]); hi=float(his[i]); x=P*hr; y=P*hi
        tc=tref(P); pc=pref(P); hcr=math.cos(pc); hci=math.sin(pc); eh=math.hypot(hr-hcr,hi-hci); et=abs(t-tc)
        rs=float(np.interp(t,rt,rr)); pp=float(np.interp(t,rt,rp)); xs=rs*math.cos(pp); ys=rs*math.sin(pp); er=abs(P-rs); exy=math.hypot(x-xs,y-ys)
        maxvals[0]=max(maxvals[0],et); maxvals[1]=max(maxvals[1],eh); maxvals[2]=max(maxvals[2],er); maxvals[3]=max(maxvals[3],exy)
        w.writerow([int(st),P,t,x,y,float(urs[i]),float(uis[i]),hr,hi,0.0,N,C,D,0.0,0.0,0.0,tc,pc,et,eh,rs,xs,ys,er,exy])
print(json.dumps({'samples':len(steps),'cross_step':int(ncross),'last_r':float(Plast),'maxvals':maxvals,'sha256':hashlib.sha256(out.read_bytes()).hexdigest()},indent=2))
