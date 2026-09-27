#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import csv, hashlib, json, math, time
from pathlib import Path
from fractions import Fraction
import numpy as np
from numba import njit

HERE = Path(__file__).resolve().parent
OUTROOT = HERE / 'results'
OUTROOT.mkdir(exist_ok=True)

STEPS_PER_ORBIT = 4000
HSTEP = 2.0 * math.pi / STEPS_PER_ORBIT
PH_RE = math.cos(HSTEP)
PH_IM = math.sin(HSTEP)
TAU0 = 1.0e4
TARGET_P = 20.0
P0 = 50.0
SAMPLE_EVERY = 1000
MAX_STEPS = 20_000_000

UR,UI,P,E,HR,HI,Q,N,C,D = range(10)
K1=10; K2=13; K3=16; K4=19
PS=22; ES=23; DV=24; PM=27; EM=28; DPHI=29; QB=30
NST=41

# G=q0=c=1 normalization from Paper 6 supplement.
G = 1.0
Q0 = 1.0
C_LIGHT = 1.0
M_OVER_Q0_EXACT = Fraction(20,3)
M_TOTAL = float(M_OVER_Q0_EXACT) * Q0
M_A = M_TOTAL / 2.0
M_B = M_TOTAL / 2.0
N0 = 0.25

# Same five source patterns as Paper 6. Self terms are derived separately for aa and bb.
SOURCE_PATTERNS = (
    ('n1_attractive', 1, +1, -1),
    ('n1_repulsive',  1, +1, +1),
    ('n3_attractive', 3, +1, -1),
    ('n3_repulsive',  3, +1, +1),
    ('n4_attractive', 4, +1, -1),
)

@njit(cache=True)
def _blend(q,c):
    x=0.0
    for j in range(11):
        x += q[j]*c[j]
    return x

@njit(cache=True)
def _rates(p,n,c,d):
    rp=math.sqrt(p); rc=math.sqrt(c)
    dp=-(4.0/3.0)*n*d*rc/rp -(64.0/5.0)*n*c*rc/(p*rp)
    dt=p*rp/rc
    return dp,0.0,dt

@njit(cache=True)
def transition(z):
    q=z[QB:QB+11]
    p=z[P]; e=z[E]; n=z[N]; c=z[C]; dd=z[D]
    k1=z[K1:K1+3]; k2=z[K2:K2+3]; k3=z[K3:K3+3]; k4=z[K4:K4+3]
    ps=z[PS]; es=z[ES]; dv=z[DV:DV+3]; pm=z[PM]; em=z[EM]; dphi=z[DPHI]
    rP,rE,rT=_rates(p,n,c,dd)
    rPsP,rPsE,rPsT=_rates(ps,n,c,dd)
    cand=np.empty((30,11),dtype=np.float64)
    for i in range(30):
        for j in range(11):
            cand[i,j]=z[i]
    cand[K1,0]=rP; cand[K1+1,0]=rE; cand[K1+2,0]=rT
    for a in range(3): cand[K1+a,10]=0.0
    cand[PS,1]=p+0.5*HSTEP*k1[0]; cand[ES,1]=e+0.5*HSTEP*k1[1]
    cand[K2,2]=rPsP; cand[K2+1,2]=rPsE; cand[K2+2,2]=rPsT
    for a in range(3): cand[K2+a,10]=0.0
    cand[PS,3]=p+0.5*HSTEP*k2[0]; cand[ES,3]=e+0.5*HSTEP*k2[1]
    cand[K3,4]=rPsP; cand[K3+1,4]=rPsE; cand[K3+2,4]=rPsT
    for a in range(3): cand[K3+a,10]=0.0
    cand[PS,5]=p+HSTEP*k3[0]; cand[ES,5]=e+HSTEP*k3[1]
    cand[K4,6]=rPsP; cand[K4+1,6]=rPsE; cand[K4+2,6]=rPsT
    for a in range(3): cand[K4+a,10]=0.0
    for a in range(3):
        cand[DV+a,7]=(HSTEP/6.0)*(k1[a]+2.0*k2[a]+2.0*k3[a]+k4[a])
        cand[DV+a,10]=0.0
    cand[PM,8]=p+0.5*dv[0]; cand[EM,8]=e+0.5*dv[1]
    cand[PM,10]=p+dv[0]; cand[EM,10]=e+dv[1]
    cand[DPHI,9]=HSTEP; cand[DPHI,10]=0.0
    ur=z[UR]; ui=z[UI]; hr=z[HR]; hi=z[HI]
    cand[UR,10]=ur*PH_RE-ui*PH_IM; cand[UI,10]=ur*PH_IM+ui*PH_RE
    cand[P,10]=p+dv[0]; cand[E,10]=e+dv[1]
    cd=math.cos(dphi); sd=math.sin(dphi)
    cand[HR,10]=hr*cd-hi*sd; cand[HI,10]=hr*sd+hi*cd
    cand[Q,10]=z[Q]*math.exp(dv[2]/TAU0)
    cand[PS,10]=p+dv[0]; cand[ES,10]=e+dv[1]
    out=np.empty(NST,dtype=np.float64)
    for i in range(30): out[i]=_blend(q,cand[i])
    out[QB]=q[10]
    for j in range(1,11): out[QB+j]=q[j-1]
    return out

@njit(cache=True)
def macrostep(z):
    for _ in range(11): z=transition(z)
    return z

@njit(cache=True)
def run_to_target(z, target_p, max_steps, sample_every):
    # sampled columns: step,P,E,Ure,Uim,Hre,Him,Q,N,C,D,t,x,y,q_index,PGW,PEM
    max_samples = max_steps // sample_every + 3
    samples=np.empty((max_samples,17),dtype=np.float64)
    ns=0
    first_bad=-1
    first_bad_p=0.0
    step=0
    while step <= max_steps:
        if step % sample_every == 0 or z[P] <= target_p:
            t = TAU0*math.log(z[Q]) if z[Q] > 0.0 and math.isfinite(z[Q]) else math.inf
            qidx=0; qmax=z[QB]
            for j in range(1,11):
                if z[QB+j]>qmax: qmax=z[QB+j]; qidx=j
            mu=z[N]  # M=1 representation retained in state map; N=nu=1/4.
            p=z[P]; cc=z[C]; dd=z[D]
            pgw=(32.0/5.0)*mu*mu*cc*cc*cc/(p**5)
            pem=(2.0/3.0)*mu*mu*dd*cc*cc/(p**4)
            samples[ns,0]=step; samples[ns,1]=p; samples[ns,2]=z[E]
            samples[ns,3]=z[UR]; samples[ns,4]=z[UI]; samples[ns,5]=z[HR]; samples[ns,6]=z[HI]
            samples[ns,7]=z[Q]; samples[ns,8]=z[N]; samples[ns,9]=cc; samples[ns,10]=dd
            samples[ns,11]=t; samples[ns,12]=p*z[HR]; samples[ns,13]=p*z[HI]; samples[ns,14]=qidx
            samples[ns,15]=pgw; samples[ns,16]=pem
            ns += 1
        if not math.isfinite(z[Q]) and first_bad < 0:
            first_bad=step; first_bad_p=z[P]
        if z[P] <= target_p:
            return z, samples[:ns], step, first_bad, first_bad_p, True
        z=macrostep(z)
        step += 1
    return z, samples[:ns], step, first_bad, first_bad_p, False


def self_relations(n:int, sign:int):
    lam=Fraction(sign*3*n,10)
    c=Fraction(1,1)-lam*lam
    d=Fraction(0,1)
    return float(lam), float(c), float(d)


def init_state(c0,d0):
    z=np.zeros(NST,dtype=np.float64)
    z[UR]=1.0; z[P]=P0; z[HR]=1.0; z[Q]=1.0
    z[N]=N0; z[C]=c0; z[D]=d0
    z[PS]=P0; z[PM]=P0; z[QB]=1.0
    return z


def write_first_macro(path,z):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f)
        w.writerow(['micro','q_index','P','N','C','D','k1p','k2p','k3p','k4p','dP','Q'])
        zm=z.copy()
        for m in range(12):
            w.writerow([m,int(np.argmax(zm[QB:QB+11])),zm[P],zm[N],zm[C],zm[D],zm[K1],zm[K2],zm[K3],zm[K4],zm[DV],zm[Q]])
            if m<11: zm=transition(zm)


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def _write_sampled_case(uid, n, max_steps, filename_prefix):
    lam,c0,d0=self_relations(n,+1)
    out=OUTROOT/uid; out.mkdir(exist_ok=True)
    z0=init_state(c0,d0)
    write_first_macro(out/'raw_first_macro_microsteps.csv',z0)
    transition(z0.copy())  # JIT warm-up, not a state mutation of z0
    t0=time.time()
    zf,samples,steps,first_bad,first_bad_p,reached=run_to_target(z0.copy(),TARGET_P,max_steps,SAMPLE_EVERY)
    elapsed=time.time()-t0
    cols=['step','P','E','U_re','U_im','H_re','H_im','Q','N','C','D','t_readout','x_readout','y_readout','q_index','P_GW','P_EM']
    np.savetxt(out/f'{filename_prefix}.csv',samples,delimiter=',',header=','.join(cols),comments='')
    np.save(out/'final_full_state.npy',zf)
    summary={
        'unique_case':uid,'n':n,'lambda_self_abs':abs(lam),'C':c0,'D':d0,
        'initial_state':{'P':P0,'N':N0,'C':c0,'D':d0},
        'target_P':TARGET_P,'reached_target':bool(reached),'macro_steps_executed':int(steps),
        'cycles_executed':float(steps/STEPS_PER_ORBIT),'final_P':float(zf[P]),
        'first_nonfinite_Q_step':None if first_bad<0 else int(first_bad),
        'first_nonfinite_Q_P':None if first_bad<0 else float(first_bad_p),
        'P_EM_identically_zero':bool(d0==0.0),'elapsed_seconds':elapsed,
        'sample_file':f'{filename_prefix}.csv'
    }
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return summary


def main():
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument('--full-n3',action='store_true',help='also execute the full n=3 self trajectory to P<=20; this is long')
    args=ap.parse_args()

    matrix=[]
    for src,n,sa,sb in SOURCE_PATTERNS:
        for side,sign in (('aa',sa),('bb',sb)):
            lam,c0,d0=self_relations(n,sign)
            matrix.append({
                'source_pattern':src,'self_channel':side,'n':n,'charge_sign':sign,
                'lambda_self':lam,'C_self':c0,'D_self':d0,
                'status':'runnable' if c0>0 else 'outside_current_real_quasicircular_domain'
            })
    (OUTROOT/'case_matrix.json').write_text(json.dumps({
        'normalization':{'G':1,'q0':1,'c':1,'M_over_q0':'20/3','M_total':M_TOTAL,'m_A':M_A,'m_B':M_B},
        'numerics':{'P0':P0,'target_P':TARGET_P,'steps_per_orbit':STEPS_PER_ORBIT,'sample_every_macrosteps':SAMPLE_EVERY},
        'cases':matrix
    },ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

    # Fresh full run for n=1.
    n1=_write_sampled_case('self_n1',1,MAX_STEPS,'sampled_trajectory')
    print('self_n1',json.dumps(n1,ensure_ascii=False),flush=True)

    # Fresh prefix run for n=3. The full n=3 raw already exists in Paper 6 with exactly the same initial state and transition.
    n3_steps=MAX_STEPS if args.full_n3 else 100_000
    n3_name='sampled_trajectory' if args.full_n3 else 'sampled_trajectory_prefix100k'
    n3=_write_sampled_case('self_n3_fresh',3,n3_steps,n3_name)
    print('self_n3_fresh',json.dumps(n3,ensure_ascii=False),flush=True)

    run_summary={
      'n1_fresh_full':n1,
      'n3_fresh':n3,
      'full_n3_requested':bool(args.full_n3),
      'n4':{'C':-0.44,'D':0.0,'status':'outside_current_real_quasicircular_domain'}
    }
    (OUTROOT/'LOCAL_RUN_SUMMARY.json').write_text(json.dumps(run_summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':
    main()
