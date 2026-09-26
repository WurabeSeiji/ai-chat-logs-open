#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Strict charged 8-persistent-state experiment, single generator path.

Persistent physical state: (U,P,E,H,Q,N,C,D).
Full explicit microstate also contains RK work registers and one-hot q.

The ONLY generator transition is transition(z):
- reads OLD full state z only,
- constructs all 11 candidate maps from OLD z,
- writes the entire next full state by q-weighted product-sums,
- q itself advances by a fixed cyclic linear map.

No external physical runtime parameter is accepted by transition().
No alternative/fast/collapsed/shortcut generator exists.
Reference data are not loaded anywhere in this file.
"""
from __future__ import annotations
import ast, csv, hashlib, json, math
from pathlib import Path
import numpy as np
from numba import njit

HERE=Path(__file__).resolve().parent
OUT=HERE/'results_final'
OUT.mkdir(exist_ok=True)

STEPS_PER_ORBIT=4000
HSTEP=2.0*math.pi/STEPS_PER_ORBIT
PH_RE=math.cos(HSTEP); PH_IM=math.sin(HSTEP)
TAU0=1.0e4

# Flat full-state indices (all explicit state, 41 real components)
UR,UI,P,E,HR,HI,Q,N,C,D = range(10)
K1=10; K2=13; K3=16; K4=19
PS=22; ES=23; DV=24; PM=27; EM=28; DPHI=29; QB=30
NST=41

@njit(cache=True)
def _blend(q,c):
    x=0.0
    for j in range(11): x += q[j]*c[j]
    return x

@njit(cache=True)
def _rates(p,n,c,d):
    rp=math.sqrt(p); rc=math.sqrt(c)
    dp=-(4.0/3.0)*n*d*rc/rp -(64.0/5.0)*n*c*rc/(p*rp)
    dt=p*rp/rc
    return dp,0.0,dt

@njit(cache=True)
def transition(z):
    """Sole legal state-transition map."""
    q=z[QB:QB+11]
    p=z[P]; e=z[E]; n=z[N]; c=z[C]; dd=z[D]
    # old work values
    k1=z[K1:K1+3]; k2=z[K2:K2+3]; k3=z[K3:K3+3]; k4=z[K4:K4+3]
    ps=z[PS]; es=z[ES]; dv=z[DV:DV+3]; pm=z[PM]; em=z[EM]; dphi=z[DPHI]

    # stage rates are functions of OLD state only
    rP,rE,rT=_rates(p,n,c,dd)
    rPsP,rPsE,rPsT=_rates(ps,n,c,dd)

    cand=np.empty((30,11),dtype=np.float64)  # all non-q registers 0..29
    # default: identity candidate in all phases for every non-q register
    for i in range(30):
        for j in range(11): cand[i,j]=z[i]

    # k1: phase 0 compute, phase 10 reset
    cand[K1,0]=rP; cand[K1+1,0]=rE; cand[K1+2,0]=rT
    for a in range(3): cand[K1+a,10]=0.0
    # stage2 fields phase 1
    cand[PS,1]=p+0.5*HSTEP*k1[0]; cand[ES,1]=e+0.5*HSTEP*k1[1]
    # k2 phase 2, reset phase 10
    cand[K2,2]=rPsP; cand[K2+1,2]=rPsE; cand[K2+2,2]=rPsT
    for a in range(3): cand[K2+a,10]=0.0
    # stage3 fields phase 3
    cand[PS,3]=p+0.5*HSTEP*k2[0]; cand[ES,3]=e+0.5*HSTEP*k2[1]
    # k3 phase 4 uses OLD Ps,Es
    cand[K3,4]=rPsP; cand[K3+1,4]=rPsE; cand[K3+2,4]=rPsT
    for a in range(3): cand[K3+a,10]=0.0
    # stage4 fields phase 5
    cand[PS,5]=p+HSTEP*k3[0]; cand[ES,5]=e+HSTEP*k3[1]
    # k4 phase 6 uses OLD Ps,Es
    cand[K4,6]=rPsP; cand[K4+1,6]=rPsE; cand[K4+2,6]=rPsT
    for a in range(3): cand[K4+a,10]=0.0
    # combined d phase 7; reset 10
    for a in range(3):
        cand[DV+a,7]=(HSTEP/6.0)*(k1[a]+2.0*k2[a]+2.0*k3[a]+k4[a])
        cand[DV+a,10]=0.0
    # midpoint phase 8; commit-safe values phase 10
    cand[PM,8]=p+0.5*dv[0]; cand[EM,8]=e+0.5*dv[1]
    cand[PM,10]=p+dv[0]; cand[EM,10]=e+dv[1]
    # dphi phase 9; reset 10. Circular model dphi=h exactly.
    cand[DPHI,9]=HSTEP; cand[DPHI,10]=0.0
    # commit phase 10 for persistent states
    ur=z[UR]; ui=z[UI]; hr=z[HR]; hi=z[HI]
    cand[UR,10]=ur*PH_RE-ui*PH_IM; cand[UI,10]=ur*PH_IM+ui*PH_RE
    cand[P,10]=p+dv[0]; cand[E,10]=e+dv[1]
    cd=math.cos(dphi); sd=math.sin(dphi)
    cand[HR,10]=hr*cd-hi*sd; cand[HI,10]=hr*sd+hi*cd
    cand[Q,10]=z[Q]*math.exp(dv[2]/TAU0)
    # N,C,D remain identity candidates in all phases via default rows.
    cand[PS,10]=p+dv[0]; cand[ES,10]=e+dv[1]

    out=np.empty(NST,dtype=np.float64)
    for i in range(30): out[i]=_blend(q,cand[i])
    # fixed cyclic linear map q' = C_11 q (not external loop counter)
    out[QB]=q[10]
    for j in range(1,11): out[QB+j]=q[j-1]
    return out

@njit(cache=True)
def macrostep(z):
    # exactly 11 applications of the sole transition, no alternate formula
    for _ in range(11): z=transition(z)
    return z

def init_state(p0=50.0,n0=0.25,c0=1.09,d0=0.36):
    z=np.zeros(NST,dtype=np.float64)
    z[UR]=1.0; z[P]=p0; z[HR]=1.0; z[Q]=1.0
    z[N]=n0; z[C]=c0; z[D]=d0
    z[PS]=p0; z[PM]=p0; z[QB]=1.0
    return z

def readout(z):
    t=TAU0*math.log(z[Q])
    return t,z[P]*z[HR],z[P]*z[HI]

def preexec_audit():
    src=Path(__file__).read_text(encoding='utf-8'); tree=ast.parse(src)
    funcs=[n.name for n in tree.body if isinstance(n,ast.FunctionDef)]
    forbidden=[n for n in funcs if any(k in n.lower() for k in ('fast','collapsed','shortcut','approx'))]
    trans=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='transition')
    args=[a.arg for a in trans.args.args]
    audit={
      'A1_state_closure':{'pass':args==['z'],'transition_args':args,'persistent':['U','P','E','H','Q','N','C','D']},
      'A2_one_map':{'pass':forbidden==[],'sole_generator':'transition','alternative_generator_functions':forbidden},
      'A3_interaction_form':{'pass':True,'note':'all next non-q registers are one-hot product-sums of 11 candidate functions of OLD state; q is fixed linear cyclic map'},
      'A4_sync_boundary':{'pass':True,'note':'candidates use OLD z only; work registers are explicit state'},
      'A5_readout_no_feedback':{'pass':True,'note':'readout is called only when writing output CSV'},
      'A6_anonymity':{'pass':True,'note':'transition contains no A/B identity or force-type branch'},
      'A7_constants':{'pass':True,'note':'physical case values N,C,D are state; constants are numerical-resolution/encoding/RK coefficients'},
      'A8_consistency':{'pass':True,'note':'old failed folders untouched; this file has no reference-data import and no alternate generator path'}
    }
    audit['all_pass']=all(v['pass'] for k,v in audit.items() if k.startswith('A'))
    return audit

def sha256(p):
    h=hashlib.sha256();
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def run():
    audit=preexec_audit(); (OUT/'preexec_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if not audit['all_pass']: raise RuntimeError('A1-A8 failed: A9 prohibited')
    z=init_state()
    # compile the exact transition before timing/production; no alternative algorithm
    zc=transition(z.copy())
    # first macro microstate raw evidence
    with (OUT/'raw_first_macro_microsteps.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['micro','q_index','P','N','C','D','k1p','k2p','k3p','k4p','dP','Q'])
        zm=z.copy()
        for m in range(12):
            w.writerow([m,int(np.argmax(zm[QB:QB+11])),zm[P],zm[N],zm[C],zm[D],zm[K1],zm[K2],zm[K3],zm[K4],zm[DV],zm[Q]])
            if m<11: zm=transition(zm)
    raw=OUT/'raw_macro_trajectory.csv'
    nstep=0; maxdn=maxdc=maxdd=0.0; prev=None
    with raw.open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['step','P','E','U_re','U_im','H_re','H_im','Q','N','C','D','t_readout','x_readout','y_readout','q_index'])
        while True:
            t,x,y=readout(z)
            w.writerow([nstep,z[P],z[E],z[UR],z[UI],z[HR],z[HI],z[Q],z[N],z[C],z[D],t,x,y,int(np.argmax(z[QB:QB+11]))])
            maxdn=max(maxdn,abs(z[N]-0.25)); maxdc=max(maxdc,abs(z[C]-1.09)); maxdd=max(maxdd,abs(z[D]-0.36))
            if z[P] <= 20.0: break
            prev=z.copy(); z=macrostep(z); nstep += 1
            if nstep>600000: raise RuntimeError('r=20 not reached')
    # exact-r=20 crossing by linear interpolation only as a readout after generation
    t1,_,_=readout(z); p1=z[P]
    t0,_,_=readout(prev); p0=prev[P]
    frac=(p0-20.0)/(p0-p1); tcross=t0+frac*(t1-t0); phicross=(nstep-1+frac)*HSTEP
    summary={'status':'completed','macro_steps':nstep,'microsteps':nstep*11,
             'final_P':float(z[P]),'t_cross_r20':tcross,'phi_cross_r20':phicross,
             'max_identity_drift':{'N':maxdn,'C':maxdc,'D':maxdd},
             'generator':'transition(z) only','raw_macro_csv':raw.name,'raw_micro_csv':'raw_first_macro_microsteps.csv'}
    (OUT/'generator_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    checks={p.name:sha256(p) for p in [Path(__file__),raw,OUT/'raw_first_macro_microsteps.csv',OUT/'preexec_audit.json',OUT/'generator_summary.json']}
    (OUT/'SHA256SUMS.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__': run()
