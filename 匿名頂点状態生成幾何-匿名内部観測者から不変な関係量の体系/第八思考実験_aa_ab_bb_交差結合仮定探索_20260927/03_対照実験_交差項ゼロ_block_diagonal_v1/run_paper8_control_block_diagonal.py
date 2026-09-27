#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 8 control experiment: 123-state full vector with strictly block-diagonal S(Psi).

Purpose
-------
Use the same long state vector Psi=(Psi_aa,Psi_ab,Psi_bb) in R^123 as Paper 8,
but set every cross-channel block of S(Psi) to exactly zero. Each 41x41 diagonal
block is constructed so that S41(z) @ z reproduces the accepted Paper-6/7
transition(z) for that channel.

No cross terms, renormalization, clipping, postprocessing feedback, or alternate
physics are introduced.
"""
from __future__ import annotations
import csv, json, math, hashlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_control'
OUT.mkdir(exist_ok=True)

STEPS_PER_ORBIT=4000
HSTEP=2.0*math.pi/STEPS_PER_ORBIT
PH_RE=math.cos(HSTEP); PH_IM=math.sin(HSTEP)
TAU0=1.0e4
UR,UI,P,E,HR,HI,Q,N,C,D=range(10)
K1=10; K2=13; K3=16; K4=19
PS=22; ES=23; DV=24; PM=27; EM=28; DPHI=29; QB=30
NST=41; NCH=3; NTOT=123
CHANNELS=('aa','ab','bb')


def rates(p,n,c,d):
    rp=math.sqrt(p); rc=math.sqrt(c)
    dp=-(4.0/3.0)*n*d*rc/rp -(64.0/5.0)*n*c*rc/(p*rp)
    dt=p*rp/rc
    return dp,0.0,dt


def transition_reference(z):
    q=z[QB:QB+11]
    p=z[P]; e=z[E]; n=z[N]; c=z[C]; dd=z[D]
    k1=z[K1:K1+3]; k2=z[K2:K2+3]; k3=z[K3:K3+3]; k4=z[K4:K4+3]
    ps=z[PS]; es=z[ES]; dv=z[DV:DV+3]; pm=z[PM]; em=z[EM]; dphi=z[DPHI]
    rP,rE,rT=rates(p,n,c,dd)
    rPsP,rPsE,rPsT=rates(ps,n,c,dd)
    cand=np.empty((30,11),dtype=np.float64)
    for i in range(30):
        for j in range(11): cand[i,j]=z[i]
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
    for i in range(30): out[i]=float(np.dot(q,cand[i]))
    out[QB]=q[10]
    for j in range(1,11): out[QB+j]=q[j-1]
    return out


def build_S41(z):
    """Construct a state-dependent 41x41 matrix satisfying S41(z) @ z = transition_reference(z)."""
    S=np.zeros((NST,NST),dtype=np.float64)
    q=z[QB:QB+11]
    j=int(np.argmax(q))
    # Old-state candidate identity is the default for all logical/work components 0..29.
    for i in range(30): S[i,i]=1.0
    # one-hot phase rotation, exactly linear
    S[QB, QB+10]=1.0
    for k in range(1,11): S[QB+k,QB+k-1]=1.0

    p=z[P]; e=z[E]; n=z[N]; c=z[C]; d=z[D]
    ps=z[PS]; dphi=z[DPHI]

    if j==0:
        dp,_,dt=rates(p,n,c,d)
        S[K1,:]=0; S[K1,P]=dp/p
        S[K1+1,:]=0
        S[K1+2,:]=0; S[K1+2,P]=dt/p
    elif j==1:
        S[PS,:]=0; S[PS,P]=1.0; S[PS,K1]=0.5*HSTEP
        S[ES,:]=0; S[ES,E]=1.0; S[ES,K1+1]=0.5*HSTEP
    elif j==2:
        dp,_,dt=rates(ps,n,c,d)
        S[K2,:]=0; S[K2,PS]=dp/ps
        S[K2+1,:]=0
        S[K2+2,:]=0; S[K2+2,PS]=dt/ps
    elif j==3:
        S[PS,:]=0; S[PS,P]=1.0; S[PS,K2]=0.5*HSTEP
        S[ES,:]=0; S[ES,E]=1.0; S[ES,K2+1]=0.5*HSTEP
    elif j==4:
        dp,_,dt=rates(ps,n,c,d)
        S[K3,:]=0; S[K3,PS]=dp/ps
        S[K3+1,:]=0
        S[K3+2,:]=0; S[K3+2,PS]=dt/ps
    elif j==5:
        S[PS,:]=0; S[PS,P]=1.0; S[PS,K3]=HSTEP
        S[ES,:]=0; S[ES,E]=1.0; S[ES,K3+1]=HSTEP
    elif j==6:
        dp,_,dt=rates(ps,n,c,d)
        S[K4,:]=0; S[K4,PS]=dp/ps
        S[K4+1,:]=0
        S[K4+2,:]=0; S[K4+2,PS]=dt/ps
    elif j==7:
        for a in range(3):
            S[DV+a,:]=0
            S[DV+a,K1+a]=HSTEP/6.0
            S[DV+a,K2+a]=HSTEP/3.0
            S[DV+a,K3+a]=HSTEP/3.0
            S[DV+a,K4+a]=HSTEP/6.0
    elif j==8:
        S[PM,:]=0; S[PM,P]=1.0; S[PM,DV]=0.5
        S[EM,:]=0; S[EM,E]=1.0; S[EM,DV+1]=0.5
    elif j==9:
        S[DPHI,:]=0; S[DPHI,P]=HSTEP/p
    elif j==10:
        # rows explicitly reset to zero at macro boundary
        for base in (K1,K2,K3,K4,DV):
            for a in range(3): S[base+a,:]=0
        S[PM,:]=0; S[PM,P]=1.0; S[PM,DV]=1.0
        S[EM,:]=0; S[EM,E]=1.0; S[EM,DV+1]=1.0
        S[DPHI,:]=0
        S[UR,:]=0; S[UR,UR]=PH_RE; S[UR,UI]=-PH_IM
        S[UI,:]=0; S[UI,UR]=PH_IM; S[UI,UI]=PH_RE
        S[P,:]=0; S[P,P]=1.0; S[P,DV]=1.0
        S[E,:]=0; S[E,E]=1.0; S[E,DV+1]=1.0
        cd=math.cos(dphi); sd=math.sin(dphi)
        S[HR,:]=0; S[HR,HR]=cd; S[HR,HI]=-sd
        S[HI,:]=0; S[HI,HR]=sd; S[HI,HI]=cd
        S[Q,:]=0; S[Q,Q]=math.exp(z[DV+2]/TAU0)
        S[PS,:]=0; S[PS,P]=1.0; S[PS,DV]=1.0
        S[ES,:]=0; S[ES,E]=1.0; S[ES,DV+1]=1.0
    return S


def build_S123(psi):
    S=np.zeros((NTOT,NTOT),dtype=np.float64)
    for ch in range(3):
        a=ch*NST; b=a+NST
        S[a:b,a:b]=build_S41(psi[a:b])
    return S


def init_channel(C0,D0):
    z=np.zeros(NST,dtype=np.float64)
    z[UR]=1.0; z[P]=50.0; z[HR]=1.0; z[Q]=1.0
    z[N]=0.25; z[C]=C0; z[D]=D0
    z[PS]=50.0; z[PM]=50.0; z[QB]=1.0
    return z


def init_psi():
    # aa, ab, bb for n=1: Paper 7 self, Paper 6 attractive, Paper 7 self.
    return np.concatenate([init_channel(0.91,0.0), init_channel(1.09,0.36), init_channel(0.91,0.0)])


def independent_microstep(psi):
    return np.concatenate([transition_reference(psi[i*NST:(i+1)*NST]) for i in range(3)])


def diagnostic_row(micro, psi):
    vals=[micro]
    for ch in range(3):
        z=psi[ch*NST:(ch+1)*NST]
        vals += [int(np.argmax(z[QB:QB+11])),z[P],z[N],z[C],z[D],z[K1],z[K2],z[K3],z[K4],z[DV],z[Q]]
    return vals


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def main():
    psi=init_psi(); ref=psi.copy()
    np.save(OUT/'initial_psi123.npy',psi)

    # Exact construction audit at each microphase of first macrostep.
    micro_rows=[]
    max_matrix_vs_ref=0.0
    max_cross_abs=0.0
    headers=['micro']
    for name in CHANNELS:
        headers += [f'{name}_q',f'{name}_P',f'{name}_N',f'{name}_C',f'{name}_D',f'{name}_k1p',f'{name}_k2p',f'{name}_k3p',f'{name}_k4p',f'{name}_dP',f'{name}_Q']
    for m in range(12):
        micro_rows.append(diagnostic_row(m,psi))
        if m<11:
            S=build_S123(psi)
            # prove off-diagonal blocks are identically zero
            for i in range(3):
                for j in range(3):
                    if i!=j:
                        blk=S[i*NST:(i+1)*NST,j*NST:(j+1)*NST]
                        max_cross_abs=max(max_cross_abs,float(np.max(np.abs(blk))))
            nxt=S@psi
            refn=independent_microstep(ref)
            max_matrix_vs_ref=max(max_matrix_vs_ref,float(np.max(np.abs(nxt-refn))))
            psi=nxt; ref=refn
    with open(OUT/'first_macro_123_diagnostics.csv','w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(headers); w.writerows(micro_rows)
    np.save(OUT/'S123_micro0.npy',build_S123(init_psi()))
    np.save(OUT/'psi123_after_one_macro.npy',psi)

    # Longer control run: 4000 macrosteps (one orbit), compare every microstep to 3 independent systems.
    psi=init_psi(); ref=psi.copy()
    max_diff=0.0; max_diff_step=(-1,-1)
    checkpoint=[]
    for macro in range(4000):
        for micro in range(11):
            nxt=build_S123(psi)@psi
            refn=independent_microstep(ref)
            d=float(np.max(np.abs(nxt-refn)))
            if d>max_diff:
                max_diff=d; max_diff_step=(macro,micro)
            psi=nxt; ref=refn
        if macro in (0,9,99,999,3999):
            checkpoint.append([macro+1,psi[P],psi[NST+P],psi[2*NST+P],max_diff])
    with open(OUT/'one_orbit_checkpoints.csv','w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['macrosteps','P_aa','P_ab','P_bb','max_abs_diff_vs_independent']); w.writerows(checkpoint)
    np.save(OUT/'psi123_after_4000_macrosteps.npy',psi)

    summary={
      'experiment':'Paper 8 control: block-diagonal S(Psi) on the same 123-state vector',
      'state_dimension':123,
      'S_dimension':[123,123],
      'channels':['aa','ab','bb'],
      'initial_relations':{'aa':{'C':0.91,'D':0.0},'ab':{'C':1.09,'D':0.36},'bb':{'C':0.91,'D':0.0}},
      'cross_blocks_exactly_zero': bool(max_cross_abs==0.0),
      'max_cross_block_abs_first_macro':max_cross_abs,
      'first_macro_max_abs_diff_matrix_vs_three_independent':max_matrix_vs_ref,
      'one_orbit_macrosteps':4000,
      'one_orbit_max_abs_diff_matrix_vs_three_independent':max_diff,
      'one_orbit_location_of_max_diff':{'macro':max_diff_step[0],'micro':max_diff_step[1]},
      'final_P_after_4000':{'aa':float(psi[P]),'ab':float(psi[NST+P]),'bb':float(psi[2*NST+P])},
      'claim_scope':'This validates the long-vector / block-diagonal-matrix implementation against three independent copies of the accepted Paper-6/7 transition. It does not add cross physics.'
    }
    (OUT/'control_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__': main()
