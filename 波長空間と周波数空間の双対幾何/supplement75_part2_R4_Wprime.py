#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺75 部2: R4検査 + 配置粒度 W' 系列 (Δ'(s), s=11/13)
import itertools, time
from collections import Counter, defaultdict
from fractions import Fraction as Fr
from math import comb
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
CAP={1:1,3:8,5:24,7:40,9:64,11:96,13:96}
def M_channel(f,used):
    tot=1
    for m,k in Counter(f).items():
        free=CAP[m]-used.get(m,0)
        if free<k: return 0
        tot*=comb(free,k)
    return tot
def seq_twin(s):
    fins=all_finals(s)
    M1={f:M_channel(f,{}) for f in fins}; T1=sum(M1.values())
    P=defaultdict(lambda: Fr(0))
    for f1 in fins:
        if M1[f1]==0: continue
        used=Counter(f1)
        M2={f:M_channel(f,used) for f in fins}; T2=sum(M2.values())
        for f2 in fins:
            if M2[f2]==0: continue
            P[tuple(sorted([f1,f2]))]+=Fr(M1[f1],T1)*Fr(M2[f2],T2)
    return dict(P)
def derived_twin(s,Wd):
    fins=all_finals(s); P={}; tot=0.0
    for i,f1 in enumerate(fins):
        for f2 in fins[i:]:
            joint=Counter(f1)+Counter(f2)
            if any(CAP[m]<k for m,k in joint.items()): continue
            w=(2 if f1!=f2 else 1)*Wd.get(f1,0)*Wd.get(f2,0)
            if w>0: P[tuple(sorted([f1,f2]))]=w; tot+=w
    return {k:v/tot for k,v in P.items()}
def TV(P,Q):
    ks=set(P)|set(Q)
    return sum(abs(float(P.get(k,0))-float(Q.get(k,0))) for k in ks)/2
def run_config(vp,s):
    """チャネル一貫 W と 配置粒度 W' を同時計算"""
    fins=all_finals(s)
    zch={f:0+0j for f in fins}; zcf={f:defaultdict(complex) for f in fins}
    for f in fins:
        for ((a,b),x,(c,d),d1) in two_step_paths(s,f):
            for ia,va in enumerate(SH[a]):
                for ib,vb in enumerate(SH[b]):
                    if a==b and not ia<ib: continue
                    e1=tphase(vp,va,vb)
                    if e1 is None: continue
                    pairs=[(va,vb),(vb,va)] if a==b else [((va if x==a else vb),(vb if x==a else va))]
                    for vx,kept in pairs:
                        for ic,vc in enumerate(SH[c]):
                            for idd,vd in enumerate(SH[d]):
                                if c==d and not ic<idd: continue
                                e2=tphase(vx,vc,vd)
                                if e2 is None: continue
                                eta=e1*e2
                                zch[f]+=eta
                                K=frozenset([kept,vc,vd])
                                if len(K)==3: zcf[f][K]+=eta
    W={f:abs(zch[f])**2 for f in fins}
    Wp={f:sum(abs(z)**2 for z in zcf[f].values()) for f in fins}
    ncf={f:len([z for z in zcf[f].values() if abs(z)>1e-9]) for f in fins}
    return W,Wp,ncf
print("==== R4a: (2,2,0,0) 3セクターの導出分布の比較 / R4b: (3,0,0,0)B の W 縮退 ====")
seq13=seq_twin(13)
t0=time.time()
res={}
for name,vp in [('(2,2,0,0)opp',(-2,2,0,0)),('(2,2,0,0)++',(2,2,0,0)),('(2,2,0,0)--',(-2,-2,0,0)),
                ('(3,0,0,0)B',(-3,0,0,0)),
                ('(2,1,1,0)A_s11',None)]:
    if vp is None: continue
    W,Wp,ncf=run_config(vp,13)
    res[name]=(W,Wp,ncf)
    nz={f:round(w,1) for f,w in W.items() if w>1e-6}
    print(f"  {name}: W={nz}")
    d=derived_twin(13,W); dp=derived_twin(13,Wp)
    print(f"    Δ(13)={TV(d,seq13):.5f}, Δ'(13)[配置粒度]={TV(dp,seq13):.5f}, W'非零配置数={ {f:n for f,n in ncf.items() if n>0} } ({time.time()-t0:.0f}s)")
# R4a: 3セクターの導出分布が一致するか
d1=derived_twin(13,res['(2,2,0,0)opp'][0]); d2=derived_twin(13,res['(2,2,0,0)++'][0]); d3=derived_twin(13,res['(2,2,0,0)--'][0])
print(f"  R4a: TV(opp,++)={TV(d1,d2):.6f}, TV(opp,--)={TV(d1,d3):.6f} → {'導出分布まで一致' if max(TV(d1,d2),TV(d1,d3))<1e-6 else 'Δのみ一致(分布は別)'}")
print()
print("==== s=11 の W'/Δ' ====")
seq11=seq_twin(11)
for name,vp in [('s11 A',(2,1,1,0)),('s11 B',(-2,1,1,0))]:
    W,Wp,ncf=run_config(vp,11)
    d=derived_twin(11,W); dp=derived_twin(11,Wp)
    print(f"  {name}: W={ {f:round(w) for f,w in W.items() if w>1e-6} }")
    print(f"    W'={ {f:round(w) for f,w in Wp.items() if w>1e-6} } 配置数={ {f:n for f,n in ncf.items() if n>0} }")
    print(f"    Δ={TV(d,seq11):.5f} / Δ'(配置粒度)={TV(dp,seq11):.5f}")
