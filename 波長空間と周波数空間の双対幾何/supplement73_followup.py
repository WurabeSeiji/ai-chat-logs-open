#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺73 追補: (1) 衝突検査の修正 (2) 結合振幅の vp2 全掃引 + インコヒーレント版
import itertools
from collections import Counter, defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
# (1) 衝突検査修正: 結合多重集合 = sorted(f1+f2)
for s in (9,11,13):
    fins=all_finals(s); jm=defaultdict(list)
    for i,f1 in enumerate(fins):
        for f2 in fins[i:]:
            jm[tuple(sorted(f1+f2))].append((f1,f2))
    coll=[v for v in jm.values() if len(v)>1]
    print(f"s={s} 結合多重集合の真の衝突: {len(coll)} 件 {coll if coll else ''}")
# (2) vp2 全掃引
def eta_paths_full(vp,parent,final):
    out=[]
    for ((a,b),x,(c,d),d1) in two_step_paths(parent,final):
        for ia,va in enumerate(SH[a]):
            for ib,vb in enumerate(SH[b]):
                if a==b and not ia<ib: continue
                e1=tphase(vp,va,vb)
                if e1 is None: continue
                for vx,kept in ([(va,vb),(vb,va)] if a==b else [((va if x==a else vb),(vb if x==a else va))]):
                    for ic,vc in enumerate(SH[c]):
                        for idd,vd in enumerate(SH[d]):
                            if c==d and not ic<idd: continue
                            e2=tphase(vx,vc,vd)
                            if e2 is None: continue
                            out.append((e1*e2, frozenset([kept,vc,vd])))
    return out
vp1=(2,1,0,0)
P1={f:eta_paths_full(vp1,9,f) for f in [(1,3,5),(3,3,3)]}
typeA=[v for v in SH[9] if sorted(map(abs,v),reverse=True)==[2,1,0,0] and v[list(map(abs,v)).index(2)]>0]
import numpy as np
vals=[]; vals_inc=[]
for vp2 in SH[9]:
    if sorted(map(abs,vp2),reverse=True)!=[2,1,0,0]: continue
    if vp2==vp1: continue
    dom=max(range(4),key=lambda j:abs(vp2[j]))
    if vp2[dom]<0: continue   # セクターA のみ
    P2={f:eta_paths_full(vp2,9,f) for f in [(1,3,5),(3,3,3)]}
    zc=defaultdict(complex); inc=defaultdict(float)
    for F1 in P1:
        for F2 in P2:
            ch=tuple(sorted([F1,F2]))
            for (n1,s1) in P1[F1]:
                if len(s1)<3: continue
                for (n2,s2) in P2[F2]:
                    if len(s2)<3: continue
                    if s1 & s2: continue
                    zc[(frozenset([s1,s2]),ch)]+=n1*n2
                    inc[ch]+=abs(n1*n2)**2
    Pj=defaultdict(float)
    for (cfg,ch),z in zc.items(): Pj[ch]+=abs(z)**2
    tot=sum(Pj.values()); ti=sum(inc.values())
    XJ=tuple(sorted([(1,3,5),(3,3,3)]))
    vals.append(Pj[XJ]/tot if tot>0 else float('nan'))
    vals_inc.append(inc[XJ]/ti if ti>0 else float('nan'))
vals=np.array(vals); vals_inc=np.array(vals_inc)
print(f"vp2 セクターA 全{len(vals)}セル掃引:")
print(f"  コヒーレント P(X): 平均 {np.nanmean(vals):.5f}, 最小 {np.nanmin(vals):.5f}, 最大 {np.nanmax(vals):.5f}, 標準偏差 {np.nanstd(vals):.5f}")
hist=Counter(round(v,3) for v in vals)
print(f"  値の分布: {dict(sorted(hist.items()))}")
print(f"  インコヒーレント版 P(X): 平均 {np.nanmean(vals_inc):.5f}, 範囲 [{np.nanmin(vals_inc):.5f},{np.nanmax(vals_inc):.5f}]")
print(f"  D3 因子化 = 0.98195 との比較")
