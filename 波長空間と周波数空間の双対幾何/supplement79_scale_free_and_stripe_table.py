#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺79: (1) 粒度スケール自由性テスト 単軸親系列 m=1,2 (予言: 鏡像で W=0 ∧ W'=W'(A)>0)
#         (2) 21軌道→縞クラス対応表 (64域、併合対の全数リスト)
import itertools
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
def run_config(vp,s):
    fins=all_finals(s)
    zch={f:0+0j for f in fins}; zcf={f:defaultdict(complex) for f in fins}
    for f in fins:
        for ((a,b),x,(c,d),d1) in two_step_paths(s,f):
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
                                eta=e1*e2; zch[f]+=eta
                                K=frozenset([kept,vc,vd])
                                if len(K)==3: zcf[f][K]+=eta
    W={f:abs(zch[f])**2 for f in fins}
    Wp={f:sum(abs(z)**2 for z in zcf[f].values()) for f in fins}
    return W,Wp
print("==== (1) スケール自由性: 単軸親系列 s(m)=m²+m+1 ====")
for m,s in ((1,3),(2,7),(3,13)):
    vA=tuple([m]+[0]*3); vB=tuple([-m]+[0]*3)
    WA,WpA=run_config(vA,s); WB,WpB=run_config(vB,s)
    nzWA={f:round(w,1) for f,w in WA.items() if w>1e-9}
    nzWB={f:round(w,1) for f,w in WB.items() if w>1e-9}
    nzpA={f:round(w,1) for f,w in WpA.items() if w>1e-9}
    nzpB={f:round(w,1) for f,w in WpB.items() if w>1e-9}
    pred = (not nzWB) and nzpB and all(abs(WpA[f]-WpB[f])<1e-6 for f in WpA)
    print(f"m={m} (s={s}): W(A)={nzWA} W(B)={nzWB if nzWB else '0(全チャネル)'}")
    print(f"          W'(A)={nzpA} W'(B)={nzpB} → 予言(W_B=0 ∧ W'_B=W'_A>0): {'成立' if pred else '不成立'}")
print()
print("==== (2) 21軌道→縞クラス対応 (64域・無順序・対角除) ====")
PERMS=list(itertools.permutations(range(4)))
SIGNS=list(itertools.product([1,-1],repeat=4))
G=[(p,s_) for p in PERMS for s_ in SIGNS]
def act(g,v):
    p,s_=g
    return tuple(s_[i]*v[p[i]] for i in range(4))
def canon_pair(u,v):
    best=None
    for g in G:
        a,b=act(g,u),act(g,v)
        for pr in ((a,b),(b,a)):
            if best is None or pr<best: best=pr
    return best
def unsigned(w):
    """線 [w]={w,-w} の正準形"""
    return max(w,tuple(-x for x in w))
def stripe_canon(u,v):
    """強度縞の可読データ: 無順序の非符号線対 {[u+v],[u-v]} の B4 正準形"""
    s1=tuple(a+b for a,b in zip(u,v)); s2=tuple(a-b for a,b in zip(u,v))
    best=None
    for g in G:
        l1=unsigned(act(g,s1)); l2=unsigned(act(g,s2))
        pr=tuple(sorted([l1,l2]))
        if best is None or pr<best: best=pr
    return best
cells64=SH[9]
orb=defaultdict(list); stripe_of_orb={}
for i,u in enumerate(cells64):
    for j in range(i+1,len(cells64)):
        v=cells64[j]
        o=canon_pair(u,v)
        orb[o].append((u,v))
for o,members in orb.items():
    u,v=members[0]
    stripe_of_orb[o]=stripe_canon(u,v)
    # 軌道内で縞が一定かの確認
    assert all(stripe_canon(a,b)==stripe_of_orb[o] for a,b in members[:5])
stripes=defaultdict(list)
for o,st in stripe_of_orb.items(): stripes[st].append(o)
print(f"軌道数 {len(orb)} → 縞クラス数 {len(stripes)}")
merged=[v for v in stripes.values() if len(v)>1]
print(f"併合のある縞クラス: {len(merged)} 組")
for i,(st,orbs) in enumerate(sorted(stripes.items(), key=lambda kv:-len(kv[1]))):
    if len(orbs)>1:
        reps=[orb[o][0] for o in orbs]
        sizes=[len(orb[o]) for o in orbs]
        print(f"  縞クラス{i}: 軌道{len(orbs)}個が併合 サイズ{sizes}")
        for r in reps: print(f"    代表対: {r[0]} , {r[1]}")
# 48域でも縞クラス数を確認 (claude.ai の8と照合)
cells48=[v for v in cells64 if tuple(sorted(map(abs,v)))==(0,0,1,2)]
orb48=set(); st48=set()
for i,u in enumerate(cells48):
    for j in range(i+1,len(cells48)):
        orb48.add(canon_pair(u,cells48[j])); st48.add(stripe_canon(u,cells48[j]))
print(f"48域: 軌道 {len(orb48)}, 縞クラス {len(st48)} (claude.ai: 13/8 と照合)")
