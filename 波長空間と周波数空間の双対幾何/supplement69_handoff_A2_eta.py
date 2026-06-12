#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺69 §A-2: η² 等分布の検証と η→±iη 対合の計算的探索 (補遺68 §A-2 引き継ぎ)
import numpy as np, itertools
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
def eta_paths(vp,parent,final):
    out=[]
    for ((a,b),x,(c,d),d1) in two_step_paths(parent,final):
        for ia,va in enumerate(SH[a]):
            for ib,vb in enumerate(SH[b]):
                if a==b and not ia<ib: continue
                e1=tphase(vp,va,vb)
                if e1 is None: continue
                vxs=[va,vb] if a==b else [va if x==a else vb]
                for vx in vxs:
                    for ic,vc in enumerate(SH[c]):
                        for idd,vd in enumerate(SH[d]):
                            if c==d and not ic<idd: continue
                            e2=tphase(vx,vc,vd)
                            if e2 is None: continue
                            out.append((e1*e2,(va,vb,vx,vc,vd,(a,b),x,(c,d))))
    return out
def z4(e):
    for p,u in enumerate([1,1j,-1,-1j]):
        if abs(e-u)<1e-9: return p
    return None
vp=(2,1,0,0)
P531=eta_paths(vp,9,(1,3,5)); P333=eta_paths(vp,9,(3,3,3))
for name,P in (("531",P531),("333",P333)):
    sq=defaultdict(int)
    for e,_ in P: sq[z4(e*e)]+=1
    print(f"{name}: 経路 {len(P)}, η²分布 {{+1:{sq[0]}, -1:{sq[2]}}} (i:{sq[1]}, -i:{sq[3]})")
# ---- 対合探索: 経路集合上の写像 σ で η→±iη ----
def keyof(t): return t[:5]
def try_sigma(P, sigma, name):
    """sigma: (va,vb,vx,vc,vd,(a,b),x,(c,d)) → 新タプル or None。η変化を集計"""
    idx={}
    for e,t in P: idx[keyof(t)]=e
    n_in=0; ratios=defaultdict(int); n_fix=0; n_invol=0
    for e,t in P:
        t2=sigma(t)
        if t2 is None: continue
        k2=keyof(t2)
        if k2 not in idx: continue
        n_in+=1
        if k2==keyof(t): n_fix+=1
        r=idx[k2]/e
        ratios[z4(r)]+=1
        t3=sigma(t2+(t[5],t[6],t[7])) if len(t2)==5 else None
    print(f"  σ={name}: 像が集合内 {n_in}/{len(P)}, 固定点 {n_fix}, η'/η 分布 {{1:{ratios[0]}, i:{ratios[1]}, -1:{ratios[2]}, -i:{ratios[3]}}}")
    return n_in,ratios
def flip_axis(v,j): return tuple(-v[i] if i==j else v[i] for i in range(4))
def first_nz(v):
    for j in range(4):
        if v[j]!=0: return j
    return None
# 候補1: vc の最初の非零軸の枝反転 (cos↔sin = ¼シフト)
def s1(t):
    va,vb,vx,vc,vd=t[:5]
    j=first_nz(vc)
    if j is None: return None
    return (va,vb,vx,flip_axis(vc,j),vd)
# 候補2: vd 側
def s2(t):
    va,vb,vx,vc,vd=t[:5]
    j=first_nz(vd)
    if j is None: return None
    return (va,vb,vx,vc,flip_axis(vd,j))
# 候補3: vx の最初の非零軸反転 (中間セルの¼シフト)
def s3(t):
    va,vb,vx,vc,vd=t[:5]
    j=first_nz(vx)
    if j is None: return None
    vx2=flip_axis(vx,j)
    # vx は va/vb のどちらかと同一: 対応する初段セルも同時に反転
    va2 = vx2 if va==vx else va
    vb2 = vx2 if vb==vx else vb
    return (va2,vb2,vx2,vc,vd)
# 候補4: 親の支配軸(axis0)で vc を反転
def s4(t):
    va,vb,vx,vc,vd=t[:5]
    if vc[0]==0: return None
    return (va,vb,vx,flip_axis(vc,0),vd)
# 候補5: 全セルの axis0 反転 (鏡映)
def s5(t):
    va,vb,vx,vc,vd=t[:5]
    return tuple(flip_axis(w,0) for w in (va,vb,vx,vc,vd))
for name,P in (("531",P531),("333",P333)):
    print(f"[{name}]")
    for s,nm in ((s1,'vc枝反転(第1非零軸)'),(s2,'vd枝反転'),(s3,'vx枝反転(初段連動)'),(s4,'vc axis0反転'),(s5,'全セルaxis0鏡映')):
        try_sigma(P,s,nm)

print()
print("==== 修正版対合 σ*: (vc≠0なら vc、原点なら vd)の第1非零軸の枝反転、娘は無順序照合 ====")
def canon_key(t):
    va,vb,vx,vc,vd=t[:5]
    return (va,vb,vx,frozenset([vc,vd]) if vc!=vd else (vc,vd))
def sigma_star(t):
    va,vb,vx,vc,vd=t[:5]
    if any(vc): j=first_nz(vc); return (va,vb,vx,flip_axis(vc,j),vd)
    j=first_nz(vd); return (va,vb,vx,vc,flip_axis(vd,j))
for name,P in (("531",P531),("333",P333)):
    idx={}
    for e,t in P: idx[canon_key(t)]=(e,t)
    n=len(P); n_in=0; n_fix=0; ratios=defaultdict(int); invol_ok=0
    for e,t in P:
        t2=sigma_star(t)
        k2=canon_key(t2)
        if k2 not in idx: continue
        n_in+=1
        e2_,t2stored=idx[k2]
        if k2==canon_key(t): n_fix+=1
        ratios[z4(e2_/e)]+=1
        t3=sigma_star(t2stored)
        if canon_key(t3)==canon_key(t): invol_ok+=1
    print(f"{name}: 完全性 {n_in}/{n}, 固定点 {n_fix}, σ²=id {invol_ok}/{n_in}, η'/η {{1:{ratios[0]}, i:{ratios[1]}, -1:{ratios[2]}, -i:{ratios[3]}}}")

print()
print("==== 区分的対合 σ**: c≠d→(vc|vd)枝反転 / c=d→反転不変選択で娘反転 / ±単軸対→保持娘反転 ====")
def abspat(v): return tuple(abs(x) for x in v)
def pick_daughter(vc,vd):
    """反転不変な娘選択: first_nz 小→その娘; 同じなら |成分| 辞書順; 完全同型(±対)なら None"""
    jc,jd=first_nz(vc),first_nz(vd)
    if jc!=jd: return 0 if jc<jd else 1
    if abspat(vc)!=abspat(vd): return 0 if abspat(vc)<abspat(vd) else 1
    return None
def sigma2(t):
    va,vb,vx,vc,vd,(a,b),x,(c,d)=t
    kept_is_b = (x==a)
    if c!=d:
        if any(vc): j=first_nz(vc); return (va,vb,vx,flip_axis(vc,j),vd),'cd'
        j=first_nz(vd); return (va,vb,vx,vc,flip_axis(vd,j)),'cd'
    p=pick_daughter(vc,vd)
    if p is not None:
        w=(vc,vd)[p]; j=first_nz(w)
        w2=flip_axis(w,j)
        return ((va,vb,vx,w2,vd) if p==0 else (va,vb,vx,vc,w2)),'cd'
    # ±単軸対: 保持娘(e1のみに入る)を反転
    if kept_is_b:
        j=first_nz(vb); return (va,flip_axis(vb,j),vx,vc,vd),'kept'
    j=first_nz(va); return (flip_axis(va,j),vb,vx,vc,vd),'kept'
def key2(t,cd_unord):
    va,vb,vx,vc,vd=t[:5]
    return (va,vb,vx,frozenset([vc,vd])) if cd_unord else (va,vb,vx,vc,vd)
for name,P in (("531",P531),("333",P333)):
    idx={}
    for e,t in P:
        c,d=t[7]
        idx[key2(t[:5],c==d)]=(e,t)
    n=len(P); n_in=0; n_fix=0; ratios=defaultdict(int); invol=0
    for e,t in P:
        c,d=t[7]
        t2,_=sigma2(t)
        k2=key2(t2,c==d)
        if k2 not in idx: continue
        n_in+=1
        e2_,t2s=idx[k2]
        if k2==key2(t[:5],c==d): n_fix+=1
        ratios[z4(e2_/e)]+=1
        t3,_=sigma2(t2s)
        if key2(t3,c==d)==key2(t[:5],c==d): invol+=1
    print(f"{name}: 完全性 {n_in}/{n}, 固定点 {n_fix}, σ²=id {invol}/{n_in}, η'/η {{1:{ratios[0]}, i:{ratios[1]}, -1:{ratios[2]}, -i:{ratios[3]}}}")
