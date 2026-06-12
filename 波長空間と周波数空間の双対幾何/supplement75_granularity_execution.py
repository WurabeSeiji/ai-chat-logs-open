#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺75: 補遺74 §4 の実行 — 合併キー掃引・関係クラス分類・A_seq・粒度系列 + R4検査
import itertools, time
from collections import Counter, defaultdict
from fractions import Fraction as Fr
from math import comb
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
def sval(v): return sum((abs(t)+0.5)**2 for t in v)
CAP={1:1,3:8,5:24,7:40,9:64,11:96,13:96}
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
FINS9=[(1,3,5),(3,3,3)]
P1={f:eta_paths_full(vp1,9,f) for f in FINS9}

# ---- §4-4(i): 衝突部分和=0 (第五の相殺候補) ----
coll=sum(n for n,s_ in P1[(3,3,3)] if len(s_)<3)
print(f"§4-4i 333 衝突経路(16本)の部分和: |Σ| = {abs(coll):.2e} → {'厳密相殺' if abs(coll)<1e-9 else '非零(位相回転)'}")

# ---- §4-4(ii): 配置粒度 W' の検証 (531: 8配置/576, 333: 2配置/40, 双子式 0.96644) ----
def config_z(paths):
    z=defaultdict(complex)
    for n,s_ in paths:
        if len(s_)<3: continue
        z[s_]+=n
    return z
z531=config_z(P1[(1,3,5)]); z333=config_z(P1[(3,3,3)])
Wp531=sum(abs(z)**2 for z in z531.values()); Wp333=sum(abs(z)**2 for z in z333.values())
tw=2*Wp531*Wp333/(2*Wp531*Wp333+Wp333**2)
print(f"§4-4ii 配置粒度: 531 配置{len(z531)} W'={Wp531:.0f} (期待8/576), 333 配置{len(z333)} W'={Wp333:.0f} (期待2/40), 双子式={tw:.5f} (期待0.96644)")

# ---- §4-1+2: 合併キー版(一回読み)掃引 + 関係クラス分類 ----
print()
print("§4-1/2 合併キー(一回読み)掃引 + 関係クラス分類 [vp2=セクターA 23セル]:")
def invariants(v2):
    dot=sum(a*b for a,b in zip(vp1,v2))
    sp=sval(tuple(a+b for a,b in zip(vp1,v2))); sm=sval(tuple(a-b for a,b in zip(vp1,v2)))
    ov=sum(1 for a,b in zip(vp1,v2) if a!=0 and b!=0)
    return (dot,round(min(sp,sm),1),round(max(sp,sm),1),ov)
classes=defaultdict(list)
for vp2 in SH[9]:
    if sorted(map(abs,vp2),reverse=True)!=[2,1,0,0] or vp2==vp1: continue
    dom=max(range(4),key=lambda j:abs(vp2[j]))
    if vp2[dom]<0: continue
    P2={f:eta_paths_full(vp2,9,f) for f in FINS9}
    # 分割キー(逐次的) と 合併キー(一回読み) を同時集計
    z_part=defaultdict(complex); z_merge=defaultdict(complex)
    for F1 in FINS9:
        for F2 in FINS9:
            for (n1,s1) in P1[F1]:
                if len(s1)<3: continue
                for (n2,s2) in P2[F2]:
                    if len(s2)<3: continue
                    if s1 & s2: continue
                    u=frozenset(s1|s2)
                    z_part[(frozenset([s1,s2]) if s1!=s2 else s1, tuple(sorted([F1,F2])))]+=n1*n2
                    z_merge[u]+=n1*n2
    def PX_from(zdict, merged):
        Pch=defaultdict(float)
        for k,z in zdict.items():
            if merged:
                shells=tuple(sorted(int(2*( sval(c) )**0)*0+round(sval(c)) for c in k))  # シェル多重集合
                X = (5 in shells)
            else:
                X = ((1,3,5) in k[1])
            Pch[X]+=abs(z)**2
        tot=sum(Pch.values())
        return (Pch[True]/tot if tot>0 else float('nan')), tot
    pxp,totp=PX_from(z_part,False)
    pxm,totm=PX_from(z_merge,True)
    classes[invariants(vp2)].append((vp2,pxp,pxm))
for inv,members in sorted(classes.items()):
    vals_p=[round(m[1],5) for m in members]; vals_m=[round(m[2],5) for m in members]
    const_p=len(set(map(str,vals_p)))==1; const_m=len(set(map(str,vals_m)))==1
    print(f"  クラス inv={inv} ({len(members)}対): 分割キーP(X)={vals_p[0] if const_p else vals_p} {'[一定]' if const_p else '[非一定!]'} / 合併キー={vals_m[0] if const_m else vals_m} {'[一定]' if const_m else '[非一定!]'}")

# ---- §4-3: A_seq (逐次振幅) ----
print()
print("§4-3 A_seq (第一段=配置粒度コヒーレント選択+追記, 第二段=排他条件付き):")
import numpy as np
vals=[]
for vp2 in SH[9]:
    if sorted(map(abs,vp2),reverse=True)!=[2,1,0,0] or vp2==vp1: continue
    dom=max(range(4),key=lambda j:abs(vp2[j]))
    if vp2[dom]<0: continue
    P2={f:eta_paths_full(vp2,9,f) for f in FINS9}
    z1={F:config_z(P1[F]) for F in FINS9}
    T1=sum(abs(z)**2 for F in FINS9 for z in z1[F].values())
    PX=0.0
    for F1 in FINS9:
        for K1,za in z1[F1].items():
            p1=abs(za)**2/T1
            z2={F:defaultdict(complex) for F in FINS9}
            for F2 in FINS9:
                for (n2,s2) in P2[F2]:
                    if len(s2)<3 or (s2 & K1): continue
                    z2[F2][s2]+=n2
            T2=sum(abs(z)**2 for F in FINS9 for z in z2[F].values())
            if T2==0: continue
            for F2 in FINS9:
                pF2=sum(abs(z)**2 for z in z2[F2].values())/T2
                if (1,3,5) in (F1,F2): PX+=p1*pF2
    vals.append(PX)
vals=np.array(vals)
print(f"  A_seq P(X): 平均 {np.nanmean(vals):.5f}, 範囲 [{np.nanmin(vals):.5f},{np.nanmax(vals):.5f}] vs D1(逐次数え上げ)=0.98263, D3=0.98195, 配置粒度床=0.96644")
