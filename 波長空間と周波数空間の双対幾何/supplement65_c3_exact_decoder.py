#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺65 C3確定版: 厳密配置読み — 上位候補の組合せを記録と厳密照合(配置読み=記録が集合を決める)
import numpy as np, random
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('[1]')[0].split('print(')[0])
def sval(v): return sum((abs(t)+0.5)**2 for t in v)
Ng=8; g=np.arange(Ng)/Ng
def axwave(k):
    if k==0: return np.ones(Ng)
    return np.sqrt(2)*(np.cos(2*np.pi*abs(k)*g) if k>0 else np.sin(2*np.pi*abs(k)*g))
def cellwave(v):
    w=axwave(v[0]).reshape(-1,1,1,1)*axwave(v[1]).reshape(1,-1,1,1)
    return w*axwave(v[2]).reshape(1,1,-1,1)*axwave(v[3]).reshape(1,1,1,-1)
ALL=[v for r in (1,3,5,7,9) for v in SH[r]]
CW={v:cellwave(v) for v in ALL}
random.seed(11)
vp=(2,1,0,0); cases=set()
for f in all_finals(9):
    for ((a,b),x,(c,d),d1) in two_step_paths(9,f):
        for va in SH[a]:
            for vb in SH[b]:
                if tphase(vp,va,vb) is None: continue
                kept = vb if x==a else va; vdec = va if x==a else vb
                for vc in SH[c]:
                    for vd in SH[d]:
                        if tphase(vdec,vc,vd) is None: continue
                        if len({kept,vc,vd})==3: cases.add((kept,vc,vd))
cases=list(cases); random.shuffle(cases); cases=cases[:60]
import itertools as it
rec=0; sadd=0; amb=0
for (k1,k2,k3) in cases:
    I=(1.0+CW[k1]+CW[k2]+CW[k3])**2
    scores=sorted(ALL,key=lambda v:-abs(float(np.mean(I*CW[v]))))[:10]
    hits=[]
    for trio in it.combinations(scores,3):
        Ip=(1.0+CW[trio[0]]+CW[trio[1]]+CW[trio[2]])**2
        if np.max(np.abs(Ip-I))<1e-9: hits.append(set(trio))
    if len(hits)==1 and hits[0]=={k1,k2,k3}:
        rec+=1
        if abs(sum(sval(v) for v in hits[0])-9.0)<1e-9: sadd+=1
    elif len(hits)>1: amb+=1
print(f"C3確定版(厳密配置読み): 60終状態 → 一意復元 {rec}/60, s_read加法 {sadd}/60, 多義 {amb}")
