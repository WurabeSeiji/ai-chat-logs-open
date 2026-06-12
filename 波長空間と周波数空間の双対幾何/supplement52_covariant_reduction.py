#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺52 検証: 位相管理連鎖の共変化 — 三案の全数検査
#  (a) η² (Z2商, 向き付け完全不変) のガウス和
#  (b) 向き平均 (Re η1 · Re η2) の和
#  (c) |z|² を保つ B4 元の個数 (位相水準ゲージ群の指数)
import itertools, numpy as np
from collections import defaultdict
SQ2=np.sqrt(2.0)
def expcoeffs(v):
    acc=[((),1.0+0j)]
    for i in range(4):
        m=v[i]
        if m==0: opts=[(0,1.0+0j)]
        elif m>0: opts=[(+m,1/SQ2+0j),(-m,1/SQ2+0j)]
        else:
            mm=-m; opts=[(+mm,-1j/SQ2),(-mm,+1j/SQ2)]
        acc=[(n+(f,),c*w) for n,c in acc for f,w in opts]
    d={}
    for n,c in acc: d[n]=d.get(n,0)+c
    return d
def cross_at(va,vb,t):
    A=expcoeffs(tuple(va)); B=expcoeffs(tuple(vb)); tot=0+0j
    for na,ca in A.items():
        nb=tuple(x-y for x,y in zip(t,na))
        if nb in B: tot+=ca*B[nb]
    return tot
def self_at(v,t): return expcoeffs(tuple(v)).get(tuple(t),0+0j)
def tphase(vp,va,vb):
    vp=tuple(vp); Cp=self_at(vp,vp)
    if abs(Cp)<1e-12: return None
    Cc=cross_at(va,vb,vp)
    if abs(Cc)<1e-12: return None
    return (Cc/abs(Cc))/(Cp/abs(Cp))
def two_step_paths(parent,final):
    paths=[]; fin=tuple(sorted(final))
    for a in range(1,parent+2,2):
        for b in range(a,parent+2,2):
            d1=a+b-parent
            if d1 not in (1,-1): continue
            for (x,spec) in ((a,b),(b,a)):
                for c in range(1,x+2,2):
                    dd=x-d1-c
                    if dd<1 or dd%2==0: continue
                    if tuple(sorted((spec,c,dd)))==fin: paths.append(((a,b),x,(c,dd),d1))
    return sorted(set(paths))
def shell_cells(m):
    out=[]; K=4
    for k in itertools.product(range(-K,K+1),repeat=4):
        if abs(sum((abs(t)+0.5)**2 for t in k)-m)<1e-9: out.append(np.array(k))
    return out
SH={m:shell_cells(float(m)) for m in (1,3,5,7,9)}
LBL={lbl:two_step_paths(9,fin) for lbl,fin in [('531',(5,3,1)),('333',(3,3,3))]}
def run(v9):
    res={'531':[],'333':[]}
    for lbl,lp in LBL.items():
        for ((a,b),x,(c,d),d1) in lp:
            for va in SH[a]:
                for vb in SH[b]:
                    if a==b and tuple(va)>=tuple(vb): continue
                    e1=tphase(v9,va,vb)
                    if e1 is None: continue
                    vx=va if x==a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            if c==d and tuple(vc)>=tuple(vd): continue
                            e2=tphase(vx,vc,vd)
                            if e2 is None: continue
                            res[lbl].append((e1,e2))
    return res
parents=[np.array(v) for v in itertools.product(range(-4,5),repeat=4)
         if sorted(map(abs,v))==[0,0,1,2] and abs(sum((abs(t)+0.5)**2 for t in v)-9)<1e-9]
za_set=set(); zb_set=set(); split_set=set()
for vp in parents:
    r=run(vp); za={}; zb={}; spl={}
    for lbl in ('531','333'):
        e2s=[(e1*e2)**2 for e1,e2 in r[lbl]]
        za[lbl]=complex(np.round(sum(e2s).real,6), np.round(sum(e2s).imag,6))
        plus=sum(1 for v in e2s if abs(v-1)<1e-9)
        spl[lbl]=(plus,len(e2s)-plus)
        zb[lbl]=round(sum(np.real(e1)*np.real(e2) for e1,e2 in r[lbl]),6)
    za_set.add((za['531'],za['333'])); zb_set.add((zb['531'],zb['333'])); split_set.add((spl['531'],spl['333']))
print(f"(a) z(η²): {za_set} / 分布 {split_set}")
print(f"(b) z(ReRe): {zb_set}")
PERMS=list(itertools.permutations(range(4)))
SIGNS=list(itertools.product([1,-1],repeat=4))
def G(p,s,v): return np.array([s[i]*v[p[i]] for i in range(4)])
v9=np.array((2,1,0,0))
r0=run(v9); base=tuple(round(abs(sum(e1*e2 for e1,e2 in r0[lbl]))**2,3) for lbl in ('531','333'))
keep=0
for p in PERMS:
    for s in SIGNS:
        vp2=G(p,s,v9)
        if sorted(map(abs,vp2))!=[0,0,1,2]: continue
        r2=run(vp2)
        val=tuple(round(abs(sum(e1*e2 for e1,e2 in r2[lbl]))**2,3) for lbl in ('531','333'))
        if val==base: keep+=1
print(f"(c) 同値 |z|² に移す B4 元: {keep}/384 → 指数 {384//keep}")
