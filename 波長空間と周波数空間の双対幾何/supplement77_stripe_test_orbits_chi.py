#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺77: (B)縞クラステスト即答 (C)21/13・11/8 定義照合 (D)χ(K)指標直交性
import itertools
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])

# ---- (C) B4×swap 軌道数の定義照合 ----
cells48=[v for v in SH[9] if tuple(sorted(map(abs,v)))==(0,0,1,2)]
PERMS=list(itertools.permutations(range(4)))
SIGNS=list(itertools.product([1,-1],repeat=4))
G=[(p,s) for p in PERMS for s in SIGNS]
def act(g,v):
    p,s=g
    return tuple(s[i]*v[p[i]] for i in range(4))
def canon_pair(u,v,unordered=True):
    best=None
    for g in G:
        a,b=act(g,u),act(g,v)
        for pair in ([ (a,b),(b,a) ] if unordered else [(a,b)]):
            if best is None or pair<best: best=pair
    return best
def count_orbits(unordered, diag):
    seen=set()
    for i,u in enumerate(cells48):
        for j,v in enumerate(cells48):
            if unordered and j<i: continue
            if not diag and i==j: continue
            seen.add(canon_pair(u,v,unordered))
    return len(seen)
for unordered in (True,False):
    for diag in (True,False):
        n=count_orbits(unordered,diag)
        print(f"(C) B4×swap軌道: {'無順序' if unordered else '順序'}・対角{'込' if diag else '除'} → {n}")

# ---- (D) χ(K) 指標直交性: (±3,0,0,0)@13 チャネル(3,5,5) ----
def config_z(vp,s,final):
    z=defaultdict(complex)
    for ((a,b),x,(c,d),d1) in two_step_paths(s,final):
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
                            K=frozenset([kept,vc,vd])
                            if len(K)==3: z[K]+=e1*e2
    return z
zA=config_z((3,0,0,0),13,(3,5,5)); zB=config_z((-3,0,0,0),13,(3,5,5))
common=[K for K in zA if abs(zA[K])>1e-9 and K in zB and abs(zB[K])>1e-9]
chis=defaultdict(int); modeq=0
for K in common:
    r=zB[K]/zA[K]
    key=(round(r.real,6),round(r.imag,6))
    chis[key]+=1
    if abs(abs(zB[K])-abs(zA[K]))<1e-9: modeq+=1
print(f"(D) (±3,0,0,0)@13 (3,5,5): 共通配置 {len(common)}, |z|一致 {modeq}, χ分布 {dict(chis)}")
sA=sum(zA.values()); sB=sum(zB.values())
print(f"    チャネル和: z^A={sA:.1f}, z^B={sB:.3f} (Bは厳密零か)")
# 指標直交性: Σ χ(K) z^A_K = z^B = 0
# (2,2,0,0) 3セクターの配置比較
z1=config_z((-2,2,0,0),13,(1,3,9)); z2=config_z((2,2,0,0),13,(1,3,9)); z3=config_z((-2,-2,0,0),13,(1,3,9))
m1=sorted(round(abs(z),3) for z in z1.values() if abs(z)>1e-9)
m2=sorted(round(abs(z),3) for z in z2.values() if abs(z)>1e-9)
m3=sorted(round(abs(z),3) for z in z3.values() if abs(z)>1e-9)
print(f"(D2) (2,2,0,0)セクターの(1,3,9)配置|z|多重集合: opp{len(m1)}個 / ++{len(m2)}個 / --{len(m3)}個")
print(f"     opp={m1}")
print(f"     ++ ={m2}")
print(f"     -- ={m3}")
