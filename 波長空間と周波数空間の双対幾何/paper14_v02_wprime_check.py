#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 論文14 v0.2 範囲限定付記の検証: s=9 両セクターの W(チャネル一貫) vs W'(配置粒度)
# 結果: W'_531 はセクター非依存 (576=576)、W'_333 は依存 (40 vs 16)
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
def config_W(vp,s):
    fins=all_finals(s)
    zcf={f:defaultdict(complex) for f in fins}; zch={f:0+0j for f in fins}
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
                                zch[f]+=e1*e2
                                K=frozenset([kept,vc,vd])
                                if len(K)==3: zcf[f][K]+=e1*e2
    return ({f:abs(z)**2 for f,z in zch.items()},
            {f:sum(abs(z)**2 for z in zk.values()) for f,zk in zcf.items()})
for vp,nm in (((2,1,0,0),'A(自明)'),((-2,1,0,0),'B(ねじれ)')):
    W,Wp=config_W(vp,9)
    print(nm, {f:round(w) for f,w in W.items() if w>1e-9}, {f:round(w) for f,w in Wp.items() if w>1e-9})
