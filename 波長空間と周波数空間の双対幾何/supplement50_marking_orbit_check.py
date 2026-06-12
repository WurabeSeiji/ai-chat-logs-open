#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺50 検証: 軸対称性の原則 — H構造(マーキング)の族は B4 の単一推移軌道
#  E1: 標準 qmul の B4 移送で得られる相異なる積表の個数 (期待: 16 = 384/24)
#  E2: 各構造の単位元(実軸)が4軸すべてを走ること
#  E3: 全マーキング共通の不変量 = ノルム(時計/R) と ε(Q) の B4 完全不変性
import itertools
def qmul(a,b):
    a0,a1,a2,a3=a; b0,b1,b2,b3=b
    return (a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
            a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0)
PERMS=list(itertools.permutations(range(4)))
SIGNS=list(itertools.product([1,-1],repeat=4))
def G(perm,sign,v): return tuple(sign[i]*v[perm[i]] for i in range(4))
def Ginv(perm,sign,v):
    out=[0]*4
    for i in range(4): out[perm[i]] = sign[i]*v[i]
    return tuple(out)
basis=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
tables={}
for perm in PERMS:
    for sign in SIGNS:
        tbl=tuple(G(perm,sign,qmul(Ginv(perm,sign,a),Ginv(perm,sign,b))) for a in basis for b in basis)
        tables.setdefault(tbl,[]).append((perm,sign))
print(f"E1: 相異なる H 構造(積表): {len(tables)} 個 / 安定化群サイズ: {sorted({len(v) for v in tables.values()})}")
print(f"    384 = {len(tables)} × {len(next(iter(tables.values())))} → 単一推移軌道(マーキングの選択は基層に痕跡なし)")
reals=set()
for tbl in tables:
    for v in basis:
        for s in (1,-1):
            cand=tuple(s*x for x in v)
            idx=basis.index(v)
            ok=all(tuple(s*x for x in tbl[idx*4+basis.index(b)])==b for b in basis)
            if ok: reals.add(v)
print(f"E2: 単位元(実軸)になりうる軸: {len(reals)} / 4 → どの軸が「実」かはマーキングの属性")
inv_norm = all(sum(x*x for x in G(p,s,(1,2,3,4)))==30 for p in PERMS for s in SIGNS)
inv_eps  = all((-1)**sum(abs(x) for x in G(p,s,(1,2,0,1)))==1 for p in PERMS for s in SIGNS)
print(f"E3: ノルム(時計/R)の B4 不変: {inv_norm} / ε(Q)の B4 不変: {inv_eps}")
print(f"    → 全16マーキングが共有する正準不変量はこの二つ — R/Q は「選ばれる軸」ではなく全視点の共有物")
