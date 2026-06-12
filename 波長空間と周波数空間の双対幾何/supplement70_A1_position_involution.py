#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺70 §2: 位置水準 A1 対合 — ベクトル一致管理連鎖(68本)上の「借入先行↔貸出先行」対合の構成
import itertools
from collections import defaultdict
def sval(v): return sum((abs(t)+0.5)**2 for t in v)
K=6
def cells_upto(K):
    return None
vp=(2,1,0,0)
RNG=range(-3,4)
# 管理連鎖: vp = v_kept + vc + vd (ベクトル一致), v_x = vc+vd
# 終状態シェル {s(kept), s(c), s(d)} = {5,3,1} or {3,3,3}, δ1 = s(kept)+s(vx) − 9 ∈ {±1}
def chains(vp, final):
    out=[]
    fin=tuple(sorted(final))
    for vk in itertools.product(RNG,repeat=4):
        sk=sval(vk)
        if sk not in final: continue
        rem=tuple(p-k for p,k in zip(vp,vk))
        for vc in itertools.product(RNG,repeat=4):
            sc=sval(vc)
            if sc not in final: continue
            vd=tuple(r-c for r,c in zip(rem,vc))
            sd=sval(vd)
            if tuple(sorted((sk,sc,sd)))!=fin: continue
            vx=tuple(c+d for c,d in zip(vc,vd))
            d1=sk+sval(vx)-9
            if d1 not in (1,-1): continue
            out.append((vk,vc,vd,int(sval(vx)),int(d1)))
    return out
for final,nexp in (((5,3,1),48),((3,3,3),20)):
    C=chains(vp,final)
    # (vk,{vc,vd}) の無順序で数える(管理連鎖の本数勘定)
    uno={}
    for (vk,vc,vd,x,d1) in C:
        key=(vk,frozenset([vc,vd]))
        uno[key]=(vk,vc,vd,x,d1)
    plus=[t for t in uno.values() if t[4]==1]; minus=[t for t in uno.values() if t[4]==-1]
    print(f"final {final}: 連鎖 {len(uno)} (期待 {nexp}), δ1=+1: {len(plus)}, −1: {len(minus)} → A1相殺の前提 {'OK' if len(plus)==len(minus) else 'NG'}")
    # 対合探索: 再括弧化 — kept と娘の役割交換 (vp = vk+vc+vd は固定、束ね直し)
    # σ_rb: (vk, vc, vd) → (vc, vk, vd) / (vd, vc, vk) のうち δ1 が反転し連鎖集合内に残るもの
    idx={(t[0],frozenset([t[1],t[2]])):t for t in uno.values()}
    def try_rebracket(t):
        vk,vc,vd,x,d1=t
        cands=[]
        for newk,n1,n2 in ((vc,vk,vd),(vd,vc,vk)):
            key=(newk,frozenset([n1,n2]))
            if key in idx and idx[key][4]==-d1:
                cands.append(key)
        return cands
    n1=0; n2=0; n0=0
    for t in uno.values():
        c=try_rebracket(t)
        if len(c)==0: n0+=1
        elif len(c)==1: n1+=1
        else: n2+=1
    print(f"  再括弧化候補: 一意 {n1}, 複数 {n2}, なし {n0}")
    # 一意なら自然な対合になっているか (σ²=id)
    if n0==0:
        ok=0
        for t in uno.values():
            c=try_rebracket(t)
            k2=c[0]   # 一意 or 第一候補
            t2=idx[k2]
            c2=try_rebracket(t2)
            if c2 and idx[c2[0]]==t: ok+=1
        print(f"  σ²=id: {ok}/{len(uno)}")
