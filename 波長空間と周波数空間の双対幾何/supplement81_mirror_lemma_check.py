#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺81: 鏡像反転補題の機械検証 (補遺80 §A 引き継ぎ)
#  L80: 共通経路で η_B/η_A = −i·(−1)^{σ0}, σ0 = 第一段娘対の軸0 sin枝数
#  C1: σ0 の配置定数性 (m=3 で定数 ⟹ χ=±i / m=2 で混在 ⟹ W'破れ)
#  C2: 第六相殺の還元 W_B=0 ⟺ S+ = S− (σパリティ符号付き和)
import itertools
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
def paths_with_e1(vp,s,final):
    out=[]
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
                            out.append(((va,vb,vx,vc,vd),e1*e2,frozenset([kept,vc,vd])))
    return out
for m,s,fin in ((2,7,(1,3,3)),(3,13,(3,5,5))):
    vA=tuple([m]+[0]*3); vB=tuple([-m]+[0]*3)
    PA=paths_with_e1(vA,s,fin); PB=paths_with_e1(vB,s,fin)
    dA={t:(e,K) for t,e,K in PA}; dB={t:(e,K) for t,e,K in PB}
    common=set(dA)&set(dB)
    onlyA=len(dA)-len(common); onlyB=len(dB)-len(common)
    ok=0; bad=0
    Sp=0+0j; Sm=0+0j
    cfg_par=defaultdict(set)
    for t in common:
        va,vb,vx,vc,vd=t
        s0=sum(1 for w in (va,vb) if w[0]<0)
        r=dB[t][0]/dA[t][0]
        pred=(-1j)*((-1)**s0)
        if abs(r-pred)<1e-9: ok+=1
        else: bad+=1
        if s0%2==0: Sp+=dA[t][0]
        else: Sm+=dA[t][0]
        K=dA[t][1]
        if len(K)==3: cfg_par[K].add(s0%2)
    mixed=sum(1 for v in cfg_par.values() if len(v)>1)
    zB_pred=(-1j)*(Sp-Sm)
    zB_actual=sum(e for _,e,_ in PB)
    zA=sum(e for _,e,_ in PA)
    print(f"m={m} (s={s}): 経路 A={len(dA)} B={len(dB)} 共通={len(common)} (A単独{onlyA}/B単独{onlyB})")
    print(f"  L80 補題 η_B/η_A=−i(−1)^σ0: {ok}/{len(common)} (反例 {bad})")
    print(f"  C1 σパリティ配置定数性: 混在配置 {mixed}/{len(cfg_par)} → {'定数(χ=±i説明)' if mixed==0 else '混在(W´破れの所在)'}")
    print(f"  C2 還元: S+−S−={Sp-Sm:.3f} → z_B予測=−i(S+−S−)={zB_pred:.3f} vs 実測={zB_actual:.3f} {'一致' if abs(zB_pred-zB_actual)<1e-6 else '不一致'}")
    print(f"     z_A=S++S−={Sp+Sm:.3f} 実測={zA:.3f}")

print()
print("==== C1改: 配置ごとのσパリティ部分和 (混在でもχ=±iの真機構) ====")
for m,s,fin in ((2,7,(1,3,3)),(3,13,(3,5,5))):
    vA=tuple([m]+[0]*3)
    PA=paths_with_e1(vA,s,fin)
    SK=defaultdict(lambda:[0+0j,0+0j])
    for (va,vb,vx,vc,vd),e,K in PA:
        if len(K)<3: continue
        s0=sum(1 for w in (va,vb) if w[0]<0)
        SK[K][s0%2]+=e
    onezero=0; both=0
    for K,(sp,sm) in SK.items():
        if min(abs(sp),abs(sm))<1e-9 and max(abs(sp),abs(sm))>1e-9: onezero+=1
        elif abs(sp)>1e-9 and abs(sm)>1e-9: both+=1
    print(f"m={m}: 非零配置 {len(SK)} → 片側部分和=0 が {onezero}, 両側非零 {both}")
