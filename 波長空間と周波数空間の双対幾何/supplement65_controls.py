#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺65 補強: C2対照(ε違反候補は存在するが輸送がゼロにする) + C3拡大
import numpy as np, itertools, random
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('[1]')[0].split('print(')[0])
def sval(v): return sum((abs(t)+0.5)**2 for t in v)
def eps(v): return (-1)**sum(abs(t) for t in v)
# C2対照: s=13、混在シェル(7,9)を含む終状態で、e フィルタ前の全セル組合せの ε 違反数 vs 寄与経路の違反数
s=13; fins=all_finals(s)
cand_viol=0; cand_tot=0; contrib_viol=0; contrib_tot=0
for vp in SH[s][:3]:
    for f in fins:
        for ((a,b),x,(c,d),d1) in two_step_paths(s,f):
            for va in SH[a]:
                for vb in SH[b]:
                    kept = vb if x==a else va; vdec = va if x==a else vb
                    e1=tphase(vp,va,vb)
                    for vc in SH[c]:
                        for vd in SH[d]:
                            ok_eps = (eps(vp)==eps(kept)*eps(vc)*eps(vd))
                            cand_tot+=1
                            if not ok_eps: cand_viol+=1
                            if e1 is None: continue
                            if tphase(vdec,vc,vd) is None: continue
                            contrib_tot+=1
                            if not ok_eps: contrib_viol+=1
print(f"C2対照 (s=13, 親3例): 全セル候補 {cand_tot} 中 ε違反 {cand_viol} ({100*cand_viol/cand_tot:.1f}%)")
print(f"           寄与経路   {contrib_tot} 中 ε違反 {contrib_viol}")
print(f"→ {'輸送(管理連鎖の三角恒等式)がε違反候補を厳密にゼロ化している' if cand_viol>0 and contrib_viol==0 else '判定要再検'}")
# C3拡大: 全シェルから無作為に60終状態
Ng=8; g=np.arange(Ng)/Ng
def axwave(k):
    if k==0: return np.ones(Ng)
    return np.sqrt(2)*(np.cos(2*np.pi*abs(k)*g) if k>0 else np.sin(2*np.pi*abs(k)*g))
def cellwave(v):
    w=axwave(v[0]).reshape(-1,1,1,1)*axwave(v[1]).reshape(1,-1,1,1)
    return w*axwave(v[2]).reshape(1,1,-1,1)*axwave(v[3]).reshape(1,1,1,-1)
ALL=[v for r in (1,3,5,7,9) for v in SH[r]]
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
rec=0; sadd=0
for (k1,k2,k3) in cases:
    I=(1.0+cellwave(k1)+cellwave(k2)+cellwave(k3))**2
    scores={v: float(np.mean(I*cellwave(v))) for v in ALL}
    top=sorted(scores,key=lambda v:-abs(scores[v]))[:3]
    if set(top)=={k1,k2,k3}:
        rec+=1
        if abs(sum(sval(v) for v in top)-9.0)<1e-9: sadd+=1
print(f"C3拡大: 異なり終状態 {len(cases)} 件 → 構成完全復元 {rec}/{len(cases)}, s_read 加法 {sadd}/{len(cases)}")
