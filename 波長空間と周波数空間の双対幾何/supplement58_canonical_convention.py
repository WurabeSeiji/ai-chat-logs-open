#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺58 検証: 正準規約の特定 — 配置読み(無順序)+排他(対角除外)+事象多重度(両崩壊分岐の和)
# 結果 (2026-06-12, ローカル実行):
#   s=9 (2,1,0,0): 531=1088 333=40 P(X)=0.98195 (補遺51 値の復元)
#   s=13 (1,1,1,2): セクター 2 [32,32] (補遺57 の人工5を解消)
#   s=13 (0,0,2,2): セクター 3 [12,6,6] (真の破れは保持)
#   s=9  (0,0,1,2): セクター 2 [24,24]
# 基盤関数は補遺57 スクリプトと同一 (exec で再利用)
import itertools, numpy as np
from collections import defaultdict
src = open('supplement57_convention_dependence.py', encoding='utf-8').read()
exec(src.split('[1]')[0].split('print(')[0])  # expcoeffs/cross_at/self_at/tphase/two_step_paths/all_finals/SH

def run_canon(vp, parent, finals, o1=+1, o2=+1):
    """正準規約: 第一段=無順序+対角除外(排他)、a==b なら崩壊分岐を両方和、第二段=無順序+対角除外"""
    z={f:0+0j for f in finals}
    for f in finals:
        for ((a,b),x,(c,d),d1) in two_step_paths(parent,f):
            for ia,va in enumerate(SH[a]):
                for ib,vb in enumerate(SH[b]):
                    if a==b and not (ia<ib): continue
                    e1=tphase(vp,va,vb,o1)
                    if e1 is None: continue
                    vxs = [va,vb] if a==b else [va if x==a else vb]
                    for vx in vxs:
                        for ic,vc in enumerate(SH[c]):
                            for idd,vd in enumerate(SH[d]):
                                if c==d and not (ic<idd): continue
                                e2=tphase(vx,vc,vd,o2)
                                if e2 is None: continue
                                z[f]+=e1*e2
    return z

if __name__=='__main__':
    fins9=all_finals(9)
    z=run_canon((2,1,0,0),9,fins9)
    W531=abs(z[(1,3,5)])**2; W333=abs(z[(3,3,3)])**2
    print(f"s=9 (2,1,0,0): 531={W531:.0f} 333={W333:.0f} P(X)={2*W531/(2*W531+W333):.5f}")
    for s,t in [(13,(1,1,1,2)),(13,(0,0,2,2)),(9,(0,0,1,2))]:
        fins=all_finals(s)
        cells=[v for v in SH[s] if tuple(sorted(map(abs,v)))==t]
        sg=defaultdict(int)
        for v in cells:
            zz=run_canon(v,s,fins)
            sg[tuple(round(abs(zz[f])**2,1) for f in fins)]+=1
        print(f"s={s} 型{t} ({len(cells)}): セクター {len(sg)} {sorted(sg.values(),reverse=True)}")
