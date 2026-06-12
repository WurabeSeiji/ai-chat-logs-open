#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺71 §1: 位置水準管理連鎖(線連続則)の再現と A1/B2 の構造
import itertools, numpy as np
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])   # two_step_paths, SH(タプル), all_finals
SHa={m:[np.array(v) for v in SH[m]] for m in (1,3,5,7,9)}
def line_ok(vp,va,vb):
    for x in (va+vb, va-vb):
        if np.array_equal(x,vp) or np.array_equal(-x,vp): return True
    return False
v9=np.array((2,1,0,0))
chains=[]
for final in [(5,3,1),(3,3,3)]:
    for ((a,b),x,(c,d),d1) in two_step_paths(9,final):
        for va in SHa[a]:
            for vb in SHa[b]:
                if a==b and tuple(va)>=tuple(vb): continue
                if not line_ok(v9,va,vb): continue
                vx = va if x==a else vb
                for vc in SHa[c]:
                    for vd in SHa[d]:
                        if c==d and tuple(vc)>=tuple(vd): continue
                        if line_ok(vx,vc,vd):
                            chains.append((final,((a,b),x,(c,d),d1),tuple(va),tuple(vb),tuple(vx),tuple(vc),tuple(vd)))
cnt=defaultdict(int); dcnt=defaultdict(lambda: defaultdict(int))
for ch in chains:
    cnt[ch[0]]+=1; dcnt[ch[0]][ch[1][3]]+=1
print(f"連鎖再現: 531終 {cnt[(5,3,1)]} (期待48), 333終 {cnt[(3,3,3)]} (期待20), 計 {len(chains)}")
print(f"A1セル水準: 531 δ1内訳 {dict(dcnt[(5,3,1)])}, 333 {dict(dcnt[(3,3,3)])}")
s531=[i for i in dcnt[(5,3,1)]]
A1_531=sum((1j)**d * n for d,n in dcnt[(5,3,1)].items())
A1_333=sum((1j)**d * n for d,n in dcnt[(3,3,3)].items())
print(f"A1セル水準和 Σi^δ1: 531={A1_531}, 333={A1_333}  (シェル水準6/6・1/1は補遺70で証明済)")
# B2: sin枝総数 mod 4 (sin枝 = 負成分の数。集計対象の変種を試す)
def nsin(*vs): return sum(sum(1 for t in v if t<0) for v in vs)
for name,sel in (("全5セル",lambda ch:(ch[2],ch[3],ch[5],ch[6],(2,1,0,0))),
                 ("娘4セル(va,vb,vc,vd)",lambda ch:(ch[2],ch[3],ch[5],ch[6])),
                 ("最終3セル(spec,vc,vd)",None)):
    hist=defaultdict(int)
    for ch in chains:
        if ch[0]!=(5,3,1): continue
        if sel is None:
            (a,b),x,(c,d),d1=ch[1]
            spec = ch[3] if x==a else ch[2]
            vs=(spec,ch[5],ch[6])
        else: vs=sel(ch)
        hist[nsin(*vs)%4]+=1
    print(f"B2 sin枝総数 mod4 ({name}): {dict(sorted(hist.items()))}")

print()
print("==== B2 の局所構造: シェル経路ごと・第一段ごとの sin枝 mod4 ====")
from collections import defaultdict as dd
per_path=dd(lambda: dd(int))
per_stage=dd(lambda: dd(int))
for ch in chains:
    if ch[0]!=(5,3,1): continue
    (ab,x,cdp,d1)=ch[1]
    n4=nsin(ch[2],ch[3],ch[5],ch[6])%4
    per_path[(ab,x,cdp,d1)][n4]+=1
    n2=nsin(ch[5],ch[6])%4   # 第二段のみ
    per_stage[(ab,x,cdp,d1,'2nd')][n2]+=1
for k in sorted(per_path):
    print(f"  path {k}: 全体mod4 {dict(sorted(per_path[k].items()))}  第二段のみ {dict(sorted(per_stage[(k[0],k[1],k[2],k[3],'2nd')].items()))}")
# 第二段の4解 (補遺38 の Z2×Z2 トーサー) の sin枝が {q,q+1,q+2,q+3} を走るか?
print()
print("==== 第二段4解の sin枝: トーサー軌道内で mod4 完全代表か ====")
grp=dd(list)
for ch in chains:
    if ch[0]!=(5,3,1): continue
    key=(ch[1],ch[2],ch[3])   # path + 第一段セル固定
    grp[key].append((ch[5],ch[6]))
allgood=True; nz_groups=0
for key,sols in grp.items():
    if len(sols)<2: continue
    nz_groups+=1
    s4=sorted(nsin(vc,vd)%4 for (vc,vd) in sols)
    if s4!=[0,1,2,3]: allgood=False; print(f"  非代表: {key[0]} 解{len(sols)} mod4={s4}")
print(f"  4解グループ数 {nz_groups}, 全グループが mod4 完全代表 {{0,1,2,3}}: {allgood}")

print()
print("==== トーサー残差補題: 4解軌道(第一段・第二段とも)の sin残差は完全代表系 {0,1,2,3} か ====")
# 第一段グループ: path + 第二段セル固定 → 第一段解の集合
g1=dd(list)
for ch in chains:
    if ch[0]!=(5,3,1): continue
    g1[(ch[1],ch[5],ch[6])].append((ch[2],ch[3]))
n4=0; bad4=0; n2=0; bad2=0
for key,sols in g1.items():
    rs=sorted(nsin(*s)%4 for s in sols)
    if len(sols)==4:
        n4+=1
        if rs!=[0,1,2,3]: bad4+=1
    elif len(sols)==2:
        n2+=1
        if rs not in ([0,2],[1,3]): bad2+=1
print(f"第一段: 4解軌道 {n4} (非代表 {bad4}), 2解軌道 {n2} (非{{q,q+2}} {bad2})")
g2=dd(list)
for ch in chains:
    if ch[0]!=(5,3,1): continue
    g2[(ch[1],ch[2],ch[3])].append((ch[5],ch[6]))
n4=0; bad4=0; pair_types=dd(int)
for key,sols in g2.items():
    rs=sorted(nsin(*s)%4 for s in sols)
    if len(sols)==4:
        n4+=1
        if rs!=[0,1,2,3]: bad4+=1
    elif len(sols)==2:
        pair_types[tuple(rs)]+=1
print(f"第二段: 4解軌道 {n4} (非代表 {bad4}), 2解軌道タイプ {dict(pair_types)}")
print()
print("==== 333終(20連鎖)の同構造 ====")
hist=dd(int)
for ch in chains:
    if ch[0]!=(3,3,3): continue
    hist[nsin(ch[2],ch[3],ch[5],ch[6])%4]+=1
print(f"333 sin枝 mod4: {dict(sorted(hist.items()))}")
