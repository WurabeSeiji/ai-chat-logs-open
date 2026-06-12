#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺69 §A-1: 3軸タイ非退化最小例 (2,2,2,0)@shell19 のセクター構造 (補遺68 §A-1 引き継ぎ)
import itertools, time
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
for m in (15,17,19):
    if m not in SH: SH[m]=shell_cells(float(m),5)
print("shell sizes:", {m:len(SH[m]) for m in sorted(SH)})
S=19; fins=all_finals(S)
print(f"finals({S}):", fins)
cells=[v for v in SH[S] if tuple(sorted(map(abs,v)))==(0,2,2,2)]
print(f"(2,2,2,0) cells: {len(cells)}")
CH2={}
def ch2(vx,c,d):
    key=(vx,c,d)
    if key in CH2: return CH2[key]
    tot=0+0j
    if c==d:
        sh=SH[c]
        for ic in range(len(sh)):
            for idd in range(ic+1,len(sh)):
                e2=tphase(vx,sh[ic],sh[idd])
                if e2 is not None: tot+=e2
    else:
        for vc in SH[c]:
            for vd in SH[d]:
                e2=tphase(vx,vc,vd)
                if e2 is not None: tot+=e2
    CH2[key]=tot; return tot
def run19(vp):
    z={f:0+0j for f in fins}
    for f in fins:
        for ((a,b),x,(c,d),d1) in two_step_paths(S,f):
            for ia,va in enumerate(SH[a]):
                for ib,vb in enumerate(SH[b]):
                    if a==b and not ia<ib: continue
                    e1=tphase(vp,va,vb)
                    if e1 is None: continue
                    vxs=[va,vb] if a==b else [va if x==a else vb]
                    for vx in vxs:
                        z[f]+=e1*ch2(vx,c,d)
    return z
def tie3_class(v):
    nz=[v[j] for j in range(4) if v[j]!=0]
    npos=sum(1 for t in nz if t>0)
    return abs(npos-(3-npos))   # 3=全同符号, 1=2-1混合
t0=time.time()
sigs=defaultdict(list)
for i,v in enumerate(cells):
    z=run19(v)
    sg=tuple(round(abs(z[f])**2,1) for f in fins)
    sigs[sg].append(v)
    if i%8==7: print(f"  {i+1}/{len(cells)} 完了 ({time.time()-t0:.0f}s)")
print(f"総時間 {time.time()-t0:.0f}s")
print(f"セクター数: {len(sigs)}")
allzero=all(all(x==0 for x in sg) for sg in sigs)
print(f"全零退化か: {allzero}")
for sg,mem in sorted(sigs.items(), key=lambda kv:-len(kv[1])):
    classes=defaultdict(int)
    for v in mem: classes[tie3_class(v)]+=1
    nzsig={f:x for f,x in zip(fins,sg) if x!=0}
    print(f"  サイズ {len(mem)}: 相対符号クラス {dict(classes)} 非零チャネル {nzsig if nzsig else '全零'}")
    print(f"    例: {mem[:3]}")
