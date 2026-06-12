#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺71 §3: (a) 単符号支持補題の s=11/13 拡張 (b) m=3 タイ則 (3,3,0,0)@shell25 (予言: 3セクター[12,6,6])
import itertools, time
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
def support_signs(va,vb,t):
    A=expcoeffs(tuple(va)); B=expcoeffs(tuple(vb))
    ax=[set() for _ in range(4)]
    for na,ca in A.items():
        nb=tuple(x-y for x,y in zip(t,na))
        if nb in B and abs(ca*B[nb])>1e-15:
            for j in range(4): ax[j].add((na[j]>0)-(na[j]<0))
    return ax
def flip_axis(v,j): return tuple(-v[i] if i==j else v[i] for i in range(4))
for s in (11,13):
    chans=sorted({(p[1],p[2][0],p[2][1]) for f in all_finals(s) for p in two_step_paths(s,f)})
    nc=0; ns=0; nf=0; npi=0; nz=0; bad=0
    tested=set()
    for (x,c,d) in chans:
        for vx in SH[x]:
            t=tuple(vx)
            for vc in SH[c]:
                for vd in SH[d]:
                    key=(vx,vc,vd)
                    if key in tested: continue
                    tested.add(key)
                    e2=tphase(vx,vc,vd)
                    if e2 is None: continue
                    nc+=1
                    if all(len(a)<=1 for a in support_signs(vc,vd,t)): ns+=1
                    for j in range(4):
                        if vc[j]==0: continue
                        nf+=1
                        e2f=tphase(vx,flip_axis(vc,j),vd)
                        if e2f is None: nz+=1
                        else:
                            r=e2f/e2
                            if abs(r-1j)<1e-9 or abs(r+1j)<1e-9: npi+=1
                            else: bad+=1
    print(f"補題拡張 s={s}: 寄与係数 {nc}, 単符号支持 {ns}/{nc}, 枝反転 {nf} → ±i {npi}, 零化 {nz}, 反例 {bad}")
# ---- m=3 タイ ----
for m in (15,17,19,21,23,25):
    if m not in SH: SH[m]=shell_cells(float(m),5)
S=25; fins=all_finals(S)
cells=[v for v in SH[S] if tuple(sorted(map(abs,v)))==(0,0,3,3)]
print(f"(3,3,0,0)@shell25 cells: {len(cells)} 予言: 3セクター[12,6,6] (n=2タイ則のm=3版)")
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
def runS(vp):
    z={f:0+0j for f in fins}
    for f in fins:
        for ((a,b),x,(c,d),d1) in two_step_paths(S,f):
            for ia,va in enumerate(SH[a]):
                for ib,vb in enumerate(SH[b]):
                    if a==b and not ia<ib: continue
                    e1=tphase(vp,va,vb)
                    if e1 is None: continue
                    for vx in ([va,vb] if a==b else [va if x==a else vb]):
                        z[f]+=e1*ch2(vx,c,d)
    return z
def cls2(v):
    nzc=[t for t in v if t!=0]
    npos=sum(1 for t in nzc if t>0)
    return f"{npos}+/{2-npos}-"
t0=time.time()
sigs=defaultdict(list)
for i,v in enumerate(cells):
    z=runS(v)
    sigs[tuple(round(abs(z[f])**2,1) for f in fins)].append(v)
print(f"m=3タイ 総時間 {time.time()-t0:.0f}s セクター数: {len(sigs)}")
for sg,mem in sorted(sigs.items(), key=lambda kv:-len(kv[1])):
    cl=defaultdict(int)
    for v in mem: cl[cls2(v)]+=1
    print(f"  サイズ {len(mem)}: 符号クラス {dict(cl)} 例 {mem[:2]}")
