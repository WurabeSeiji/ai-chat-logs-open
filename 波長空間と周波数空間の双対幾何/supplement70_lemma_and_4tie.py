#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺70 検証:
#  L: 単軸枝反転補題 — (i) 寄与する輸送係数の交差和は各軸で単符号支持 (ii) 枝反転は厳密に±iを掛けるか0にする
#  T: 4軸タイ (2,2,2,2)@shell25 — 統一則の事前予言「5セクター [6,4,4,1,1]、2-2クラスのみゲージ退化」のテスト
import itertools, time
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])

# ---------- L: 補題の悉皆検証 (s=9 の全チャネル) ----------
def support_signs(va,vb,t):
    """cross_at の寄与項の na_j 符号集合を軸ごとに返す"""
    A=expcoeffs(tuple(va)); B=expcoeffs(tuple(vb))
    axsigns=[set() for _ in range(4)]
    for na,ca in A.items():
        nb=tuple(x-y for x,y in zip(t,na))
        if nb in B and abs(ca*B[nb])>1e-15:
            for j in range(4): axsigns[j].add((na[j]>0)-(na[j]<0))
    return axsigns
def flip_axis(v,j): return tuple(-v[i] if i==j else v[i] for i in range(4))
n_combo=0; n_single=0; n_flip=0; n_pm_i=0; n_zero=0; bad=0
chans=[]   # (vx, vc, vd) 第二段 + (vp, va, vb) 第一段 を同列に検査
for f in all_finals(9):
    for ((a,b),x,(c,d),d1) in two_step_paths(9,f):
        chans.append((x,c,d))
chans=sorted(set(chans))
vps=[v for v in SH[9] if tuple(sorted(map(abs,v)))==(0,0,1,2)][:4]
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
                n_combo+=1
                ss=support_signs(vc,vd,t)
                if all(len(s)<=1 for s in ss): n_single+=1
                for j in range(4):
                    if vc[j]==0: continue
                    n_flip+=1
                    e2f=tphase(vx,flip_axis(vc,j),vd)
                    if e2f is None: n_zero+=1
                    else:
                        r=e2f/e2
                        if abs(r-1j)<1e-9 or abs(r+1j)<1e-9: n_pm_i+=1
                        else: bad+=1
print(f"L(第二段, s=9 全チャネル): 寄与係数 {n_combo}, 単符号支持 {n_single}/{n_combo}")
print(f"  枝反転 {n_flip} 件 → ±i {n_pm_i}, 零化 {n_zero}, 反例 {bad}")
# 第一段も
n1=0; s1=0; nf=0; pi=0; z0=0; bad1=0
for vp in vps:
    t=tuple(vp)
    for (a,b) in sorted({p[0] for f in all_finals(9) for p in two_step_paths(9,f)}):
        for va in SH[a]:
            for vb in SH[b]:
                e1=tphase(vp,va,vb)
                if e1 is None: continue
                n1+=1
                ss=support_signs(va,vb,t)
                if all(len(s)<=1 for s in ss): s1+=1
                for j in range(4):
                    if va[j]==0: continue
                    nf+=1
                    e1f=tphase(vp,flip_axis(va,j),vb)
                    if e1f is None: z0+=1
                    else:
                        r=e1f/e1
                        if abs(r-1j)<1e-9 or abs(r+1j)<1e-9: pi+=1
                        else: bad1+=1
print(f"L(第一段, 代表親4): 寄与係数 {n1}, 単符号支持 {s1}/{n1}")
print(f"  枝反転 {nf} 件 → ±i {pi}, 零化 {z0}, 反例 {bad1}")

# ---------- T: 4軸タイ (2,2,2,2)@shell25 ----------
for m in (15,17,19,21,23,25):
    if m not in SH: SH[m]=shell_cells(float(m),5)
S=25; fins=all_finals(S)
cells=[v for v in SH[S] if tuple(sorted(map(abs,v)))==(2,2,2,2)]
print(f"shell25 cells: {len(SH[25])}, (2,2,2,2): {len(cells)}, finals: {len(fins)}")
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
                    vxs=[va,vb] if a==b else [va if x==a else vb]
                    for vx in vxs:
                        z[f]+=e1*ch2(vx,c,d)
    return z
def cls4(v):
    npos=sum(1 for t in v if t>0)
    return f"{npos}+/{4-npos}-"
t0=time.time()
sigs=defaultdict(list)
for i,v in enumerate(cells):
    z=runS(v)
    sg=tuple(round(abs(z[f])**2,1) for f in fins)
    sigs[sg].append(v)
    print(f"  cell {i+1}/{len(cells)} ({time.time()-t0:.0f}s)")
print(f"総時間 {time.time()-t0:.0f}s")
print(f"セクター数: {len(sigs)} (予言: 5 [6,4,4,1,1])")
for sg,mem in sorted(sigs.items(), key=lambda kv:-len(kv[1])):
    cl=defaultdict(int)
    for v in mem: cl[cls4(v)]+=1
    nz=sum(1 for x in sg if x!=0)
    print(f"  サイズ {len(mem)}: 符号クラス {dict(cl)} 非零チャネル数 {nz}  例 {mem[:2]}")
