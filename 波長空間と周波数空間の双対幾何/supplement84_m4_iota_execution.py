#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺84: 論文16査読#5の実行 — (1) m=4@s21 パリティ仮説判定 (2) ι構成(隠れvxの軸0枝反転)の検証
import itertools, time
from collections import defaultdict
src=open('supplement57_convention_dependence.py',encoding='utf-8').read()
exec(src.split('# 現行(字句順序)規約')[0])
for m in (15,17,19,21):
    if m not in SH: SH[m]=shell_cells(float(m),5)

def paths_full(vp,s,final):
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
                            out.append((va,vb,vx,kept,vc,vd,e1*e2))
    return out

# ---- (2) ι 検証 (m=2,3: 軽い方を先に) ----
print("==== ι構成の検証 (隠れvxの軸0枝反転) ====")
def flip0(v): return (-v[0],)+tuple(v[1:])
for m,s,fin in ((2,7,(1,3,3)),(3,13,(3,5,5))):
    vA=tuple([m]+[0]*3)
    P=paths_full(vA,s,fin)
    idx={}
    for (va,vb,vx,kept,vc,vd,eta) in P:
        idx[(va,vb,vx,vc,vd)]=eta
    n_dom=0; n_in=0; n_flip=0; mult=defaultdict(int); n_fix=0
    for (va,vb,vx,kept,vc,vd,eta) in P:
        if vx[0]==0: continue
        n_dom+=1
        vx2=flip0(vx)
        va2 = vx2 if va==vx else va
        vb2 = vx2 if vb==vx else vb
        key=(va2,vb2,vx2,vc,vd)
        if key not in idx: continue
        n_in+=1
        s0a=sum(1 for w in (va,vb) if w[0]<0); s0b=sum(1 for w in (va2,vb2) if w[0]<0)
        if (s0a-s0b)%2==1: n_flip+=1
        r=idx[key]/eta
        mult[(round(r.real,6),round(r.imag,6))]+=1
    nz_dom0=sum(1 for (va,vb,vx,kept,vc,vd,eta) in P if vx[0]==0)
    print(f"m={m}: 経路{len(P)} 定義域(vx0≠0) {n_dom} (域外 {nz_dom0}), 像が集合内 {n_in}/{n_dom}, σパリティ反転 {n_flip}/{n_in}")
    print(f"   η乗数の分布: {dict(mult)}")

# ---- (1) m=4@s21 ----
print()
print("==== m=4 (s=21) パリティ仮説判定: 予言=偶mなので相殺せず (W_B≠0) ====")
t0=time.time()
fins=all_finals(21)
res={}
for vp,nm in (((4,0,0,0),'A'),((-4,0,0,0),'B')):
    zch=defaultdict(complex); zcf=defaultdict(lambda: defaultdict(complex))
    Sp=0+0j; Sm=0+0j
    for f in fins:
        for (va,vb,vx,kept,vc,vd,eta) in paths_full(vp,21,f):
            zch[f]+=eta
            K=frozenset([kept,vc,vd])
            if len(K)==3: zcf[f][K]+=eta
            if nm=='A':
                s0=sum(1 for w in (va,vb) if w[0]<0)
                if s0%2==0: Sp+=eta
                else: Sm+=eta
    W={f:abs(z)**2 for f,z in zch.items() if abs(z)>1e-9}
    Wp={f:sum(abs(z)**2 for z in zk.values()) for f,zk in zcf.items()}
    Wp={f:w for f,w in Wp.items() if w>1e-9}
    res[nm]=(W,Wp)
    print(f"  {nm}: W={ {f:round(w,1) for f,w in W.items()} }")
    print(f"     W'={ {f:round(w,1) for f,w in Wp.items()} } ({time.time()-t0:.0f}s)")
    if nm=='A': print(f"     符号付き和 S+−S− = {Sp-Sm:.3f} (0なら相殺=仮説反例, 非0なら偶m非相殺=仮説支持)")
