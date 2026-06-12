#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺54 検証: 強制ビットの同定テスト (T1: 矢反転 / T2: 部分向き反転 / T3: 親平均)
# 実行結果 (2026-06-12):
#   基準 (o1=+,o2=+): (1088,40) クラスA
#   T1 全段共役(矢反転): (1088,40) 不変 → ビット≠矢
#   T2a 第1段のみ反転: (320,32) / T2b 第2段のみ反転: (320,32) / T2c 両段反転: (1088,40)
#     → クラス = リンク向きの相対パリティ (Z2 ホロノミー)
#   鏡像親 (-2,1,0,0): (320,32) → 親の支配枝は鎖の境界値として同じ Z2 に入る
#   T3 親等重・非干渉平均: (704,36) → 分岐比 0.95135, 双子 P(X)=0.97507
import itertools, numpy as np
SQ2=np.sqrt(2.0)
def expcoeffs(v):
    acc=[((),1.0+0j)]
    for i in range(4):
        m=v[i]
        if m==0: opts=[(0,1.0+0j)]
        elif m>0: opts=[(+m,1/SQ2+0j),(-m,1/SQ2+0j)]
        else:
            mm=-m; opts=[(+mm,-1j/SQ2),(-mm,+1j/SQ2)]
        acc=[(n+(f,),c*w) for n,c in acc for f,w in opts]
    d={}
    for n,c in acc: d[n]=d.get(n,0)+c
    return d
def cross_at(va,vb,t):
    A=expcoeffs(tuple(va)); B=expcoeffs(tuple(vb)); tot=0+0j
    for na,ca in A.items():
        nb=tuple(x-y for x,y in zip(t,na))
        if nb in B: tot+=ca*B[nb]
    return tot
def self_at(v,t): return expcoeffs(tuple(v)).get(tuple(t),0+0j)
def tphase(vp,va,vb,orient=+1):
    tgt=tuple(orient*x for x in vp)
    Cp=self_at(vp,tgt)
    if abs(Cp)<1e-12: return None
    Cc=cross_at(va,vb,tgt)
    if abs(Cc)<1e-12: return None
    return (Cc/abs(Cc))/(Cp/abs(Cp))
def two_step_paths(parent,final):
    paths=[]; fin=tuple(sorted(final))
    for a in range(1,parent+2,2):
        for b in range(a,parent+2,2):
            d1=a+b-parent
            if d1 not in (1,-1): continue
            for (x,spec) in ((a,b),(b,a)):
                for c in range(1,x+2,2):
                    dd=x-d1-c
                    if dd<1 or dd%2==0: continue
                    if tuple(sorted((spec,c,dd)))==fin: paths.append(((a,b),x,(c,dd),d1))
    return sorted(set(paths))
def shell_cells(m):
    out=[]; K=4
    for k in itertools.product(range(-K,K+1),repeat=4):
        if abs(sum((abs(t)+0.5)**2 for t in k)-m)<1e-9: out.append(np.array(k))
    return out
SH={m:shell_cells(float(m)) for m in (1,3,5,7,9)}
LBL={lbl:two_step_paths(9,fin) for lbl,fin in [('531',(5,3,1)),('333',(3,3,3))]}
def run(v9,o1=+1,o2=+1,conj_all=False):
    z={'531':0+0j,'333':0+0j}
    for lbl,lp in LBL.items():
        for ((a,b),x,(c,d),d1) in lp:
            for va in SH[a]:
                for vb in SH[b]:
                    if a==b and tuple(va)>=tuple(vb): continue
                    e1=tphase(v9,va,vb,o1)
                    if e1 is None: continue
                    vx=va if x==a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            if c==d and tuple(vc)>=tuple(vd): continue
                            e2=tphase(vx,vc,vd,o2)
                            if e2 is None: continue
                            eta=e1*e2
                            if conj_all: eta=np.conj(eta)
                            z[lbl]+=eta
    return z
if __name__=='__main__':
    v9=np.array((2,1,0,0))
    for name,kw in [("基準",{}),("T1 矢反転",{'conj_all':True}),("T2a o1=-1",{'o1':-1}),
                    ("T2b o2=-1",{'o2':-1}),("T2c 両反転",{'o1':-1,'o2':-1})]:
        z=run(v9,**kw)
        print(name, {l: round(abs(z[l])**2) for l in z})
    z=run(np.array((-2,1,0,0)))
    print("鏡像親", {l: round(abs(z[l])**2) for l in z})
