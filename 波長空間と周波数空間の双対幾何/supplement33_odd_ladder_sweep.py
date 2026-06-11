#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺33 追加2: 奇数倍音のみの帯域でのアイ開口 + 頑健性チェック (ppu 依存)
import numpy as np, itertools, math

OUT=[]
def log(s):
    print(s); OUT.append(s)

def centers4(R):
    K=int(math.floor(R))+1
    out=[]
    for k in itertools.product(range(-K,K+1),repeat=4):
        if sum((abs(t)+0.5)**2 for t in k)<=R*R+1e-9: out.append(k)
    return out

def run_case(s, ppu_target=8.0, odd_cuts=(1,3,5)):
    R=math.sqrt(s); P=2*R
    G=int(round(ppu_target*P))
    cs=centers4(R)
    K=max(abs(t) for k in cs for t in k)
    L=2*K+3
    occ=np.zeros((L,)*4,dtype=bool)
    for k in cs:
        occ[tuple(t+K+1 for t in k)]=True
    x=-R+(np.arange(G)+0.5)*(P/G)
    col=np.clip(np.rint(x).astype(int),-(K+1),K+1)+K+1
    o=occ[col[:,None,None,None],col[None,:,None,None],col[None,None,:,None],col[None,None,None,:]].astype(np.float64)
    F=np.fft.fftn(o)
    idx=(np.fft.fftfreq(G)*G).astype(int)
    slots=list(itertools.product(range(-(K+1),K+2),repeat=4))
    occset=set(cs)
    def near_idx(v):
        p=((v+R)%P)-R
        return int(np.argmin(np.abs(x-p)))
    axidx={k:near_idx(k) for k in range(-(K+1),K+2)}
    occ_pts=[tuple(axidx[t] for t in k) for k in cs]
    emp_pts=[tuple(axidx[t] for t in k) for k in slots if k not in occset]
    res={}
    for nc in odd_cuts:
        # 奇数倍音のみ: |n| 奇数かつ <=nc、または n=0 (DC)
        m1=((np.abs(idx)%2==1)&(np.abs(idx)<=nc))|(idx==0)
        M=m1[:,None,None,None]&m1[None,:,None,None]&m1[None,None,:,None]&m1[None,None,None,:]
        ff=np.fft.ifftn(np.where(M,F,0)).real
        vo=np.array([ff[i] for i in occ_pts]); ve=np.array([ff[i] for i in emp_pts])
        res[nc]=(vo.min(),ve.max(),vo.min()-ve.max())
        del M,ff
    return dict(s=s,R=R,G=G,n=len(cs),res=res)

log("奇数倍音のみの帯域 (DC + |n| 奇数 <= nc) のアイ開口  [ppu=8]")
log("  s    R     セル数   nc=1 eye      nc=3 eye      nc=5 eye")
for s in [1,3,5,7,9,11,13]:
    r=run_case(s)
    e=[r['res'][nc][2] for nc in (1,3,5)]
    mk=lambda v:"O" if v>0 else "x"
    log(f"  {s:3d} {r['R']:.3f} {r['n']:6d}  {e[0]:+.4f} {mk(e[0])}   {e[1]:+.4f} {mk(e[1])}   {e[2]:+.4f} {mk(e[2])}")

log("")
log("頑健性: nc=1 の eye を ppu=6/8/10 で比較")
log("  s    ppu=6        ppu=8        ppu=10")
for s in [3,5,7,9,11,13]:
    vals=[]
    for ppu in (6.0,8.0,10.0):
        r=run_case(s,ppu_target=ppu,odd_cuts=(1,))
        vals.append(r['res'][1][2])
    log(f"  {s:3d}  {vals[0]:+.4f}     {vals[1]:+.4f}     {vals[2]:+.4f}")

with open('/tmp/closure_sweep2_results.txt','w') as f: f.write("\n".join(OUT))
log("done")
