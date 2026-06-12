#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Paper 13 検証付録: 4軸積波の全数復号 (補遺64 の 2軸→4軸昇格)
#  予言: 対角商 = 偶数枚符号反転群 Z2^3 (位数8) → 256配置/8 = 32クラス、読み出しは厳密にその商
#  パリティ混在(奇活性軸セル追加)で 256/256 完全復号
import numpy as np, itertools
from collections import defaultdict
N=8
g=np.arange(N)/N
ax=[g.reshape(-1,1,1,1),g.reshape(1,-1,1,1),g.reshape(1,1,-1,1),g.reshape(1,1,1,-1)]
def prod4(ds):
    w=1.0
    for j in range(4): w=w*np.sqrt(2)*np.cos(2*np.pi*(ax[j]-ds[j]/4.0))
    return w
def line4(I,fv):
    e=np.exp(-2j*np.pi*sum(f*x for f,x in zip(fv,ax)))
    return 2*np.mean(I*np.conj(e))
def qc(z):
    p=(np.angle(z)/(2*np.pi))%1.0
    return int(np.round(p*4))%4, abs(p*4-np.round(p*4))
SIGNS=[s for s in itertools.product([1,-1],repeat=4) if s[0]==1]  # (1,±1,±1,±1) 8本
# T1: 状態水準の商 = 偶数枚反転群
def key_state(ds):
    orb=set()
    for flips in itertools.product([0,1],repeat=4):
        if sum(flips)%2: continue
        orb.add(tuple((d+2*f)%4 for d,f in zip(ds,flips)))
    return frozenset(orb)
classes=defaultdict(list); sig2cfg=defaultdict(set)
ok_q=0
for ds in itertools.product(range(4),repeat=4):
    I=(1.0+prod4(ds))**2
    sig=[]
    for sv in SIGNS:
        r,e=qc(line4(I,sv)); sig.append(r)
        assert e<1e-9
    sig=tuple(sig)
    classes[key_state(ds)].append(ds)
    sig2cfg[sig].add(key_state(ds))
n_state=len(classes)
sig_classes=len(sig2cfg)
quot_ok = all(len(v)==1 for v in sig2cfg.values()) and sig_classes==n_state
print(f"T1/T2 4軸: 状態クラス(偶数枚反転商) {n_state} (予言32), 読み出しシグネチャ {sig_classes}, 商と一対一: {quot_ok}")
# T5: パリティ混在 — 1軸セル(軸1,2,3)を追加: 全ての非自明偶フリップがいずれかの奇活性セルと奇交差 → Z2^3 完全破れ
# 復号: d1,d2,d3 = 単独軸線の直読み、d4 = (Σd線 − d1−d2−d3) mod 4。交差項は使用線に衝突しない(確認済み設計)
ok=0
for ds in itertools.product(range(4),repeat=4):
    Psi=1.0+prod4(ds)
    for j in range(3):
        Psi=Psi+np.sqrt(2)*np.cos(2*np.pi*(ax[j]-ds[j]/4.0))
    I=Psi**2
    es=[]
    d=[0]*4
    for j,fv in enumerate([(1,0,0,0),(0,1,0,0),(0,0,1,0)]):
        d[j],e=qc(line4(I,fv)); es.append(e)
    r4,e=qc(line4(I,(1,1,1,1))); es.append(e)
    d[3]=(r4-d[0]-d[1]-d[2])%4
    if max(es)<1e-9 and tuple(d)==ds: ok+=1
print(f"T5 4軸パリティ混在(1軸セル×3): 全256配置 一意復号 {ok}/256")
