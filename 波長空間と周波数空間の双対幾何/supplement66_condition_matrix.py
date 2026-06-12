#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺66 検証: 射影・逆射影の条件行列 — 特殊条件依存性の排除
import numpy as np, random, itertools
N=2160; xs=np.arange(N)/N; MMAX=81
def comb(k,a,jit=None):
    w=np.zeros(N); m=1; i=0
    while m*k<=MMAX:
        amp=(4/np.pi)*((-1)**((m-1)//2))/m
        if jit is not None: amp*=(1.0+jit[i])
        w+=amp*np.cos(2*np.pi*m*k*(xs-a)); m+=2; i+=1
    return w
def line(I,f):
    c=2*np.mean(I*np.cos(2*np.pi*f*xs)); s=2*np.mean(I*np.sin(2*np.pi*f*xs))
    return complex(c,s)
def phase(z): return (np.angle(z)/(2*np.pi))%1.0
def qc(z):
    p=phase(z); return int(np.round(p*4))%4, abs(p*4-np.round(p*4))
random.seed(5); np.random.seed(5)

# --- A: 同種3断片(基本波1)、12スロットから3占有 → 厳密配置読み(仮説照合) ---
# 同種断片は同一の線を共有(線値=3位相子の和)。復号=C(12,3)=220仮説の予測記録と厳密照合(配置読み)
import itertools as _it
HYP=list(_it.combinations(range(12),3))
PRED={}
for hyp in HYP:
    Ih=(1.0+sum(comb(1,t/12.0) for t in hyp))**2
    PRED[hyp]=np.array([line(Ih,f) for f in (1,3,5,7)])
okA=0; ambA=0; TA=100
for _ in range(TA):
    slots=tuple(sorted(random.sample(range(12),3)))
    I=(1.0+sum(comb(1,t/12.0) for t in slots))**2
    sig=np.array([line(I,f) for f in (1,3,5,7)])
    hits=[h for h in HYP if np.max(np.abs(PRED[h]-sig))<1e-9]
    if len(hits)==1 and hits[0]==slots: okA+=1
    elif len(hits)>1: ambA+=1
print(f"A 同種3断片・厳密配置読み(12スロット, {TA}試行): {okA}/{TA} (多義 {ambA})")
# --- B: 異種3断片(基本波1,3,5=セル|k|=0,1,2)、ピーリング連鎖、全64配置 ---
okB=0
for d1,d3,d5 in itertools.product(range(4),repeat=3):
    a1,a3,a5=d1/4.0,d3/12.0,d5/20.0
    I=(1.0+comb(1,a1)+comb(3,a3)+comb(5,a5))**2
    r1,e1=qc(line(I,1)); A1=r1/4.0
    z3=line(I,3)-2*(4/np.pi)*(-1)/3*np.exp(2j*np.pi*3*A1)
    r3,e3=qc(z3)
    z5=line(I,5)-2*(4/np.pi)*(+1)/5*np.exp(2j*np.pi*5*A1)
    r5,e5=qc(z5)
    if (r1,r3,r5)==(d1,d3,d5) and max(e1,e3,e5)<1e-6: okB+=1
print(f"B 異種3断片(1/3/5)・ピーリング、全64配置: {okB}/64")

# --- C: 連続位置3断片、復元精度(丸めなしの位相直読み) ---
okC=0; TC=200; werr=0
for _ in range(TC):
    a1=random.random(); a3=random.random()/3; a5=random.random()/5
    I=(1.0+comb(1,a1)+comb(3,a3)+comb(5,a5))**2
    h1=phase(line(I,1))
    z3=line(I,3)-2*(4/np.pi)*(-1)/3*np.exp(2j*np.pi*3*h1)
    h3=phase(z3)/3
    z5=line(I,5)-2*(4/np.pi)*(+1)/5*np.exp(2j*np.pi*5*h1)
    h5=phase(z5)/5
    err=max(abs(h1-a1),abs(h3-a3),abs(h5-a5))
    werr=max(werr,err)
    if err<1e-9: okC+=1
print(f"C 連続位置3断片({TC}試行): 位相直読み厳密復元(誤差<1e-9) {okC}/{TC}, 最大誤差 {werr:.2e}")

# --- D: B + 高調波振幅±20%ジッタ + 付加ガウスノイズ(σ=0.5)、全64×3draw ---
okD=0; TD=0
for d1,d3,d5 in itertools.product(range(4),repeat=3):
    for _ in range(3):
        TD+=1
        j1=np.random.uniform(-0.2,0.2,64); j3=np.random.uniform(-0.2,0.2,64); j5=np.random.uniform(-0.2,0.2,64)
        Psi=1.0+comb(1,d1/4.0,j1)+comb(3,d3/12.0,j3)+comb(5,d5/20.0,j5)
        I=Psi**2+np.random.normal(0,0.5,N)
        r1,_=qc(line(I,1)); A1=r1/4.0
        z3=line(I,3)-2*(4/np.pi)*(-1)/3*np.exp(2j*np.pi*3*A1)
        z5=line(I,5)-2*(4/np.pi)*(+1)/5*np.exp(2j*np.pi*5*A1)
        r3,_=qc(z3); r5,_=qc(z5)
        if (r1,r3,r5)==(d1,d3,d5): okD+=1
print(f"D 摂動下(±20%ジッタ+σ0.5ノイズ)、64配置×3: {okD}/{TD}")
