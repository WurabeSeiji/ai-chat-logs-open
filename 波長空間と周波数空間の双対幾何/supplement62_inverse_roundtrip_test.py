#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺62 追補検証: 逆写像の合成往復
#   連続 x → [逆写像: 貪欲丸め→桁] → [状態の実構成: 実コム3階層] → [順写像: ピーリング復号→桁→座標 x̂]
#   判定: 桁の完全一致 かつ |x−x̂| ≤ (1/4)·3^{-L}·(3/2) (深さL=2, 上界1/24)
import numpy as np, random
N=2160; xs=np.arange(N)/N; MMAX=81
def comb(k,a):
    w=np.zeros(N); m=1
    while m*k<=MMAX:
        w+=(4/np.pi)*((-1)**((m-1)//2))*np.cos(2*np.pi*m*k*(xs-a))/m; m+=2
    return w
def line(I,f):
    c=2*np.mean(I*np.cos(2*np.pi*f*xs)); s=2*np.mean(I*np.sin(2*np.pi*f*xs))
    return complex(c,s)
def expected(f,kl,a):
    if f%kl: return 0j
    m=f//kl
    if m%2==0 or m*kl>MMAX: return 0j
    return 2*(4/np.pi)*((-1)**((m-1)//2))/m*np.exp(2j*np.pi*f*a)
def qc(z):
    phi=(np.angle(z)/(2*np.pi))%1.0
    return int(np.round(phi*4))%4, abs(phi*4-np.round(phi*4))
def inverse_map(x):
    """逆写像 (補遺62 D4): 貪欲丸め → 桁 (d, c1, c2)"""
    d=int(np.floor(4*x))%4; rem=4*x-d
    c1=int(np.floor(3*rem)); c1=min(c1,2); rem=3*rem-c1
    c2=int(np.floor(3*rem)); c2=min(c2,2)
    return d,c1,c2
def construct(d,c1,c2):
    """桁 → 状態の実構成 (実コム3階層)"""
    a0=d/4.0; a1=a0+c1/12.0; a2=a1+c2/36.0
    return (1.0+comb(1,a0)+comb(3,a1)+comb(9,a2))**2, (a0,a1,a2)
def forward_map(I):
    """順写像 (補遺63 ピーリング復号): 記録 → 桁 → 座標"""
    r0,e0=qc(line(I,1)); A0=r0/4.0
    r1,e1=qc(line(I,3)-expected(3,1,A0)); c1=(r1-3*r0)%4; A1=A0+c1/12.0
    r2,e2=qc(line(I,9)-expected(9,1,A0)-expected(9,3,A1)); c2=(r2-3*r1)%4
    assert max(e0,e1,e2)<1e-9
    xhat=(r0 + c1/3.0 + c2/9.0)/4.0
    return (r0,c1,c2), xhat
random.seed(7)
BOUND=(1.0/4)*(1.0/9)*1.5
ok=0; T=500; werr=0
for _ in range(T):
    x=random.random()
    dg = inverse_map(x)
    I,_ = construct(*dg)
    dg2, xhat = forward_map(I)
    err=abs(x-xhat)
    werr=max(werr,err)
    if dg2==dg and err<=BOUND: ok+=1
    else: print(f"FAIL x={x:.6f} 桁{dg}→{dg2} 誤差{err:.5f}")
print(f"合成往復 (連続x→逆写像→状態実構成→順写像→x̂): {ok}/{T} PASS")
print(f"最大誤差 {werr:.5f} ≤ 上界 {BOUND:.5f} : {'PASS' if werr<=BOUND else 'FAIL'}")
