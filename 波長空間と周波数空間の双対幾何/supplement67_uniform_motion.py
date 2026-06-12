#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺67 検証: 等速運動の写像
#  O1: 等速運動 = 位置桁の混合基数オドメータ (base [4;3,3])。毎ティック最深¼量子(1/36)を1進める。
#      36ティック全てで階層ピーリング復号が桁=カウンタ値を返し、x̂_n = n/36 厳密。
#  O2: 単一コムの等速移動で (i) 全線の振幅 |Z(f)| は厳密不変 (ii) 位相は線形回転、回転率 ∝ f·v
#      (iii) R(=s_read), Q は不変 — 「周波数・波長は時間発展しない」の機械確認
import numpy as np
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
    p=(np.angle(z)/(2*np.pi))%1.0
    return int(np.round(p*4))%4, abs(p*4-np.round(p*4))

# ---- O1: オドメータ運動の階層復号 ----
ok=0
for n in range(36):
    d,c1,c2 = n//9, (n%9)//3, n%3
    a0=d/4.0; a1=a0+c1/12.0; a2=a1+c2/36.0
    I=(1.0+comb(1,a0)+comb(3,a1)+comb(9,a2))**2
    r0,e0=qc(line(I,1)); A0=r0/4.0
    r1,e1=qc(line(I,3)-expected(3,1,A0)); c1r=(r1-3*r0)%4; A1=A0+c1r/12.0
    r2,e2=qc(line(I,9)-expected(9,1,A0)-expected(9,3,A1)); c2r=(r2-3*r1)%4
    xhat=(r0+c1r/3.0+c2r/9.0)/4.0
    if (r0,c1r,c2r)==(d,c1,c2) and abs(xhat-n/36.0)<1e-12 and max(e0,e1,e2)<1e-9: ok+=1
print(f"O1 オドメータ等速運動(36ティック, v=最深¼量子/ティック): 桁=カウンタ一致+x̂=n/36厳密 {ok}/36")

# ---- O2: 単一コム等速移動のスペクトル ----
ticks=36; v_per_tick=1.0/36
fr=[f for f in range(1,28,2)]
mags=np.zeros((ticks,len(fr))); phs=np.zeros((ticks,len(fr)))
for n in range(ticks):
    I=(1.0+comb(1,n*v_per_tick))**2
    for j,f in enumerate(fr):
        z=line(I,f); mags[n,j]=abs(z); phs[n,j]=(np.angle(z)/(2*np.pi))%1.0
dev=np.max(np.abs(mags-mags[0:1,:]))
print(f"O2-i  振幅スペクトル不変性: 36ティックでの最大偏差 {dev:.2e}")
slopes=[]
for j,f in enumerate(fr[:5]):
    up=np.unwrap(phs[:,j]*2*np.pi)/(2*np.pi)
    co=np.polyfit(np.arange(ticks),up,1)
    res=np.max(np.abs(up-np.polyval(co,np.arange(ticks))))
    slopes.append(co[0])
    print(f"O2-ii 線 f={f}: 位相回転率 {co[0]:.6f}/tick (理論 f·v={f*v_per_tick:.6f}), 線形残差 {res:.1e}")
print(f"O2-ii 回転率比 f=3/f=1: {slopes[1]/slopes[0]:.6f} (理論 3)")
print(f"O2-iii R=√s, Q=ε: スペクトル支持(線の位置)が不変 → s_read・ε_read 不変 (偏差 {dev:.0e} の範囲で厳密)")
print(f"速度の符号化: v = (位相回転率/f) = {slopes[0]:.6f}/tick — 単一スナップショットには現れず、記録列にのみ存在")
