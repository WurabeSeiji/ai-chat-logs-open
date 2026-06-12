#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺63 追補: 実コム(全奇数高調波=矩形波)での周波数衝突線のピーリング復号
# 衝突: 子コム基本波(周波数3)は親コムの第3高調波と同一線。孫(9)は親m=9・子m=3と同一線。
# 鍵1: 奇数×奇数の干渉項(コム×コム)は全て偶数周波数に落ちる → 奇数の読み出し線は転写項のみ。
# 鍵2: 浅層の寄与は復号済みの位置から厳密に既知 → 線から剥がす(ピーリング) → 残りが当該レベルの桁。
import numpy as np
N=2160; xs=np.arange(N)/N
MMAX=81   # 高調波打ち切り
def comb(k, a, Mmax=MMAX):
    """矩形波コム: sign波 = (4/π) Σ_{m odd} (-1)^{(m-1)/2} cos(2π m k (x-a))/m"""
    w=np.zeros(N)
    m=1
    while m*k<=Mmax:
        w += (4/np.pi)*((-1)**((m-1)//2))*np.cos(2*np.pi*m*k*(xs-a))/m
        m+=2
    return w
def line(I,f):
    c=2*np.mean(I*np.cos(2*np.pi*f*xs)); s=2*np.mean(I*np.sin(2*np.pi*f*xs))
    return complex(c,s)    # 位相規約: cos(2πf(x-a)) → amp·e^{+2πifa} (位相 +fa を読む)
def expected(f, klevel, a):
    """レベル klevel=3^ℓ のコム(位置 a)が周波数 f に置く転写寄与 (I の交差項 2·1·comb)"""
    if f % klevel != 0: return 0j
    m = f // klevel
    if m % 2 == 0 or m*klevel > MMAX: return 0j
    amp = 2*(4/np.pi)*((-1)**((m-1)//2))/m
    return amp*np.exp(+2j*np.pi*f*a)
def qclass(z):
    phi=(np.angle(z)/(2*np.pi))%1.0
    r=int(np.round(phi*4))%4
    return r, abs(phi*4-np.round(phi*4))
ok=0; tot=0; worst=0
for d in range(4):
    for c1 in range(3):
        for c2 in range(3):
            tot+=1
            a0=d/4.0; a1=a0+c1/12.0; a2=a1+c2/36.0
            Psi = 1.0 + comb(1,a0) + comb(3,a1) + comb(9,a2)
            I = Psi**2
            # 周波数1: 親のみ
            r0,e0 = qclass(line(I,1))
            d_r = r0
            A0 = d_r/4.0
            # 周波数3: 親m=3 を剥がす → 子m=1
            z3 = line(I,3) - expected(3,1,A0)
            r1,e1 = qclass(z3)
            c1_r = (r1 - 3*d_r) % 4
            A1 = A0 + c1_r/12.0
            # 周波数9: 親m=9・子m=3 を剥がす → 孫m=1
            z9 = line(I,9) - expected(9,1,A0) - expected(9,3,A1)
            r2,e2 = qclass(z9)
            c2_r = (r2 - 3*r1) % 4
            worst=max(worst,e0,e1,e2)
            if (d_r,c1_r,c2_r)==(d,c1,c2): ok+=1
            elif tot<6: print(f"  FAIL d={d},c1={c1},c2={c2} → {d_r},{c1_r},{c2_r} (誤差 {e0:.3f},{e1:.3f},{e2:.3f})")
print(f"実コム(高調波≤{MMAX})・3階層・全36配置: 一意復号 {ok}/{tot}, 最大¼量子化誤差 {worst:.2e}")
# 偶数落ちの確認: コム×コム干渉が奇数線に漏れないこと
a0,a1,a2=1/4,1/4+1/12,1/4+1/12+2/36
X = (comb(1,a0)+comb(3,a1)+comb(9,a2))**2 + 2*0  # コム同士の積のみ(参照なし) → 自己+交差
leak=max(abs(line(X,f)) for f in (1,3,5,7,9,11,13))
print(f"コム×コム項の奇数線への漏れ(最大|係数|): {leak:.2e} → {'ゼロ(偶数落ち確認)' if leak<1e-10 else '非ゼロ!'}")
