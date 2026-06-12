#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺63 検証: 階層¼桁の読み出しプロトコル (順写像の完成テスト)
# 構成(1D): 参照=原点断片(DC, 転写線) + 親コム(基本波1, 位置 a0=d/4) + 子コム(基本波k, 位置 a0+c/(4k))
# 読み出し: 周波数1の線の位相クラス → d、周波数kの線の位相クラス r=(k·d+c) mod 4 → c=(r−k·d) mod 4
# 判定: k=3 は全12配置で一意復号のはず / k=5 は c∈{0,4} の縮退(Z4情報不足)が出るはず
import numpy as np
N=240; xs=np.arange(N)/N
def proj_phase(I, f):
    """周波数 f の線の cos/sin 振幅 → 位相 φ∈[0,1) (周期1の割合) と振幅"""
    c = 2*np.mean(I*np.cos(2*np.pi*f*xs))
    s = 2*np.mean(I*np.sin(2*np.pi*f*xs))
    amp = np.hypot(c,s)
    phi = (np.arctan2(s,c)/(2*np.pi)) % 1.0
    return phi, amp
def quarter_class(phi):
    """位相→Z4 クラス (最近接¼) と量子化誤差"""
    q = int(np.round(phi*4)) % 4
    err = abs(phi*4 - np.round(phi*4))
    return q, err
for k in (3,5):
    ok=0; amb=0; tot=0
    fails=[]
    for d in range(4):
        for c in range(k):
            tot+=1
            a0 = d/4.0
            ac = a0 + c/(4.0*k)
            Psi = 1.0 + np.sqrt(2)*np.cos(2*np.pi*1*(xs-a0)) + np.sqrt(2)*np.cos(2*np.pi*k*(xs-ac))
            I = Psi**2
            phi1,_ = proj_phase(I,1)
            phiK,_ = proj_phase(I,k)
            d_read, e1 = quarter_class(phi1)
            r, e2 = quarter_class(phiK)
            c_read = (r - k*d_read) % 4
            if e1>0.01 or e2>0.01:
                fails.append((d,c,'非¼位相'))
                continue
            if d_read==d and c_read==c:
                ok+=1
            elif d_read==d and (c_read%4)==(c%4) and c>=4:
                amb+=1; fails.append((d,c,f'縮退: c={c} が c={c%4} と区別不能'))
            else:
                fails.append((d,c,f'誤復号 d={d_read} c={c_read}'))
    print(f"k={k}: 全{tot}配置 → 一意復号 {ok}, 縮退 {amb}, その他 {len(fails)-amb}")
    for f in fails[:6]: print(f"   {f}")
# 一般再帰 (3進・2階層: d, c1, c2) の復号
print()
print("3進・3レベル (d, c1, c2) の再帰復号テスト:")
k=3; ok=0; tot=0
for d in range(4):
    for c1 in range(3):
        for c2 in range(3):
            tot+=1
            a0=d/4.0; a1=a0+c1/12.0; a2=a1+c2/36.0
            Psi = (1.0 + np.sqrt(2)*np.cos(2*np.pi*(xs-a0))
                       + np.sqrt(2)*np.cos(2*np.pi*3*(xs-a1))
                       + np.sqrt(2)*np.cos(2*np.pi*9*(xs-a2)))
            I=Psi**2
            d_r,_  = quarter_class(proj_phase(I,1)[0])
            r1,_   = quarter_class(proj_phase(I,3)[0])
            r2,_   = quarter_class(proj_phase(I,9)[0])
            c1_r = (r1 - 3*d_r) % 4
            c2_r = (r2 - 3*r1) % 4          # 再帰: r2 = 9a2 = 3·(3a1) + c2 (mod 4 で読む)
            if (d_r,c1_r,c2_r)==(d,c1,c2): ok+=1
print(f"  全{tot}配置 → 一意復号 {ok}/{tot}")
