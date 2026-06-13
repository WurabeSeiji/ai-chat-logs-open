#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
論文0 検算プログラム:正曲率定曲率空間 S^d(R) 上の辺長 a=1 測地正則超立方体の歪み。
- d=2: 頂点角 θ(R)=2 arcsin(1/(√2 cos(1/2R)))、面積 A=R²(4θ−2π)（ガウス・ボネ）
- d≥2 統一: 中心射影で構成した測地胞体の体積を厳密ヤコビアン積分で評価
        V_d(R)=∫_{[-t,t]^d} R^d w /(|y|²+w²)^{(d+1)/2} dy,
        t=R sin(1/2R), w=R√(1−d sin²(1/2R))
- 存在閾値 R ≥ 1/(2 arcsin(1/√d))
- d=2 で「ガウス・ボネ面積」と「体積積分」が一致することを検証（手法の妥当性）
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

# ---------- d=2 の厳密閉形式（ガウス・ボネ） ----------
def theta_d2(R, a=1.0):
    arg = 1.0/(np.sqrt(2.0)*np.cos(a/(2.0*R)))
    if arg > 1.0+1e-15:   # 存在しない
        return None
    return 2.0*np.arcsin(min(arg,1.0))
def area_d2_gaussbonnet(R, a=1.0):
    th = theta_d2(R,a)
    if th is None: return None
    return R*R*(4.0*th - 2.0*np.pi)

# ---------- 統一: 中心射影体積の Gauss-Legendre 積分 ----------
def threshold(d):
    return 1.0/(2.0*np.arcsin(1.0/np.sqrt(d)))
def t_w(R, d, a=1.0):
    s = np.sin(a/(2.0*R))
    if d*s*s > 1.0:        # w² < 0 → 構成不能
        return None, None
    t = R*s
    w = R*np.sqrt(1.0 - d*s*s)
    return t, w
def volume_box_integral(R, d, npts, a=1.0):
    t, w = t_w(R, d, a)
    if t is None: return None
    # [-t,t] 上の GL ノード（対称性で 1 octant にせず全域で安全に）
    x, wt = leggauss(npts)
    x = t*x; wt = t*wt           # スケール [-1,1]->[-t,t]
    # テンソル積を逐次縮約（d 次元、メモリ節約）
    # 積分核 g(y)=R^d w /(|y|²+w²)^{(d+1)/2}
    # |y|² = Σ x_i²。各軸を順に畳み込む。
    grids = np.meshgrid(*([x]*d), indexing='ij')
    r2 = sum(g*g for g in grids) + w*w
    kern = (R**d)*w / r2**((d+1)/2.0)
    W = wt
    for _ in range(d-1):
        W = np.multiply.outer(W, wt)
    return float(np.sum(kern*W))

def main():
    print("="*72)
    print("検証1: 存在閾値 R* = 1/(2 arcsin(1/√d))")
    for d in (2,3,4,5):
        print(f"  d={d}: R* = {threshold(d):.6f}")
    print()
    print("="*72)
    print("検証2: d=2 で『ガウス・ボネ面積』と『体積積分』が一致するか")
    for R in (0.7, 1.0, 2.0, 3.0, 5.0, 10.0):
        A_gb = area_d2_gaussbonnet(R)
        A_int = volume_box_integral(R, 2, 64)
        if A_gb is None:
            print(f"  R={R}: 存在せず（閾値 {threshold(2):.4f} 未満）")
        else:
            print(f"  R={R:5.1f}: GB={A_gb:.10f}  積分={A_int:.10f}  差={abs(A_gb-A_int):.2e}")
    print()
    print("="*72)
    print("検証3: 小角係数 c_d の確定（V_d(R)=1+c_d/R²+...）")
    # 解析予測: c_2 = 1/6
    for d in (2,3,4,5):
        npts = {2:80,3:48,4:32,5:24}[d]
        vals=[]
        for R in (1000.0, 3000.0, 10000.0):
            V = volume_box_integral(R, d, npts)
            vals.append((R, (V-1.0)*R*R))
        # Richardson 風: 大 R 極限が c_d
        cd = vals[-1][1]
        extras = "  (解析値 1/6=%.6f)"%(1/6) if d==2 else ""
        print(f"  d={d}: c_d ≈ {cd:.6f}  [R=1000:{vals[0][1]:.6f}, 3000:{vals[1][1]:.6f}, 10000:{vals[2][1]:.6f}]{extras}")
    print()
    print("="*72)
    print("検証4: 極限検算 R=100,1000,10000 で平坦値への 1/R² 減衰")
    for d in (2,3,4):
        npts = {2:80,3:48,4:32}[d]
        row=[]
        for R in (100.0,1000.0,10000.0):
            V=volume_box_integral(R,d,npts)
            row.append(V-1.0)
        print(f"  d={d}: V-1 = [{row[0]:.3e}, {row[1]:.3e}, {row[2]:.3e}]  (比 ~100² ごと)")
    print()
    print("="*72)
    print("数表（θ[d=2,rad] / A=V₂ / V₃ / V₄ / V₅）— n/a=存在閾値未満")
    print(f"{'R':>8} | {'theta':>10} | {'V2(A)':>11} | {'V3':>11} | {'V4':>11} | {'V5':>11}")
    Rlist = [round(0.5+0.5*i,1) for i in range(20)] + [100.0,1000.0,10000.0]
    npts_d = {2:80,3:48,4:32,5:24}
    table={}
    for R in Rlist:
        th = theta_d2(R)
        ths = f"{th:.7f}" if th is not None else "   n/a   "
        cells={}
        for d in (2,3,4,5):
            V = volume_box_integral(R,d,npts_d[d])
            cells[d] = f"{V:.8f}" if V is not None else "   n/a    "
        table[R]=(th,)+tuple(volume_box_integral(R,d,npts_d[d]) for d in (2,3,4,5))
        print(f"{R:>8} | {ths:>10} | {cells[2]:>11} | {cells[3]:>11} | {cells[4]:>11} | {cells[5]:>11}")
    print()
    print("="*72)
    print("検証5: R のバラツキ ±½（中心 R と区間端 R±½ の体積、d=4 例）")
    for R in (1.0,2.0,3.0,5.0,10.0):
        vc = volume_box_integral(R,4,32)
        vm = volume_box_integral(R-0.5,4,32) if R-0.5>threshold(4) else None
        vp = volume_box_integral(R+0.5,4,32)
        vms = f"{vm:.7f}" if vm is not None else "n/a"
        print(f"  R={R:4.1f}: V4(R-½)={vms}  V4(R)={vc:.7f}  V4(R+½)={vp:.7f}")

if __name__=="__main__":
    main()

def verify_cd_closedform():
    print()
    print("="*72)
    print("検証6: 係数の閉形式 c_d = d(d-1)/12 = C(d,2)/6 の確認")
    print("（解釈: d次元箱の体積超過 = C(d,2) 個の座標2平面それぞれの面積超過1/6 の和）")
    for d in (2,3,4,5,6):
        npts={2:80,3:48,4:32,5:24,6:18}[d]
        V=volume_box_integral(10000.0,d,npts)
        cd_num=(V-1.0)*1e8
        cd_cf=d*(d-1)/12.0
        print(f"  d={d}: 数値 c_d={cd_num:.6f}  閉形式 d(d-1)/12={cd_cf:.6f}  一致={abs(cd_num-cd_cf)<1e-3}")
if __name__=="__main__":
    verify_cd_closedform()
