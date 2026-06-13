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

# ============ 観点別の4表（一辺/角度/面積/体積）============
def exists(R,d,a=1.0):
    return d*np.sin(a/(2*R))**2 <= 1.0
def face_angle_deg(R):       # 頂点角[度]（次元独立: cosθ=-tan²(1/2R)）
    return np.degrees(np.arccos(-np.tan(1/(2*R))**2))
def viewpoint_tables():
    Rlist=[round(0.5+0.5*i,1) for i in range(20)]+[100.0,1000.0,10000.0]
    npts_d={2:80,3:48,4:32,5:24}
    def cell(v): return f"{v:.7f}" if v is not None else "  n/a   "
    print("\n"+"="*72)
    print("【表1】一辺の長さ（1-content）— 全 d・全 R で 1.000000（不変）")
    print("【表2】頂点角 θ[度]（直角=90°の歪み）— d=2..5 で同値、n/a=存在閾値未満")
    print(f"{'R':>8} | {'d=1':>8} | {'d=2':>9} | {'d=3':>9} | {'d=4':>9} | {'d=5':>9}")
    for R in Rlist:
        th=face_angle_deg(R)
        cols=[]
        for d in (2,3,4,5):
            cols.append(f"{th:.5f}" if exists(R,d) else "  n/a  ")
        print(f"{R:>8} | {'—(辺のみ)':>8} | {cols[0]:>9} | {cols[1]:>9} | {cols[2]:>9} | {cols[3]:>9}")
    print("\n  直角からの超過 Δθ=θ-90°[度]（d非依存・正曲率で正）:")
    for R in (1.0,1.5,2.0,3.0,5.0,10.0,100.0):
        if exists(R,2): print(f"    R={R:>7}: Δθ={face_angle_deg(R)-90:.6f}°")
    print("\n"+"="*72)
    print("【表3】2-面の面積 A（2-content）— d=2..5 で同値（等質性）、n/a=存在閾値未満")
    print(f"{'R':>8} | {'d=1':>8} | {'d=2':>11} | {'d=3':>11} | {'d=4':>11} | {'d=5':>11}")
    for R in Rlist:
        A=volume_box_integral(R,2,80) if exists(R,2) else None
        cols=[cell(A if exists(R,d) else None) for d in (2,3,4,5)]
        print(f"{R:>8} | {'—':>8} | {cols[0]:>11} | {cols[1]:>11} | {cols[2]:>11} | {cols[3]:>11}")
    print("\n"+"="*72)
    print("【表4】体積 V_d（d-content）— 次元で真に異なる唯一の量、c_d=d(d-1)/12")
    print(f"{'R':>8} | {'d=1':>6} | {'d=2':>11} | {'d=3':>11} | {'d=4':>11} | {'d=5':>11}")
    for R in Rlist:
        cols=[]
        for d in (2,3,4,5):
            V=volume_box_integral(R,d,npts_d[d]) if exists(R,d) else None
            cols.append(cell(V))
        print(f"{R:>8} | {'1.000000':>6} | {cols[0]:>11} | {cols[1]:>11} | {cols[2]:>11} | {cols[3]:>11}")
viewpoint_tables()

# ============ §4.5 角度からの曲率逆算（内部観測者の曲率計）============
def curvature_from_angle(theta_deg):
    """頂点角θ[度]から曲率Kと半径Rを符号つきで逆算（次元非依存・a=1）。"""
    c=np.cos(np.radians(theta_deg))
    if abs(c)<1e-15: return 0.0, np.inf, "flat"
    if c<0:   # θ>90°: 正曲率 cosθ=-tan²(1/2R)
        h=np.arctan(np.sqrt(-c)); R=1/(2*h); return +1/R**2, R, "positive(spherical)"
    else:     # θ<90°: 負曲率 cosθ=+tanh²(1/2R)
        h=np.arctanh(np.sqrt(c)); R=1/(2*h); return -1/R**2, R, "negative(hyperbolic)"
def angle_from_R(R, sign):
    if sign>0: return np.degrees(np.arccos(-np.tan(1/(2*R))**2))
    else:      return np.degrees(np.arccos(+np.tanh(1/(2*R))**2))
def inversion_demo():
    print("\n"+"="*72)
    print("§4.5 曲率計: 測定した頂点角 θ → 曲率 K（符号つき）と半径 R")
    print(f"{'θ[度]':>10} | {'符号':>6} | {'K':>12} | {'R':>10} | 種別")
    for th in (84.0,88.0,89.5,90.0,90.5,92.0,96.0,107.36431):
        K,R,kind=curvature_from_angle(th)
        sgn="負" if K<0 else ("零" if K==0 else "正")
        Rs=f"{R:.5f}" if np.isfinite(R) else "∞"
        print(f"{th:>10} | {sgn:>6} | {K:>+12.6f} | {Rs:>10} | {kind}")
    print("\n往復検算（R→θ→逆算R, 正/負 両曲率）:")
    ok=True
    for sign in (+1,-1):
        for R in (1.5,3.0,10.0,100.0):
            th=angle_from_R(R,sign); K,Rb,_=curvature_from_angle(th)
            ok &= abs(Rb-R)<1e-6
    print(f"  全往復一致: {ok}")
inversion_demo()

# ============ §4.6 次元の曖昧さと次元の共役量 ============
def dmax(R):
    if R < 1/np.pi: return None
    return int(np.floor(1.0/np.sin(1/(2*R))**2))
def dimension_ambiguity():
    print("\n"+"="*72)
    print("§4.6 次元の天井 d_max(R)=⌊1/sin²(1/2R)⌋ と ±½ による次元の曖昧さ")
    print(f"{'R':>6} | {'κ=sin²(1/2R)':>12} | {'d_max(R-½)':>10} | {'d_max(R)':>9} | {'d_max(R+½)':>10} | {'相対幅Δd/d':>10}")
    for R in [0.5,0.7,1.0,1.5,2.0,3.0,5.0,10.0]:
        k=np.sin(1/(2*R))**2; dm=dmax(R)
        dl=dmax(R-0.5); dp=dmax(R+0.5)
        rel=(dp-(dl if dl else 1))/dm
        print(f"{R:>6} | {k:>12.6f} | {str(dl) if dl else '—':>10} | {dm:>9} | {dp:>10} | {rel:>10.3f}")
    print("\n  容量共役 d·κ≤1（予算1）、飽和 d·κ=1 が天井 R*_d。残容量 1−dκ=(w/R)²。κ≈K/4。")
    for d in (2,3,4,5):
        Rs=1/(2*np.arcsin(1/np.sqrt(d)))
        print(f"    d={d}: R*={Rs:.5f} で d·κ={d*np.sin(1/(2*Rs))**2:.8f}（飽和）")
dimension_ambiguity()
