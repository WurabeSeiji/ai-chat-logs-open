#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺46 §3 決定実験: 12/2経路への Z4 位相割当とガウス和 → 測度三候補との照合
# 割当規則の候補を明示的に列挙し、各規則の下で z=Σ i^p を計算、|z|^2 比を判定表に掛ける。
import itertools, numpy as np
from fractions import Fraction as F
from collections import defaultdict

# ---------- ラベル水準の2段経路 (補遺39 E2 と同一) ----------
def two_step_paths(parent, final):
    paths=[]
    fin=tuple(sorted(final))
    for a in range(1,parent+2,2):
        for b in range(a,parent+2,2):
            d1=a+b-parent
            if d1 not in (1,-1): continue
            for (x,spec) in ((a,b),(b,a)):
                for c in range(1,x+2,2):
                    dd=x-d1-c
                    if dd<1 or dd%2==0: continue
                    if tuple(sorted((spec,c,dd)))==fin:
                        paths.append(((a,b),x,(c,dd),d1))
    return sorted(set(paths))

P531 = two_step_paths(9,(5,3,1))
P333 = two_step_paths(9,(3,3,3))
print("="*72)
print(f"ラベル水準経路: 9→(5,3,1): {len(P531)} 本 / 9→(3,3,3): {len(P333)} 本")
for p in P531: print("   531:", p)
for p in P333: print("   333:", p)

# ---------- 規則A: 不可視方向の符号 (ラベル水準) ----------
print()
print("="*72)
print("規則A: 位相 = 各仮想頂点で渡る不可視方向の符号")
print("="*72)
def gauss(zs):
    z=sum(zs); return z, abs(z)**2
# A1: i^{δ1} (最初のR取引の符号)
zA1_531,_ = gauss([1j**p[3] for p in P531]); zA1_333,_ = gauss([1j**p[3] for p in P333])
# A2: i^{δ1}·i^{δ2} = i^{δ1-δ1} = 1 (取引の総和: 正味ゼロ)
zA2_531,_ = gauss([1j**0 for p in P531]); zA2_333,_ = gauss([1j**0 for p in P333])
n531_plus = sum(1 for p in P531 if p[3]>0); n333_plus = sum(1 for p in P333 if p[3]>0)
print(f"  δ1 の符号分布: 531 → +1:{n531_plus} / -1:{len(P531)-n531_plus},  333 → +1:{n333_plus} / -1:{len(P333)-n333_plus}")
print(f"  A1 (i^δ1):       z531={zA1_531}, |z|²={abs(zA1_531)**2:.0f} / z333={zA1_333}, |z|²={abs(zA1_333)**2:.0f}")
print(f"  A2 (i^(δ1+δ2)):  z531={zA2_531}, |z|²={abs(zA2_531)**2:.0f} / z333={zA2_333}, |z|²={abs(zA2_333)**2:.0f}  ← 全位相同符号=古典退化")

# ---------- セル水準: 管理連鎖経路 (補遺39 E3) + 枝XOR位相 ----------
print()
print("="*72)
print("規則B: 位相 = 娘の相対枝データ (セル水準・管理連鎖つき全経路)")
print("="*72)
def shell_cells(m):
    out=[]; K=4
    for k in itertools.product(range(-K,K+1),repeat=4):
        if abs(sum((abs(t)+0.5)**2 for t in k)-m)<1e-9: out.append(np.array(k))
    return out
SH={m: shell_cells(float(m)) for m in (1,3,5,7,9)}
def line_ok(vp,va,vb):
    for x in (va+vb,va-vb):
        if np.array_equal(x,vp) or np.array_equal(-x,vp): return True
    return False

def branch_bits(v):
    """軸ごとの枝ビット: k_i>0=cos(0), k_i<0=sin(1), k_i=0=不活性"""
    return [(1 if t<0 else 0) if t!=0 else None for t in v]

def step_phase_B1(vc, vd):
    """B1: 娘対の枝XOR — 枝ビットが異なる活性軸の本数 (mod 4)"""
    bc, bd = branch_bits(vc), branch_bits(vd)
    q=0
    for i in range(4):
        if bc[i] is not None and bd[i] is not None and bc[i]!=bd[i]: q+=1
    return q%4

def step_phase_B2(vc, vd):
    """B2: 娘対の sin 枝総数 (¼オフセットの総和, mod 4)"""
    q=0
    for v in (vc,vd):
        for i,t in enumerate(v):
            if t<0: q+=1
    return q%4

def eps(v): return (-1)**int(np.abs(v).sum())
def step_phase_B3(vx, vc, vd):
    """B3: 仮想中間のε簿記: 分裂者と娘対のε積 (Z2⊂Z4: 0 or 2)"""
    return 0 if eps(vx)*eps(vc)*eps(vd)>0 else 2

def cell_paths(v9):
    """管理連鎖が両ステップ成立する全セル経路を列挙 (補遺39 E3 + 経路データ保持)"""
    out=[]
    for final,label in [((5,3,1),'531'),((3,3,3),'333')]:
        for ((a,b),x,(c,d),d1) in two_step_paths(9,final):
            for va in SH[a]:
                for vb in SH[b]:
                    if a==b and tuple(va)>=tuple(vb): continue
                    if not line_ok(v9,va,vb): continue
                    vx = va if x==a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            if c==d and tuple(vc)>=tuple(vd): continue
                            if line_ok(vx,vc,vd):
                                out.append((label,(va,vb),vx,(vc,vd),d1))
    return out

v9 = np.array((2,1,0,0))
paths = cell_paths(v9)
n531 = sum(1 for p in paths if p[0]=='531'); n333 = len(paths)-n531
print(f"  親 v9=(2,1,0,0): 管理連鎖セル経路 = {len(paths)} 本 (531終: {n531} / 333終: {n333})")

for name, rule in [("B1 枝XOR", 'B1'), ("B2 sin総数", 'B2'), ("B3 ε簿記", 'B3')]:
    zz=defaultdict(complex); cnt=defaultdict(lambda: defaultdict(int))
    for (label,(va,vb),vx,(vc,vd),d1) in paths:
        if rule=='B1':
            p = (step_phase_B1(va,vb)+step_phase_B1(vc,vd))%4
        elif rule=='B2':
            p = (step_phase_B2(va,vb)+step_phase_B2(vc,vd))%4
        else:
            p = (step_phase_B3(v9,va,vb)+step_phase_B3(vx,vc,vd))%4
        zz[label]+= 1j**p; cnt[label][p]+=1
    z1,z2 = zz['531'], zz['333']
    W1,W2 = abs(z1)**2, abs(z2)**2
    print(f"  [{name}] 位相分布 531:{dict(cnt['531'])} 333:{dict(cnt['333'])}")
    print(f"           z531={z1:.0f} |z|²={W1:.0f} / z333={z2:.0f} |z|²={W2:.0f}", end="")
    if W1+W2>0:
        f1=W1/(W1+W2)
        print(f"  → 分岐比 {f1:.4f}:{1-f1:.4f}")
        if W2>0:
            PX = 2*W1*W2/(2*W1*W2+W2*W2)
            print(f"           双子二段 P(X)=2W₁W₂/(2W₁W₂+W₂²)={PX:.5f}")
    else:
        print("  → 完全相殺")

# ---------- 照合表 ----------
print()
print("="*72)
print("測度候補との照合 (論文12 §5.3 / 補遺36 §5)")
print("="*72)
print(f"  一括(配置等重): P(X)=60/61={60/61:.5f}   単一崩壊比 192:56 → {192/248:.4f}")
print(f"  逐次(条件付き): P(X)=396/403={396/403:.5f}")
print(f"  経路等重(配置): P(X)=0.96000")
print(f"  R頂点経路等重:  単一崩壊比 12:2 → {12/14:.4f}")
print(f"  管理連鎖等重:   単一崩壊比 {n531}:{n333} → {n531/(n531+n333):.4f}")
