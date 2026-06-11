#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺39 検証: R認可つき頂点
# E1: mod 8 選択定理 / E2: ±1 R取引つき2体分解 / E3: 段階的同一性連続の全数検査
import itertools, numpy as np
from collections import defaultdict

print("="*70)
print("E1: mod 8 選択定理 — 主観セクターが4軸でなければならない算術的理由")
print("="*70)
for d in range(1,9):
    # 4s = d個の奇数平方の和 ≡ d (mod 8)
    integer = (d%4==0)
    parity = "奇数" if d==4 else ("偶数" if d==8 else "-")
    print(f"  d={d}: 4s≡{d%8} (mod 8) → ラベル整数: {'YES' if integer else 'no '}"
          f"  整数時のラベル偶奇: {parity}")
print("  → 整数ラベル: d≡0 mod 4 のみ。奇数整数(単一状態の分類定理): d=4 が唯一。")
print()
print("="*70)
print("E2: ±1 R取引つき2体分解 (ラベル水準の全経路列挙)")
print("="*70)
def two_step_paths(parent, final):
    """parent → (a,b) [δ1] → spectator + (c,d) [δ2], δ∈{+1,-1}, 正味0"""
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
    # 重複除去
    return sorted(set(paths))
for final in [(5,3,1),(3,3,3)]:
    ps=two_step_paths(9,final)
    print(f"  9 → {final}: 経路 {len(ps)} 本")
    for ((a,b),x,(c,d),d1) in ps:
        t1="借入+1" if d1>0 else "貸出-1"
        print(f"    9→({a},{b}) [{t1}] → {x}→({c},{d}) [返済{-d1:+d}]")
    print(f"    各ステップは2体分裂 = 可視帳簿では奇パリティ違反(禁止)。R取引が3次頂点を認可。")
print()
print("="*70)
print("E3: 段階的同一性連続の全数検査 (ベクトル則 v_parent ∈ {±(v_a±v_b)} を毎ステップ)")
print("="*70)
def shell_cells(m):
    out=[]
    K=4
    for k in itertools.product(range(-K,K+1),repeat=4):
        if abs(sum((abs(t)+0.5)**2 for t in k)-m)<1e-9: out.append(np.array(k))
    return out
SH={m: shell_cells(float(m)) for m in (1,3,5,7,9)}
def line_ok(vp, va, vb):
    for x in (va+vb, va-vb):
        if np.array_equal(x,vp) or np.array_equal(-x,vp): return True
    return False
def count_full_paths(v9):
    """全ラベル経路×全セル割当で、両ステップの線連続が成立する経路数"""
    total=0; detail=defaultdict(int)
    for final in [(5,3,1),(3,3,3)]:
        for ((a,b),x,(c,d),d1) in two_step_paths(9,final):
            spec = b if x==a else a
            for va in SH[a]:
                for vb in SH[b]:
                    if a==b and tuple(va)>=tuple(vb): continue
                    if not line_ok(v9, va, vb): continue
                    vx = va if x==a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            if c==d and tuple(vc)>=tuple(vd): continue
                            if line_ok(vx, vc, vd):
                                total+=1
                                detail[((a,b),x,(c,d))]+=1
    return total, detail
for v9, name in [(np.array((2,1,0,0)),"(2,1,0,0)型"), (np.array((1,1,1,1)),"(1,1,1,1)型")]:
    tot, det = count_full_paths(v9)
    print(f"  親 v9={tuple(v9)} [{name}]: 完全連続経路 = {tot}")
    for k,v in sorted(det.items()):
        print(f"    経路 {k}: {v} 通り")
print()
print("  注: s=1 の断片は原点ベクトル0。x→(c,1) の線則は v_c=±v_x (ノルム一致が必要)。")
