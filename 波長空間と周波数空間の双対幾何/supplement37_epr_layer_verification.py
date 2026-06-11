#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺37 検証: (1) 関係限定読み出し定理 (自己スペクトルは偶数ベクトル → 奇数セクター不可視)
#             (2) 条件付き統計シフト (Aのイベントが Bの統計を瞬時に変える定量表)
import numpy as np, itertools
from fractions import Fraction as F
from math import comb

print("="*70)
print("(1) 親記録の線の分類: 自己線 / 交差線 × 奇数アルファベット可視性")
print("="*70)
def classify(a,b):
    a=np.array(a); b=np.array(b)
    selfA=set(); selfB=set(); cross=set()
    # 自己線: 成分 in {0, ±2k_i} (少なくとも1成分非零)
    def selfs(k):
        out=set()
        for mask in itertools.product([0,1],repeat=4):
            if not any(mask): continue
            for sg in itertools.product([1,-1],repeat=4):
                v=tuple(2*k[i]*mask[i]*sg[i] for i in range(4))
                if not any(v): continue
                t=v if v>=tuple(-x for x in v) else tuple(-x for x in v)
                out.add(t)
        return out
    selfA=selfs(a); selfB=selfs(b)
    for sg in (1,-1):
        v=tuple(a[i]+sg*b[i] for i in range(4))
        if any(v):
            t=v if v>=tuple(-x for x in v) else tuple(-x for x in v)
            cross.add(t)
    def odd_visible(v):  # 全成分が {0, 奇数}
        return all(x==0 or abs(x)%2==1 for x in v)
    return selfA, selfB, cross, odd_visible

for (a,b,name) in [((2,1,0,0),(0,-1,2,0),"同型対 (2,1,0,0)x(2,1,0,0), |Δ|²=10"),
                   ((2,1,0,0),(1,1,1,1),"混合対 (2,1,0,0)x(1,1,1,1), |Δ|²=3"),
                   ((1,1,1,1),(1,1,1,-1),"同型対 (1,1,1,1)x(1,1,1,1), |Δ|²=4")]:
    sA,sB,cr,ov=classify(a,b)
    sA_vis=[v for v in sA if ov(v)]; sB_vis=[v for v in sB if ov(v)]
    cr_vis=[v for v in cr if ov(v)]
    print(f"  {name}")
    print(f"   自己線A: {len(sA)}本 → 奇数セクター可視 {len(sA_vis)}本")
    print(f"   自己線B: {len(sB)}本 → 奇数セクター可視 {len(sB_vis)}本")
    print(f"   交差線 : {len(cr)}本 {sorted(cr)} → 奇数セクター可視 {len(cr_vis)}本 {sorted(cr_vis)}")
print()
print("  定理(一行証明): 単一断片の自己スペクトルは成分 {0,±2k_i} = 全成分偶数。")
print("  非零の偶成分は {0,奇数} に属さない → 奇数(論理波)セクターでは自己線は恒等的に不可視。")
print("  → 論理波復号器が読めるのは DC(計数) と 奇数可視の交差線(関係) のみ。個体は読めない。")
print()
print("="*70)
print("(2) 条件付き統計シフト: Aのイベントは Bの分岐統計を瞬時に書き換える")
print("="*70)
# B 単独(Aが未崩壊): shell5=24 shell3=8 origin=1 全部空き
W531=F(24*8*1); W333=F(comb(8,3))
P0=(W531/(W531+W333), W333/(W531+W333))
# A が (5,3,1) を取った後: origin 占有 → B: 531=0, 333=C(7,3)
P1=(F(0), F(1))
# A が (3,3,3) を取った後: shell3 残5 → B: 531=24*5*1, 333=C(5,3)
W531b=F(24*5); W333b=F(comb(5,3))
P2=(W531b/(W531b+W333b), W333b/(W531b+W333b))
print(f"  Bの (531, 333) 確率:")
print(f"   A未崩壊      : ({float(P0[0]):.4f}, {float(P0[1]):.4f})")
print(f"   A=(5,3,1)後  : ({float(P1[0]):.4f}, {float(P1[1]):.4f})   全変動距離 = {float(abs(P1[0]-P0[0])):.4f}")
print(f"   A=(3,3,3)後  : ({float(P2[0]):.4f}, {float(P2[1]):.4f})   全変動距離 = {float(abs(P2[0]-P0[0])):.4f}")
print(f"  シフトは |Δ|(双子の格子分離) に一切依存しない (大域制約のため) — 距離非依存の瞬時相関")
print()
print("(3) パリティ保護: S=18 は偶数 → 単一状態スロットなし (分類定理)")
print("  対の合成ラベル 18 = 偶数: もつれた対は上のレイヤで「一個の粒子」になれない")
print("  (粒子=奇数ラベル単一状態。対は偶数 → 別種の対象として算術的に区別される)")
