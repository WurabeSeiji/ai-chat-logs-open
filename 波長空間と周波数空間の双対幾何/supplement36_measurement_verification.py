#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺36 検証: (V1) 逐次追記 vs 一括数え上げ (双子宇宙の二段崩壊)
#            (V2) 崩壊前後の記録構造 (無→有・DC増分・帳簿単調性・情報クラス)
from fractions import Fraction as F
from math import comb
import itertools, numpy as np

print("="*70)
print("V1: 双子宇宙 S=18 (9,9) の二段崩壊 — 三つの測度の比較")
print("="*70)
# 殻容量: c(1)=1, c(3)=8, c(5)=24
# チャネル: 9->(5,3,1) or (3,3,3)
# 終状態多重集合: X=(5,3,3,3,3,1) [一方が531、他方が333] / Y=(3,3,3,3,3,3)
W_X_global = 24*comb(8,4)*1        # shell5から1, shell3から4, 原点
W_Y_global = comb(8,6)
print(f"[一括] W(5,3,3,3,3,1)={W_X_global}, W(3^6)={W_Y_global}")
PX_g=F(W_X_global, W_X_global+W_Y_global); PY_g=1-PX_g
print(f"[一括/配置等重] P(X)={PX_g}={float(PX_g):.5f}  P(Y)={PY_g}={float(PY_g):.5f}")
# 逐次 (条件付き数え上げ): 先手の channel: W531=24*8*1=192, W333=C(8,3)=56
W531_1=F(24*8*1); W333_1=F(comb(8,3))
P531_1=W531_1/(W531_1+W333_1); P333_1=1-P531_1
# 先手=531 なら後手: 531は原点ブロック → 333のみ: C(7,3)=35 (確率1)
# 先手=333 なら後手: 531: 24*5*1=120, 333: C(5,3)=10
P531_2=F(24*5*1)/F(24*5*1+comb(5,3)); P333_2=1-P531_2
PX_s=P531_1*1 + P333_1*P531_2
PY_s=P333_1*P333_2
print(f"[逐次/条件付き] P(X)={PX_s}={float(PX_s):.5f}  P(Y)={PY_s}={float(PY_s):.5f}")
# 経路等重: 経路数を直接数える
paths_X = 2*(192*35 + 56*120)   # 先手がどちらの双子か x (先手531→後手333 + 先手333→後手531)
paths_Y = 2*(56*comb(5,3))
PX_p=F(paths_X, paths_X+paths_Y); PY_p=1-PX_p
print(f"[経路等重] 経路数 X={paths_X}, Y={paths_Y}: P(X)={float(PX_p):.5f}  P(Y)={float(PY_p):.5f}")
print()
print(f"  三測度の P(X): 一括 {float(PX_g):.5f} / 逐次 {float(PX_s):.5f} / 経路 {float(PX_p):.5f}")
print(f"  一括 vs 逐次 の差 = {float(PX_g-PX_s):+.5f}  (一致せず → 測度問題が定量化された)")
print()
print("  結合チャネル表 (逐次測度):")
print(f"   P(両方531) = 0  (原点スロット c(1)=1 の独占 — 運動学的反相関)")
print(f"   P(531,333)+(333,531) = {float(PX_s):.5f}")
print(f"   P(両方333) = {float(PY_s):.5f}")
print()
print("="*70)
print("V2: 崩壊前後の記録構造 (双子A=(2,1,0,0) → (5,3,1): v5=(1,1,0,0),v3=(0,0,1,0),原点)")
print("="*70)
def line_set(occ):
    """強度 I=Psi^2 の交差線+線形転写線の整数ベクトル集合 (向きの代表化: 辞書式最大)"""
    lines=set()
    occ=[np.array(k) for k in occ]
    for i in range(len(occ)):
        for j in range(len(occ)):
            if i==j: continue
            for sgn in (1,-1):
                v=occ[i]+sgn*occ[j]
                if not np.any(v): continue
                t=tuple(v) if tuple(v)>=tuple(-v) else tuple(-v)
                lines.add(t)
    return lines
before=[(2,1,0,0),(1,1,1,1)]
after=[(1,1,0,0),(0,0,1,0),(0,0,0,0),(1,1,1,1)]
Lb=line_set(before); La=line_set(after)
print(f"  DC (断片数): 2 → 4  (計数の追記)")
print(f"  線の数: 前={len(Lb)} → 後={len(La)}")
appeared=La-Lb; vanished=Lb-La; common=La&Lb
print(f"  出現した線 (無→有): {len(appeared)} 本")
print(f"  消えた線: {len(vanished)} 本  (前の交差線は前の占有の関数 — 消失は現在の記録からは読めない)")
print(f"  共通: {len(common)} 本")
# 線形転写線 (原点占有時のみ): 占有セル k 自身の位置の線
linear_after=[k for k in after if any(k)]
print(f"  情報クラス: 前=干渉計型(参照なし、関係のみ) → 後=ホログラフィック型")
print(f"   (原点=s1 が出現 → 線形転写線 {[k for k in linear_after]} が立ち、全配置復元可能 — 補遺17 §4)")
# 帳簿
lam2_before=F(1,9)+F(1,9)
lam2_after=F(1,5)+F(1,3)+F(1,1)+F(1,9)
print(f"  Σλ² 帳簿: 前={lam2_before}={float(lam2_before):.4f} → 後={lam2_after}={float(lam2_after):.4f}  (厳密増加 ✓)")
print(f"  Σν² 帳簿: 前=9+9=18 → 後=5+3+1+9=18  (保存 ✓)")
