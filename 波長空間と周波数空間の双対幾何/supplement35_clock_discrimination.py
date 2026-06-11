#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 時計帳簿の判別計算
# (1) 3進カスケード (補遺4, k=3, m=14) を ノルム時計 τ=1/√s' と 面時計 τ=2√s' で再計算
# (2) (9,5) 対と双子系のサイクル比: ノルム/面/裸線 の三候補
import numpy as np, math

print("="*70)
print("(1) カスケード時計の整合性 (k=3, m=14, S=3^14)")
k=3; m=14; S=float(k**m)
js=np.arange(0,m+1)
sp=np.array([float(k**(m-j)) for j in js])     # s'_j
# a (主項): a=(8pi^2/3 * S/s')^(1/4)
a=(8*math.pi**2/3*S/sp)**0.25
# ノルム時計: 周期 1/√s'
t_norm=np.cumsum(1/np.sqrt(sp))
# 面時計: 周期 2√s'
t_face=np.cumsum(2*np.sqrt(sp))
t_face_inf=2*math.sqrt(S)/(1-1/math.sqrt(k))
print(f"  面時計の飽和値 t_inf = 2√S/(1-k^(-1/2)) = {t_face_inf:.1f}")
print("  段j   s'        a       t_norm     p_norm    t_face     残り t_inf-t   p_face")
for j in range(1,m+1):
    pn=(math.log(a[j])-math.log(a[j-1]))/(math.log(t_norm[j])-math.log(t_norm[j-1]))
    dpf=math.log(t_face[j])-math.log(t_face[j-1])
    pf=(math.log(a[j])-math.log(a[j-1]))/dpf if dpf>1e-15 else float('inf')
    print(f"  {j:3d} {sp[j]:9.0f} {a[j]:8.2f} {t_norm[j]:10.4f} {pn:8.3f}  {t_face[j]:9.1f}  {t_face_inf-t_face[j]:10.3f} {pf:9.2f}")
# 中央窓フィット (j=4..9) ノルム時計
w=slice(4,10)
pn_fit=np.polyfit(np.log(t_norm[w]),np.log(a[w]),1)[0]
print(f"  ノルム時計 中央窓(j=4..9) フィット指数 p = {pn_fit:.4f}  (補遺4 の 0.507 と比較)")
print(f"  面時計: t は {t_face_inf:.1f} に飽和 → 冪則なし、有限 t で a→∞ (発散読み)")
print()
print("="*70)
print("(2) サイクル比(無次元・t 不要): 「Aが1周する間にBは何周するか」")
print()
print("  候補時計の定義:")
print("   ノルム: ν=√s (dressed エネルギー)")
print("   面    : ν=1/(2√s) (容器基本波)")
print("   裸線  : ν=|k| (親記録の実在線; ω²=s−‖k‖₁−1)")
print()
import itertools
def bare(v): return math.sqrt(sum(x*x for x in v))
v9a=(1,1,1,1); v9b=(2,1,0,0); v5=(1,1,0,0); v3=(1,0,0,0)
print("  --- (9,5) 不等対 (S=14, 2体チャネル) ---")
print(f"   ノルム: ν9:ν5 = 3:√5 = {3/math.sqrt(5):.4f}  (重い方が速い)")
print(f"   面    : ν9:ν5 = √5:3 = {math.sqrt(5)/3:.4f}  (軽い方が速い)")
print(f"   裸線  : s=9(1,1,1,1)型: |k|=2,  s=9(2,1,0,0)型: |k|=√5={math.sqrt(5):.4f}, s=5: |k|=√2={math.sqrt(2):.4f}")
print(f"           ν9:ν5 = 2:√2 = {2/math.sqrt(2):.4f}  または √5:√2 = {math.sqrt(5)/math.sqrt(2):.4f}  (軌道型に依存!)")
print()
print("  --- 双子系 (9,9) 異軌道対 ((1,1,1,1)型 × (2,1,0,0)型) ---")
print(f"   ノルム: 3:3 = 1 (縮退、ビートなし)")
print(f"   面    : 1 (縮退、ビートなし)")
print(f"   裸線  : 2:√5 = {2/math.sqrt(5):.4f} ≠ 1 → 同エネルギーでもビートする!")
print()
print("="*70)
print("(3) 構造チェック")
print(f"   ½共役: ν_norm × ν_face = √s × 1/(2√s) = 1/2 (全種) — 二候補は零点を介した双対")
print(f"   大域反転: 全時計の帳簿を一斉に入替えるとサイクル比は逆数になるだけ")
print(f"     (9,5): {3/math.sqrt(5):.4f} ↔ {math.sqrt(5)/3:.4f} (積=1)")
print(f"   → 1対の比単独では『どちらが速いか』のラベル替えと区別できない(ゲージ)")
print(f"   裸線時計のみ: 同エネルギー異軌道で比≠1 を予言 → ゲージ的反転で消えない実差")
print(f"   ドレッシング検算: s=ω²+‖k‖₁+1: 9=4+4+1 ✓ / 9=5+3+1 ✓ / 5=2+2+1 ✓ / 3=1+1+1 ✓")
