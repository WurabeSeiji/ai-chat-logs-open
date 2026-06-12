#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺64 検証: 写像の4D昇格 — 積波の読み出しと対角縮退の有無
# T1: 状態水準 — 積波 Φ=Πφ(x_j−a_j) は全軸同時の半波シフト (d_j→d_j+2 ∀j) で不変か?
# T2: 読み出し水準 — 和差線 (1,±1,...) の¼クラスから桁タプルは「対角Z2商」まで一意復号できるか?
# T3: 商の構造 — 縮退クラスの対構造は厳密に (d)~(d+2,...,+2) のみか?
import numpy as np, itertools
N=32; xs=np.arange(N)/N
X = [xs.reshape(-1,1), xs.reshape(1,-1)]   # 2軸で原理判定(対角構造は軸数に依らない)
def wave2(a1,a2):
    return 2*np.cos(2*np.pi*(X[0]-a1))*np.cos(2*np.pi*(X[1]-a2))
def line2(I,f1,f2):
    e = np.exp(-2j*np.pi*(f1*X[0]+f2*X[1]))
    return 2*np.mean(I*np.conj(e))   # e^{+2πi(f·x)} 係数 → 位相 +f·a
def qc(z):
    phi=(np.angle(z)/(2*np.pi))%1.0
    return int(np.round(phi*4))%4, abs(phi*4-np.round(phi*4))
# T1: 状態水準の同一性
pairs_same=0; pairs_diff=0
for d1,d2 in itertools.product(range(4),repeat=2):
    w  = wave2(d1/4,d2/4)
    w2 = wave2((d1+2)%4/4,(d2+2)%4/4)
    if np.allclose(w,w2,atol=1e-12): pairs_same+=1
    else: pairs_diff+=1
print(f"T1: Φ(d1,d2) ≡ Φ(d1+2,d2+2) : {pairs_same}/16 (全16なら対角半波シフトは状態の恒等変換)")
# 片軸のみ+2 は別状態であることの対照
w  = wave2(0,0); w3 = wave2(2/4,0)
print(f"    対照: 片軸のみ+2 は別状態: {not np.allclose(w,w3)}")
# T2: 読み出し — 和差線から復号、商クラスの分離能
sigs={}
for d1,d2 in itertools.product(range(4),repeat=2):
    I=(1.0+wave2(d1/4,d2/4))**2
    rp,e1 = qc(line2(I,1, 1))   # 位相 a1+a2 → (d1+d2) mod 4
    rm,e2 = qc(line2(I,1,-1))   # 位相 a1−a2 → (d1−d2) mod 4
    # 積波自身の単独線 (1,0) は存在しない(積のスペクトルは(±1,±1)のみ)ことも確認
    z10=line2(I,1,0)
    sigs.setdefault((rp,rm),[]).append((d1,d2))
    assert abs(z10)<1e-10, "純軸線が存在?!"
    assert e1<1e-9 and e2<1e-9
nclass=len(sigs)
ok_quot=True
for sig,members in sigs.items():
    if len(members)!=2: ok_quot=False
    else:
        (a,b),(c,d)=members
        if not ((c-a)%4==2 and (d-b)%4==2): ok_quot=False
print(f"T2: 読み出しシグネチャ(和線,差線)のクラス数: {nclass}/16 配置 (8なら対角Z2商と一致)")
print(f"T3: 全クラスが厳密に (d1,d2)~(d1+2,d2+2) の対: {'PASS' if ok_quot else 'FAIL'}")
print(f"    純軸線 (1,0) の振幅: ゼロ(積波に単独軸線なし) 確認済み")
# T4: 階層(2レベル)でも商構造が保たれるか — レベルごとに独立な対角Z2?
sigs2={}
for d1,d2,c1,c2 in itertools.product(range(4),range(4),range(3),range(3)):
    a1=d1/4+c1/12; a2=d2/4+c2/12
    Psi=1.0+wave2(d1/4,d2/4)+2*np.cos(2*np.pi*3*(X[0]-a1))*np.cos(2*np.pi*3*(X[1]-a2))
    I=Psi**2
    key=[]
    for (f1,f2) in [(1,1),(1,-1),(3,3),(3,-3)]:
        r,e=qc(line2(I,f1,f2)); key.append(r)
        assert e<1e-9
    sigs2.setdefault(tuple(key),[]).append((d1,d2,c1,c2))
sizes=sorted(set(len(v) for v in sigs2.values()))
print(f"T4: 2レベル4桁 全{4*4*3*3}配置 → クラス数 {len(sigs2)}, クラスサイズの種類 {sizes}")
# クラス内の対の構造を標本表示
sample=[v for v in sigs2.values() if len(v)>1][:3]
for v in sample: print(f"    クラス例: {v}")

# T5: 複合体(セル和・パリティ混在)では縮退が破れるか
# Φ_A: m=(1,1) (活性軸数 偶) + Φ_B: m=(1,0) (活性軸数 奇) → 半波同時シフトで A は不変・B は符号反転
print()
sigs5={}
for d1,d2 in itertools.product(range(4),repeat=2):
    a1,a2=d1/4.0,d2/4.0
    Psi = 1.0 + 2*np.cos(2*np.pi*(X[0]-a1))*np.cos(2*np.pi*(X[1]-a2)) + np.sqrt(2)*np.cos(2*np.pi*(X[0]-a1))
    I=Psi**2
    r10,e0=qc(line2(I,1,0))    # B の単独軸線 → d1 直読み
    r11,e1=qc(line2(I,1,1))    # A の和線 → (d1+d2) mod 4
    assert e0<1e-9 and e1<1e-9
    d1_r=r10; d2_r=(r11-r10)%4
    sigs5.setdefault((r10,r11),[]).append((d1,d2))
    assert (d1_r,d2_r)==(d1,d2), f"誤復号 {d1},{d2}→{d1_r},{d2_r}"
print(f"T5: パリティ混在複合体 全16配置 → シグネチャ {len(sigs5)} 種, 全配置一意復号 PASS")
print("    機構: 活性軸数が奇のセル(Φ_B)が単独軸線を供給し、対角縮退を完全に破る")
