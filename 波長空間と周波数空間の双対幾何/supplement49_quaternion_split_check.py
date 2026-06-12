#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺49 検証: (1,3)集約 = C⊂H の選択、と振幅平面の同時出現
#  E1: B4 のうち四元数環自己同型であるものの数 (格子の量子化が H 構造をどこまで残すか)
#  E2: Lipschitz 単数群 Q8 と、時間方向 u の選択 → Z4=⟨u⟩=ガウス単数群の出現
#  E3: ノルム分裂の恒等式: 4s = N(w1)+N(w2) (w=奇ガウス整数), w=(1+i)v → 2s = N(v1)+N(v2) (N(v) 奇数)
#      — 時間段(4平方)が (1+i)² を介して振幅段(2平方,奇ノルム)の対に厳密分解する
#  E4: 極分解 q=|q|e^{uθ}: 位相 θ が ¼ 回転(Z4)に量子化されるのは q∈C_u のときに限ることの数え上げ
#  T1: 正準性の必要条件テスト: 規則 B1/B3 のガウス和を48親セル全数で計算し、|z|² の親非依存性を判定
import itertools, math, numpy as np
from collections import defaultdict

# ---------- 四元数 (基底 1,i,j,k ↔ 軸 0,1,2,3) ----------
def qmul(a, b):
    a0,a1,a2,a3 = a; b0,b1,b2,b3 = b
    return (a0*b0 - a1*b1 - a2*b2 - a3*b3,
            a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1,
            a0*b3 + a1*b2 - a2*b1 + a3*b0)

PERMS = list(itertools.permutations(range(4)))
SIGNS = list(itertools.product([1,-1], repeat=4))
def apply_g(perm, sign, v):
    return tuple(sign[i]*v[perm[i]] for i in range(4))

print("="*72)
print("E1: B4 (384) のうち四元数環自己同型 (積を保つ) の個数")
print("="*72)
basis = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)]
autos = []
for perm in PERMS:
    for sign in SIGNS:
        ok = True
        for a in basis:
            for b in basis:
                lhs = apply_g(perm, sign, qmul(a,b))
                rhs = qmul(apply_g(perm,sign,a), apply_g(perm,sign,b))
                if lhs != rhs: ok=False; break
            if not ok: break
        if ok: autos.append((perm,sign))
print(f"  環自己同型: {len(autos)} / 384")
fix_real = [g for g in autos if apply_g(g[0],g[1],(1,0,0,0))==(1,0,0,0)]
print(f"  うち実軸(中心)を固定: {len(fix_real)}  ← H 構造の選択は B4 を 384→{len(autos)} に破る")

print()
print("="*72)
print("E2: Lipschitz 単数群と時間方向 u の選択")
print("="*72)
units = [q for q in itertools.product([-1,0,1],repeat=4) if sum(x*x for x in q)==1]
# 群閉性チェック
S=set(units); closed = all(qmul(a,b) in S for a in units for b in units)
print(f"  単数: {len(units)} 個 (Q8), 群として閉じる: {closed}")
for u,name in [((0,1,0,0),'i'),((0,0,1,0),'j'),((0,0,0,1),'k')]:
    sub={(1,0,0,0)}; x=u
    while x not in sub: sub.add(x); x=qmul(x,u)
    print(f"  時間方向 u={name}: ⟨u⟩ = 位数{len(sub)} の巡回群 (Z4 = C_u のガウス単数群)")

print()
print("="*72)
print("E3: ノルム分裂の恒等式 (全セル検査, |k_i|<=4)")
print("="*72)
def gauss_norm(z): return z[0]*z[0]+z[1]*z[1]
def div_1pi(w):
    """w/(1+i) = w(1-i)/2 がガウス整数なら返す"""
    a,b = w
    p,q = (a+b), (b-a)
    if p%2 or q%2: return None
    return (p//2, q//2)
n_cells=0; okA=True; okB=True; okC=True
for k in itertools.product(range(-4,5),repeat=4):
    n_cells+=1
    # 裸ノルム: u=i のペアリング q=z1+z2 j, z1=(k0,k1), z2=(k2,k3)
    if gauss_norm((k[0],k[1]))+gauss_norm((k[2],k[3])) != sum(x*x for x in k): okA=False
    # ドレス版: w1=(2|k0|+1, 2|k1|+1), w2=(2|k2|+1, 2|k3|+1) — 奇ガウス整数 (ドレスは |k| の梯子)
    w1=(2*abs(k[0])+1, 2*abs(k[1])+1); w2=(2*abs(k[2])+1, 2*abs(k[3])+1)
    s4 = sum((2*abs(x)+1)**2 for x in k)   # 4s (ドレスラベル)
    if gauss_norm(w1)+gauss_norm(w2) != s4: okB=False
    v1=div_1pi(w1); v2=div_1pi(w2)
    if v1 is None or v2 is None: okC=False; continue
    if gauss_norm(v1)%2==0 or gauss_norm(v2)%2==0: okC=False
    if s4 != 2*(gauss_norm(v1)+gauss_norm(v2)):  # 2s = N(v1)+N(v2) ⇔ 4s = 2(N1+N2)
        okC=False
print(f"  検査セル数: {n_cells}")
print(f"  (a) 裸: |k|² = |z1|²+|z2|² : {'PASS' if okA else 'FAIL'}")
print(f"  (b) ドレス: 4s = N(w1)+N(w2), w=奇ガウス整数: {'PASS' if okB else 'FAIL'}")
print(f"  (c) w=(1+i)v の厳密因数分解, N(v) 奇数, 2s = N(v1)+N(v2): {'PASS' if okC else 'FAIL'}")
print(f"  → 時間段(4平方)のドレスラベルは、記録アルファベットの基底 (1+i) を1回ずつ")
print(f"    くくり出すと、奇ノルムの振幅段(2平方)成分の対 2s=N(v1)+N(v2) に厳密分解する。")

print()
print("="*72)
print("E4: 極分解 q=|q|e^{uθ} と位相の Z4 量子化 (u=i, |k_i|<=3)")
print("="*72)
in_Cu=0; tot=0; z4_exact=0
for k in itertools.product(range(-3,4),repeat=4):
    if k==(0,0,0,0): continue
    tot+=1
    if k[2]==0 and k[3]==0:
        in_Cu+=1
        # C_u 内: θ = arg(k0+k1 i)。Z4 (¼回転) に乗るのは軸上のみ
        if k[0]==0 or k[1]==0: z4_exact+=1
print(f"  全セル {tot} のうち C_u=⟨1,i⟩ 平面内: {in_Cu}")
print(f"  うち位相が厳密に Z4 (¼回転格子) 上: {z4_exact}")
print(f"  → 連続位相 θ は C_u 内の格子方向の階層極限としてのみ稠密化 (補遺43 §3 と同型)")

print()
print("="*72)
print("T1: 正準性テスト — B1/B3 のガウス和を48親セル全数で (親非依存なら必要条件PASS)")
print("="*72)
def two_step_paths(parent, final):
    paths=[]; fin=tuple(sorted(final))
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
def branch_xor(vc,vd):
    q=0
    for i in range(4):
        c,d = vc[i], vd[i]
        if c!=0 and d!=0 and ((c<0)!=(d<0)): q+=1
    return q%4
def eps(v): return (-1)**int(np.abs(v).sum())

LBL = {lbl: two_step_paths(9, fin) for lbl,fin in [('531',(5,3,1)),('333',(3,3,3))]}
results=defaultdict(set); npaths=set()
parents=[np.array(v) for v in itertools.product(range(-4,5),repeat=4)
         if sorted(map(abs,v))==[0,0,1,2] and abs(sum((abs(t)+0.5)**2 for t in v)-9)<1e-9]
assert len(parents)==48
for v9 in parents:
    z={'B1':defaultdict(complex),'B3':defaultdict(complex)}
    cnt=0
    for lbl, lpaths in LBL.items():
        for ((a,b),x,(c,d),d1) in lpaths:
            for va in SH[a]:
                for vb in SH[b]:
                    if a==b and tuple(va)>=tuple(vb): continue
                    if not line_ok(v9,va,vb): continue
                    vx = va if x==a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            if c==d and tuple(vc)>=tuple(vd): continue
                            if line_ok(vx,vc,vd):
                                cnt+=1
                                p1 = (branch_xor(va,vb)+branch_xor(vc,vd))%4
                                z['B1'][lbl] += 1j**p1
                                p3 = 0 if eps(v9)*eps(va)*eps(vb)>0 else 2
                                p3 = (p3 + (0 if eps(vx)*eps(vc)*eps(vd)>0 else 2))%4
                                z['B3'][lbl] += 1j**p3
    npaths.add(cnt)
    for rule in ('B1','B3'):
        results[rule].add((round(abs(z[rule]['531'])**2,6), round(abs(z[rule]['333'])**2,6)))
print(f"  管理連鎖経路数の集合 (48親): {npaths}")
for rule in ('B1','B3'):
    vals = results[rule]
    print(f"  規則{rule}: (|z531|²,|z333|²) の値の集合 = {vals}")
    print(f"    → 親非依存性: {'PASS (必要条件クリア)' if len(vals)==1 else 'FAIL (ゲージ依存→非正準)'}")
