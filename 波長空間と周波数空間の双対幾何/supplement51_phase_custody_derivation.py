#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺51 検証: 位相つき管理連鎖 — 接続則の導出
# 原理: 記録連続則の強い形 (補遺38 §6.1 で未実施と登録されたもの)。
#   親の同一性線 v_p は娘対の交差スペクトルに「位置」だけでなく「位相」ごと存続する。
#   娘交差線の複素指数係数 C_cross(v_p) と親自身の係数 C_par(v_p) の比が、
#   そのステップの輸送位相 η を【辞書の三角恒等式から一意に】与える — 割当でなく導出。
# 検査:
#   D1: η が常に Z4 ({1,i,-1,-i}) に乗ること (導出位相のアルファベット確認)
#   D2: 導出規則 D のガウス和 (経路位相 = η1·η2) → B1/B3/新規との照合
#   D3: 強い連続 (η=1) の経路数
#   D4: 等価性: 48親セル全数で |z|² が親非依存か / B4 移送での位相多重集合の不変性
import itertools, numpy as np, cmath, random
from collections import defaultdict

SQ2 = np.sqrt(2.0)

def expcoeffs(v):
    """断片セル v の波 Φ_v の複素指数係数 {周波数ベクトル n: 係数}
       φ_0=1; φ_m (m>0) = √2cos(2πmx) = (e^{+}+e^{-})/√2
       φ_m (m<0) = √2sin(2π|m|x) = (e^{+}-e^{-})/(i√2) = -i/√2 e^{+} + i/√2 e^{-}"""
    coeffs = {(): 1.0+0j}
    out = {((),): None}
    acc = [((), 1.0+0j)]
    for i in range(4):
        m = v[i]
        new = []
        if m == 0:
            opts = [(0, 1.0+0j)]
        elif m > 0:
            opts = [(+m, 1/SQ2+0j), (-m, 1/SQ2+0j)]
        else:
            mm = -m
            opts = [(+mm, -1j/SQ2), (-mm, +1j/SQ2)]
        for n, c in acc:
            for f, w in opts:
                new.append((n+(f,), c*w))
        acc = new
    d = {}
    for n, c in acc:
        d[n] = d.get(n, 0) + c
    return d

def cross_coeff_at(va, vb, target):
    """Φ_va·Φ_vb の周波数 target における複素指数係数"""
    A = expcoeffs(tuple(va)); B = expcoeffs(tuple(vb))
    tot = 0+0j
    for na, ca in A.items():
        nb = tuple(t-x for t, x in zip(target, na))
        if nb in B:
            tot += ca * B[nb]
    return tot

def self_coeff_at(v, target):
    return expcoeffs(tuple(v)).get(tuple(target), 0+0j)

def transport_phase(vp, va, vb):
    """親線 vp が娘対 (va,vb) の交差線に存続するときの輸送位相 η。
       s=1 娘 (vb=0) は線形転写: Φ_va·Φ_0 = Φ_va。存続しなければ None。"""
    vp = tuple(vp)
    C_par = self_coeff_at(vp, vp)
    if abs(C_par) < 1e-12:
        return None
    C_cross = cross_coeff_at(va, vb, vp)
    if abs(C_cross) < 1e-12:
        return None
    eta = (C_cross/abs(C_cross)) / (C_par/abs(C_par))
    return eta

# ---------- 経路機械 (補遺39/47/49 と同一) ----------
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
LBL = {lbl: two_step_paths(9, fin) for lbl,fin in [('531',(5,3,1)),('333',(3,3,3))]}

def z4_round(eta):
    """ηが Z4 元に一致するか判定し丸める"""
    for p,u in enumerate([1,1j,-1,-1j]):
        if abs(eta-u)<1e-9: return p
    return None

def run_parent(v9):
    """位相つき管理連鎖: 各経路の導出位相 η=η1·η2 を計算"""
    res = {'531': [], '333': []}
    nonz4 = 0
    for lbl, lpaths in LBL.items():
        for ((a,b),x,(c,d),d1) in lpaths:
            for va in SH[a]:
                for vb in SH[b]:
                    if a==b and tuple(va)>=tuple(vb): continue
                    e1 = transport_phase(v9, va, vb)
                    if e1 is None: continue
                    vx = va if x==a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            if c==d and tuple(vc)>=tuple(vd): continue
                            e2 = transport_phase(vx, vc, vd)
                            if e2 is None: continue
                            eta = e1*e2
                            p = z4_round(eta)
                            if p is None: nonz4 += 1
                            res[lbl].append(eta)
    return res, nonz4

print("="*72)
print("D1/D2/D3: 代表親 v9=(2,1,0,0) での導出位相")
print("="*72)
v9 = np.array((2,1,0,0))
res, nonz4 = run_parent(v9)
for lbl in ('531','333'):
    etas = res[lbl]
    dist = defaultdict(int)
    for e in etas:
        p = z4_round(e)
        dist['i^%d'%p if p is not None else 'nonZ4'] += 1
    z = sum(etas)
    print(f"  {lbl}: 経路 {len(etas)} 本, 位相分布 {dict(dist)}")
    print(f"       z = {z:.4f}, |z|² = {abs(z)**2:.4f}")
print(f"  Z4 外の位相: {nonz4} 件 → {'D1 PASS (位相は Z4 に量子化)' if nonz4==0 else 'D1 FAIL'}")
W1 = abs(sum(res['531']))**2; W2 = abs(sum(res['333']))**2
if W1+W2>0:
    print(f"  分岐比 {W1/(W1+W2):.4f} : {W2/(W1+W2):.4f}")
    if W2>0:
        print(f"  双子二段 P(X) = {2*W1*W2/(2*W1*W2+W2*W2):.5f}")
strong = sum(1 for lbl in res for e in res[lbl] if abs(e-1)<1e-9)
print(f"  D3: 強い連続 (η=1) の経路数: {strong}")
print(f"  照合: B1=(1152,208)→0.8471 / B3=(2304,400)→0.8521 / 位置のみ管理連鎖 68本(48:20)")

print()
print("="*72)
print("D4: 等価性 — 48親セル全数")
print("="*72)
parents=[np.array(v) for v in itertools.product(range(-4,5),repeat=4)
         if sorted(map(abs,v))==[0,0,1,2] and abs(sum((abs(t)+0.5)**2 for t in v)-9)<1e-9]
vals=set(); counts=set(); dists=set()
for vp in parents:
    r,_ = run_parent(vp)
    z1=sum(r['531']); z2=sum(r['333'])
    vals.add((round(abs(z1)**2,6), round(abs(z2)**2,6)))
    counts.add((len(r['531']), len(r['333'])))
    d1=tuple(sorted(defaultdict(int, {z4_round(e):1 for e in []}).items()))  # placeholder
for vp in parents[:6]:
    r,_ = run_parent(vp)
    dd=defaultdict(int)
    for e in r['531']: dd[z4_round(e)]+=1
    dists.add(tuple(sorted(dd.items())))
print(f"  経路数 (531,333) の集合: {counts}")
print(f"  (|z531|²,|z333|²) の集合: {vals}")
print(f"  位相分布(531, 先頭6親): {dists}")
print(f"  → 親非依存性: {'PASS' if len(vals)==1 else 'FAIL'}")

# B4 移送チェック (位相多重集合の不変性)
print()
random.seed(11)
PERMS=list(itertools.permutations(range(4)))
SIGNS=list(itertools.product([1,-1],repeat=4))
def apply_g(perm,sign,v): return np.array([sign[i]*v[perm[i]] for i in range(4)])
base_r,_ = run_parent(np.array((2,1,0,0)))
base_ms = {lbl: tuple(sorted(z4_round(e) for e in base_r[lbl])) for lbl in base_r}
ok_g=True
for _ in range(3):
    g=(random.choice(PERMS), random.choice(SIGNS))
    vp2 = apply_g(*g, np.array((2,1,0,0)))
    r2,_ = run_parent(vp2)
    ms2 = {lbl: tuple(sorted(z4_round(e) for e in r2[lbl])) for lbl in r2}
    if ms2 != base_ms: ok_g=False
print(f"  B4 移送での位相多重集合の不変性 (3例): {'PASS' if ok_g else 'FAIL'}")
