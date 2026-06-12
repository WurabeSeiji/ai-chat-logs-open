# -*- coding: utf-8 -*-
# 補遺56 検証: s=11,13 頑健性検査。補遺51/52/54 の機械をそのまま拡張(箱 K=5)。
# 実行結果(2026-06-12, claude.ai):
#   アンカー s=9 親(2,1,0,0): (1,1,7)->0, (1,3,5)->1088, (3,3,3)->40, P(X)=0.98195
#   TEST A (Z4): s=9 2816/0違反, s=11 6144/0, s=13 17744/0  → η は常に Z4
#   TEST C (セクター): s=11 型(2,1,1,0) 2セクター(48/48, 支配軸符号)
#                      s=13 型(3,0,0,0) 2, 型(2,2,0,0) 3, 型(2,1,1,1) 5(単一B4軌道)
#   → 「ちょうど2/1ビット」は支配軸一意のとき限定。縮退で破れる。
import itertools, numpy as np
from collections import defaultdict
SQ2 = np.sqrt(2.0)

_EXP = {}
def expcoeffs(v):
    v = tuple(v)
    if v in _EXP: return _EXP[v]
    acc = [((), 1.0+0j)]
    for i in range(4):
        m = v[i]
        if m == 0: opts = [(0, 1.0+0j)]
        elif m > 0: opts = [(+m, 1/SQ2+0j), (-m, 1/SQ2+0j)]
        else:
            mm = -m; opts = [(+mm, -1j/SQ2), (-mm, +1j/SQ2)]
        acc = [(n+(f,), c*w) for n, c in acc for f, w in opts]
    d = {}
    for n, c in acc: d[n] = d.get(n, 0) + c
    _EXP[v] = d; return d

def cross_at(va, vb, t):
    A = expcoeffs(va); B = expcoeffs(vb); tot = 0+0j
    for na, ca in A.items():
        nb = tuple(x-y for x, y in zip(t, na))
        if nb in B: tot += ca*B[nb]
    return tot
def self_at(v, t): return expcoeffs(tuple(v)).get(tuple(t), 0+0j)

def tphase(vp, va, vb, orient=+1):
    tgt = tuple(orient*x for x in vp)
    Cp = self_at(vp, tgt)
    if abs(Cp) < 1e-12: return None
    Cc = cross_at(tuple(va), tuple(vb), tgt)
    if abs(Cc) < 1e-12: return None
    return (Cc/abs(Cc)) / (Cp/abs(Cp))

def z4_idx(e):
    for p, u in enumerate([1, 1j, -1, -1j]):
        if abs(e-u) < 1e-9: return p
    return None

def shell_cells(m, K=5):
    return [tuple(k) for k in itertools.product(range(-K, K+1), repeat=4)
            if abs(sum((abs(t)+0.5)**2 for t in k) - m) < 1e-9]
SH = {m: shell_cells(float(m)) for m in (1,3,5,7,9,11,13)}

def two_step_paths(parent, final):
    paths = []; fin = tuple(sorted(final))
    for a in range(1, parent+2, 2):
        for b in range(a, parent+2, 2):
            d1 = a+b-parent
            if d1 not in (1, -1): continue
            for (x, spec) in ((a, b), (b, a)):
                for c in range(1, x+2, 2):
                    dd = x-d1-c
                    if dd < 1 or dd % 2 == 0: continue
                    if tuple(sorted((spec, c, dd))) == fin:
                        paths.append(((a, b), x, (c, dd), d1))
    return sorted(set(paths))

def all_finals(parent):
    fins = set()
    for a in range(1, parent+2, 2):
        for b in range(a, parent+2, 2):
            d1 = a+b-parent
            if d1 not in (1, -1): continue
            for (x, spec) in ((a, b), (b, a)):
                for c in range(1, x+2, 2):
                    dd = x-d1-c
                    if dd < 1 or dd % 2 == 0: continue
                    fins.add(tuple(sorted((spec, c, dd))))
    return sorted(fins)

def run_parent(vp, parent, finals, o1=+1, o2=+1):
    z = {f: 0+0j for f in finals}
    for f in finals:
        for ((a, b), x, (c, d), d1) in two_step_paths(parent, f):
            for va in SH[a]:
                for vb in SH[b]:
                    if a == b and va >= vb: continue
                    e1 = tphase(vp, va, vb, o1)
                    if e1 is None: continue
                    vx = va if x == a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            if c == d and vc >= vd: continue
                            e2 = tphase(vx, vc, vd, o2)
                            if e2 is None: continue
                            z[f] += e1*e2
    return z

if __name__ == "__main__":
    print("[ANCHOR] s=9 parent (2,1,0,0):")
    z9 = run_parent((2,1,0,0), 9, all_finals(9))
    for f in all_finals(9):
        print(f"  {f}: |z|^2={abs(z9[f])**2:.0f}")
    W1, W2 = abs(z9[(1,3,5)])**2, abs(z9[(3,3,3)])**2
    print(f"  P(X)=2W531/(2W531+W333)={2*W1/(2*W1+W2):.5f} (expect 0.98195)")

    print("\n[TEST A] Z4 quantization (single-step, exhaustive):")
    def daughter_pairs(s):
        return [(a, b) for a in range(1, s+2, 2) for b in range(a, s+2, 2) if a+b-s in (1, -1)]
    for s in (9, 11, 13):
        tot = viol = 0; dist = defaultdict(int)
        for vp in SH[s]:
            for (a, b) in daughter_pairs(s):
                for va in SH[a]:
                    for vb in SH[b]:
                        if a == b and va >= vb: continue
                        e = tphase(vp, va, vb, +1)
                        if e is None: continue
                        tot += 1; p = z4_idx(e)
                        if p is None: viol += 1
                        else: dist['i^%d' % p] += 1
        print(f"  s={s}: total {tot}, non-Z4 {viol}  dist {dict(dist)}")

    print("\n[TEST C] sector structure (all cells, all types):")
    def dom_sign(v):
        mx = max(abs(x) for x in v)
        ax = [j for j in range(4) if abs(v[j]) == mx]
        return None if len(ax) > 1 else (1 if v[ax[0]] > 0 else -1)
    for s in (11, 13):
        fins = all_finals(s)
        types = defaultdict(list)
        for v in SH[s]: types[tuple(sorted(map(abs, v)))].append(v)
        for t in sorted(types):
            sigs = defaultdict(list); bysign = defaultdict(set)
            for vp in types[t]:
                z = run_parent(vp, s, fins)
                sg = tuple(round(abs(z[f])**2, 1) for f in fins)
                sigs[sg].append(vp); bysign[dom_sign(vp)].add(sg)
            note = "/".join(f"{('+' if k>0 else '-') if k is not None else 'tie'}:{len(v)}sec"
                            for k, v in sorted(bysign.items(), key=lambda kv:(kv[0] is None, kv[0])))
            print(f"  s={s} type {t} ({len(types[t])} cells): {len(sigs)} sectors  [{note}]")
