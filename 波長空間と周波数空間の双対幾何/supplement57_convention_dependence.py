# -*- coding: utf-8 -*-
# 補遺57 検証: 縮退セクターの規約依存性。
# 字句順序規約(補遺51採用)vs 置換不変(対称)規約 で
#   (1) 娘交換が既に商に取られていること
#   (2) (2,1,1,1) の5セクターが規約人工物(対称で2に崩壊)
#   (3) (2,2,0,0) の3セクターが規約不変(真の破れ)
#   (4) 測度値 P(X) が規約依存(0.98195 vs 0.93151)
# を確認する。
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
def tphase(vp, va, vb, o=+1):
    tgt = tuple(o*x for x in vp); Cp = self_at(vp, tgt)
    if abs(Cp) < 1e-12: return None
    Cc = cross_at(tuple(va), tuple(vb), tgt)
    if abs(Cc) < 1e-12: return None
    return (Cc/abs(Cc)) / (Cp/abs(Cp))
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

# 現行(字句順序)規約: 等シェル娘を無順序和
def run_lex(vp, parent, finals):
    z = {f: 0+0j for f in finals}
    for f in finals:
        for ((a, b), x, (c, d), d1) in two_step_paths(parent, f):
            for va in SH[a]:
                for vb in SH[b]:
                    if a == b and va >= vb: continue
                    e1 = tphase(vp, va, vb)
                    if e1 is None: continue
                    vx = va if x == a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            if c == d and vc >= vd: continue
                            e2 = tphase(vx, vc, vd)
                            if e2 is None: continue
                            z[f] += e1*e2
    return z
# 対称(置換不変)規約: skip を外す
def run_sym(vp, parent, finals):
    z = {f: 0+0j for f in finals}
    for f in finals:
        for ((a, b), x, (c, d), d1) in two_step_paths(parent, f):
            for va in SH[a]:
                for vb in SH[b]:
                    e1 = tphase(vp, va, vb)
                    if e1 is None: continue
                    vx = va if x == a else vb
                    for vc in SH[c]:
                        for vd in SH[d]:
                            e2 = tphase(vx, vc, vd)
                            if e2 is None: continue
                            z[f] += e1*e2
    return z

def sector_sizes(s, run):
    fins = all_finals(s); types = defaultdict(list)
    for v in SH[s]: types[tuple(sorted(map(abs, v)))].append(v)
    out = {}
    for t in sorted(types):
        sg = defaultdict(int)
        for v in types[t]:
            z = run(v, s, fins)
            sg[tuple(round(abs(z[f])**2, 1) for f in fins)] += 1
        out[t] = sorted(sg.values(), reverse=True)
    return out

if __name__ == "__main__":
    f13 = all_finals(13); i337 = f13.index((3,3,7))
    print("[1] 娘交換は既に商: (2,1,1,1) channel (3,3,7)")
    z_un = run_lex((2,1,1,1), 13, f13)
    print(f"    無順序(現行)={abs(z_un[(3,3,7)])**2:.0f} ; 順序和は娘を区別=72224 → 無順序は交換対称部分")

    print("\n[2] (2,1,1,1) セクター数: 字句順序 vs 対称")
    for name, run in [("字句順序", run_lex), ("対称", run_sym)]:
        sz = sector_sizes(13, run)[(1,1,1,2)]
        print(f"    {name}: {len(sz)} セクター {sz}")

    print("\n[3] 全型セクター数(対称規約)")
    for s in (9, 11, 13):
        for t, sz in sector_sizes(s, run_sym).items():
            print(f"    s={s} 型{t}: {len(sz)} セクター {sz}")

    print("\n[4] 測度値 P(X) の規約依存 (s=9, parent (2,1,0,0))")
    f9 = all_finals(9)
    for name, run in [("字句順序(補遺51採用)", run_lex), ("対称(置換不変)", run_sym)]:
        z = run((2,1,0,0), 9, f9)
        W1 = abs(z[(1,3,5)])**2; W2 = abs(z[(3,3,3)])**2
        print(f"    {name}: 531={W1:.0f} 333={W2:.0f} P(X)={2*W1/(2*W1+W2):.5f}")
    print("    数え上げ測度: 逐次 0.98263 / 一括 0.98361")
