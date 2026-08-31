#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス13（木原の方法論の実装）：k は導出しない。当てはまる {k_e} を数え上げる。
等波長状態（λ_e 全て等しい）の場合、L_e = k_e·λ/2 なので長さ=整数 k_e。
大域整合条件＝三角不等式（縮退込み判定）＋ Cayley–Menger 半正定値。
N=3（M=3）と N=4（M=6）について k_e ∈ {1..K_MAX} の全組合せを検定し、
全体スケール（k 全体の公倍）と頂点ラベルの置換で同値類にまとめて列挙する。"""
import itertools, math, os, csv
from math import gcd
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
K_MAX = 4

def edges(N): return [(i, j) for i in range(N) for j in range(i+1, N)]

def admissible(N, E, k, strict=True):
    D = np.zeros((N, N))
    for ke, (i, j) in zip(k, E): D[i, j] = D[j, i] = float(ke)**2
    for i in range(N):
        for j in range(N):
            for m in range(N):
                if len({i, j, m}) < 3: continue
                if math.sqrt(D[i, j]) > math.sqrt(D[i, m])+math.sqrt(D[m, j])+1e-12:
                    return None
    J = np.eye(N)-np.ones((N, N))/N; B = -0.5*J@D@J
    ev = np.linalg.eigvalsh(B); scale = max(abs(ev).max(), 1e-300)
    if ev.min()/scale < -1e-9: return None
    rank = int((ev/scale > 1e-9).sum())
    return rank

def canonical(N, E, k):
    """頂点置換とスケール（gcd で割る）で正規形に。"""
    g = 0
    for x in k: g = gcd(g, x)
    k = tuple(x//g for x in k)
    best = None
    idx = {e: a for a, e in enumerate(E)}
    for perm in itertools.permutations(range(N)):
        kk = tuple(k[idx[tuple(sorted((perm[i], perm[j])))]] for (i, j) in E)
        if best is None or kk < best: best = kk
    return best

rows = []
for N in (3, 4):
    E = edges(N); M = len(E)
    seen = {}
    for k in itertools.product(range(1, K_MAX+1), repeat=M):
        c = canonical(N, E, k)
        if c in seen: continue
        r = admissible(N, E, c)
        seen[c] = r
    ok = {c: r for c, r in seen.items() if r is not None}
    full = {c: r for c, r in ok.items() if r == N-1}
    degen = {c: r for c, r in ok.items() if r < N-1}
    print(f"N={N} (M={M}, k≤{K_MAX}): 同値類 {len(seen)} 中、整合 {len(ok)}（非縮退 rank={N-1}: {len(full)}、縮退: {len(degen)}）")
    for c, r in sorted(ok.items()):
        rows.append(dict(N=N, k_pattern=str(c), rank=r, degenerate=(r < N-1)))
md = ["# {k_e} の数え上げ（パス13）：等波長状態で大域整合する整数長パターンの全列挙", "",
      f"方法：k_e ∈ {{1..{K_MAX}}}、頂点置換と全体スケールで同値類化、三角不等式＋CM 半正定値で判定。", ""]
for N in (3, 4):
    sel = [r for r in rows if r["N"] == N]
    md.append(f"## N={N}：整合パターン {len(sel)} 種")
    md.append("")
    md.append("| {k_e}（正規形） | rank | 縮退 |")
    md.append("|---|---|---|")
    for r in sel:
        md.append(f"| {r['k_pattern']} | {r['rank']} | {'縮退' if r['degenerate'] else '—'} |")
    md.append("")
with open(os.path.join(ROOT, "results", "k_enumeration.md"), "w") as f:
    f.write("\n".join(md)+"\n")
with open(os.path.join(ROOT, "results", "k_enumeration.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["N", "k_pattern", "rank", "degenerate"]); w.writeheader(); w.writerows(rows)
print("PASS13 OK")
