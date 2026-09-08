#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三角数分解の全列挙（2026-09-08、木原指示）。M(N)=N(N-1)/2 が他の M の和に分解できるパターン。

(A) 等分分解 M(N)=m×M(k)（K_N を K_k コピーへ辺分割する必要条件も判定: (k-1)|(N-1) と Fisher m>=N）
(B) 2部分分解 M(N)=M(a)+M(b)（混合含む）
(C) 3部分分解 M(N)=M(a)+M(b)+M(c)（部品 M(4) 以上・全同一除外・N<=24）
(D) 部品 M(4) 以上でどの分解も持たない「素」な N（N<=16）
純粋な数え上げ（決定論・乱数なし）。実験21（N9=N7+N6 混合和の加法法則検証）の選定根拠。
"""
from itertools import combinations_with_replacement

def M(n):
    return n * (n - 1) // 2

NMAX = 40
ks = list(range(3, NMAX + 1))
tri = {k: M(k) for k in ks}

out = []
out.append('=== (A) 等分分解 M(N) = m × M(k)  (m>=2, k>=3, N<=40) ===')
out.append('  [辺分割の必要条件] (k-1)|(N-1) かつ Fisher m>=N')
for n in ks:
    for k in range(3, n):
        if tri[n] % tri[k] == 0:
            m = tri[n] // tri[k]
            if m < 2:
                continue
            div_ok = (n - 1) % (k - 1) == 0
            fisher_ok = m >= n
            if div_ok and fisher_ok:
                tag = '  <= 辺分割の必要条件 全充足'
            elif div_ok:
                tag = '  (次数条件OK, Fisher m>=N 違反 → 分割不能)'
            else:
                tag = '  (次数条件 (k-1)|(N-1) 違反 → 分割不能)'
            out.append(f'N={n:2d} M={tri[n]:3d} = {m:3d} × M({k})={tri[k]:3d}{tag}')

out.append('')
out.append('=== (B) 2部分分解 M(N) = M(a) + M(b)  (a>=b>=3, N<=40) — 混合含む ===')
c2 = 0
for n in ks:
    for a, b in combinations_with_replacement(ks, 2):
        if tri[a] + tri[b] == tri[n]:
            c2 += 1
            mark = ' *' if a != b else ''
            out.append(f'N={n:2d} M={tri[n]:3d} = M({a})+M({b}) = {tri[a]}+{tri[b]}{mark}')
out.append(f'計 {c2} 件（* = 異種混合）')

out.append('')
out.append('=== (C) 3部分分解 M(N)=M(a)+M(b)+M(c)  (a>=b>=c>=4; 全同一除外) — N<=24 ===')
c3 = 0
for n in ks:
    if n > 24:
        break
    for combo in combinations_with_replacement([k for k in ks if k >= 4], 3):
        a, b, c = sorted(combo, reverse=True)
        if len(set(combo)) == 1:
            continue
        if tri[a] + tri[b] + tri[c] == tri[n]:
            c3 += 1
            out.append(f'N={n:2d} M={tri[n]:3d} = M({a})+M({b})+M({c}) = {tri[a]}+{tri[b]}+{tri[c]}')
out.append(f'計 {c3} 件（部品 M(4)以上・全同一除外・N<=24）')

out.append('')
out.append('=== (D) 部品 M(4) 以上でどの分解も持たない N（N<=16）===')
for n in ks:
    if n > 16:
        break
    has = False
    for a, b in combinations_with_replacement([k for k in ks if k >= 4], 2):
        if tri[a] + tri[b] == tri[n]:
            has = True
    for combo in combinations_with_replacement([k for k in ks if k >= 4], 3):
        if sum(tri[x] for x in combo) == tri[n]:
            has = True
    if not has:
        out.append(f'N={n} M={tri[n]}: 分解なし（素）')

text = '\n'.join(out) + '\n'
import os
BASE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(BASE, '三角数分解一覧_N3_N40.txt'), 'w', encoding='utf-8') as f:
    f.write(text)
print(text)
print('ALL DONE', flush=True)
