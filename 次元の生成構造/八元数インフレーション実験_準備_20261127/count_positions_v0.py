#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ダンプを読み、全部の行・全部の塊について、向きが取る違う位置の数を数える。2026-11-27。

図は一つの行しか描けないので、「行 2 で見えたことが、ほかの行でも同じか」を数で確かめるためのもの。
数え方は図化プログラムと同じ: n >= 3 の回について、塊の θn (虚部の全成分。ダンプに書かれた 10 桁) を並べ、違うものの数を数える。
向きが決まらない回 (全成分 0、または向き不定) は数えない。
出すもの (標準出力):
  - 違う位置の数ごとの、(行, 塊) の個数
  - 行ごとに「位置の数がいちばん多い塊」がどれか。その塊が、行の番号と同じ番号の成分を含む塊と、その反転の相手
    (成分の番号が N/2 だけずれた塊) の二つになっているかどうか
使い方: python3 count_positions_v0.py ダンプのファイル
"""
import collections
import re
import sys

lines = open(sys.argv[1], encoding="utf-8").read().split("\n")
m0 = re.match(r"^N = (\d+)。.*。回数 (\d+)。(?:読み方 1(\d+) \(.*\)。)?実部と読む成分: (\d+) 番目", lines[0])
N, N_MAX = int(m0.group(1)), int(m0.group(2))
BLK = (int(m0.group(3)) + 1) if m0.group(3) else N
NB = N // BLK
assert NB >= 2, "塊に区切って読むダンプだけを対象にする"
pos = collections.defaultdict(set)          # (行, 塊) -> 違う位置の集まり
cnt = collections.Counter()                 # (行, 塊) -> 向きが決まった回の数
n = None
k = -1
for line in lines:
    mh = re.match(r"^n = (\d+)    \(H を", line)
    if mh:
        n = int(mh.group(1))
        k = -1
        continue
    if n is None or n < 3:
        continue
    if line.startswith("[行の指数部"):
        k = 0
        continue
    if k >= 0 and k < N * NB and line.startswith("  行"):
        i, q = divmod(k, NB)
        k += 1
        mm = re.search(r"θn = \((.*)\)$", line)
        if mm:
            pos[(i, q)].add(mm.group(1))
            cnt[(i, q)] += 1
print("ダンプ: %s" % sys.argv[1])
print(lines[0])
print("行 %d 本 × 塊 %d 個 = %d 個。n = 3〜%d の %d 回" % (N, NB, N * NB, N_MAX, N_MAX - 2))
dist = collections.Counter(len(pos[(i, q)]) for i in range(N) for q in range(NB))
print("違う位置の数ごとの個数:")
for v, c in sorted(dist.items()):
    print("  違う位置 %2d 個: %d 個" % (v, c))
und = sum(1 for i in range(N) for q in range(NB) if cnt[(i, q)] < N_MAX - 2)
print("向きが決まらない回を持つもの: %d 個" % und)
top = max(dist)
ok = 0
bad = []
for i in range(N):
    most = sorted(q for q in range(NB) if len(pos[(i, q)]) == top)
    expect = sorted({i // BLK, ((i + N // 2) % N) // BLK})
    if most == expect:
        ok += 1
    else:
        bad.append((i + 1, [q + 1 for q in most], [q + 1 for q in expect]))
print("行ごとに、位置の数が最大 (%d 個) の塊が「行の番号の成分を含む塊と、その反転の相手の塊」の二つか: %d 行で成り立つ / %d 行" % (top, ok, N))
for r, m_, e_ in bad[:10]:
    print("  成り立たない行 %d: 最大の塊 %s、見込み %s" % (r, m_, e_))
