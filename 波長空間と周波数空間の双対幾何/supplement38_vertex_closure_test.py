#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 頂点閉鎖テスト: 親の線ベクトル v9 が娘系の交差スペクトルに存続する条件
# (記録連続則の最小形)。全数検査。
import itertools, numpy as np
from collections import defaultdict

def shell_cells(m):
    out=[]
    for k in itertools.product(range(-3,4),repeat=4):
        if sum((abs(t)+0.5)**2 for t in k)==m: out.append(np.array(k))
    return out
S9=shell_cells(9.0); S5=shell_cells(5.0); S3=shell_cells(3.0)
print(f"shell9={len(S9)}, shell5={len(S5)}, shell3={len(S3)}")

def lines_of_pair(a,b):
    return [tuple(a+b), tuple(a-b), tuple(b-a), tuple(-(a+b))]

def persists_531(v9):
    """(5,3,1): 娘 {v5, v3, 0}。交差線: v5±v3, v5±0, v3±0。v9 に一致する配置を数える"""
    sols=[]
    t9=tuple(v9)
    for v5 in S5:
        for v3 in S3:
            cand=set()
            for x in (v5+v3, v5-v3): cand.add(tuple(x)); cand.add(tuple(-x))
            # 線形転写線 (原点との交差): v5, v3 自身 — ノルムが違うので v9 には届かないが念のため
            for x in (v5, v3): cand.add(tuple(x)); cand.add(tuple(-x))
            if t9 in cand: sols.append((tuple(v5),tuple(v3)))
    return sols

def persists_333(v9):
    sols=[]
    t9=tuple(v9); tm9=tuple(-v9)
    for trio in itertools.combinations(range(len(S3)),3):
        cells=[S3[i] for i in trio]
        cand=set()
        for i,j in itertools.combinations(range(3),2):
            for x in (cells[i]+cells[j], cells[i]-cells[j]):
                cand.add(tuple(x)); cand.add(tuple(-x))
        if t9 in cand or tm9 in cand: sols.append(tuple(map(tuple,cells)))
    return sols

def persists_55(v9):
    """仮想チャネル (5,5) — s=10=9+1、エネルギー1単位超過"""
    sols=[]
    t9=tuple(v9); tm9=tuple(-v9)
    for i,j in itertools.combinations(range(len(S5)),2):
        a,b=S5[i],S5[j]
        cand=set()
        for x in (a+b, a-b): cand.add(tuple(x)); cand.add(tuple(-x))
        if t9 in cand or tm9 in cand: sols.append((tuple(a),tuple(b)))
    return sols

typ=lambda v: tuple(sorted(np.abs(v)))
res=defaultdict(lambda: [0,0,0,0])   # orbit -> [n_cells, sum531, sum333, sum55]
examples={}
for v9 in S9:
    t=typ(v9)
    s531=persists_531(v9); s333=persists_333(v9); s55=persists_55(v9)
    r=res[t]; r[0]+=1; r[1]+=len(s531); r[2]+=len(s333); r[3]+=len(s55)
    if t not in examples: examples[t]=(tuple(v9), s531[:4], s55[:4])

print()
print("親軌道型ごとの存続解の総数 (全セル合計):")
print("  軌道型          セル数  (5,3,1)解  (3,3,3)解  [仮想(5,5)解]")
for t,r in sorted(res.items()):
    print(f"  {t}  {r[0]:4d}   {r[1]:6d}    {r[2]:6d}    {r[3]:8d}")
for t,(v9ex,s531ex,s55ex) in examples.items():
    print(f"  例 v9={v9ex} (型{t}):")
    if s531ex:
        print(f"    (5,3,1) 解: {s531ex}  → ベクトル保存 v5±v3=v9 (運動量様)")
    else:
        print(f"    (5,3,1) 解: なし")
        if s55ex: print(f"    仮想(5,5) 解の例: {s55ex}  → s=10 が必要 (9+1: 1量子不足)")
print()
print("支持の数え上げ: 娘対の非零成分の和は最大 2+1=3 軸 < 4 軸")
print("→ (1,1,1,1) 型 (4軸) の親線は許容チャネルでは構成不能 — 不足はちょうど1軸/1量子")
