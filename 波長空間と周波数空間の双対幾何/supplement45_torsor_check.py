#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺45 §3: 残余 Z2×Z2 の torsor 検査
# 48 親セル × 4 解への「反転 I: (v5,v3)→(-v5,-v3)」×「娘符号 D: (v5,v3)→(v5,-v3)」の
# 作用が単純推移的 (Z2×Z2 torsor) かを全数検査する。
import itertools, numpy as np

def shell_cells(m):
    out=[]; K=4
    for k in itertools.product(range(-K,K+1),repeat=4):
        if abs(sum((abs(t)+0.5)**2 for t in k)-m)<1e-9: out.append(tuple(k))
    return out
SH5 = shell_cells(5.0); SH3 = shell_cells(3.0); SH9 = shell_cells(9.0)
parents = [v for v in SH9 if sorted(map(abs,v))==[0,0,1,2]]
print(f"(2,1,0,0)型 親セル: {len(parents)} 個 / shell5: {len(SH5)} / shell3: {len(SH3)}")

def solutions(v9):
    v9 = np.array(v9); sols=[]
    for v5 in SH5:
        for v3 in SH3:
            a5,a3 = np.array(v5), np.array(v3)
            if any(np.array_equal(x,v9) or np.array_equal(-x,v9) for x in (a5+a3, a5-a3)):
                sols.append((v5,v3))
    return sols

def I_act(s): return (tuple(-x for x in s[0]), tuple(-x for x in s[1]))
def D_act(s): return (s[0], tuple(-x for x in s[1]))

ok_all=True; sizes=set()
for v9 in parents:
    sols = solutions(v9); S=set(sols); sizes.add(len(sols))
    # 群作用の閉性
    closed = all(I_act(s) in S and D_act(s) in S for s in sols)
    # 可換性 (定義上自明だが念のため)
    comm = all(I_act(D_act(s))==D_act(I_act(s)) for s in sols)
    # 固定点なし
    free = all(I_act(s)!=s and D_act(s)!=s and I_act(D_act(s))!=s for s in sols)
    # 単純推移性: 1点の軌道が全体
    s0 = sols[0]
    orbit = {s0, I_act(s0), D_act(s0), I_act(D_act(s0))}
    transitive = (orbit==S)
    if not (len(sols)==4 and closed and comm and free and transitive):
        ok_all=False
        print(f"  NG: v9={v9} 解{len(sols)} closed={closed} free={free} trans={transitive}")
print(f"解の個数の集合: {sizes}")
print(f"全48セル: {'PASS — 4解は Z2×Z2 (反転×娘符号) の単純推移的軌道 (torsor)' if ok_all else 'FAIL'}")

# 追加検査: (1,1,1,1)型 (解0のはず) と、他の自己同型候補 (v5側符号 D5) での同値性
parents2 = [v for v in SH9 if sorted(map(abs,v))==[1,1,1,1]]
n0 = sum(1 for v9 in parents2 if len(solutions(v9))==0)
print(f"(1,1,1,1)型 {len(parents2)} セル: 解0 のセル数 = {n0} (全0なら補遺38と整合)")

def D5_act(s): return (tuple(-x for x in s[0]), s[1])
v9=parents[0]; S=set(solutions(v9))
same = all(D5_act(s) in S for s in S)
# D5 = I∘D なので独立な第3生成元ではないことの確認
dep = all(D5_act(s)==I_act(D_act(s)) for s in S)
print(f"v5側符号 D5 の閉性: {same} / D5=I∘D (独立でない): {dep}")
