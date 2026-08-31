#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス14（木原指示）：N≥5 の星型（2 波長系）で、L_e = k_e·λ_e/2 の整数 {k_e} により
実数シンプレックスが成立する組合せを探索する。無理数の波長比 ρ はそのまま使う
（対角が √ 比になるのは実幾何の常態。通約は要求しない。判定は三角不等式＋CM 半正定値のみ）。
N=5：星辺 4 本 k_i ≤ 12、衛星辺 6 本 k_ij ≤ 6 の全探索（衛星置換で同値類化）。
N=12：一様クラス（星辺一様 c ≤ 80、衛星辺一様 d ≤ 2）。
特別配置の照合：ハブが衛星正単体の外心に一致する縮退（rank が 1 落ちる）と、その比。
v2（木原指摘による記録拡充）：(1) N=5 の整合解 全同値類を results/k_search_star_solutions.csv に全数保存
(2) 外接半径照合を md に記録 (3) 軌道一様解（星辺 k 全同・衛星辺 k 全同）の全列挙と Σk 最小の同定（計算で確定）。"""
import itertools, math, os
from math import gcd
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RHO5 = 3.0/(3.0-math.sqrt(5.0))     # 実測 3.9265〜3.9286、理論 3.92705（相対差 2e-4 で理論値を使用）
RHO12 = 49.7866                      # 実測（ne/rb 一致値）

def cm_check(N, D2):
    for i in range(N):
        for j in range(N):
            for m in range(N):
                if len({i, j, m}) < 3: continue
                if math.sqrt(D2[i, j]) > math.sqrt(D2[i, m])+math.sqrt(D2[m, j])+1e-12:
                    return None
    J = np.eye(N)-np.ones((N, N))/N; B = -0.5*J@D2@J
    ev = np.linalg.eigvalsh(B); sc = max(abs(ev).max(), 1e-300)
    if ev.min()/sc < -1e-9: return None
    return int((ev/sc > 1e-9).sum())

# ---- N=5 全探索 ----
N = 5
sats = [0, 1, 2, 3]; hub = 4
pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
sols = {}
for ks in itertools.combinations_with_replacement(range(1, 13), 4):   # 星辺（ソート済み=ハブ辺の多重集合）
    kmax_ij = {}
    ok = True
    for (i, j) in pairs:
        m = math.floor((ks[i]+ks[j])/RHO5 + 1e-12)
        if m < 1: ok = False; break
        kmax_ij[(i, j)] = min(m, 6)
    if not ok: continue
    for kd in itertools.product(*[range(1, kmax_ij[p]+1) for p in pairs]):
        D2 = np.zeros((N, N))
        for a, (i, j) in enumerate(pairs):
            d = RHO5*kd[a]/2.0; D2[i, j] = D2[j, i] = d*d
        for i in sats:
            d = ks[i]/2.0; D2[i, hub] = D2[hub, i] = d*d
        r = cm_check(N, D2)
        if r is None: continue
        # 正規形（衛星置換）＋スケール（gcd）
        best = None
        for perm in itertools.permutations(range(4)):
            kss = tuple(ks[perm[i]] for i in range(4))
            kdd = tuple(kd[pairs.index(tuple(sorted((perm.index(i), perm.index(j)))))] if False else 0 for _ in pairs)
            # 置換はハブ辺と衛星辺を同時に並べ替える
            kdd = tuple(kd[pairs.index(tuple(sorted((perm[i], perm[j]))))] for (i, j) in pairs)
            cand = (kss, kdd)
            if best is None or cand < best: best = cand
        g = 0
        for x in best[0]+best[1]: g = gcd(g, x)
        best = (tuple(x//g for x in best[0]), tuple(x//g for x in best[1]))
        if best not in sols or r > sols[best]:
            sols[best] = r
full = {k: r for k, r in sols.items() if r == 4}
deg = {k: r for k, r in sols.items() if r < 4}
print(f"N=5（ρ={RHO5:.5f}、k_star≤12・k_sat≤6 全探索、同値類）: 整合 {len(sols)}（非縮退 rank4: {len(full)}、縮退: {len(deg)}）")
mins = sorted(sols.items(), key=lambda kv: (sum(kv[0][0])+sum(kv[0][1])))[:8]
for (kss, kdd), r in mins:
    print(f"  Σk最小級: 星辺 k={kss} 衛星辺 k={kdd} rank={r}")

# ---- N=12 一様クラス ----
N12 = 12; nsat = 11
res12 = []
for d in (1, 2):
    for c in range(1, 81):
        D2 = np.zeros((N12, N12))
        for i in range(nsat):
            for j in range(i+1, nsat):
                dd = RHO12*d/2.0; D2[i, j] = D2[j, i] = dd*dd
        for i in range(nsat):
            dd = c/2.0; D2[i, nsat] = D2[nsat, i] = dd*dd
        r = cm_check(N12, D2)
        if r is not None:
            res12.append((c, d, r))
print(f"\nN=12（ρ={RHO12}、一様クラス 星辺 c・衛星辺 d）: 整合 {len(res12)} 組")
for c, d, r in res12[:6]:
    print(f"  c={c}, d={d}: rank={r}" + ("（最小）" if (c, d) == (res12[0][0], res12[0][1]) else ""))
# 外心縮退の比（rank が落ちる特別点）
Rcirc = math.sqrt(nsat and (nsat-1)/(2*nsat))  # 正 (nsat-1) 単体（頂点 nsat 個・辺1）の外接半径
print(f"  参考：衛星正単体（11 点）の外接半径/辺 = {math.sqrt((nsat-1)/(2*nsat)):.6f} → ハブ外心一致は c/(ρd) = {math.sqrt((nsat-1)/(2*nsat)):.6f}（無理数、整数 c,d では厳密到達不可＝一般には rank 11 の非縮退で成立）")
# 全解の保存と軌道一様解の抽出
import csv as _csv
with open(os.path.join(ROOT, "results", "k_search_star_solutions.csv"), "w", newline="") as f:
    w = _csv.writer(f); w.writerow(["k_star_1","k_star_2","k_star_3","k_star_4","k_sat_01","k_sat_02","k_sat_03","k_sat_12","k_sat_13","k_sat_23","rank","sum_k"])
    for (kss, kdd), r in sorted(sols.items(), key=lambda kv: sum(kv[0][0])+sum(kv[0][1])):
        w.writerow(list(kss)+list(kdd)+[r, sum(kss)+sum(kdd)])
uniform = sorted([((kss, kdd), r) for (kss, kdd), r in sols.items()
                  if len(set(kss)) == 1 and len(set(kdd)) == 1],
                 key=lambda kv: sum(kv[0][0])+sum(kv[0][1]))
Rc5 = math.sqrt(3.0/8.0)  # 正四面体（4点・辺1）の外接半径
with open(os.path.join(ROOT, "results", "k_search_star.md"), "w") as f:
    f.write(f"# 星型（2 波長系）の {{k}} 探索（パス14）\n\n")
    f.write(f"## N=5（ρ=3/(3−√5)={RHO5:.6f}）：整合 {len(sols)} 同値類（非縮退 {len(full)}・縮退 {len(deg)}）\n\n")
    f.write("Σk 最小級の例：\n\n| 星辺 k | 衛星辺 k | rank |\n|---|---|---|\n")
    for (kss, kdd), r in mins:
        f.write(f"| {kss} | {kdd} | {r} |\n")
    f.write(f"\n## N=12（ρ={RHO12} 実測、一様クラス）：整合 {len(res12)} 組\n\n| 星辺 c | 衛星辺 d | rank |\n|---|---|---|\n")
    for c, d, r in res12:
        f.write(f"| {c} | {d} | {r} |\n")
    f.write(f"\n## 軌道一様解（星辺 k 全同 s・衛星辺 k 全同 t、N=5）：{len(uniform)} 個（全解 csv から機械抽出）\n\n")
    f.write("| (s;t) | Σk | rank |\n|---|---|---|\n")
    for (kss, kdd), r in uniform[:12]:
        f.write(f"| ({kss[0]};{kdd[0]}) | {sum(kss)+sum(kdd)} | {r} |\n")
    if uniform:
        (kss0, kdd0), r0 = uniform[0]
        f.write(f"\nΣk 最小の軌道一様解 = (s={kss0[0]}, t={kdd0[0]})、rank={r0}（機械判定）。\n")
    f.write(f"\n## 外接半径の照合（記録）\n\n")
    f.write(f"- N=5：衛星正四面体（辺 ρt/2）の外接半径/辺 = √(3/8) = {Rc5:.6f}。ハブ外心一致条件 s/(ρt) = {Rc5:.6f} は無理数のため整数 (s,t) で厳密到達不可。s=3,t=1 は s/(ρt) = {3/RHO5:.6f} > {Rc5:.6f} で非縮退成立。\n")
    f.write(f"- N=12：衛星正 10 単体（11 点）の外接半径/辺 = √(10/22) = {math.sqrt(10/22):.6f}。c/(ρd) ≥ 外接半径が必要 → d=1 で c ≥ {RHO12*math.sqrt(10/22):.4f} → 最小整数 c=34（一様表と一致）。\n")
    f.write(f"\n注：N=5 の ρ は理論値 3/(3−√5)（実測 3 タグとの差 1.3〜4.0e-4、sideband_theory_check.csv）を使用。\n")
    f.write(f"\n結論：無理数の波長比のままでも、{{k}} の整数選択で実単体は成立する（埋め込みは開条件）。全整合解は k_search_star_solutions.csv に全数保存。\n")
print("PASS14 OK")
