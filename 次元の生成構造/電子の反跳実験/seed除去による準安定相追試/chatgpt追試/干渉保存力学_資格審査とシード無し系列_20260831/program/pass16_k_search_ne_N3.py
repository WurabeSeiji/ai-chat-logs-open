#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス16：ne_N3（3 波長系）の {k} 探索。実測の波長比（wavelength_table.csv）をそのまま使い、
k_e ∈ {1..12} の全探索で三角不等式＋CM を判定。全整合解を保存。"""
import os, csv, math, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
wt = [r for r in csv.DictReader(open(os.path.join(ROOT, "results", "wavelength_table.csv"))) if r["tag"] == "ne_N3"]
lam = {r["edge"]: float(r["lam_norm"]) for r in wt}
E = [(0, 1), (0, 2), (1, 2)]
lams = [lam[str(e)] for e in E]
sols = []
for k in itertools.product(range(1, 13), repeat=3):
    L = [k[a]*lams[a]/2.0 for a in range(3)]
    a, b, c = sorted(L)
    if a+b < c-1e-12: continue
    s = (a+b+c)/2
    area2 = s*(s-a)*(s-b)*(s-c)
    rank = 2 if area2 > 1e-9*(max(L)**4) else 1
    sols.append((k, rank, sum(k)))
sols.sort(key=lambda x: x[2])
with open(os.path.join(ROOT, "results", "k_search_ne_N3.md"), "w") as f:
    f.write("# ne_N3（3 波長系）の {k} 全探索（パス16）\n\n")
    f.write(f"実測波長（最短=1）：辺(0,1)={lams[0]:.4f}、辺(0,2)={lams[1]:.4f}、辺(1,2)={lams[2]:.4f}（wavelength_table.csv）\n\n")
    f.write(f"k ≤ 12 全 1728 組中、整合 {len(sols)} 組（非縮退 rank2: {sum(1 for s in sols if s[1]==2)}）\n\n")
    f.write("| k=(k01,k02,k12) | Σk | rank | 辺長 L=(k·λ/2) |\n|---|---|---|---|\n")
    for k, r, sk in sols[:20]:
        L = [k[a]*lams[a]/2 for a in range(3)]
        f.write(f"| {k} | {sk} | {r} | ({L[0]:.3f}, {L[1]:.3f}, {L[2]:.3f}) |\n")
with open(os.path.join(ROOT, "results", "k_search_ne_N3.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["k01", "k02", "k12", "rank", "sum_k"])
    for k, r, sk in sols: w.writerow(list(k)+[r, sk])
print(f"ne_N3: λ=({lams[0]:.4f},{lams[1]:.4f},{lams[2]:.4f})  整合 {len(sols)}/1728（非縮退 {sum(1 for s in sols if s[1]==2)}）")
print("Σk 最小の解:", sols[:5])
print("PASS16 OK")
