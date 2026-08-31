#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス17（木原質問「hm 系列では全て成立したのか」への全数回答）：
hm_N3〜N16 の全 14 終状態について、辺別波長を測定し、{k} の整数選択で実数シンプレックスが
成立するかを判定する。方法：
  1. 終窓 16384 step の辺別支配振動数（整合フィルタ）→ λ_e = 1/|ν_e|（最短=1 規格化）
  2. λ を 0.5% クラスタで群化（報告用）。{k} 探索には 2% クラスタを使う（未収束状態では
     星辺群自体が ~1% 幅に割れ、0.5% 群化だとアンザッツが構造的に合わないため。hm_N6 で判明）
  3. {k} 探索：群一様アンザッツ（同じ λ 群は同じ k）。単一群なら k≡1（等長正単体）で自明成立。
     複数群なら 2 段：まず全衛星群に共通 d ∈ {1,2,3}（高速経路。衛星 λ は互いに数 % 以内なので自然）、
     見つからなければ群別 k_g ∈ {1..3} の全組合せ。c は下限 floor(0.5·max(k_g λ_g)) から 100 まで。
  4. 記録：λ 群構成・成立可否・最小解・rank。アンザッツ内で解が無い場合は「群一様の範囲では不成立」と報告。"""
import os, csv, math, itertools
import numpy as np
from scipy.optimize import minimize_scalar
import sys
sys.path.insert(0, os.path.dirname(__file__))
from common import edges
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
LW = 16384; PAD = 8

def dominant(sig):
    n = len(sig); h = np.hanning(n)
    F = np.fft.fft(sig*h, n*PAD); P = np.abs(F)**2
    fr = np.fft.fftfreq(n*PAD)*2*math.pi
    om0 = fr[int(np.argmax(P))]
    t = np.arange(n); bw = 2*math.pi/n
    f = lambda om: -abs(np.sum(sig*h*np.exp(-1j*om*t)))
    return float(minimize_scalar(f, bounds=(om0-bw, om0+bw), method="bounded",
                                 options=dict(xatol=1e-9)).x)

def cm_ok(N, D2):
    Dm = np.sqrt(D2)
    # 三角不等式（ベクトル化）：Dm[i,j] <= min_k (Dm[i,k]+Dm[k,j])。k=i,j は等号になり無害
    mn = (Dm[:, None, :] + Dm[None, :, :]).min(axis=2)
    if (Dm > mn + 1e-12).any(): return None
    J = np.eye(N)-np.ones((N, N))/N; B = -0.5*J@D2@J
    ev = np.linalg.eigvalsh(B); sc = max(abs(ev).max(), 1e-300)
    if ev.min()/sc < -1e-9: return None
    return int((ev/sc > 1e-9).sum())

md = ["# hm 系列全 14 状態の {k} 成立判定（パス17）", "",
      "| N | λ 群（本数） | 判定 | 最小解（群 k） | rank |", "|---|---|---|---|---|"]
rows = []
for N in range(3, 17):
    tag = f"hm_N{N}"; E = edges(N); M = len(E)
    Z = np.load(os.path.join(DATA, tag, "states_treatment.npz"))["Z"][-LW:]
    nus = [dominant(Z[:, e]) for e in range(M)]
    lam = [1.0/abs(v) for v in nus]
    lmin = min(lam); lamn = [x/lmin for x in lam]
    groups = []
    for e in range(M):
        for g in groups:
            if abs(lamn[e]-g[0])/g[0] < 0.005:
                g[1].append(e); g[0] = (g[0]*(len(g[1])-1)+lamn[e])/len(g[1]); break
        else:
            groups.append([lamn[e], [e]])
    groups.sort()
    gdesc = " / ".join(f"{g[0]:.3f}×{len(g[1])}" for g in groups)
    # 探索用の粗い群化（2%）
    sgroups = []
    for e in range(M):
        for g in sgroups:
            if abs(lamn[e]-g[0])/g[0] < 0.02:
                g[1].append(e); g[0] = (g[0]*(len(g[1])-1)+lamn[e])/len(g[1]); break
        else:
            sgroups.append([lamn[e], [e]])
    sgroups.sort()
    groups = sgroups
    if len(groups) == 1:
        verdict, sol, rank = "成立（k≡1、等長正単体）", "k≡1", N-1
    else:
        idx = {}
        for gi, g in enumerate(groups):
            for e in g[1]: idx[e] = gi
        found = None
        lam_arr = np.array([groups[idx[e]][0] for e in range(M)])
        def try_kg(kg):
            L = np.array([kg[idx[e]]*groups[idx[e]][0]/2.0 for e in range(M)])
            D2 = np.zeros((N, N))
            for e, (i, j) in enumerate(E):
                D2[i, j] = D2[j, i] = L[e]*L[e]
            return cm_ok(N, D2)
        # 高速経路：全衛星群 共通 d
        for d in (1, 2, 3):
            cmin = max(1, int(0.5*max(g[0] for g in groups[1:])*d))
            for c in range(cmin, 101):
                kg = (c,)+(d,)*(len(groups)-1)
                r = try_kg(kg)
                if r is not None:
                    found = (kg, r); break
            if found: break
        # 代替経路：群別 k_g（群数が多い場合は共通 d で見つからなかったときのみ）
        if not found:
            for c in range(1, 101):
                for ks in itertools.product(range(1, 4), repeat=len(groups)-1):
                    r = try_kg((c,)+ks)
                    if r is not None:
                        found = ((c,)+ks, r); break
                if found: break
        if found:
            verdict, sol, rank = "成立（群一様）", str(found[0]), found[1]
        else:
            verdict, sol, rank = "群一様の範囲では不成立", "—", None
    md.append(f"| {N} | {gdesc} | {verdict} | {sol} | {rank if rank else '—'} |")
    rows.append(dict(N=N, groups=gdesc, verdict=verdict, min_solution=sol, rank=rank))
    print(f"hm_N{N}: 群[{gdesc}] → {verdict} {sol} rank={rank}")
with open(os.path.join(ROOT, "results", "hm_series_k.md"), "w") as f:
    f.write("\n".join(md)+"\n\n注：群一様アンザッツ（同 λ 群は同 k、最短群 c≤100・他群 ≤3）。単一群は k≡1 で自明成立。\n")
with open(os.path.join(ROOT, "results", "hm_series_k.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["N", "groups", "verdict", "min_solution", "rank"]); w.writeheader(); w.writerows(rows)
print("PASS17 OK")
