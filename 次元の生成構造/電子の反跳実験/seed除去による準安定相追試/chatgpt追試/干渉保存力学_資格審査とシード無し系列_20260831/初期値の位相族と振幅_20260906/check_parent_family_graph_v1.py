#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""現行の親（N=3..40）の 0°族／90°族が辺集合としてどんなグラフか（読み出しのみ）。
木原（2026-09-06）「どちらの構造かは今のデータを見て判断すべき」→ 検討28で全Nが厳密90°二軸と判明済み。
残る問いは辺への割り当て構造。各 N について:
 (1) 族A（0°）の頂点次数列: 正則か（設計B=(N−1)/2 正則）／次数の範囲
 (2) 各辺の「反対族の隣接辺本数」c_e（隣接辺は 2(N−2) 本）: 等振幅なら和則 c_e = N−1 が固定点条件。
     現行親は非等振幅なので重み付き和則 (Σ_f r_f sin²φ_ef)/r_e = σ を併記（自己無撞着の確認）
 (3) 族A と族B の次数列が一致するか（自己補グラフの必要条件）
辺順はエンジンの build_edges(N) と同一（親 npz の Z0 の順）。
出力: parent_family_graph_v1.csv / check_parent_family_graph_v1.json"""
import csv, hashlib, json, os, sys
import numpy as np
BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
ROOT = os.path.abspath(os.path.join(BASE, '..', '..', '..', '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, '時間軸Q軸とフェルミオンの生成構造', '検証_対照実験', '第5論文原本_自発的分裂予備実験_v1'))
from run_n_scaling_lowrank_v1 import build_edges
led = {l.split()[1]: l.split()[0] for l in open(os.path.join(PKG, 'SHA256SUMS.txt')) if len(l.split()) == 2}
rows, out = [], {}
for N in range(3, 41):
    rel = f'parents/parent_static_N{N:05d}_makeparent_20260905.npz'
    path = os.path.join(PKG, rel)
    assert led[rel] == hashlib.sha256(open(path, 'rb').read()).hexdigest(), rel
    z = np.asarray(np.load(path)['Z0'], dtype=np.complex128); M = z.size
    ea, eb = build_edges(N); ea = np.asarray(ea); eb = np.asarray(eb)
    u = np.angle(z) % np.pi
    phi = (np.angle(np.mean(np.exp(4j * u))) / 4.0) % (np.pi / 2)
    fam = np.round((u - phi) / (np.pi / 2)).astype(int) % 2          # 0=族A, 1=族B
    r2 = np.abs(z) ** 2
    # 次数列
    degA = np.zeros(N, int); degB = np.zeros(N, int)
    for e in range(M):
        (degA if fam[e] == 0 else degB)[ea[e]] += 1
        (degA if fam[e] == 0 else degB)[eb[e]] += 1
    # 各辺の反対族隣接数と重み付き和則
    inc = [[] for _ in range(N)]
    for e in range(M): inc[ea[e]].append(e); inc[eb[e]].append(e)
    c_opp = np.zeros(M, int); w_sum = np.zeros(M)
    theta = np.angle(z)
    for e in range(M):
        nb = [f for f in inc[ea[e]] + inc[eb[e]] if f != e]
        c_opp[e] = sum(1 for f in nb if fam[f] != fam[e])
        # 固有関係 KZ=iσZ ⇒ Σ_f r_f sin²φ_ef = σ r_e（r は振幅、r² ではない）。不変量は左辺/r_e
        w_sum[e] = sum(np.sqrt(r2[f]) * np.sin(theta[f] - theta[e]) ** 2 for f in nb) / np.sqrt(r2[e])
    regA = degA.min() == degA.max(); regB = degB.min() == degB.max()
    rec = dict(N=N, M=M, n_A=int((fam == 0).sum()), n_B=int((fam == 1).sum()),
               degA_min=int(degA.min()), degA_max=int(degA.max()), degA_regular=bool(regA),
               degB_min=int(degB.min()), degB_max=int(degB.max()), degB_regular=bool(regB),
               deg_seqs_equal=bool(sorted(degA) == sorted(degB)),
               c_opp_min=int(c_opp.min()), c_opp_max=int(c_opp.max()), c_opp_target_Nm1=N - 1,
               frac_edges_c_opp_eq_Nm1=float(np.mean(c_opp == N - 1)),
               weighted_sumrule_cv=float(w_sum.std() / w_sum.mean()),
               weighted_sumrule_mean=float(w_sum.mean()))
    rows.append(rec); out[N] = rec
    print(f"N={N:2d} M={M:3d} A/B={rec['n_A']}/{rec['n_B']} degA=[{degA.min()},{degA.max()}]{'正則' if regA else '　　'} "
          f"degB=[{degB.min()},{degB.max()}]{'正則' if regB else '　　'} 次数列一致={rec['deg_seqs_equal']!s:5s} "
          f"c_opp=[{c_opp.min()},{c_opp.max()}] (N−1={N-1}) 一致率={rec['frac_edges_c_opp_eq_Nm1']:.2f} 重み付き和則CV={rec['weighted_sumrule_cv']:.1e}")
with open(os.path.join(BASE, 'parent_family_graph_v1.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
json.dump(out, open(os.path.join(BASE, 'check_parent_family_graph_v1.json'), 'w'), indent=2)
print('ALL DONE')
