#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""最も対称性の高い初期値（N=3..40）の生成・検査・一覧・図化（新規走行なし、初期値の作成のみ）。

根拠（2026-09-06、検討28・29・実測 check_parent_family_graph_v1）:
 - 現行の親は全 N で「全波が直交2軸のどちらかに厳密に乗る」90°二族構造（(i)型）。
   ただし族の辺割当は乱数的（次数列が不揃い、正則グラフでない）で、閉塞は
   振幅の不均衡（P_A=P_B）で達成している。
 - 木原の要件: 全て等振幅・相対位相差 90°。等振幅なら閉塞 Σz²=r²(n_A−n_B)=0 ⇔ n_A=n_B。
 - 等振幅の固定点条件（和則）: 各辺の隣接 2(N−2) 本のうち反対族がちょうど N−1 本
   ⇔ 0°族が (N−1)/2 正則グラフ ⇔ N≡1 (mod 4)。

規則（N mod 4 で3通り。辺順はエンジン build_edges(N) と同一）:
 [B] N≡1 (mod 4): 頂点を円周に並べ、巡回距離 d≤(N−1)/4 の辺を 0°、残りを 90°。
     等振幅・閉塞厳密・(N−1)/2 正則 ⇒ 固定点。D_N 対称。
 [SC] N≡0 (mod 4): 自己補グラフ（P4 再帰構成: 4個ずつのブロック、ブロック内の道 P4 と
     中央2頂点から先行ブロック全頂点への辺を 0°族）を 0°、補グラフを 90°。
     n_A=n_B=M/2 で閉塞厳密。正則でないので固定点ではない（置いた瞬間から動く）。
 [T] M 奇数（N≡2,3 mod 4）: 二族のみでは等振幅閉塞が不可能（計数）。最小の逸脱として
     三角形 {0,1,2} の3辺を (0°,60°,120°) の閉塞三つ組（z² が 0°,120°,240° で零和）とし、
     残り M−3 本（偶数）を (距離, i, j) 順に 0°/90° を交互に割当（厳密に半々）。
     固定点ではない。N=3 は三つ組のみ（=N=3 結晶と同じ配置）。
 振幅: 全辺 r=1/√M（‖Z‖=1 規約、現行親と同じ。スケールは物理に無関係）。

出力:
 parents_symmetric/parent_symmetric_N{N:05d}_v1.npz  … 現行親と同じキー（v,g,Z0,sigma,residual,
   n,seed,delta,tol,iters）＋ family, theta, rule。走行プログラムは Z0 のみ読む。
 listing/parent_symmetric_N{N:05d}_v1.csv … 全 M 辺の (idx,i,j,d,family,theta_deg,a,b,r)
 figures/fig_symmetric_parent_N{N}.png（個別）／fig_symmetric_parents_grid_N3_N40.png（8×5）
 summary_symmetric_parents_v1.json / 最も対称性の高い初期値の作り方.md（根拠＋全一覧）
検査（各 N）: ‖Z‖、|Σz²|、n_A/n_B、和則 c_opp（反対族隣接数）の範囲と N−1 一致率、
 自己無撞着残差 ‖KZ−λZ‖/‖Z‖（λ=Rayleigh 商; 固定点なら機械零）、σ=|λ|、
 D_N（回転・鏡映 2N 個）のうち割当を保つ置換の数（厳密／大域位相を除く）。"""
import csv
import hashlib
import json
import math
import os
import sys
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(BASE, '..', '..', '..', '..', '..', '..'))
ENGINE = os.path.join(ROOT, '時間軸Q軸とフェルミオンの生成構造', '検証_対照実験', '第5論文原本_自発的分裂予備実験_v1')
sys.path.insert(0, ENGINE)
from run_n_scaling_lowrank_v1 import LowRankSystem, build_edges  # 辺順の正本

for sub in ('parents_symmetric', 'listing', 'figures'):
    os.makedirs(os.path.join(BASE, sub), exist_ok=True)

NS = list(range(3, 41))


def cyc_dist(i, j, N):
    d = abs(int(i) - int(j))
    return min(d, N - d)


def assign(N, ea, eb):
    M = len(ea)
    fam = np.full(M, -1, dtype=int)
    theta = np.zeros(M)
    d = np.array([cyc_dist(ea[e], eb[e], N) for e in range(M)])
    if N % 4 == 1:
        k = (N - 1) // 4
        fam = np.where(d <= k, 0, 1)
        theta = np.where(fam == 0, 0.0, math.pi / 2)
        rule = 'B'
    elif N % 4 == 0:
        K = N // 4
        inA = set()
        for k in range(K):
            v0, v1, v2, v3 = 4 * k, 4 * k + 1, 4 * k + 2, 4 * k + 3
            inA |= {(v0, v1), (v1, v2), (v2, v3)}
            for u in range(4 * k):
                inA.add((min(u, v1), max(u, v1)))
                inA.add((min(u, v2), max(u, v2)))
        for e in range(M):
            key = (min(ea[e], eb[e]), max(ea[e], eb[e]))
            fam[e] = 0 if key in inA else 1
        theta = np.where(fam == 0, 0.0, math.pi / 2)
        rule = 'SC'
    else:
        tri = {(0, 1): 0.0, (1, 2): math.pi / 3, (0, 2): 2 * math.pi / 3}
        rest = []
        for e in range(M):
            key = (min(ea[e], eb[e]), max(ea[e], eb[e]))
            if key in tri:
                fam[e] = 2
                theta[e] = tri[key]
            else:
                rest.append(e)
        order = sorted(rest, key=lambda e: (int(d[e]), min(ea[e], eb[e]), max(ea[e], eb[e])))
        for idx, e in enumerate(order):
            fam[e] = idx % 2
            theta[e] = 0.0 if idx % 2 == 0 else math.pi / 2
        rule = 'T'
    return fam, theta, d, rule


def dense_K(sys_lr, M):
    I = np.eye(M)
    return np.column_stack([sys_lr.kmatvec(I[:, j]) for j in range(M)])


def dihedral_perms(N):
    perms = []
    for s in range(N):
        perms.append(('rot', s, np.array([(i + s) % N for i in range(N)])))
        perms.append(('ref', s, np.array([(-i + s) % N for i in range(N)])))
    return perms


def symmetry_count(N, ea, eb, theta):
    M = len(ea)
    index = {(min(ea[e], eb[e]), max(ea[e], eb[e])): e for e in range(M)}
    strict = 0
    gauge = 0
    for kind, s, p in dihedral_perms(N):
        img = np.array([index[(min(p[ea[e]], p[eb[e]]), max(p[ea[e]], p[eb[e]]))] for e in range(M)])
        diff = (theta[img] - theta + math.pi) % (2 * math.pi) - math.pi
        if np.all(np.abs(diff) < 1e-12):
            strict += 1
        if np.all(np.abs(diff - diff[0]) < 1e-12) or np.all(np.abs(((diff - diff[0] + math.pi) % (2 * math.pi)) - math.pi) < 1e-12):
            gauge += 1
    return strict, gauge


def draw_panel(ax, z, title, small=False):
    """plot_complex_plane_N3_N40_stage123_v1.py の per-panel 描画コード verbatim 流用。"""
    for w in z:
        ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.7, alpha=0.6)
    ax.plot(z.real, z.imag, 'o', ms=2.5, color='tab:red', alpha=0.85, linestyle='none')
    cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
    for (a, b), c in cnt.items():
        if c > 1:
            ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3),
                        fontsize=5 if small else 7, color='black')
    r = float(np.abs(z).max())
    lim = r * 1.15 if r > 0 else 1.0
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.grid(alpha=.25)
    ax.tick_params(labelsize=6 if small else 8)
    ax.ticklabel_format(style='sci', scilimits=(-2, 3))
    ax.set_title(title, fontsize=9 if small else 11)
    if not small:
        ax.set_xlabel('Re z (a)', fontsize=9); ax.set_ylabel('Im z (b)', fontsize=9)


RULE_TEXT = {
    'B': '[B] 巡回距離 d≤(N−1)/4 → 0°、残り → 90°（(N−1)/2 正則・固定点）',
    'SC': '[SC] 自己補グラフ（P4再帰）→ 0°、補グラフ → 90°（半々・固定点でない）',
    'T': '[T] 三角形{0,1,2}=(0°,60°,120°) 三つ組 ＋ 残り交互 0°/90°（半々・固定点でない）',
}

summary = {}
Zs = {}
for N in NS:
    ea, eb = build_edges(N)
    ea = np.asarray(ea, dtype=int); eb = np.asarray(eb, dtype=int)
    M = len(ea)
    fam, theta, d, rule = assign(N, ea, eb)
    r = 1.0 / math.sqrt(M)
    Z = r * np.exp(1j * theta)
    # 検査
    norm = float(np.linalg.norm(Z))
    closure = float(abs(np.sum(Z ** 2)))
    nA, nB, nT = int(np.sum(fam == 0)), int(np.sum(fam == 1)), int(np.sum(fam == 2))
    sys_lr = LowRankSystem(N)
    sys_lr.set_theta(np.angle(Z))
    K = dense_K(sys_lr, M)
    KZ = K @ Z
    lam = complex(np.vdot(Z, KZ) / np.vdot(Z, Z))
    resid = float(np.linalg.norm(KZ - lam * Z) / np.linalg.norm(Z))
    sigma = float(abs(lam))
    inc = [[] for _ in range(N)]
    for e in range(M):
        inc[ea[e]].append(e); inc[eb[e]].append(e)
    c_opp = np.zeros(M, int)
    for e in range(M):
        nb = [f for f in inc[ea[e]] + inc[eb[e]] if f != e]
        c_opp[e] = sum(1 for f in nb if fam[f] != fam[e])
    two = fam < 2
    strict, gauge = symmetry_count(N, ea, eb, theta)
    degA = np.zeros(N, int)
    for e in range(M):
        if fam[e] == 0:
            degA[ea[e]] += 1; degA[eb[e]] += 1
    rec = dict(N=N, M=M, rule=rule, rule_text=RULE_TEXT[rule], r=r,
               n_A=nA, n_B=nB, n_triple=nT, norm=norm, closure_abs=closure,
               degA_min=int(degA.min()), degA_max=int(degA.max()),
               c_opp_min=(int(c_opp[two].min()) if two.any() else -1), c_opp_max=(int(c_opp[two].max()) if two.any() else -1),
               frac_c_opp_eq_Nm1=(float(np.mean(c_opp[two] == N - 1)) if two.any() else float('nan')),
               selfconsistency_residual=resid, sigma=sigma,
               is_fixed_point=bool(resid < 1e-10),
               dihedral_symmetries_strict=strict, dihedral_symmetries_gauge=gauge, dihedral_order=2 * N)
    summary[N] = rec
    Zs[N] = Z
    # npz（現行親と同じキー＋付加）
    np.savez_compressed(os.path.join(BASE, 'parents_symmetric', f'parent_symmetric_N{N:05d}_v1.npz'),
                        v=Z, g=np.zeros(M, dtype=np.complex128), Z0=Z, sigma=np.array([sigma]),
                        residual=np.array(resid), n=np.array(N), seed=np.array(-1), delta=np.array(0.0),
                        tol=np.array(0.0), iters=np.array(0),
                        family=fam, theta=theta, rule=np.array(rule))
    # 一覧 CSV
    with open(os.path.join(BASE, 'listing', f'parent_symmetric_N{N:05d}_v1.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['idx', 'i', 'j', 'd', 'family', 'theta_deg', 'a', 'b', 'r'])
        for e in range(M):
            w.writerow([e, int(ea[e]), int(eb[e]), int(d[e]), int(fam[e]), f'{math.degrees(theta[e]):.1f}',
                        f'{Z[e].real:.16e}', f'{Z[e].imag:.16e}', f'{r:.16e}'])
    # 個別図
    fig, ax = plt.subplots(figsize=(6, 6))
    draw_panel(ax, Z, f'N={N} (M={M}) rule {rule}: r=1/sqrt(M)={r:.4f}, |sum z^2|={closure:.1e}, resid={resid:.1e}')
    fig.tight_layout(); fig.savefig(os.path.join(BASE, 'figures', f'fig_symmetric_parent_N{N}.png'), dpi=150); plt.close(fig)
    print(f"N={N:2d} M={M:3d} rule={rule:2s} A/B/T={nA}/{nB}/{nT} |Σz²|={closure:.1e} degA=[{degA.min()},{degA.max()}] "
          f"c_opp=[{rec['c_opp_min']},{rec['c_opp_max']}] (N−1={N-1}) 固定点残差={resid:.1e} σ={sigma:.4f} D_N対称={strict}/{gauge}/{2*N}", flush=True)

# 俯瞰グリッド 8×5
fig, axs = plt.subplots(8, 5, figsize=(20, 24)); axs = axs.ravel()
for k, N in enumerate(NS):
    draw_panel(axs[k], Zs[N], f'N={N} (M={N*(N-1)//2}) rule {summary[N]["rule"]}', small=True)
for k in range(len(NS), 40):
    axs[k].axis('off')
fig.suptitle('Most symmetric initial states, N=3..40: equal amplitude r=1/sqrt(M), phases on 0/90 (rule B/SC) or 0/60/120 triple + 0/90 (rule T)', y=.998)
fig.tight_layout(); fig.savefig(os.path.join(BASE, 'figures', 'fig_symmetric_parents_grid_N3_N40.png'), dpi=160); plt.close(fig)

with open(os.path.join(BASE, 'summary_symmetric_parents_v1.json'), 'w') as fh:
    json.dump({'edge_order_source': 'run_n_scaling_lowrank_v1.build_edges(N)', 'rules': RULE_TEXT,
               'summary': summary}, fh, indent=2, ensure_ascii=False)

# ---- 文書 ----
L = []
L.append('# 最も対称性の高い初期値の作り方（N=3〜40）\n')
L.append('作成 2026-09-06。生成プログラム `make_symmetric_parents_N3_N40_v1.py`（本フォルダ）。新規走行なし。\n')
L.append('## 0. 根拠（現行データの判定）\n')
L.append('- 現行の親（make_parent 静的親、N=3..40）は**全 N で全波が直交2軸のどちらかに厳密に乗る**（検討28: 4回対称秩序変数 1.0000、軸ずれ 0.00°）。したがって現行系は「90°二族」型であり、本書はその**理想化**を作る。')
L.append('- 現行の親の族割当は乱数的（`check_parent_family_graph_v1`: 次数列が不揃い、正則グラフは N=4 のみ）で、閉塞は振幅の不均衡（族別パワー均衡 P_A=P_B）で達成している。')
L.append('- 木原要件: **全て等振幅・相対位相差 90°**。等振幅 r なら Σz² = r²(n_A − n_B) ⇒ 閉塞 ⇔ **n_A = n_B**。')
L.append('- 等振幅の固定点（相対平衡）条件は和則「各辺の隣接 2(N−2) 本のうち反対族がちょうど N−1 本」⇔ 0°族が (N−1)/2 正則グラフ ⇔ **N≡1 (mod 4)**。')
L.append('- 絶対位相は U(1) ゲージ（検討29）。0°/90° を正準ゲージとして採用。振幅の絶対値は物理に無関係、規約 ‖Z‖=1 ⇒ **r = 1/√M**。\n')
L.append('## 1. 規則（N mod 4 で3通り）\n')
L.append('| 規則 | 適用 N | 位相 | 閉塞 | 固定点 |')
L.append('|---|---|---|---|---|')
L.append('| **[B]** | N≡1 mod 4（5,9,…,37） | 巡回距離 d≤(N−1)/4 の辺 0°、残り 90° | n_A=n_B で厳密 | **成立**（(N−1)/2 正則、D_N 対称） |')
L.append('| **[SC]** | N≡0 mod 4（4,8,…,40） | 自己補グラフ（P4 再帰構成）0°、補グラフ 90° | n_A=n_B で厳密 | 不成立（正則でない）→ 置いた瞬間から動く |')
L.append('| **[T]** | M 奇数（N≡2,3 mod 4） | 三角形{0,1,2}を (0°,60°,120°) の三つ組、残り M−3 本を (d,i,j) 順に 0°/90° 交互 | 三つ組 z² が 0°,120°,240° で零和＋残り半々 で厳密 | 不成立。N=3 は三つ組のみ |\n')
L.append('辺の番号順・(i,j) の定義はエンジン `run_n_scaling_lowrank_v1.build_edges(N)` と同一（npz の Z0 をそのまま初期値に使える）。\n')
L.append('## 2. 一覧（要約）\n')
L.append('| N | M | 規則 | r | n_A/n_B/三つ組 | ‖Z‖ | \\|Σz²\\| | c_opp 範囲（目標 N−1） | 固定点残差 | σ | D_N 対称（厳密/ゲージ/位数） |')
L.append('|---|---|---|---|---|---|---|---|---|---|---|')
for N in NS:
    s = summary[N]
    L.append(f"| {N} | {s['M']} | {s['rule']} | {s['r']:.5f} | {s['n_A']}/{s['n_B']}/{s['n_triple']} | {s['norm']:.12f} | {s['closure_abs']:.1e} | [{s['c_opp_min']},{s['c_opp_max']}]（{N-1}） | {s['selfconsistency_residual']:.1e} | {s['sigma']:.4f} | {s['dihedral_symmetries_strict']}/{s['dihedral_symmetries_gauge']}/{s['dihedral_order']} |")
L.append('\n固定点残差 = ‖KZ − λZ‖/‖Z‖（λ は Rayleigh 商、K はエンジンの生成子）。機械零なら相対平衡。\n')
L.append('## 3. 各 N の全一覧（M 個の複素数 z = a + ib）\n')
L.append('列: idx（エンジン辺番号）, (i,j) 体の対, d 巡回距離, family（0=0°族, 1=90°族, 2=三つ組）, θ, a, b。r=1/√M。CSV は `listing/`。\n')
for N in NS:
    s = summary[N]
    ea, eb = build_edges(N); ea = np.asarray(ea); eb = np.asarray(eb)
    Z = Zs[N]
    dd = np.load(os.path.join(BASE, 'parents_symmetric', f'parent_symmetric_N{N:05d}_v1.npz'))
    fam, theta = dd['family'], dd['theta']
    L.append(f"### N={N}（M={s['M']}、規則 {s['rule']}、r=1/√{s['M']}={s['r']:.6f}）\n")
    L.append(f"根拠: {s['rule_text']}。閉塞 \\|Σz²\\|={s['closure_abs']:.1e}、固定点残差 {s['selfconsistency_residual']:.1e}、図 `figures/fig_symmetric_parent_N{N}.png`。\n")
    L.append('| idx | (i,j) | d | family | θ | a | b |')
    L.append('|---|---|---|---|---|---|---|')
    for e in range(s['M']):
        L.append(f"| {e} | ({int(ea[e])},{int(eb[e])}) | {cyc_dist(ea[e], eb[e], N)} | {int(fam[e])} | {math.degrees(theta[e]):.0f}° | {Z[e].real:+.6f} | {Z[e].imag:+.6f} |")
    L.append('')
with open(os.path.join(BASE, '最も対称性の高い初期値の作り方.md'), 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(L))
print('ALL DONE')
