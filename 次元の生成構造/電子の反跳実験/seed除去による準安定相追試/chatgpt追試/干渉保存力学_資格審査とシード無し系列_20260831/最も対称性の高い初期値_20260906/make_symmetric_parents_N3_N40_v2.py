#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""最も対称性の高い初期値 v2（N=3..40）——木原仕様（2026-09-06 合意）の生成・検査・一覧・図化。
v1（規則T=三つ組で 60°/120° を導入）は仕様違反のため廃止（git 履歴 cb1ae755 に記録のみ残す）。

木原仕様（3項）:
 1) 位相は全 M 本が 0° か 90° のみ（例外なし）。
 2) 本数が等しくできない場合（M 奇数）は、振幅を族ごとに統一した2値 r_A, r_B に分け、
    パワーバランス n_A r_A² = n_B r_B² で閉塞 Σz²=0 を厳密に達成する。
 3) 全 N で合計パワーを同一に: n_A r_A² + n_B r_B² = P_tot = 1（旧親の ‖Z‖=1 と同じ）。
 ⇒ 閉形式: r_A = √(1/(2 n_A))、r_B = √(1/(2 n_B))。M 偶数なら r_A=r_B=1/√M（完全等振幅）。

合意済みの割当規則（どの辺を 0° 族にするか）:
 辺を（巡回距離 d, i, j）の辞書順に並べ、先頭から n_A = ⌈M/2⌉ 本を 0° 族、残り n_B = ⌊M/2⌋ 本を
 90° 族とする。距離クラス単位で詰まり、端数が出るクラスだけ (i,j) 順で分割される（決定論）。
 「多い方を 0° に置く」のは全体 90° 回転のゲージで物理的に無差別。
 N≡1 (mod 4) では M/2 = N·(N−1)/4 がクラス N 本の整数倍なので分割クラスなし
 ＝ v1 規則B（d≤(N−1)/4 正則グラフ）と厳密に一致し、固定点（σ=N−1）が保たれる。

出力: parents_symmetric/parent_symmetric_N{N:05d}_v2.npz（現行親と同キー＋family,theta,rule）、
 listing/parent_symmetric_N{N:05d}_v2.csv、figures/fig_symmetric_parent_N{N}.png（個別）、
 figures/fig_symmetric_parents_grid_N3_N40.png、summary_symmetric_parents_v2.json、
 最も対称性の高い初期値の作り方.md（仕様＋全一覧、v2 で全面改訂）。
検査: ‖Z‖・|Σz²|・n_A/n_B・振幅2値（値と比 √(n_A/n_B) 一致）・分割クラス・
 固定点残差 ‖KZ−λZ‖/‖Z‖・σ・D_N 対称数・c_opp（反対族隣接数）範囲。"""
import csv
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
    dd = abs(int(i) - int(j))
    return min(dd, N - dd)


def assign_v2(N, ea, eb):
    """（d, i, j）辞書順で先頭 ⌈M/2⌉ 本を 0° 族に。返り値: family(0/1), d, 分割クラス情報。"""
    M = len(ea)
    d = np.array([cyc_dist(ea[e], eb[e], N) for e in range(M)])
    order = sorted(range(M), key=lambda e: (int(d[e]), int(min(ea[e], eb[e])), int(max(ea[e], eb[e]))))
    nA = (M + 1) // 2
    fam = np.ones(M, dtype=int)
    for e in order[:nA]:
        fam[e] = 0
    # 分割クラスの検出（族が混在する距離クラス）
    split = []
    for dv in sorted(set(int(x) for x in d)):
        cls = fam[d == dv]
        if cls.min() != cls.max():
            split.append({'d': dv, 'n_in_A': int(np.sum(cls == 0)), 'n_in_B': int(np.sum(cls == 1))})
    return fam, d, split


def dense_K(sys_lr, M):
    I = np.eye(M)
    return np.column_stack([sys_lr.kmatvec(I[:, j]) for j in range(M)])


def dihedral_perms(N):
    perms = []
    for s in range(N):
        perms.append(np.array([(i + s) % N for i in range(N)]))
        perms.append(np.array([(-i + s) % N for i in range(N)]))
    return perms


def symmetry_count(N, ea, eb, fam):
    """D_N の 2N 置換のうち、族割当（=位相パターン）を厳密に保つものの数。"""
    M = len(ea)
    index = {(min(ea[e], eb[e]), max(ea[e], eb[e])): e for e in range(M)}
    keep = 0
    for p in dihedral_perms(N):
        img = np.array([index[(min(p[ea[e]], p[eb[e]]), max(p[ea[e]], p[eb[e]]))] for e in range(M)])
        if np.array_equal(fam[img], fam):
            keep += 1
    return keep


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


summary = {}
Zs = {}
for N in NS:
    ea, eb = build_edges(N)
    ea = np.asarray(ea, dtype=int); eb = np.asarray(eb, dtype=int)
    M = len(ea)
    fam, d, split = assign_v2(N, ea, eb)
    nA = int(np.sum(fam == 0)); nB = int(np.sum(fam == 1))
    rA = math.sqrt(1.0 / (2 * nA)); rB = math.sqrt(1.0 / (2 * nB))
    r = np.where(fam == 0, rA, rB)
    theta = np.where(fam == 0, 0.0, math.pi / 2)
    Z = r * np.exp(1j * theta)
    # 検査
    norm = float(np.linalg.norm(Z))
    closure = float(abs(np.sum(Z ** 2)))
    n_amp = len(set(np.round(np.abs(Z), 15)))
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
    keep = symmetry_count(N, ea, eb, fam)
    rec = dict(N=N, M=M, M_odd=bool(M % 2), n_A=nA, n_B=nB,
               r_A=rA, r_B=rB, amp_ratio_rB_over_rA=rB / rA,
               amp_ratio_theory=math.sqrt(nA / nB), n_distinct_amp=n_amp,
               P_A=nA * rA * rA, P_B=nB * rB * rB, P_tot=norm * norm,
               norm=norm, closure_abs=closure,
               split_classes=split, n_split_classes=len(split),
               c_opp_min=int(c_opp.min()), c_opp_max=int(c_opp.max()), c_opp_target=N - 1,
               selfconsistency_residual=resid, sigma=sigma,
               is_fixed_point=bool(resid < 1e-10),
               dihedral_symmetries=keep, dihedral_order=2 * N,
               matches_v1_ruleB=bool(N % 4 == 1 and len(split) == 0))
    summary[N] = rec
    Zs[N] = Z
    np.savez_compressed(os.path.join(BASE, 'parents_symmetric', f'parent_symmetric_N{N:05d}_v2.npz'),
                        v=Z, g=np.zeros(M, dtype=np.complex128), Z0=Z, sigma=np.array([sigma]),
                        residual=np.array(resid), n=np.array(N), seed=np.array(-1), delta=np.array(0.0),
                        tol=np.array(0.0), iters=np.array(0),
                        family=fam, theta=theta, rule=np.array('v2_two_amp'))
    with open(os.path.join(BASE, 'listing', f'parent_symmetric_N{N:05d}_v2.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['idx', 'i', 'j', 'd', 'family', 'theta_deg', 'a', 'b', 'r'])
        for e in range(M):
            w.writerow([e, int(ea[e]), int(eb[e]), int(d[e]), int(fam[e]),
                        f'{math.degrees(theta[e]):.0f}',
                        f'{Z[e].real:.16e}', f'{Z[e].imag:.16e}', f'{float(r[e]):.16e}'])
    fig, ax = plt.subplots(figsize=(6, 6))
    draw_panel(ax, Z, f'N={N} (M={M}): nA={nA}@0d rA={rA:.4f}, nB={nB}@90d rB={rB:.4f}, |sum z^2|={closure:.1e}')
    fig.tight_layout(); fig.savefig(os.path.join(BASE, 'figures', f'fig_symmetric_parent_N{N}.png'), dpi=150); plt.close(fig)
    print(f"N={N:2d} M={M:3d}{'奇' if M % 2 else '偶'} nA/nB={nA}/{nB} rA={rA:.5f} rB={rB:.5f} 比={rB/rA:.5f}(理論√(nA/nB)={math.sqrt(nA/nB):.5f}) "
          f"P_A={rec['P_A']:.15f} |Σz²|={closure:.1e} ‖Z‖²={norm*norm:.15f} 分割クラス={len(split)} "
          f"固定点残差={resid:.1e} σ={sigma:.4f} D_N={keep}/{2*N}", flush=True)

# 俯瞰グリッド 8×5
fig, axs = plt.subplots(8, 5, figsize=(20, 24)); axs = axs.ravel()
for k, N in enumerate(NS):
    s = summary[N]
    axs[k].axis('on')
    draw_panel(axs[k], Zs[N], f"N={N} (M={s['M']}) nA/nB={s['n_A']}/{s['n_B']}", small=True)
for k in range(len(NS), 40):
    axs[k].axis('off')
fig.suptitle('Most symmetric initial states v2, N=3..40: phases 0/90 only; two amplitudes rA=sqrt(1/2nA), rB=sqrt(1/2nB); P_A=P_B=1/2, ||Z||=1', y=.998)
fig.tight_layout(); fig.savefig(os.path.join(BASE, 'figures', 'fig_symmetric_parents_grid_N3_N40.png'), dpi=160); plt.close(fig)

with open(os.path.join(BASE, 'summary_symmetric_parents_v2.json'), 'w') as fh:
    json.dump({'edge_order_source': 'run_n_scaling_lowrank_v1.build_edges(N)',
               'spec': ['phases 0/90 only', 'two family amplitudes, P_A=P_B', 'P_tot=1 for all N',
                        'assignment: (d,i,j) lexicographic, first ceil(M/2) edges -> 0 deg'],
               'summary': summary}, fh, indent=2, ensure_ascii=False)

# ---- 仕様書（全面改訂）----
L = []
L.append('# 最も対称性の高い初期値の作り方（N=3〜40）v2\n')
L.append('作成 2026-09-06（v2 全面改訂）。生成プログラム `make_symmetric_parents_N3_N40_v2.py`（本フォルダ）。新規走行なし。')
L.append('v1（M 奇数で 60°/120° の三つ組を導入）は木原仕様違反のため**廃止**（git 履歴 cb1ae755 のみ）。\n')
L.append('## 0. 仕様（木原 2026-09-06 指示・合意）\n')
L.append('1. **位相は全 M 本が 0° か 90° のみ**（例外なし）。')
L.append('2. 本数を等しくできない場合（M 奇数）は、**振幅を族ごとに統一した2値** r_A, r_B に分け、')
L.append('   **パワーバランス n_A r_A² = n_B r_B²** で閉塞 Σz² = 0 を厳密に達成する。')
L.append('3. **全 N で合計パワーを同一**にする: n_A r_A² + n_B r_B² = P_tot = 1（旧親の ‖Z‖=1 と同じ規約）。\n')
L.append('閉形式（この3項から一意）:\n')
L.append('- n_A = ⌈M/2⌉（0°族）、n_B = ⌊M/2⌋（90°族）')
L.append('- **r_A = √(1/(2 n_A))、r_B = √(1/(2 n_B))**、振幅比 r_B/r_A = √(n_A/n_B)')
L.append('- M 偶数 ⇒ r_A = r_B = 1/√M（完全等振幅）。M 奇数 ⇒ 2値（比は 1+O(1/M)）')
L.append('- 閉塞: Σz² = n_A r_A² − n_B r_B² = 1/2 − 1/2 = 0（恒等的に厳密）\n')
L.append('## 1. 割当規則（どの辺を 0° 族にするか・合意済み）\n')
L.append('辺を（巡回距離 d, i, j）の**辞書順**に並べ、先頭から n_A 本を 0°族、残りを 90°族とする。')
L.append('距離クラス単位で詰まり、端数が出るクラスだけ (i,j) 順で分割（決定論・最大限 D_N を保つ）。')
L.append('多い方をどちらの軸に置くかは全体 90° 回転のゲージで物理的に無差別。')
L.append('辺の番号順・(i,j) はエンジン `run_n_scaling_lowrank_v1.build_edges(N)` と同一（npz の Z0 をそのまま初期値に使える）。')
L.append('**N≡1 (mod 4)** では M/2 がクラスサイズ N の整数倍なので分割クラスが出ず、')
L.append('0°族 = 巡回距離 d≤(N−1)/4 の (N−1)/2 正則グラフとなって**厳密固定点（σ=N−1）**。\n')
L.append('## 2. 一覧（要約・実測検査値）\n')
L.append('| N | M | n_A/n_B | r_A | r_B | r_B/r_A | 分割クラス | ‖Z‖² | \\|Σz²\\| | 固定点残差 | σ | D_N 対称/位数 |')
L.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
for N in NS:
    s = summary[N]
    L.append(f"| {N} | {s['M']} | {s['n_A']}/{s['n_B']} | {s['r_A']:.6f} | {s['r_B']:.6f} | {s['amp_ratio_rB_over_rA']:.6f} | "
             f"{s['n_split_classes']} | {s['P_tot']:.12f} | {s['closure_abs']:.1e} | {s['selfconsistency_residual']:.1e} | {s['sigma']:.4f} | {s['dihedral_symmetries']}/{s['dihedral_order']} |")
L.append('\n固定点残差 = ‖KZ − λZ‖/‖Z‖（λ は Rayleigh 商、K はエンジン生成子）。機械零なら相対平衡（置いて静止）。')
L.append('固定点は N≡1 mod 4（5,9,13,17,21,25,29,33,37）のみ。他の N は仕様上の最対称状態だが固定点でなく、走行では点火の前に**緩和**が観測される（読み分けを事前登録）。\n')
L.append('## 3. 各 N の全一覧（M 個の複素数 z = a + ib）\n')
L.append('列: idx（エンジン辺番号）, (i,j) 体の対, d 巡回距離, family（0=0°族, 1=90°族）, θ, a, b, r。CSV は `listing/`。\n')
for N in NS:
    s = summary[N]
    ea, eb = build_edges(N); ea = np.asarray(ea); eb = np.asarray(eb)
    Z = Zs[N]
    dd = np.load(os.path.join(BASE, 'parents_symmetric', f'parent_symmetric_N{N:05d}_v2.npz'))
    fam, theta = dd['family'], dd['theta']
    L.append(f"### N={N}（M={s['M']}、n_A/n_B={s['n_A']}/{s['n_B']}、r_A=√(1/{2*s['n_A']})={s['r_A']:.6f}、r_B=√(1/{2*s['n_B']})={s['r_B']:.6f}）\n")
    L.append(f"閉塞 \\|Σz²\\|={s['closure_abs']:.1e}、‖Z‖²={s['P_tot']:.12f}、固定点残差 {s['selfconsistency_residual']:.1e}、σ={s['sigma']:.4f}、図 `figures/fig_symmetric_parent_N{N}.png`。分割クラス: {s['split_classes'] if s['split_classes'] else 'なし'}\n")
    L.append('| idx | (i,j) | d | family | θ | a | b | r |')
    L.append('|---|---|---|---|---|---|---|---|')
    for e in range(s['M']):
        L.append(f"| {e} | ({int(ea[e])},{int(eb[e])}) | {cyc_dist(ea[e], eb[e], N)} | {int(fam[e])} | {math.degrees(theta[e]):.0f}° | {Z[e].real:+.6f} | {Z[e].imag:+.6f} | {abs(Z[e]):.6f} |")
    L.append('')
with open(os.path.join(BASE, '最も対称性の高い初期値の作り方.md'), 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(L))
print('ALL DONE')
