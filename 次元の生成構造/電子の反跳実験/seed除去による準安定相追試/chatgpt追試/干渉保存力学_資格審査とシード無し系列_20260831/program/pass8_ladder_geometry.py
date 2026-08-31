#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス8（読出し専用・実験目的の正本 §8 に直結）。状態は states_treatment.npz（全 step 保存）から 25 時点を抽出：
(A) z² 角度スペクトル（倍音梯子）の時間発展：各 snapshot で z_e² の角度を |z_e|² 重みでクラスタリングし、
    占有ライン数（重み≥1% のクラスタ数）と角度エントロピーを測る。親=全梯子（q 本）→ 星=1 本 の潰れ方と、
    非自明中間パターン（途中の本数の踊り場）の有無を見る。
(B) 最終状態の実読み：共通軸へ回転（γ=arg(Σz²)/2）後、d²=(z e^{−iγ})² から B=−½JD J を作り、
    虚部残差・最小固有値（PSD か）・rank・相異なる辺長の数を測る＝どんな幾何が生まれたか。
(C) 代表走行（hm_N12）の角度×重みの絵を snapshot 4 時点で描く。
出力：results/ladder_lines_timeseries.csv、results/final_geometry.{csv,md}、figures/fig9_ladder_lines_grid.png、figures/fig10_ladder_snapshots_hm_N12.png"""
import os, csv, json
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams['font.family']=['Hiragino Sans','Arial Unicode MS','sans-serif']
import sys
sys.path.insert(0, os.path.dirname(__file__))
from common import edges
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data"); FIG = os.path.join(ROOT, "figures")

def angle_lines(Z, wmin=0.01, gap=2*np.pi*0.02):
    """z² の角度を重みでクラスタリング。(ライン数[重み≥wmin], エントロピー, 最大ライン重み)"""
    w = np.abs(Z)**2; W = w.sum()
    if W <= 0: return 0, 0.0, 0.0
    ang = np.angle(Z**2)
    idx = np.argsort(ang); a = ang[idx]; ww = w[idx]/W
    # 円周上のギャップでクラスタ分割
    cl = [0]*len(a); c = 0
    for i in range(1, len(a)):
        if a[i]-a[i-1] > gap: c += 1
        cl[i] = c
    # 先頭と末尾が円周で近ければ結合
    if len(a) > 1 and (a[0]+2*np.pi)-a[-1] <= gap and c > 0:
        cl = [0 if x == c else x for x in cl]
    weights = {}
    for i, ci in enumerate(cl): weights[ci] = weights.get(ci, 0.0)+ww[i]
    wl = sorted(weights.values(), reverse=True)
    nlines = sum(1 for x in wl if x >= wmin)
    p = np.array([x for x in wl if x > 0])
    ent = float(-(p*np.log(p)).sum())
    return nlines, ent, float(wl[0])

rows_ts = []; geo = []
for tag in sorted(os.listdir(DATA)):
    if tag == "reference": continue
    sp = os.path.join(DATA, tag, "states_treatment.npz")
    if not os.path.exists(sp): continue
    allZ = np.load(sp)["Z"]
    steps = np.array(sorted(set([0,25,50,75,100,125,150,200,300,500,750,1000,1500,2000,
                                 3000,4000,5000,7500,10000,15000,20000,25000,30000,35000,40000])))
    steps = steps[steps < len(allZ)]
    ZZ = allZ[steps]
    N = json.load(open(os.path.join(DATA, tag, "parent_checks.json")))["N"]
    for t, Z in zip(steps, ZZ):
        nl, ent, wmax = angle_lines(Z)
        rows_ts.append(dict(tag=tag, N=N, method=tag.split("_")[0], step=int(t),
                            n_lines=nl, angle_entropy=ent, max_line_weight=wmax))
    # (B) 最終状態の実読み
    Z = ZZ[-1]; M = len(Z); E = edges(N)
    gamma = np.angle(np.sum(Z*Z))/2.0
    zr = Z*np.exp(-1j*gamma)
    d2 = zr*zr
    D = np.zeros((N, N), complex)
    for val, (i, j) in zip(d2, E): D[i, j] = D[j, i] = val
    J = np.eye(N)-np.ones((N, N))/N; B = -0.5*J@D@J
    im_res = float(np.linalg.norm(B.imag)/max(np.linalg.norm(B.real), 1e-300))
    ev = np.linalg.eigvalsh(B.real); scale = max(abs(ev).max(), 1e-300)
    min_eig_rel = float(ev.min()/scale)
    rank = int((ev/scale > 1e-8).sum())
    amps = np.abs(Z); aw = amps/amps.max() if amps.max() > 0 else amps
    sig = np.sort(aw)[::-1]
    dist_lengths = int(len(set(np.round(sig[sig > 1e-6], 6))))
    nl_fin, ent_fin, wmax_fin = angle_lines(Z)
    nl_par, ent_par, _ = angle_lines(ZZ[0])
    geo.append(dict(tag=tag, N=N, method=tag.split("_")[0],
                    lines_parent=nl_par, lines_final=nl_fin,
                    z2_align_final=wmax_fin,
                    imB_over_reB=im_res, min_eig_rel=min_eig_rel,
                    PSD=bool(min_eig_rel > -1e-6), rank=rank,
                    n_distinct_lengths=dist_lengths))
with open(os.path.join(ROOT, "results", "ladder_lines_timeseries.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows_ts[0].keys())); w.writeheader(); w.writerows(rows_ts)
with open(os.path.join(ROOT, "results", "final_geometry.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(geo[0].keys())); w.writeheader(); w.writerows(geo)
lines = ["# 最終状態の実読み（γ 回転後の d²=z² から B、既存 snapshots のみ）", "",
         "| tag | 梯子:親→終 | z²整列(終) | ‖ImB‖/‖ReB‖ | min eig/scale | PSD | rank | 相異長数 |", "|---|---|---|---|---|---|---|---|"]
for r in geo:
    lines.append(f"| {r['tag']} | {r['lines_parent']}→{r['lines_final']} | {r['z2_align_final']:.4f} | "
                 f"{r['imB_over_reB']:.1e} | {r['min_eig_rel']:+.1e} | {'Y' if r['PSD'] else 'N'} | {r['rank']} | {r['n_distinct_lengths']} |")
with open(os.path.join(ROOT, "results", "final_geometry.md"), "w") as f:
    f.write("\n".join(lines)+"\n")

# fig9: hm の占有ライン数 vs step（4×4 グリッド、fig3 と同レイアウト）
by = {}
for r in rows_ts: by.setdefault((r["method"], r["N"]), []).append(r)
col = {'mp': 'k', 'hm': '#1f77b4', 'ne': '#ff7f0e', 'rb': '#2ca02c'}
fg, axs = plt.subplots(4, 4, figsize=(17, 14)); axs = axs.ravel()
for i, N in enumerate(range(3, 17)):
    ax = axs[i]
    for m in ('mp', 'hm', 'ne', 'rb'):
        rs = sorted(by.get((m, N), []), key=lambda r: r["step"])
        if not rs: continue
        ax.plot([r["step"] for r in rs], [r["n_lines"] for r in rs], color=col[m], lw=1.2, marker='.', ms=3, label=m)
    ax.set_title(f'N={N}', fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6)
    ax.set_ylim(0, max(16, N))
    if i % 4 == 0: ax.set_ylabel('z² 角度ラインの占有数（重み≥1%）')
    if i >= 10: ax.set_xlabel('step')
for j in range(14, 16): axs[j].axis('off')
fg.suptitle('倍音梯子の占有ライン数の時間発展（snapshots 25 時点）：全梯子（親）→ 1 本（星）への潰れ方', fontsize=11)
plt.tight_layout(rect=(0, 0, 1, 0.97)); plt.savefig(os.path.join(FIG, 'fig9_ladder_lines_grid.png'), dpi=130); plt.close()

# fig10: hm_N12 の角度×重み 4 時点
allZ = np.load(os.path.join(DATA, 'hm_N12', 'states_treatment.npz'))["Z"]
steps = np.array(sorted(set([0,25,50,75,100,125,150,200,300,500,750,1000,1500,2000,3000,4000,5000,7500,10000,15000,20000,25000,30000,35000,40000])))
steps = steps[steps < len(allZ)]; ZZ = allZ[steps]
pick = [0, np.searchsorted(steps, 10000), np.searchsorted(steps, 20000), len(steps)-1]
fg, axs = plt.subplots(1, 4, figsize=(16, 3.6))
for ax, ip in zip(axs, pick):
    Z = ZZ[ip]; w = np.abs(Z)**2; w = w/w.sum()
    ang = (np.angle(Z**2)-np.angle(np.sum(Z*Z))+np.pi)%(2*np.pi)-np.pi
    ax.vlines(np.degrees(ang), 0, w, color='#1f77b4', lw=2)
    ax.set_xlim(-180, 180); ax.set_ylim(0, 1.0); ax.grid(alpha=.3)
    ax.set_title(f'step {steps[ip]}', fontsize=9); ax.set_xlabel('arg z² − 共通軸 [deg]')
axs[0].set_ylabel('重み |z|²/H')
fg.suptitle('hm_N12：z² 角度スペクトル（重み付き）の時間発展——12 本の梯子 → 1 本への凝縮', fontsize=10)
plt.tight_layout(rect=(0, 0, 1, 0.93)); plt.savefig(os.path.join(FIG, 'fig10_ladder_snapshots_hm_N12.png'), dpi=150); plt.close()
print("PASS8 OK")
# コンソール要約
import collections
tr = collections.defaultdict(list)
for r in rows_ts:
    if r["method"] == "hm": tr[r["N"]].append((r["step"], r["n_lines"]))
for N in sorted(tr):
    seq = [x[1] for x in sorted(tr[N])]
    print(f"hm_N{N}: lines {seq}")
