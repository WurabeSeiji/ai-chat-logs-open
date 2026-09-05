#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文7 図3（5色占有 log・比較）を、段1+2+3 の 10,000 歩データで再現（読み出しのみ）。

問い（木原 2026-09-06）: 種なし系論文の「第3の方向まで含めた三方向構造の定着」は
その論文（旧エンジン＋5色読出し）のアーティファクトか。同じ読出しを今回のデータ
（stage1+2+3、N=5, 40、den=N、10,000歩）に適用して判定する。

方法（系列規約: 過去論文依拠はコピー→対照テスト→import）:
- 基底構成（parent_plane_split_exact / gram_reduce / dominant_plane / make_parent /
  zero_closure_kernel_seed）は第5論文原本エンジンから**そのまま import**。
- 5色読出し（occ / s4_new_dirs / align_2d / E_d3・E_d4・残余・核の式）は
  run_paper7_5color_timeseries.py（正本）からの逐語コピー。
- 時間発展だけを今回の状態 npz の読込みに置換（力学は一切走らせない）。
- 対照ゲート: (a) 入力 npz を SHA 台帳と照合、(b) make_parent の v・種 g・初期 Z が
  今回の静的親 npz および states[0] と bit 一致すること。
出力: fig3_compare_stage123_N5_N40.png/.svg、bands CSV、meta JSON。"""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[5]                     # ai-chat-logs-open/
ENGINE = ROOT / "時間軸Q軸とフェルミオンの生成構造" / "検証_対照実験" / "第5論文原本_自発的分裂予備実験_v1"
V2 = ENGINE / "exact_lowN_eigenspectrum_v2"
sys.path.insert(0, str(ENGINE))
sys.path.insert(0, str(V2 / "code"))
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed
from run_plane_flow_exact_v1 import parent_plane_split_exact
from run_n300_dimension_saturation_v2 import gram_reduce, dominant_plane

SWEEP = BASE.parent / "N3_N40_long10000_20260905"
PARENTS = BASE.parent / "N3_N40_stage123_sweep_20260905" / "parents"
DELTA = 1e-15
XMAX = 10000
SAMPLE = 25
NS = [5, 40]
COLORS = ["#4C78A8", "#E45756", "#F58518", "#B0B0B0", "#54A24B"]  # P1, d3, d4, other残, 核
LABELS = ["P1 (dominant plane)", "direction 3", "direction 4", "remaining other-rotation", "kernel"]

ledger = {}
with open(SWEEP / "SHA256SUMS.txt") as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

def gate_npz(rel):
    h = hashlib.sha256(open(SWEEP / rel, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'

# ---- run_paper7_5color_timeseries.py 42-45, 66-80 行の逐語コピー ----
def occ(B, Z):
    if B is None or (hasattr(B, "shape") and B.shape[1] == 0):
        return 0.0
    return float(np.sum((B.T @ Z.real) ** 2) + np.sum((B.T @ Z.imag) ** 2))

def s4_new_dirs(B0, Bdom):
    """S4=orthonormalize[B0|Bdom] の B0 直交補2方向 e3,e4 を返す。"""
    R = Bdom - B0 @ (B0.T @ Bdom)
    Qr, _ = np.linalg.qr(R)
    return Qr[:, :2]

def align_2d(f_prev, f_new):
    """f_new(M×2) を前時刻 f_prev へ 2×2 回転で整列（連続基底固定・色反転防止）。"""
    if f_prev is None:
        return f_new
    Ov = f_prev.T @ f_new                # 2×2
    U, _, Vt = np.linalg.svd(Ov)
    Rot = U @ Vt                          # 直交 2×2
    return f_new @ Rot.T
# ---- コピーここまで ----

def run(n):
    rel = f'results/hm_N{n}_den_{n}_states_10000.npz'
    gate_npz(rel)
    Zs = np.asarray(np.load(SWEEP / rel)['Z'], dtype=np.complex128)
    # 基底構成（paper7 build() と同一手順）
    sys_lr = LowRankSystem(n); M = sys_lr.m
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    p1s, B_p1, B_rot, spectrum = parent_plane_split_exact(sys_lr, v)
    gr0 = gram_reduce(sys_lr, v)
    _, B0, _, _, _ = dominant_plane(sys_lr, gr0)
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z0 = v + DELTA * g; Z0 = Z0 / np.linalg.norm(Z0)
    # 対照ゲート(b): 静的親 npz と states[0] に bit 一致
    par = np.load(PARENTS / f'parent_static_N{n:05d}_makeparent_20260905.npz')
    assert np.array_equal(np.asarray(par['v']), v), f'PARENT GATE FAIL v N={n}'
    assert np.array_equal(np.asarray(par['Z0']), Z0), f'PARENT GATE FAIL Z0 N={n}'
    assert np.array_equal(Zs[0], Z0), f'STATE0 GATE FAIL N={n}'
    print(f'[gate] N={n}: npz SHA / v / Z0 / states[0] すべて一致', flush=True)

    ts, bands, fs = [], [], []
    f_prev = None
    crossing = None
    for t in range(0, XMAX + 1):
        Zr = Zs[t]
        totZ = float(np.real(np.conj(Zr) @ Zr))
        E_P1 = occ(B_p1, Zr)
        f = 1.0 - E_P1 / totZ
        if crossing is None and f > 0.05:
            crossing = t
        if t % SAMPLE != 0 and t != XMAX:
            continue
        E_other = occ(B_rot, Zr)
        E_ker = totZ - E_P1 - E_other
        gr = gram_reduce(sys_lr, Zr)
        _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
        e34 = s4_new_dirs(B0, Bdom)
        proj = B_rot @ (B_rot.T @ e34)
        fq, _ = np.linalg.qr(proj)
        f34 = fq[:, :2] if fq.shape[1] >= 2 else fq
        f34 = align_2d(f_prev, f34); f_prev = f34
        E_d3 = occ(f34[:, [0]], Zr)
        E_d4 = occ(f34[:, [1]], Zr) if f34.shape[1] > 1 else 0.0
        E_rem = max(0.0, E_other - E_d3 - E_d4)
        ts.append(t)
        bands.append([E_P1 / totZ, E_d3 / totZ, E_d4 / totZ, E_rem / totZ, E_ker / totZ])
        fs.append(f)
    return {'t': np.array(ts), 'bands': np.array(bands).T, 'f': np.array(fs),
            'crossing': crossing, 'M': M}

results = {}
meta = {}
for n in NS:
    d = run(n)
    results[n] = d
    meta[str(n)] = {'M': d['M'], 'crossing': d['crossing'],
                    'final_bands_P1_d3_d4_rem_ker': [float(x) for x in d['bands'][:, -1]],
                    'final_f': float(d['f'][-1])}
    with open(BASE / f'bands5_stage123_N{n:05d}.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['step', 'P1', 'dir3', 'dir4', 'remaining_other', 'kernel', 'f'])
        for i, t in enumerate(d['t']):
            w.writerow([int(t)] + [f'{x:.10e}' for x in d['bands'][:, i]] + [f'{d["f"][i]:.10e}'])
    print(f"N={n}: crossing={d['crossing']} 終端 P1/d3/d4/rem/ker="
          f"{[round(float(x),4) for x in d['bands'][:, -1]]} f={d['f'][-1]:.4f}", flush=True)

fig, axes = plt.subplots(len(NS), 1, figsize=(11, 3.6 * len(NS)), sharex=True, squeeze=False)
for ax, n in zip(axes[:, 0], NS):
    d = results[n]
    for band, c in zip(d['bands'], COLORS):
        ax.semilogy(d['t'], np.clip(band, 1e-6, None), lw=0.8, color=c)
    ax.semilogy(d['t'], np.clip(d['f'], 1e-6, None), 'k-', lw=0.9)
    if d['crossing'] is not None:
        ax.axvline(d['crossing'], color='k', ls=':', lw=0.8)
    ax.set_xlim(0, XMAX); ax.set_xticks(np.arange(0, XMAX + 1, 1000))
    ax.set_ylabel(f'N={n}')
axes[0, 0].legend(LABELS + ['f'], fontsize=6, loc='center right')
axes[-1, 0].set_xlabel('step (absolute)')
fig.suptitle('Figure3 compare (5-color, log) — common axis  [stage1+2+3, dt=2pi/N, 10000 steps]')
fig.savefig(BASE / 'fig3_compare_stage123_N5_N40.png', dpi=130)
fig.savefig(BASE / 'fig3_compare_stage123_N5_N40.svg')
plt.close(fig)
with open(BASE / 'make_fig3_5color_stage123_v1_meta.json', 'w') as fh:
    json.dump(meta, fh, indent=2)
print('ALL DONE')
