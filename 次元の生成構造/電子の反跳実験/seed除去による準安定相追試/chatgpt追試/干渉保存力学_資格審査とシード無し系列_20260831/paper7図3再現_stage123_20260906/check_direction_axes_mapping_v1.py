#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""状態空間の方向（P1/新平面/核）と複素平面の位相軸（60°3軸）の対応検定（読み出しのみ）。
仮説（木原 2026-09-06）: 5色読出しの「第3・第4方向」（新平面）が、検討7の
3次元空間読出し（60°3軸構造）の第3軸を運んでいるのではないか。
方法: 終状態 Z(10000) を Z = Z_P1 + Z_new + Z_rem + Z_ker に直交分解し、
各成分の波ごとの位相軸（arg(z_e²)/2、|z_e| 重み付き・微小成分除外）を別々に測る。
全体の3軸と、各成分の軸の対応（P1 が2軸を、新平面が第3軸を運ぶか）を機械判定。
対象: N=3, 4（結晶）＋N=40（ガラス参考）。基底は実験8と同一（原本 import・ゲート済み）。
出力: check_direction_axes_mapping_v1.json"""
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[5]
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

ledger = {}
with open(SWEEP / "SHA256SUMS.txt") as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

def s4_new_dirs(B0, Bdom):
    R = Bdom - B0 @ (B0.T @ Bdom)
    Qr, _ = np.linalg.qr(R)
    return Qr[:, :2]

def axes_of(z, amp_floor_rel=1e-6):
    """位相軸（deg, mod 180）とその重み（Σ|z|²比）。微小成分は除外。"""
    a = np.abs(z)
    if a.max() == 0:
        return []
    keep = a > amp_floor_rel * a.max()
    zz = z[keep]
    ph2 = np.angle(zz ** 2)
    order = np.argsort(ph2)
    ps = ph2[order]
    gaps = np.diff(np.concatenate([ps, [ps[0] + 2 * math.pi]]))
    cut = np.flatnonzero(gaps > 0.15)
    groups = []
    if cut.size == 0:
        groups = [list(range(len(zz)))]
    else:
        start = (cut[-1] + 1) % len(ps)
        idx = list(order[start:]) + list(order[:start])
        bounds = sorted(((c - start) % len(ps)) for c in cut)
        prev = 0
        for b in bounds:
            groups.append(idx[prev:b + 1])
            prev = b + 1
    tot = float(np.sum(np.abs(zz) ** 2))
    out = []
    for mem in groups:
        c = float(np.angle(np.sum(np.exp(1j * ph2[mem]))))
        w = float(np.sum(np.abs(zz[mem]) ** 2)) / tot
        out.append((round(math.degrees(c) / 2.0 % 180.0, 2), round(w, 4), len(mem)))
    return sorted(out)

out = {}
for n in (3, 4, 40):
    rel = f'results/hm_N{n}_den_{n}_states_10000.npz'
    h = hashlib.sha256(open(SWEEP / rel, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    Z = np.asarray(np.load(SWEEP / rel)['Z'][10000], dtype=np.complex128)
    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    par = np.load(PARENTS / f'parent_static_N{n:05d}_makeparent_20260905.npz')
    assert np.array_equal(np.asarray(par['v']), v), f'PARENT GATE FAIL N={n}'
    p1s, B_p1, B_rot, spectrum = parent_plane_split_exact(sys_lr, v)
    gr0 = gram_reduce(sys_lr, v)
    _, B0, _, _, _ = dominant_plane(sys_lr, gr0)
    Z_P1 = B_p1 @ (B_p1.T @ Z)
    if B_rot is not None and B_rot.shape[1] > 0:
        gr = gram_reduce(sys_lr, Z)
        _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
        e34 = s4_new_dirs(B0, Bdom)
        proj = B_rot @ (B_rot.T @ e34)
        fq, _ = np.linalg.qr(proj)
        f34 = fq[:, :2] if fq.shape[1] >= 2 else fq
        Z_other = B_rot @ (B_rot.T @ Z)
        Z_new = f34 @ (f34.T @ Z)
        Z_rem = Z_other - Z_new
    else:
        # N=3: other 回転空間は空（M=3 = P1 2次元＋核1次元）。垂直成長は核に住む
        Z_other = np.zeros_like(Z); Z_new = np.zeros_like(Z); Z_rem = np.zeros_like(Z)
    Z_ker = Z - Z_P1 - Z_other
    comp = {'total': Z, 'P1': Z_P1, 'new_plane_d3d4': Z_new, 'remaining_other': Z_rem,
            'kernel': Z_ker}
    res = {}
    for name, zc in comp.items():
        e = float(np.sum(np.abs(zc) ** 2))
        res[name] = {'energy_frac': round(e / float(np.sum(np.abs(Z) ** 2)), 6),
                     'axes_deg_weight_count': axes_of(zc)}
    out[str(n)] = res
    print(f"=== N={n} ===")
    for name in ('total', 'P1', 'new_plane_d3d4', 'kernel', 'remaining_other'):
        r = res[name]
        print(f"  {name}: E={r['energy_frac']:.4f} axes={r['axes_deg_weight_count'][:6]}")

with open(BASE / 'check_direction_axes_mapping_v1.json', 'w') as fh:
    json.dump(out, fh, indent=2, default=str)
print('ALL DONE')
