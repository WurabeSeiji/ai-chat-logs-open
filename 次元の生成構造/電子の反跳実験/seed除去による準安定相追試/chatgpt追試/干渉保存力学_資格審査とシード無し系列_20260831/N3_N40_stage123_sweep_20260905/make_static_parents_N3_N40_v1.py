#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=3..40 の静的親データ生成（2026-09-05）。

本フォルダの検証済みエンジン run_n_scaling_lowrank_v1.py（7月正本の bit 同一コピー）から
make_parent / zero_closure_kernel_seed を import し、7月 N=40 走行と同一の手順・
同一の rng 式・同一の引数で各 N の親 v・零閉塞種 g・正規化初期状態 Z0 を生成する。

    rng = default_rng(40260722 + 1000*N + 0)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z0 = (v + 1e-15*g) / ||v + 1e-15*g||

保存先: parents/parent_static_N{N:05d}_makeparent_20260905.npz（v, g, Z0, sigma, residual, 条件）
既存ファイルがあれば上書きしない（スキップして検証のみ）。

合格ゲート:
  (1) N=40 の v・g・Z0 が既存正本
      自発的分裂予備実験_v1_N40対照実験系_20260904/largeN_splitting_result_v1/
      parent_static_N40_makeparent_20260904.npz と bit 一致（不一致なら exit 1）
  (2) 全 N の親残差を parents/parents_summary.csv に記録。residual >= 1e-8 の N は
      NOT_CONVERGED と明示（隠さない）。
"""
import csv
import os
import sys

import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed

PARENT_DIR = os.path.join(BASE_DIR, "parents")
REF40 = '/Users/kiharahanakira/Library/CloudStorage/GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/ai-chat-logs-open/次元の生成構造/自発的分裂予備実験_v1_N40対照実験系_20260904/largeN_splitting_result_v1/parent_static_N40_makeparent_20260904.npz'

SEED = 0
DELTA = 1e-15
TOL = 1e-12
ITERS = 1200

def main():
    os.makedirs(PARENT_DIR, exist_ok=True)
    rows = []
    gate_ok = True
    for N in range(3, 41):
        out_path = os.path.join(PARENT_DIR, f"parent_static_N{N:05d}_makeparent_20260905.npz")
        sys_lr = LowRankSystem(N)
        rng = np.random.default_rng(40260722 + 1000 * N + SEED)
        v, residual, sig = make_parent(sys_lr, rng, iters=ITERS, tol=TOL)
        g = zero_closure_kernel_seed(sys_lr, rng)
        Z = v + DELTA * g
        Z = Z / np.linalg.norm(Z)
        converged = bool(residual < 1e-8)
        status = "ok" if converged else "NOT_CONVERGED"

        if N == 40:
            ref = np.load(REF40)
            same = all(np.array_equal(x, ref[k]) for x, k in ((v, 'v'), (g, 'g'), (Z, 'Z0')))
            print(f"GATE N=40 v/g/Z0 bit-identical to canonical static parent: {same}")
            if not same:
                gate_ok = False
                status = "GATE_FAIL"

        if os.path.exists(out_path):
            prev = np.load(out_path)
            same_prev = all(np.array_equal(x, prev[k]) for x, k in ((v, 'v'), (g, 'g'), (Z, 'Z0')))
            print(f"N={N}: 既存ファイルあり（上書きせず検証のみ）一致={same_prev} residual={residual:.3e} {status}")
            if not same_prev:
                gate_ok = False
                status = "EXISTING_MISMATCH"
        else:
            np.savez_compressed(out_path, v=v, g=g, Z0=Z,
                                sigma=sig, residual=np.float64(residual),
                                n=np.int64(N), seed=np.int64(SEED),
                                delta=np.float64(DELTA), tol=np.float64(TOL),
                                iters=np.int64(ITERS))
            print(f"N={N}: saved M={sys_lr.m} residual={residual:.3e} planes={len(sig)} {status}")
        rows.append((N, sys_lr.m, f"{residual:.6e}", len(sig), status))

    with open(os.path.join(PARENT_DIR, "parents_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["N", "M", "parent_residual", "rank_planes", "status"])
        w.writerows(rows)
    n_bad = sum(1 for r in rows if r[4] != "ok")
    print(f"summary: {len(rows)} parents, {n_bad} non-ok")
    if not gate_ok:
        print("GATE FAIL")
        sys.exit(1)
    print("STATIC PARENTS DONE")

if __name__ == "__main__":
    main()
