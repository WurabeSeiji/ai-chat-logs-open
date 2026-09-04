#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=40 静的親データ生成（2026-09-04）。

本フォルダの検証済みエンジン run_n_scaling_lowrank_v1.py から make_parent /
zero_closure_kernel_seed を import し、正本走行 run_spontaneous_splitting_largeN_v1.py
（引数 40 1e-15 --after=1500 --tol=1e-12）と同一の手順・同一の rng 消費順で
親 v・零閉塞種 g・正規化初期状態 Z0 を生成し、別名の静的ファイルに保存する。

保存先（既存データは一切上書きしない）:
    largeN_splitting_result_v1/parent_static_N40_makeparent_20260904.npz
        v        : 親（make_parent 出力そのまま、正規化前）
        g        : 零閉塞核種（zero_closure_kernel_seed 出力そのまま）
        Z0       : (v + delta*g)/‖v + delta*g‖ （正本走行の step0 状態）
        sigma    : 親のσスペクトル、residual : 親残差
        n, seed, delta, tol, iters : 生成条件

合格ゲート（不合格なら exit 1）:
    Z0 が既存 largeN_splitting_result_v1/states_N00040_delta1e-15_seed0.npz の
    Z0 と bit 一致すること（=7月走行が使った初期値そのものである証明）。
"""
import os
import sys

import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent, zero_closure_kernel_seed

RESULT_DIR = os.path.join(BASE_DIR, "largeN_splitting_result_v1")
OUT_PATH = os.path.join(RESULT_DIR, "parent_static_N40_makeparent_20260904.npz")
REF_PATH = os.path.join(RESULT_DIR, "states_N00040_delta1e-15_seed0.npz")

N = 40
SEED = 0
DELTA = 1e-15
TOL = 1e-12
ITERS = 1200

def main():
    # 正本走行 run() 冒頭と同一の手順・同一の rng 消費順
    sys_lr = LowRankSystem(N)
    rng = np.random.default_rng(40260722 + 1000 * N + SEED)
    v, residual, sig = make_parent(sys_lr, rng, iters=ITERS, tol=TOL)
    g = zero_closure_kernel_seed(sys_lr, rng)
    Z = v + DELTA * g
    Z = Z / np.linalg.norm(Z)

    # 合格ゲート: 既存 states npz の Z0 と bit 一致
    ref = np.load(REF_PATH)
    same = bool(np.array_equal(Z, ref["Z0"]))
    print(f"parent residual = {residual:.15e}")
    print(f"GATE Z0 bit-identical to {os.path.basename(REF_PATH)}: {same}")
    if not same:
        print("GATE FAIL: 保存せず終了")
        sys.exit(1)

    np.savez_compressed(OUT_PATH, v=v, g=g, Z0=Z,
                        sigma=sig, residual=np.float64(residual),
                        n=np.int64(N), seed=np.int64(SEED),
                        delta=np.float64(DELTA), tol=np.float64(TOL),
                        iters=np.int64(ITERS))
    print(f"saved: {OUT_PATH}")
    print("STATIC PARENT DONE")

if __name__ == "__main__":
    main()
