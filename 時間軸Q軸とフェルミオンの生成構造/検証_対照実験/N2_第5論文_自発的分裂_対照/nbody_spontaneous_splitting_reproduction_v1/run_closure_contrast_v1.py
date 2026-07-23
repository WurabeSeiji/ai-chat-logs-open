#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
実験O5: 閉鎖規約の対照——複素双線形閉鎖 Z^T Z の値は帰結クラスを変えるか

問い（問題1の検証）:
  力学（凍結解除・実直交Cayley）はエルミートノルム ‖Z‖² と複素双線形 Z^T Z を
  ともに厳密保存する。初期条件の Z^T Z が 0 か非零かは、帰結クラス
  （拡大／有界混合／回帰）を変えるのか。

対照4族（他は全て同一、分類計量はO4と同一）:
  genericZ0 : 一般状態、Z^T Z = 0 厳密（Z=X+iY, |X|=|Y|, X⊥Y）
  genericNZ : 一般状態、Z^T Z ≠ 0（複素ガウス、|Z^T Z| ~ O(1/√M)）
  parentZ0  : 固有モード親 + 零閉鎖核種 δ（g=(u+iw)/√2, g^T g=0）
  parentNZ  : 固有モード親 + 複素ガウス種 δ（g^T g ≠ 0、旧O4のparent構成）

判定:
  各 (N, 族) の帰結クラス分布を比較。genericZ0 と genericNZ、
  parentZ0 と parentNZ の間でクラスが系統的に分かれれば「閉鎖規約は物理を変える」、
  分かれなければ「Z^T Z の初期値は帰結クラスに影響しない（表示の違い）」。
"""

import json
import math
import os
import sys

import numpy as np

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE_DIR)
import run_spontaneous_splitting_preliminary_v1 as base
import run_outcome_classification_preliminary_v1 as oc

RESULT_DIR = os.path.join(BASE_DIR, "closure_contrast_result_v1")
N_LIST = [3, 4, 5, 6]
SEEDS = [0, 1, 2]
DELTA = 1e-3
STEPS = 20000
RECORD_EVERY = 2


def generic_z0(rng, m):
    """Z^T Z = 0 厳密の一般状態。"""
    X = rng.normal(size=m)
    Y = rng.normal(size=m)
    Y = Y - (X @ Y) / (X @ X) * X
    Y = Y * (np.linalg.norm(X) / np.linalg.norm(Y))
    Z = X + 1j * Y
    return Z / np.linalg.norm(Z)


def generic_nz(rng, m):
    """Z^T Z ≠ 0 の一般状態（旧O4 generic と同構成）。"""
    Z = rng.normal(size=m) + 1j * rng.normal(size=m)
    return Z / np.linalg.norm(Z)


def parent_nz(rng, A, m, delta):
    """固有モード親 + 複素ガウス種（g^T g ≠ 0、旧O4 parent と同構成）。"""
    return oc.prepare_parent_perturbed(rng, A, m, delta)


def parent_z0(rng, A, m, delta):
    """固有モード親 + 零閉鎖核種（g=(u+iw)/√2, g^T g = 0 厳密）。"""
    v, residual = base.prepare_initial_state(rng, A, m, 0.0)
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)

    def project_off(g):
        g = g - (g @ p) * p - (g @ q) * q
        return g

    u = project_off(rng.normal(size=m))
    u = u / np.linalg.norm(u)
    w = project_off(rng.normal(size=m))
    w = w - (w @ u) * u
    w = w / np.linalg.norm(w)
    g = (u + 1j * w) / math.sqrt(2.0)
    Z = v + delta * g
    return Z / np.linalg.norm(Z), residual


def main():
    os.makedirs(RESULT_DIR, exist_ok=True)
    results = []
    print("=== 実験O5: 閉鎖規約（Z^T Z）の対照 ===", flush=True)
    for n in N_LIST:
        A = base.line_graph_adjacency(n)
        m = A.shape[0]
        for family in ["genericZ0", "genericNZ", "parentZ0", "parentNZ"]:
            for seed in SEEDS:
                rng = np.random.default_rng(11220722 + 1000 * n + seed)
                residual = None
                if family == "genericZ0":
                    Z0 = generic_z0(rng, m)
                elif family == "genericNZ":
                    Z0 = generic_nz(rng, m)
                elif family == "parentZ0":
                    Z0, residual = parent_z0(rng, A, m, DELTA)
                else:
                    Z0, residual = parent_nz(rng, A, m, DELTA)
                ztz0 = abs(complex(Z0 @ Z0))
                f, s2, pr, dev_c = oc.run_and_measure(Z0, A, STEPS, RECORD_EVERY)
                metrics = oc.classify(f, s2, family)
                entry = {
                    "n": n, "m": m, "family": family, "seed": seed,
                    "abs_ztz_initial": ztz0,
                    "max_dev_closure": dev_c,
                    "parent_residual": residual,
                    **metrics,
                }
                results.append(entry)
                print(f"  N={n} {family:>9} seed={seed}: |Z^TZ|₀={ztz0:.2e}"
                      f" → {metrics['class']:>12}"
                      f" (f_min={metrics['f_late_min']:.3f},"
                      f" σ2act={metrics['sigma2_active_frac']:.2f})", flush=True)

    with open(os.path.join(RESULT_DIR, "summary_v1.json"), "w") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False, default=float)

    print("\n=== 集計: (N, 族) ごとのクラス分布 ===")
    hdr = f"{'N':>2} {'family':>9} | {'expansion':>9} {'intermediate':>12} {'recurrent':>9} | |Z^TZ|₀(med)"
    print(hdr)
    print("-" * len(hdr))
    for n in N_LIST:
        for family in ["genericZ0", "genericNZ", "parentZ0", "parentNZ"]:
            runs = [r for r in results if r["n"] == n and r["family"] == family]
            cnt = {"expansion": 0, "intermediate": 0, "recurrent": 0}
            for r in runs:
                cnt[r["class"]] += 1
            zmed = np.median([r["abs_ztz_initial"] for r in runs])
            print(f"{n:>2} {family:>9} | {cnt['expansion']:>9}"
                  f" {cnt['intermediate']:>12} {cnt['recurrent']:>9} | {zmed:.2e}")
    dev = max(r["max_dev_closure"] for r in results)
    print(f"\n閉鎖保存の最大偏差（全run）: {dev:.2e}")
    print(f"出力: {RESULT_DIR}")


if __name__ == "__main__":
    main()
