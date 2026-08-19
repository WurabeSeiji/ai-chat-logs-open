#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL0 真空走行 v1 — 正単体定理の統一エンジン上の検証（判定 J1〜J8・実行前固定）

導出: DL0_導出ノート.md（真空正単体定理）。判定はノート §6 の J1〜J8。
力学: unified_interaction_v1（F は v1——実装決定 2026-08-19、真空湧きゼロを厳密判定するため。
      v2 分岐は対照走行として奇数帯湧きを記録する）。
読出し: x_e = Σ_{k,η} C2[e,k,η]（ノート §7 実装仕様）。

規格化の扱い: 親構成の Σ|x_e|² = c は走行で厳密保存される（単一スライス＋直交 Cayley）。
判定は c で正規化した量（λ·2M/c 等）に適用する（定理は規格化に対して同次）。

出力: result_dl0_vacuum_v1.json（判定表・時系列要約）・dl0_series_v1.npz（系列）
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNI = HERE.parent.parent / "統一万能関数_v1"

N = 16
M = N * (N - 1) // 2
T_STEPS = 20000
SAMPLE = 10
LATE_FRAC = 0.25  # 後期窓 = 末尾25%


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def analyze(x):
    """関係振幅 x (M,) → (c, bilin, λ降順, trB, top3occ, n_neg, V3)."""
    c = float(np.sum(np.abs(x) ** 2))
    bilin = complex(np.sum(x * x))
    d2 = np.abs(x) ** 2
    # 距離行列（M 本の辺 → N×N）
    D2 = np.zeros((N, N))
    idx = 0
    for i in range(N):
        for j in range(i + 1, N):
            D2[i, j] = D2[j, i] = d2[idx]
            idx += 1
    J = np.eye(N) - np.ones((N, N)) / N
    B = -0.5 * J @ D2 @ J
    lam, V = np.linalg.eigh(B)
    # 自明零＝固有ベクトルが定数ベクトル 1/√N に最も近いもの（降順末尾とは限らない。
    # 一般状態では負固有値が末尾に来るため、v1 初版の lam[:-1] はバグだった）
    ones = np.ones(N) / np.sqrt(N)
    i_triv = int(np.argmax(np.abs(V.T @ ones)))
    keep = np.ones(N, bool); keep[i_triv] = False
    lam, V = lam[keep], V[:, keep]
    order = np.argsort(lam)[::-1]
    lam_nt, V = lam[order], V[:, order]
    trB = float(np.sum(lam_nt))
    top3 = float(np.sum(lam_nt[:3]) / trB) if trB > 0 else float("nan")
    n_neg = int(np.sum(lam_nt < -1e-12 * max(trB, 1e-300)))
    return c, bilin, lam_nt, trB, top3, n_neg, V[:, :3]


def main():
    t0 = time.time()
    u1 = _load("uni_dl0_v1", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, 0.0)

    frames = {"tau": [], "c": [], "bilin_abs": [], "trB": [], "top3": [],
              "n_neg": [], "lam_min": [], "lam_max": [], "lam_med": [],
              "V3": []}
    LAG_FR = 12  # ラグ 12 フレーム = 120 步（v4 のラグ124步に対応する実装固定）
    for tau in range(T_STEPS):
        eng.step()
        if tau % SAMPLE:
            continue
        x = eng.C2().sum(axis=(1, 2))
        c, bilin, lam, trB, top3, n_neg, V3 = analyze(x)
        frames["tau"].append(tau)
        frames["c"].append(c)
        frames["bilin_abs"].append(abs(bilin))
        frames["trB"].append(trB)
        frames["top3"].append(top3)
        frames["n_neg"].append(n_neg)
        frames["lam_min"].append(float(lam.min()))
        frames["lam_max"].append(float(lam.max()))
        frames["lam_med"].append(float(np.median(lam)))
        frames["V3"].append(V3.copy())

    V3s = frames.pop("V3")
    overlaps = [float(np.sum((V3s[i].T @ V3s[i + LAG_FR]) ** 2) / 3.0)
                for i in range(len(V3s) - LAG_FR)]
    frames["overlap"] = overlaps + [float("nan")] * LAG_FR
    F = {k: np.array(v) for k, v in frames.items()}
    nf = len(F["tau"])
    late = slice(int(nf * (1 - LATE_FRAC)), nf)
    c_late = float(np.mean(F["c"][late]))

    # ---- 判定 J1〜J8（正規化: 理論 λ = c/(2M), trB = c/N）----
    lam_norm_med = F["lam_med"][late] * (2 * M) / F["c"][late]
    lam_ratio = (F["lam_max"][late] / F["lam_min"][late]) - 1.0
    J1 = {"lam_maxmin_minus1_late_mean": float(np.mean(lam_ratio)),
          "lam_med_relerr_late_mean": float(np.mean(np.abs(lam_norm_med - 1.0))),
          "pass": bool(np.mean(lam_ratio) < 1e-2
                       and np.mean(np.abs(lam_norm_med - 1.0)) < 1e-2)}
    trB_rel = np.abs(F["trB"] * N / F["c"] - 1.0)
    J2 = {"trB_relerr_max": float(trB_rel.max()), "pass": bool(trB_rel.max() < 1e-6)}
    t3l = F["top3"][late]
    J3 = {"top3_late_mean": float(t3l.mean()),
          "pass": bool(abs(t3l.mean() - 0.200) < 0.005)}
    J4 = {"n_neg_late_max": int(F["n_neg"][late].max()),
          "pass": bool(F["n_neg"][late].max() == 0)}
    bil_rel = F["bilin_abs"] / F["c"]
    J5 = {"bilin_over_c_max": float(bil_rel.max()), "pass": bool(bil_rel.max() < 1e-12)}
    # r_rms² = trB/N = c/N² → r_rms·N/√c = 1
    rrms_rel = np.abs(np.sqrt(F["trB"] / N) * N / np.sqrt(F["c"]) - 1.0)
    J6 = {"rrms_relerr_max": float(rrms_rel.max()), "pass": bool(rrms_rel.max() < 1e-6)}
    J7 = {"top3_first": float(F["top3"][0]), "top3_late_mean": float(t3l.mean()),
          "pass": bool(abs(t3l.mean() - 0.200) < 0.005)}
    ovl = F["overlap"][late]
    ovl = ovl[~np.isnan(ovl)]
    J8 = {"overlap_late_mean": float(ovl.mean()),
          "pass": bool(abs(ovl.mean() - 0.200) < 0.05)}

    # ---- v2 対照: 真空湧き（奇数帯パワー）----
    u2 = _load("uni_dl0_v2", UNI / "unified_interaction_v2.py")
    eng2, _, _ = u2.build_standard_universe(N, 0.0)
    odd_leak_v2 = 0.0
    odd_leak_v1 = 0.0
    eng1b, _, _ = u1.build_standard_universe(N, 0.0)
    for _ in range(2000):
        eng2.step()
        eng1b.step()
        P2 = np.abs(eng2.C2()) ** 2
        odd_leak_v2 = max(odd_leak_v2, float(P2[:, eng2.odd_k, :].sum()))
        P1 = np.abs(eng1b.C2()) ** 2
        odd_leak_v1 = max(odd_leak_v1, float(P1[:, eng1b.odd_k, :].sum()))

    res = {
        "config": {"N": N, "M": M, "T": T_STEPS, "SAMPLE": SAMPLE,
                   "engine": "unified_interaction_v1（実装決定）",
                   "late_frac": LATE_FRAC, "c_late": c_late},
        "J1_degeneracy": J1, "J2_trB": J2, "J3_top3": J3, "J4_imag_dirs": J4,
        "J5_bilinear": J5, "J6_rrms": J6, "J7_relaxation": J7, "J8_no_frame": J8,
        "v1_odd_leak_max_2000": odd_leak_v1,
        "v2_odd_leak_max_2000": odd_leak_v2,
        "elapsed_sec": time.time() - t0,
    }
    passes = [res[k]["pass"] for k in
              ("J1_degeneracy", "J2_trB", "J3_top3", "J4_imag_dirs",
               "J5_bilinear", "J6_rrms", "J7_relaxation", "J8_no_frame")]
    res["all_pass"] = bool(all(passes))

    np.savez_compressed(HERE / "dl0_series_v1.npz", **F)
    (HERE / "result_dl0_vacuum_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))

    print(f"DL0 真空走行 N={N} T={T_STEPS}  (v1 engine)  c={c_late:.6f}")
    for k in ("J1_degeneracy", "J2_trB", "J3_top3", "J4_imag_dirs",
              "J5_bilinear", "J6_rrms", "J7_relaxation", "J8_no_frame"):
        print(f"  {k}: {'PASS' if res[k]['pass'] else 'FAIL'}  {res[k]}")
    print(f"  v1 奇数帯湧き(2000步max) = {odd_leak_v1:.3e}   "
          f"v2 = {odd_leak_v2:.3e}")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'}  "
          f"({res['elapsed_sec']:.1f}s)")


if __name__ == "__main__":
    main()
