#!/usr/bin/env python3
"""make_parent（倍音対応）単体テスト v2 — v1 の反証を受けた訂正版。

v1 からの変更（v1 の結果 test_result_v1.json は反証記録として保存済み）:
  (1) T8 の訂正: v1 は「固有値枝の反転＝カイラリティ反転（μ 符号反転）」を
      予想したが反証された（両枝とも μ<0。枝は収束先の不動点族を変えるだけ）。
      鏡像の正しい定義は複素共役 v→v̄ であり、代数的に μ(v̄)=μ(v)（不変）。
      T8v2 = 共役状態が全構造検査を同値通過し μ が不変であることの検証。
  (2) 収束堅牢化: restarts=10, iters=2000（v1 は 3/1200 で一部盆地が非収束）。
      設計変更ではなくパラメータの増強。
  (3) 追加観測: 段ごとの σ₁ を不動点族ラベルとして記録
      （v1 実測: N=5 で √14≈3.7417 と 4.000=N−1、N=40 で 38.867/38.836/39.000=N−1
       ——自己無撞着閉包は離散族を成す）。
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from make_parent_harmonic_v1 import (make_parent_harmonic, self_consistent_mode,
                                       _eigenmode_mu_residual, eng, ENGINE_SHA256)

HERE = Path(__file__).resolve().parent

ITERS = 2000
RESTARTS = 10
TOL = 1e-12
TOL_RES = 1e-10
TOL_CLO = 1e-12
TOL_NRM = 1e-12
TOL_CIR = 1e-10
MIN_AMP = 1e-3
TOL_OCC = 1e-8
TOL_MU_INV = 1e-9


def check_case(n, H, seed):
    Z, info = make_parent_harmonic(n, H, seed, iters=ITERS, restarts=RESTARTS, tol=TOL)
    lv = info["levels"]
    t = {}
    t["T1_residual"] = all(l["residual"] < TOL_RES for l in lv)
    t["T2_level_closure"] = all(l["closure_abs"] < TOL_CLO for l in lv)
    t["T3_total_closure"] = info["total_closure_abs"] < TOL_CLO
    t["T4_norms"] = (abs(info["frobenius_norm"] - 1.0) < TOL_NRM
                     and all(abs(a - 1.0 / np.sqrt(H)) < TOL_NRM for a in info["level_amp"]))
    t["T5_circular"] = all(l["rank_re_im"] == 2
                            and abs(l["norm_re"] - l["norm_im"]) < TOL_CIR
                            and abs(l["re_dot_im"]) < TOL_CIR for l in lv)
    t["T6_all_nonzero"] = all(l["min_amp"] > MIN_AMP for l in lv)
    t["T7_occupancy"] = all(l["occupancy_top_plane"] > 1.0 - TOL_OCC for l in lv)
    return Z, info, t


def mu_of(n, v):
    sys_lr = eng.LowRankSystem(n)
    sys_lr.set_theta(np.angle(v))
    mu, res = _eigenmode_mu_residual(sys_lr, v)
    return mu, res


def main():
    t0 = time.time()
    print(f"make_parent（倍音対応）単体テスト v2  engine sha256={ENGINE_SHA256[:16]}…")
    print(f"  パラメータ: iters={ITERS} restarts={RESTARTS} tol={TOL}")
    results = {"engine_sha256": ENGINE_SHA256,
               "params": {"iters": ITERS, "restarts": RESTARTS, "tol": TOL},
               "v1_falsification_note": ("v1 T8（枝反転=μ符号反転）は反証。枝は不動点族の"
                                          "選択子でありカイラリティのつまみではない。鏡像=複素共役、μは共役不変")}
    all_pass = True

    for (n, H, seed) in ((5, 8, 40260801), (40, 4, 40260802)):
        label = f"N{n}_H{H}"
        Z, info, t = check_case(n, H, seed)
        print(f"\n[{label}] seed={seed}")
        for l in info["levels"]:
            print(f"  段n={l['n']}: branch={l['branch']} σ₁={l['sigma1']:.9f} μ={l['mu']:+.6f} "
                  f"残差={l['residual']:.2e} |vᵀv|={l['closure_abs']:.2e} 占有率={l['occupancy_top_plane']:.10f} "
                  f"再start={l['restarts_used']}")
        print(f"  総閉塞|ΣZ²|={info['total_closure_abs']:.2e} ‖Z‖={info['frobenius_norm']:.15f}")
        print(f"  不動点族（σ₁の値）: {sorted(set(round(l['sigma1'], 6) for l in info['levels']))}")
        for k, v in t.items():
            print(f"  {k}: {'PASS' if v else 'FAIL'}")
        all_pass &= all(t.values())

        # T8v2 鏡像=共役: conj(Z) の構造検査＋μ 共役不変
        t8_items = []
        for k, l in enumerate(info["levels"]):
            v = Z[:, k] * np.sqrt(H)
            vc = np.conj(v)
            mu_c, res_c = mu_of(n, vc)
            t8_items.append(abs(mu_c - l["mu"]) < TOL_MU_INV and res_c < TOL_RES
                            and abs(complex(vc @ vc)) < TOL_CLO)
        zc = np.conj(Z)
        t8 = all(t8_items) and abs(complex(np.sum(zc * zc))) < TOL_CLO
        print(f"  T8v2_mirror_conjugation（μ共役不変・共役状態全検査通過）: {'PASS' if t8 else 'FAIL'}")
        all_pass &= t8

        # T9 再現性
        Z2, _, _ = check_case(n, H, seed)
        t9_same = bool(np.array_equal(Z, Z2))
        _, info3, t3 = check_case(n, H, seed + 999)
        t9_diff = all(t3.values())
        t9 = t9_same and t9_diff
        print(f"  T9_reproducibility: {'PASS' if t9 else 'FAIL'} "
              f"(同seed bit一致={t9_same}, 異seed全通過={t9_diff})")
        all_pass &= t9

        results[label] = {"info": info, "tests": t, "T8v2": bool(t8),
                           "T9": {"pass": bool(t9), "same_seed_bitwise": t9_same,
                                   "diff_seed_all_pass": t9_diff},
                           "sigma_families": sorted(set(round(l["sigma1"], 6)
                                                         for l in info["levels"]))}

    results["all_pass"] = bool(all_pass)
    results["runtime_sec"] = time.time() - t0
    (HERE / "test_result_v2.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\n総合判定: {'ALL PASS' if all_pass else 'FAIL あり'}  ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
