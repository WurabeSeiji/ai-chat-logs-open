#!/usr/bin/env python3
"""make_parent（倍音対応）単体テスト v1 — 設計書 §5 の T1〜T9 を機械判定。

対象: (N=5, H=8) 主検査、(N=40, H=4) 実機スケール確認。
全結果を test_result_v1.json に保存。インフレーション系への組込みは行わない。
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from make_parent_harmonic_v1 import make_parent_harmonic, ENGINE_SHA256

HERE = Path(__file__).resolve().parent

TOL_RES = 1e-10       # T1
TOL_CLO = 1e-12       # T2, T3
TOL_NRM = 1e-12       # T4
TOL_CIR = 1e-10       # T5 等ノルム直交
MIN_AMP = 1e-3        # T6
TOL_OCC = 1e-8        # T7: 占有率 > 1 − 1e-8


def check_case(n, H, seed):
    Z, info = make_parent_harmonic(n, H, seed)
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


def main():
    t0 = time.time()
    print(f"make_parent（倍音対応）単体テスト v1  engine sha256={ENGINE_SHA256[:16]}…")
    results = {"engine_sha256": ENGINE_SHA256,
               "thresholds": {"T1": TOL_RES, "T2_T3": TOL_CLO, "T4": TOL_NRM,
                               "T5": TOL_CIR, "T6": MIN_AMP, "T7": TOL_OCC}}
    all_pass = True

    for (n, H, seed) in ((5, 8, 40260801), (40, 4, 40260802)):
        label = f"N{n}_H{H}"
        Z, info, t = check_case(n, H, seed)
        print(f"\n[{label}] seed={seed}")
        for l in info["levels"]:
            print(f"  段n={l['n']}: branch={l['branch']} μ={l['mu']:+.6f} σ₁={l['sigma1']:.6f} "
                  f"残差={l['residual']:.2e} |vᵀv|={l['closure_abs']:.2e} "
                  f"min|v|={l['min_amp']:.4f} rank={l['rank_re_im']} 占有率={l['occupancy_top_plane']:.12f}")
        print(f"  総閉塞|ΣZ²|={info['total_closure_abs']:.2e} ‖Z‖={info['frobenius_norm']:.15f}")
        for k, v in t.items():
            print(f"  {k}: {'PASS' if v else 'FAIL'}")
        all_pass &= all(t.values())

        # T8 鏡像: 全段 '−' と全段 '+' で μ 符号が段ごとに反転し、他は同水準
        Zm, im_ = make_parent_harmonic(n, H, seed, force_branch="-")
        Zp, ip_ = make_parent_harmonic(n, H, seed, force_branch="+")
        mu_m = [l["mu"] for l in im_["levels"]]
        mu_p = [l["mu"] for l in ip_["levels"]]
        t8_sign = all(a < 0 for a in mu_m) != all(a < 0 for a in mu_p) and \
                  (all(a < 0 for a in mu_m) or all(a > 0 for a in mu_m)) and \
                  (all(a < 0 for a in mu_p) or all(a > 0 for a in mu_p))
        t8_qual = (all(l["residual"] < TOL_RES and l["closure_abs"] < TOL_CLO
                       for l in im_["levels"] + ip_["levels"])
                   and im_["total_closure_abs"] < TOL_CLO
                   and ip_["total_closure_abs"] < TOL_CLO)
        t8 = t8_sign and t8_qual
        print(f"  T8_mirror: {'PASS' if t8 else 'FAIL'} "
              f"(μ['-']平均={np.mean(mu_m):+.4f}, μ['+']平均={np.mean(mu_p):+.4f})")
        all_pass &= t8

        # T9 再現性: 同 seed 再実行で bit 一致、異 seed で T1-T7 通過
        Z2, _, _ = check_case(n, H, seed)
        t9_same = bool(np.array_equal(Z, Z2))
        _, info3, t3 = check_case(n, H, seed + 999)
        t9_diff = all(t3.values())
        t9 = t9_same and t9_diff
        print(f"  T9_reproducibility: {'PASS' if t9 else 'FAIL'} "
              f"(同seed bit一致={t9_same}, 異seed全通過={t9_diff})")
        all_pass &= t9

        results[label] = {"info": info, "tests": t,
                           "T8": {"pass": bool(t8), "mu_minus": mu_m, "mu_plus": mu_p},
                           "T9": {"pass": bool(t9), "same_seed_bitwise": t9_same,
                                   "diff_seed_all_pass": t9_diff}}

    results["all_pass"] = bool(all_pass)
    results["runtime_sec"] = time.time() - t0
    (HERE / "test_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\n総合判定: {'ALL PASS' if all_pass else 'FAIL あり'}  ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
