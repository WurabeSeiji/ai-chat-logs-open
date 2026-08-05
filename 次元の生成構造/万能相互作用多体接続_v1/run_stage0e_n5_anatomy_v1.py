#!/usr/bin/env python3
"""GATE-0e: N=5安定解の解剖——黄金比構造とKAM共鳴保護仮説の検定

目的（2026-08-05 木原氏指示「なぜN=5が特殊か」の検証）:
    仮説: 5は「素数×二次無理角（badly approximable）」を満たす唯一の位数であり、
    点火（凝縮体の不安定化）は内部モード間の共鳴現象、五重レジスタは
    黄金比の位相関係によりKAM的に共鳴から最大限保護される→唯一の安定真空。
    さらに五角形の黄金自己相似＝スケール流の固定点（走らない真空）。

予言（実行前固定・事後変更禁止）:
    P-A（黄金比構造）: K₅の辺は D₅ の下で二軌道（五角形辺: 距離1の対、
        対角辺: 距離2の対）に分かれ、正五角形の対角/辺比は φ=(1+√5)/2。
        安定解の振幅は二クラス構造を持ち、クラス比が φ または 1/φ の
        2%以内に入る。burst解では入らない。
    P-B（共鳴距離）: 接線写像の固有周波数の対比 ω_i/ω_j について、
        最近接有理数（分母≤10）距離の中央値が、安定解 > burst解。
    P-C（全N横断）: N=3..12 の burst親の λ_max は、固有周波数スペクトルの
        共鳴近接度（最近接有理数距離の最小値の逆数的指標）と正相関する。
    探索的記録: 位相間隔の 2π/5 構造、振幅階層の φ 冪、固有値の分布。

方法: gen3.make_parent(5, seed) を GATE-0d の分類済みシードで再構成
    （安定: 2,12,20 / burst: 0dのJSONから取得）。接線写像はONS-2と同一の
    中心差分（h=1e-7、位相整合、wp基点固定）。N=3..12は各Nの最初の収束
    burstシード1本。全てO(小)——数分で完了。

使い方: python3 run_stage0e_n5_anatomy_v1.py
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SERIES = HERE.parent
ABL = SERIES / "第8論文_二段階seed除去による準安定相の因果分離" / "code" / "run_preliminary_seed_ablation_v1.py"
GEN3 = SERIES / "make_parent_white_managed_v1" / "make_parent_white_harmonics_n_only_v3.py"
PHI = (1 + math.sqrt(5)) / 2
H_FD = 1e-7

spec = importlib.util.spec_from_file_location("abl_e", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("gen_e", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def make_step(n, wp_base):
    sys_lr = abl.LowRankSystem(n)

    def F(Z):
        wp = wp_base.copy()
        sys_lr.set_theta(np.angle(Z))
        se, _ = sys_lr.sigma_max_power(wp)
        return sys_lr.cayley_step(Z.copy(), se)

    return F


def jacobian_eigs(F, Z0):
    m = Z0.shape[0]
    G0 = F(Z0)
    ip = np.conj(Z0) @ G0
    ph = ip / abs(ip)
    cols = [np.eye(m)[i] + 0j for i in range(m)] + [1j * np.eye(m)[i] for i in range(m)]
    J = np.zeros((2 * m, 2 * m))
    for idx, d in enumerate(cols):
        Gp = np.conj(ph) * F(Z0 + H_FD * d)
        Gm = np.conj(ph) * F(Z0 - H_FD * d)
        gcol = (Gp - Gm) / (2 * H_FD)
        J[:, idx] = np.concatenate([gcol.real, gcol.imag])
    return np.linalg.eigvals(J)


def nearest_rational_dist(x, qmax=10):
    best = 1e9
    for q in range(1, qmax + 1):
        p = round(x * q)
        best = min(best, abs(x - p / q))
    return best


def edge_classes_n5():
    ia, ib = np.triu_indices(5, k=1)
    side, diag = [], []
    for idx, (a, b) in enumerate(zip(ia, ib)):
        d = min((b - a) % 5, (a - b) % 5)
        (side if d == 1 else diag).append(idx)
    return side, diag


def analyze_amplitudes_n5(v):
    side, diag = edge_classes_n5()
    amps = np.abs(v)
    ms, md = float(np.mean(amps[side])), float(np.mean(amps[diag]))
    ratio = max(ms, md) / min(ms, md) if min(ms, md) > 0 else np.inf
    return {"mean_side": ms, "mean_diag": md, "class_ratio": ratio,
            "phi_hit": bool(abs(ratio - PHI) < 0.02 * PHI or abs(ratio - 1) < 1e-9 and False),
            "dist_to_phi": float(abs(ratio - PHI)),
            "amps_sorted": sorted(float(a) for a in amps),
            "cv_side": float(np.std(amps[side]) / ms) if ms > 0 else None,
            "cv_diag": float(np.std(amps[diag]) / md) if md > 0 else None}


def freq_ratio_stats(eigs, unstable_floor=1e-3):
    mods = np.abs(eigs)
    lam_max = float(np.log(mods.max()))
    # 中立殻（|μ|~1）の回転周波数
    neutral = eigs[np.abs(np.log(np.clip(mods, 1e-12, None))) < 1e-4]
    ws = sorted({round(abs(float(np.angle(e))), 10) for e in neutral if abs(np.angle(e)) > 1e-6})
    dists = []
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            r = ws[i] / ws[j]
            dists.append(nearest_rational_dist(r))
    return {"lambda_max": lam_max,
            "n_unstable": int(np.sum(np.log(mods) > unstable_floor)),
            "n_neutral_freqs": len(ws),
            "freqs": [float(w) for w in ws[:12]],
            "ratio_dists_median": float(np.median(dists)) if dists else None,
            "ratio_dists_min": float(np.min(dists)) if dists else None}


def main() -> None:
    t0 = time.time()
    d0d = json.load(open(HERE / "gate0d_parent_seed_sweep_v1.json"))
    n5 = d0d["5"]
    stable_seeds = [s for s, c, _ in n5 if c == "stable"]
    burst_seeds = [s for s, c, _ in n5 if c == "burst"]
    print(f"GATE-0e N=5解剖  安定={stable_seeds} burst={burst_seeds}")

    results = {"N5": {}, "crossN": {}}
    m5 = 10
    for seed in stable_seeds + burst_seeds:
        r = gen3.make_parent(5, seed=seed)
        v = r.parent_vector / np.linalg.norm(r.parent_vector)
        wp = np.random.default_rng(91000).normal(size=m5)
        F = make_step(5, wp)
        eigs = jacobian_eigs(F, v)
        amp = analyze_amplitudes_n5(v)
        fr = freq_ratio_stats(eigs)
        kind = "stable" if seed in stable_seeds else "burst"
        results["N5"][str(seed)] = {"kind": kind, "amp": amp, "spec": fr}
        print(f"  seed={seed:2d} [{kind}] クラス比={amp['class_ratio']:.4f} "
              f"(φとの距離 {amp['dist_to_phi']:.4f}) λ_max={fr['lambda_max']:.5f} "
              f"不安定数={fr['n_unstable']} 比距離中央値={fr['ratio_dists_median']}")

    # P-A判定
    sa = [results["N5"][str(s)]["amp"] for s in stable_seeds]
    ba = [results["N5"][str(s)]["amp"] for s in burst_seeds]
    pa_stable = [a["dist_to_phi"] < 0.02 * PHI for a in sa]
    pa_burst = [a["dist_to_phi"] < 0.02 * PHI for a in ba]
    # P-B判定
    sb = [results["N5"][str(s)]["spec"]["ratio_dists_median"] for s in stable_seeds]
    bb = [results["N5"][str(s)]["spec"]["ratio_dists_median"] for s in burst_seeds]
    sb = [x for x in sb if x is not None]
    bb = [x for x in bb if x is not None]
    pb = bool(sb and bb and np.median(sb) > np.median(bb))

    # P-C: N=3..12 の burst 親 1 本ずつ
    for n in range(3, 13):
        rows = d0d[str(n)]
        seed = next((s for s, c, _ in rows if c == "burst"), None)
        if seed is None:
            continue
        r = gen3.make_parent(n, seed=seed)
        v = r.parent_vector / np.linalg.norm(r.parent_vector)
        m = n * (n - 1) // 2
        wp = np.random.default_rng(91000).normal(size=m)
        F = make_step(n, wp)
        eigs = jacobian_eigs(F, v)
        fr = freq_ratio_stats(eigs)
        results["crossN"][str(n)] = {"seed": seed, **fr}
        print(f"  N={n:2d} burst(seed={seed}): λ_max={fr['lambda_max']:.5f} "
              f"比距離min={fr['ratio_dists_min']}")

    xs = [(v["ratio_dists_min"], v["lambda_max"]) for v in results["crossN"].values()
          if v["ratio_dists_min"] is not None]
    corr = None
    if len(xs) >= 3:
        a = np.array(xs)
        # 共鳴近接度 = 1/比距離min。相関は順位相関（簡易）
        prox = 1.0 / np.clip(a[:, 0], 1e-9, None)
        lam = a[:, 1]
        rp = np.argsort(np.argsort(prox)).astype(float)
        rl = np.argsort(np.argsort(lam)).astype(float)
        corr = float(np.corrcoef(rp, rl)[0, 1])

    verdict = {"P_A_stable_phi_hits": pa_stable, "P_A_burst_phi_hits": pa_burst,
               "P_A": bool(all(pa_stable) and not any(pa_burst)) if sa else None,
               "P_B_stable_median": float(np.median(sb)) if sb else None,
               "P_B_burst_median": float(np.median(bb)) if bb else None,
               "P_B": pb,
               "P_C_rank_corr": corr,
               "P_C": bool(corr is not None and corr > 0)}
    print(f"\n判定: P-A={verdict['P_A']} (stable hits={pa_stable}, burst hits={pa_burst})")
    print(f"      P-B={verdict['P_B']} (stable med={verdict['P_B_stable_median']}, "
          f"burst med={verdict['P_B_burst_median']})")
    print(f"      P-C={verdict['P_C']} (順位相関={corr})")

    out = {"PHI": PHI, "criteria": {"P_A": "stable class ratio within 2% of phi, burst not",
                                     "P_B": "median nearest-rational distance stable > burst",
                                     "P_C": "rank corr(1/min_dist, lambda_max) > 0"},
           "verdict": verdict, "results": results,
           "runtime_sec": time.time() - t0}
    (HERE / "gate0e_n5_anatomy_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
