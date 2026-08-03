#!/usr/bin/env python3
"""E-M9 補足：N=40 / N=300（過去実験と同一条件）での倍音海初期値テスト v1

過去条件との対応:
    - N=40: E-M3 と同一（種なし条件A、abl.build_init(40, False)、crossing 実測 2011）
    - N=300: E-M3 補足と同一（abl.build_init(300, False)、crossing 実測 4849）
    - XMAX=12000（E-M1/E-M9 主測定と同一の観測窓）、SAMPLE は abl.SAMPLE 準拠
    - plane split は build_init 74-77行と同一の分岐（n≤40 exact / n>40 approx）

固定予言:
    P0-40:  対照 crossing = 2011 厳密再現（E-M3 実測）
    P0-300: 対照 crossing = 4849 厳密再現（E-M3 補足実測）
    P1/P2（開いた問い・E-M9 N=5 の族二分の追試）:
        破れ族（σ₁<N−1）→ crossing 発生・rank_Q 対照一致
        等振幅族（σ₁=N−1）→ crossing なし（f≈0 のまま）
    倍音海: N=40 は seed 40260802（破れ族2段＋等振幅族2段）、
            N=300 は seed 40260803（全4段が等振幅族＝単体実測）。
            N=300 では破れ族の追試は不可能（この種の盆地に出現しないため）で、
            等振幅族の安定性テストのみになる——これ自体が記録対象。

既存ファイルは無改変（read-only import、abl.run 96-173行の測定部を初期値引数化で複製）。
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
MPH = REPO / "次元の生成構造" / "make_parent_harmonic_unit_v1" / "make_parent_harmonic_v1.py"

spec = importlib.util.spec_from_file_location("abl_m9s", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)

spec2 = importlib.util.spec_from_file_location("mph_m9s", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)

XMAX = 12000
CASES = {40: {"H": 4, "seed": 40260802, "expected_control_crossing": 2011},
         300: {"H": 4, "seed": 40260803, "expected_control_crossing": 4849}}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run_injected(n, v0, wp, label):
    """abl.run()（96-173行）条件A測定部の複製（初期値引数化・plane split 分岐は build_init 74-77行と同一）。"""
    sample_ev = abl.SAMPLE[n]
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(v0))
    if n <= 40:
        p1s, B_p1, B_rot, spectrum = abl.parent_plane_split_exact(sys_lr, v0)
    else:
        p1s, B_p1, B_rot, smax, thr = abl.parent_plane_split_approx(sys_lr, v0, abl.SIG_REL)
    gr0 = abl.gram_reduce(sys_lr, v0)
    _, B0, _, _, _ = abl.dominant_plane(sys_lr, gr0)
    p = v0.real / np.linalg.norm(v0.real)
    q = v0.imag - (v0.imag @ p) * p
    q = q / np.linalg.norm(q)

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    Z = v0.copy()
    crossing = None
    ranks, fs_sampled = [], []
    t = 0
    while True:
        f = fval(Z)
        if crossing is None and f > 0.05:
            crossing = t
        if t % sample_ev == 0 or t == XMAX:
            gr = abl.gram_reduce(sys_lr, Z)
            _, Bdom, _, _, _ = abl.dominant_plane(sys_lr, gr)
            qs = abl.qsv4(B0, Bdom)
            rankQ = int(np.sum(qs > abl.Q_REL_TAU * qs[0]))
            ranks.append((t, rankQ)); fs_sampled.append((t, f))
        if t >= XMAX:
            break
        Z, wp = abl.evolve(sys_lr, Z, wp); t += 1

    meta_start = (crossing + abl.GUARD) if crossing is not None else None
    if meta_start is not None:
        meta_ranks = [r for (tt, r) in ranks if tt >= meta_start]
        rank_meta = Counter(meta_ranks).most_common(1)[0][0] if meta_ranks else None
    else:
        rank_meta = None
    f_late = float(np.mean([f for (tt, f) in fs_sampled if tt >= XMAX - 2000]))
    print(f"  [{label}] crossing={crossing} 準安定開始={meta_start} "
          f"rank_Q(準安定最頻値)={rank_meta} f(終盤平均)={f_late:.4f} "
          f"|Z·Z|={abs(complex(Z @ Z)):.1e}", flush=True)
    return {"crossing": crossing, "metastable_start": meta_start,
            "rank_Q_metastable_mode": rank_meta, "f_late_mean": f_late,
            "final_zero_square_abs": abs(complex(Z @ Z)),
            "rank_timeline": ranks}


def main() -> None:
    t0 = time.time()
    print("E-M9 補足（N=40 / N=300、過去実験と同一条件）実行", flush=True)
    print(f"  import: ABL sha256={sha256(ABL)[:16]}…  MPH sha256={sha256(MPH)[:16]}…", flush=True)
    results = {"imports": {"abl": sha256(ABL), "mph": sha256(MPH), "engine": mph.ENGINE_SHA256},
               "params": {"XMAX": XMAX, "cases": {str(k): v for k, v in CASES.items()}}}

    for n, cfg in CASES.items():
        print(f"\n===== N={n} =====", flush=True)
        print(f"[対照] abl.build_init({n}, initial_seed=False)", flush=True)
        sys_lr, v, B_p1, B_rot, B0, p, q, Z0, wp0 = abl.build_init(n, False)
        ctrl = run_injected(n, Z0, wp0.copy(), "control")
        p0 = ctrl["crossing"] == cfg["expected_control_crossing"]
        print(f"  P0 駆動検証（crossing={cfg['expected_control_crossing']} 再現）: "
              f"{'PASS' if p0 else 'FAIL'}", flush=True)

        print(f"[倍音海] make_parent_harmonic(N={n}, H={cfg['H']}, seed={cfg['seed']})", flush=True)
        Z, info = mph.make_parent_harmonic(n, cfg["H"], cfg["seed"],
                                            iters=2000, restarts=10, tol=1e-12)
        per_level = {}
        for h in range(cfg["H"]):
            lv = info["levels"][h]
            v0 = Z[:, h] * np.sqrt(cfg["H"])
            wp = np.random.default_rng(90000 + h).normal(size=len(v0))
            fam = "N-1" if abs(lv["sigma1"] - (n - 1)) < 1e-9 else "broken"
            r = run_injected(n, v0, wp, f"段n={h+1} σ₁={lv['sigma1']:.6f} 族={fam}")
            r["sigma1"] = lv["sigma1"]; r["family"] = fam
            per_level[f"n{h+1}"] = r

        fam_summary = {}
        for k, r in per_level.items():
            fam_summary.setdefault(r["family"], []).append(
                (r["crossing"], r["rank_Q_metastable_mode"]))
        print(f"  族別まとめ: { {f: v for f, v in fam_summary.items()} }", flush=True)
        results[f"N{n}"] = {"control": ctrl, "P0": bool(p0), "per_level": per_level,
                             "family_summary": {f: [list(x) for x in v]
                                                 for f, v in fam_summary.items()}}

    results["runtime_sec"] = time.time() - t0
    (HERE / "paper8_em9_N40_N300_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({results['runtime_sec']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
