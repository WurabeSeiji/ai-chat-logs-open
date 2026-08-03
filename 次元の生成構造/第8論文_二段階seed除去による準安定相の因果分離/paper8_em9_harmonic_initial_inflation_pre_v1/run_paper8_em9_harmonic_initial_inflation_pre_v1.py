#!/usr/bin/env python3
"""E-M9：倍音対応 make_parent の初期値によるインフレーション再現テスト v1

目的（2026-08-04 木原氏指示）:
    倍音対応 make_parent（単体テスト ALL PASS 済み）が生成する各段の閉包を、
    三方向が生まれた既存のインフレーションプログラム（第8論文 ablation 系）の
    初期値として注入し、種なしインフレーション（crossing）と方向数（rank_Q）が
    再現されるかをテストする。

方針（既存プログラムは無改変）:
    - abl = code/run_preliminary_seed_ablation_v1.py を read-only import。
      本スクリプトは abl.run()（96-196行）の測定部を「初期状態を引数化」した
      形で複製する。エンジン・測定関数（evolve / parent_plane_split_exact /
      gram_reduce / dominant_plane / occ / qsv4 / GUARD / Q_REL_TAU / SAMPLE）は
      全て abl 名前空間から呼び、いかなる既存ファイルも変更しない。
    - 初期値は make_parent_harmonic（次元の生成構造/make_parent_harmonic_unit_v1）
      の各段 v⁽ⁿ⁾（単位ノルム円偏波閉包）。段は現時点で独立に走らせる
      （段間結合の力学は未設計・範囲外）。周波数 n·ω₀ はメタデータ。

条件: 種なし（条件A相当）。XMAX=12000（E-M1 と同一の観測窓。abl.XMAX=55000 の
    短縮であり、crossing≈1166・準安定開始≈4166 を覆う）。N=5。

固定予言:
    P0（駆動検証・厳格）: 対照＝abl.build_init(5, False) の初期値で、
        crossing = 1166（E-M1 実測）を厳密再現する。不一致なら本駆動は無効。
    P1（開いた問い）: 倍音海の各段（√14族・4.0=N−1族の両方を含む）でも
        種なし crossing が起こるか。族ごとに記録する。
    P2（開いた問い）: 準安定窓（crossing+GUARD 以降）の rank_Q（方向数）が
        対照と同数になるか。族ごとに記録する。

再現性: 乱数は make_parent_harmonic(seed=40260801) と wp 用 rng(90000+段) のみ。
    import 元 SHA-256 を JSON に記録。
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

spec = importlib.util.spec_from_file_location("abl_m9", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)

spec2 = importlib.util.spec_from_file_location("mph_m9", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)

N = 5
H = 8
SEED = 40260801            # 単体テスト v2 と同一（√14族4段・4.0族4段が出る構成）
XMAX = 12000               # E-M1 と同一の観測窓
SAMPLE_EV = abl.SAMPLE[N]  # 25


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run_injected(v0, wp, label):
    """abl.run()（96-173行）の条件A測定部の複製。初期状態 v0 を引数化した以外は同一。"""
    sys_lr = abl.LowRankSystem(N)
    sys_lr.set_theta(np.angle(v0))
    p1s, B_p1, B_rot, spectrum = abl.parent_plane_split_exact(sys_lr, v0)   # build_init 75行と同一
    gr0 = abl.gram_reduce(sys_lr, v0)
    _, B0, _, _, _ = abl.dominant_plane(sys_lr, gr0)                        # build_init 78-79行と同一
    p = v0.real / np.linalg.norm(v0.real)
    q = v0.imag - (v0.imag @ p) * p
    q = q / np.linalg.norm(q)                                               # build_init 80-81行と同一

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
        if t % SAMPLE_EV == 0 or t == XMAX:
            gr = abl.gram_reduce(sys_lr, Z)
            _, Bdom, _, _, _ = abl.dominant_plane(sys_lr, gr)               # run 152行と同一
            qs = abl.qsv4(B0, Bdom)
            rankQ = int(np.sum(qs > abl.Q_REL_TAU * qs[0]))                 # run 154-155行と同一
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
          f"|Z·Z|={abs(complex(Z @ Z)):.1e} ‖Z‖²誤差={abs(float(np.real(np.conj(Z) @ Z)) - 1):.1e}")
    return {"crossing": crossing, "metastable_start": meta_start,
            "rank_Q_metastable_mode": rank_meta,
            "rank_timeline": ranks, "f_late_mean": f_late,
            "final_zero_square_abs": abs(complex(Z @ Z)),
            "final_norm_error": abs(float(np.real(np.conj(Z) @ Z)) - 1)}


def main() -> None:
    t0 = time.time()
    print("E-M9 倍音海初期値インフレーションテスト 実行")
    print(f"  import: ABL sha256={sha256(ABL)[:16]}…")
    print(f"  import: MPH sha256={sha256(MPH)[:16]}…")

    results = {"imports": {"abl": sha256(ABL), "mph": sha256(MPH),
                            "engine": mph.ENGINE_SHA256},
               "params": {"N": N, "H": H, "SEED": SEED, "XMAX": XMAX,
                           "sample_every": SAMPLE_EV}}

    # ---- P0: 対照（既存 build_init の初期値・条件A）----
    print("\n[対照] abl.build_init(5, initial_seed=False)（E-M1 と同一初期値）")
    sys_lr, v, B_p1, B_rot, B0, p, q, Z0, wp0 = abl.build_init(N, False)
    ctrl = run_injected(Z0, wp0.copy(), "control")
    p0 = ctrl["crossing"] == 1166
    print(f"  P0 駆動検証（crossing=1166 厳密再現）: {'PASS' if p0 else 'FAIL'}")
    results["control"] = ctrl
    results["P0_driver_validated"] = bool(p0)

    # ---- 倍音海の各段 ----
    print(f"\n[倍音海] make_parent_harmonic(N={N}, H={H}, seed={SEED})")
    Z, info = mph.make_parent_harmonic(N, H, SEED, iters=2000, restarts=10, tol=1e-12)
    per_level = {}
    for h in range(H):
        lv = info["levels"][h]
        v0 = Z[:, h] * np.sqrt(H)
        wp = np.random.default_rng(90000 + h).normal(size=len(v0))
        r = run_injected(v0, wp, f"段n={h+1} σ₁={lv['sigma1']:.6f}")
        r["sigma1"] = lv["sigma1"]
        r["family"] = "N-1" if abs(lv["sigma1"] - (N - 1)) < 1e-9 else "sqrt14"
        per_level[f"n{h+1}"] = r

    # ---- 判定 ----
    fam = {}
    for k, r in per_level.items():
        fam.setdefault(r["family"], []).append(r)
    print("\n==== 族別まとめ ====")
    for fname, rs in fam.items():
        cr = [r["crossing"] for r in rs]
        rk = [r["rank_Q_metastable_mode"] for r in rs]
        print(f"  族 {fname}（{len(rs)}段）: crossing={cr} rank_Q準安定={rk}")
    p1 = all(r["crossing"] is not None for r in per_level.values())
    p2 = all(r["rank_Q_metastable_mode"] == ctrl["rank_Q_metastable_mode"]
             for r in per_level.values() if r["rank_Q_metastable_mode"] is not None) \
         and all(r["rank_Q_metastable_mode"] is not None for r in per_level.values())
    print(f"\nP0 駆動検証: {'PASS' if p0 else 'FAIL'}")
    print(f"P1 全段で種なし crossing 発生: {'YES' if p1 else 'NO'}")
    print(f"P2 全段の準安定 rank_Q = 対照({ctrl['rank_Q_metastable_mode']}): {'YES' if p2 else 'NO'}")

    results["per_level"] = per_level
    results["verdicts"] = {"P0": bool(p0), "P1_all_crossing": bool(p1),
                            "P2_rank_match": bool(p2)}
    results["runtime_sec"] = time.time() - t0
    (HERE / "paper8_em9_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
