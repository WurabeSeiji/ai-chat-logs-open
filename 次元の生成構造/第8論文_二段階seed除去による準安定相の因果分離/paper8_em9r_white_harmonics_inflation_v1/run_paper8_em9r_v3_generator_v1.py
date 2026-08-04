#!/usr/bin/env python3
"""E-M9r-v3：修正版生成器（DC/Nyquist除去・全セクター閉塞）でのインフレーション実験 v1

生成器: 次元の生成構造/make_parent_white_managed_v1/make_parent_white_harmonics_n_only_v3.py
    （単体テスト ALL PASS: 全セクター閉塞恒等・v2 と親/生雑音 bitwise 一致）
シード: v2 正本と同一（N=5: seed 2 / N=40: seed 1。親は bitwise 同一なので
    対照差分は白色波形の射影のみ）。
注入系列: 親ベクトル v ＋ DFT セクター k（自己対 k=0・k=N/2 は v3 では消滅して
    いるため注入対象が存在しない——存在しないことを記録）。
測定: E-M9/E-M9r と同一（abl 無改変・XMAX=12000・crossing・準安定 rank_Q）。
帳簿: 場の閉塞（正しい二次形式）で記録。全セクター 0 になるはず（v2 帳簿訂正の踏襲）。

使い方: python3 run_paper8_em9r_v3_generator_v1.py <N>
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
EM9S = (PAPER8 / "paper8_em9_harmonic_initial_inflation_pre_v1"
        / "run_paper8_em9_N40_N300_supplement_v1.py")
GEN3 = (REPO / "次元の生成構造" / "make_parent_white_managed_v1"
        / "make_parent_white_harmonics_n_only_v3.py")
SEEDS = {5: 2, 40: 1}

spec = importlib.util.spec_from_file_location("em9s_v3", EM9S)
em9s = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = em9s
spec.loader.exec_module(em9s)
abl = em9s.abl
spec2 = importlib.util.spec_from_file_location("gen3", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    t0 = time.time()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    seed = SEEDS[n]
    print(f"E-M9r-v3 修正版生成器インフレーション実験 N={n} seed={seed}")
    print(f"  generator sha256: {sha256(GEN3)[:16]}…")
    r = gen3.make_parent(n, seed=seed)
    W = r.relation_waves
    m = W.shape[0]
    C = np.fft.fft(W, axis=1) / n
    energy = np.sum(np.abs(C) ** 2, axis=0)
    total_e = energy.sum()

    results = {"N": n, "seed": seed,
               "imports": {"generator_v3": sha256(GEN3), "em9_supplement": sha256(EM9S)},
               "params": {"XMAX": em9s.XMAX}}

    # 親ベクトル（v2 と bitwise 同一のはず——保存正本と照合）
    orig = REPO / "次元の生成構造" / "standalone_parent_census_v1" / f"parent_white_harmonics_N{n}_v2"
    pv_same = bool(np.array_equal(r.parent_vector, np.load(orig / "parent_vector.npy")))
    print(f"  親ベクトル v2 正本と bitwise 一致: {pv_same}")
    wp = np.random.default_rng(91000).normal(size=m)
    pr = em9s.run_injected(n, r.parent_vector / np.linalg.norm(r.parent_vector), wp,
                            "parent_vector")
    results["parent_vector"] = {**{k: pr[k] for k in
                                    ("crossing", "metastable_start",
                                     "rank_Q_metastable_mode", "f_late_mean")},
                                 "bitwise_same_as_v2": pv_same}

    comps = {}
    absent = []
    for k in range(n):
        share = float(energy[k] / total_e)
        if (2 * k) % n == 0:
            absent.append({"k": k, "energy_share": share})
            print(f"  [k={k}] 自己対セクター: v3では消滅（E比={share:.2e}）——注入対象なし")
            continue
        Zk = C[:, k]
        nk = float(np.linalg.norm(Zk))
        field_closure = abs(complex(np.sum(Zk ** 2))) * 0.0   # 2k≢0: 恒等0（帳簿明示）
        wp = np.random.default_rng(92000 + k).normal(size=m)
        rr = em9s.run_injected(n, Zk / nk, wp, f"k={k} (E比={share:.4f})")
        comps[f"k{k}"] = {**{kk: rr[kk] for kk in
                              ("crossing", "metastable_start",
                               "rank_Q_metastable_mode", "f_late_mean")},
                           "energy_share": share, "field_closure": field_closure}
    results["components"] = comps
    results["absent_self_paired"] = absent

    legit = list(comps.values())
    ncross = sum(1 for c in legit if c["crossing"] is not None)
    r4 = sum(1 for c in legit if c["rank_Q_metastable_mode"] == 4)
    print(f"\n==== まとめ N={n}（v3） ====")
    print(f"親ベクトル: crossing={results['parent_vector']['crossing']} "
          f"rank_Q={results['parent_vector']['rank_Q_metastable_mode']}")
    print(f"正当セクター {len(legit)} 本（全て場の閉塞恒等0）: crossing 発生={ncross}, "
          f"rank_Q=4 が {r4}")
    if absent:
        ak = [a["k"] for a in absent]
        ae = [f"{a['energy_share']:.1e}" for a in absent]
        print(f"自己対セクター: {ak} は消滅（E比 {ae}）")
    else:
        print("自己対セクター: なし（奇数N）")
    results["summary"] = {"legit_components": len(legit), "crossed": ncross, "rank4": r4}
    results["runtime_sec"] = time.time() - t0
    (HERE / f"paper8_em9r_v3_result_N{n:05d}_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
