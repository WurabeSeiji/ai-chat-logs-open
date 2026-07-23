#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N4 H12（セクター底の絶対安定性・寿命）検証実験（対照テスト済み第5論文コードを import）

H11 の修正版（超選択量＝活性平面数）を踏まえ、H12 を正面から検定する。

【H12 の予言（測定前固定）】
  (F1) 各セクター（保存量 c でラベル）内に「底」＝最小活性平面数の配置が存在し、
       それに達すると系はそこから出られない（絶対安定・吸収）。
  (F2) 底でない配置の寿命は有限で、下位に通約可能なチャネルの多寡に依存する。
  (F3) 底の活性平面数はセクター c に依存する（c=0 と c>0 で底が異なる）。

  合否:
    PASS(F1): 各逐次ランで n(t) が最小値 n_min に達した後、再上昇せず n_min に留まる
      割合（吸収率）が高い。→ 底の絶対安定。
    もし n_min から頻繁に再上昇するなら F1 は不成立＝H12 を見直す。

  本実験は長時間ラン（40000ステップ）で吸収の有無を直接見る。

実行: python3 sector_floor_lifetime_test_v1.py
"""

import importlib.util
import json
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "nbody_spontaneous_splitting_reproduction_v1",
                    "run_spontaneous_splitting_preliminary_v1.py")
_spec = importlib.util.spec_from_file_location("paper5_verified", _SRC)
p5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p5)

sine_generator = p5.sine_generator
cayley = p5.cayley
prepare_initial_state = p5.prepare_initial_state
line_graph_adjacency = p5.line_graph_adjacency


def active_planes(K, thr=0.05):
    return int(np.sum(np.linalg.eigvalsh(1j * K) > thr))


def run_floor(N, delta, seed, steps=40000, sub=20, thr=0.05):
    A = line_graph_adjacency(N)
    m = A.shape[0]
    rng = np.random.default_rng(20260721 + seed)
    Z, _ = prepare_initial_state(rng, A, m, delta)
    ztz0 = Z @ Z
    hn0 = float(np.real(np.conj(Z) @ Z))
    c0 = abs(ztz0) / hn0
    n_series = []
    for t in range(steps + 1):
        K = sine_generator(np.angle(Z), A)
        if t % sub == 0:
            n_series.append(active_planes(K, thr))
        if t < steps:
            Z = cayley(K) @ Z
    n_series = np.array(n_series)
    n_min = int(n_series.min())
    n_max = int(n_series.max())
    # 吸収検定: 最初に n_min に達した後、n_min より上へ戻る回数
    first_at_min = int(np.argmax(n_series == n_min))
    after = n_series[first_at_min:]
    reentries = int(np.sum(np.diff((after > n_min).astype(int)) == 1))
    frac_at_min_after = float(np.mean(after == n_min))
    # 最終域（後半）での平面数分布
    tail = n_series[len(n_series) // 2:]
    return {
        "N": N, "delta": delta, "seed": seed, "c0": c0,
        "n_initial": int(n_series[0]), "n_min": n_min, "n_max": n_max,
        "first_reach_min_step": first_at_min * sub,
        "reentries_above_min_after_first": reentries,
        "frac_time_at_min_after_first": frac_at_min_after,
        "tail_mean_n": float(tail.mean()),
        "tail_min": int(tail.min()), "tail_max": int(tail.max()),
        "absorbed": bool(reentries == 0),
        "_n_series": n_series.tolist(), "_sub": sub,
    }


if __name__ == "__main__":
    print(f"import 元: {_SRC}\n")
    print("=== H12(F1): セクター底の吸収検定（40000ステップ・逐次）===")
    print(f"{'N':>2} {'seed':>4} {'c0':>7} {'n初期':>5} {'n最小':>5} {'n最大':>5} "
          f"{'底到達step':>9} {'底以降再上昇':>10} {'底滞在率':>8} {'吸収?':>5}")
    results = []
    for N in (5, 6):
        for seed in range(4):
            r = run_floor(N, 1e-3, seed)
            results.append(r)
            print(f"{r['N']:>2} {r['seed']:>4} {r['c0']:>7.3f} {r['n_initial']:>5} "
                  f"{r['n_min']:>5} {r['n_max']:>5} {r['first_reach_min_step']:>9} "
                  f"{r['reentries_above_min_after_first']:>10} "
                  f"{r['frac_time_at_min_after_first']:>8.3f} "
                  f"{'YES' if r['absorbed'] else 'no':>5}")

    n_absorbed = sum(1 for r in results if r["absorbed"])
    print("\n=== 判定 ===")
    print(f"(F1) 底吸収したラン: {n_absorbed}/{len(results)}")
    print(f"     底到達後の平均再上昇回数: "
          f"{np.mean([r['reentries_above_min_after_first'] for r in results]):.1f}")
    print(f"     底滞在率の中央値: "
          f"{np.median([r['frac_time_at_min_after_first'] for r in results]):.3f}")
    if n_absorbed == len(results):
        print("⇒ F1 PASS（底は絶対安定・吸収）")
    elif n_absorbed == 0:
        print("⇒ F1 FAIL（底から頻繁に再上昇＝絶対安定な吸収底は無い。H12 を見直す）")
    else:
        print("⇒ F1 部分的（一部のランのみ吸収。セクター依存の可能性）")

    with open(os.path.join(_HERE, "sector_floor_result_v1.json"), "w") as f:
        slim = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
        json.dump({"import_source": _SRC, "runs": slim,
                   "n_absorbed": n_absorbed, "n_total": len(results)},
                  f, indent=2, ensure_ascii=False)
    print("\n結果を sector_floor_result_v1.json に保存")

    # --- 図: セクター底の吸収反証（永久揺らぎ）---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    outdir = os.path.join(_HERE, "figures")
    os.makedirs(outdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.0, 4.2))
    r = next(x for x in results if x["N"] == 6 and x["seed"] == 0)
    t = np.arange(len(r["_n_series"])) * r["_sub"]
    ax.step(t, r["_n_series"], where="post", color="tab:purple", lw=1.0)
    ax.axhline(r["n_min"], color="tab:red", ls="--", lw=1.5,
               label=f"minimum n = {r['n_min']} (NOT absorbing)")
    ax.set_xlabel("step  t (40000 total)")
    ax.set_ylabel("active plane count  n(t)")
    ax.set_yticks(range(0, 8))
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_title("Experiment 3 (REFUTED): minimum n is NOT absorbing; system settles above it "
                 f"(N=6; {r['reentries_above_min_after_first']} re-entries; "
                 f"{r['frac_time_at_min_after_first']*100:.0f}% time at min)")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "N4_fig3_no_absorbing_floor.svg"))
    fig.savefig(os.path.join(outdir, "N4_fig3_no_absorbing_floor.png"), dpi=150)
    plt.close(fig)
    print("図を生成: figures/N4_fig3_no_absorbing_floor.(svg|png)")
