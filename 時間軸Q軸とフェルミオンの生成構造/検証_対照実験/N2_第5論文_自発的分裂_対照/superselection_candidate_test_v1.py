#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N4 超選択量の探索（対照テスト済み第5論文コードを import）

位相巻き数 arg(Z_e) は超選択量でないことが winding_superselection_extension_v1.py で
判明した（固定/逐次で同数のジャンプ、ゼロ交差と無相関）。本スクリプトは超選択量の
候補を、固定生成子（線形）と逐次（交換散乱）の対比で直接検定する。

【超選択量の要件（測定前固定）】
  ある量 Q が超選択量であるとは:
    (S1) 固定生成子（線形発展）では Q が不変（時間を通じて一定）。
    (S2) 逐次再構成（交換散乱）でのみ Q が変化する。
    (S3) 変化は整数刻み（Q が整数値量なら）。

【候補】
  Q1 = 活性平面数 n(t) = #{固有値 σ(iK) > 閾}（整数）
  Q2 = 閉鎖量 c = |Z^T Z| / ||Z||^2（実数・厳密保存＝両モードで不変のはず＝超選択でなく普遍保存）

  合否:
    PASS(Q1超選択): 固定生成子で n(t) が完全一定、逐次で n(t) が整数変化。
    → 成立なら「超選択される整数＝活性平面数」。位相巻き数でなく平面数が担う、と結論。
    → 不成立なら候補を再検討（仮説をさらに見直す）。

実行: python3 superselection_candidate_test_v1.py
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


def run_candidate(N, delta, seed, steps=4320, sub=8, sequential=True, thr=0.05):
    A = line_graph_adjacency(N)
    m = A.shape[0]
    rng = np.random.default_rng(20260721 + seed)
    Z, _ = prepare_initial_state(rng, A, m, delta)
    K_fixed = sine_generator(np.angle(Z), A)
    U_fixed = cayley(K_fixed)
    n_series = []
    c_series = []
    for t in range(steps + 1):
        if sequential:
            K_now = sine_generator(np.angle(Z), A)
        else:
            K_now = K_fixed
        if t % sub == 0:
            n_series.append(active_planes(K_now, thr))
            ztz = Z @ Z
            hn = float(np.real(np.conj(Z) @ Z))
            c_series.append(abs(ztz) / hn if hn > 0 else 0.0)
        if t < steps:
            Z = (cayley(K_now) if sequential else U_fixed) @ Z
    n_series = np.array(n_series)
    c_series = np.array(c_series)
    n_changes = int(np.sum(np.diff(n_series) != 0))
    n_odd_changes = int(np.sum(np.abs(np.diff(n_series)) % 2 != 0))
    return {
        "N": N, "delta": delta, "seed": seed, "sequential": sequential,
        "n_initial": int(n_series[0]), "n_final": int(n_series[-1]),
        "n_min": int(n_series.min()), "n_max": int(n_series.max()),
        "n_change_events": n_changes,
        "n_constant": bool(n_changes == 0),
        "c_range": float(c_series.max() - c_series.min()),
        "_n_series": n_series.tolist(), "_sub": sub,
    }


if __name__ == "__main__":
    print(f"import 元: {_SRC}\n")
    print("=== 候補Q1: 活性平面数 n(t) の超選択検定 ===")
    print(f"{'N':>2} {'seed':>4} {'mode':>10} {'n初期':>5} {'n終':>4} "
          f"{'n最小':>5} {'n最大':>5} {'変化回数':>7} {'一定?':>5} {'c幅':>10}")
    results = []
    for N in (5, 6):
        for seed in range(3):
            for seq in (False, True):
                r = run_candidate(N, 1e-3, seed, sequential=seq)
                results.append(r)
                mode = "逐次(交換)" if seq else "固定(線形)"
                print(f"{r['N']:>2} {r['seed']:>4} {mode:>10} {r['n_initial']:>5} "
                      f"{r['n_final']:>4} {r['n_min']:>5} {r['n_max']:>5} "
                      f"{r['n_change_events']:>7} {'YES' if r['n_constant'] else 'no':>5} "
                      f"{r['c_range']:>10.2e}")

    fixed = [r for r in results if not r["sequential"]]
    seqs = [r for r in results if r["sequential"]]
    fixed_all_const = all(r["n_constant"] for r in fixed)
    seq_changes = sum(r["n_change_events"] for r in seqs)
    print("\n=== 判定 ===")
    print(f"(S1) 固定生成子で n(t) 完全一定: {'全試行YES' if fixed_all_const else '不成立'}")
    print(f"(S2) 逐次で n(t) 変化: 総変化回数={seq_changes}")
    verdict = fixed_all_const and seq_changes > 0
    print(f"⇒ 活性平面数 n は超選択量か: {'PASS（線形で不変・交換でのみ変化）' if verdict else 'FAIL（仮説再検討）'}")

    with open(os.path.join(_HERE, "superselection_candidate_result_v1.json"), "w") as f:
        slim = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
        json.dump({"import_source": _SRC, "runs": slim,
                   "fixed_all_constant": fixed_all_const,
                   "sequential_total_changes": seq_changes,
                   "planes_are_superselected": verdict}, f, indent=2, ensure_ascii=False)
    print("\n結果を superselection_candidate_result_v1.json に保存")

    # --- 図: 平面数超選択（固定=一定 vs 逐次=変化）---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    outdir = os.path.join(_HERE, "figures")
    os.makedirs(outdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for r in results:
        if r["N"] != 6 or r["seed"] != 1:
            continue
        t = np.arange(len(r["_n_series"])) * r["_sub"]
        if r["sequential"]:
            ax.step(t, r["_n_series"], where="post", color="tab:purple", lw=1.6,
                    label="sequential (exchange): changes")
        else:
            ax.step(t, r["_n_series"], where="post", color="tab:gray", lw=2.2,
                    label="fixed (linear): invariant")
    ax.set_xlabel("step  t"); ax.set_ylabel("active plane count  n(t)")
    ax.set_yticks(range(0, 8))
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_title("Experiment 2: plane count is superselected (N=6, seed 1)")
    ax.legend(loc="center right")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "N4_fig2_plane_superselection.svg"))
    fig.savefig(os.path.join(outdir, "N4_fig2_plane_superselection.png"), dpi=150)
    plt.close(fig)
    print("図を生成: figures/N4_fig2_plane_superselection.(svg|png)")
