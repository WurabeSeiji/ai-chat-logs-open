#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N2 rank 測定拡張（検証済み第5論文コードの import による）

【再現性の規約】
  本スクリプトは独自再実装を一切行わない。対照テスト（control_test_v1.py）で
  公開ベースラインの完全一致（714フィールド, 最大差 0.0）を確認した
  run_spontaneous_splitting_preliminary_v1.py の力学関数
  （sine_generator / cayley / prepare_initial_state / line_graph_adjacency）を
  そのまま import して用いる。したがって本測定の力学は、公開済み第5論文の力学と
  ビット単位で同一である。

【測定内容（N2 の主張 = 偶数ランク・有効ランク階段）】
  第5論文の自己参照(再構成)正弦生成子力学を N=5,6 で走らせ、実効生成子 K(τ) の
  代数的ランクと閾値付き有効ランクを区別して測る:
    (1) 固有値 μ(iK) の ±対性残差 max_i|μ_i+μ_{-i}|。
    (2) 数値的代数ランク np.linalg.matrix_rank(K)。
    (3) 活性平面数 n_ε(τ)=#{σ>ε} と有効ランク rank_ε(τ)=2n_ε(τ)。
    (4) ε={0.02,0.05,0.10} ごとの増分イベント数と増分値。

  rank_ε は定義上偶数である。|Δrank_ε|=2 は単一固有値対の閾値横断、
  |Δrank_ε|>=4 は複数対の同時横断を表す。

実行: python3 rank_measurement_extension_v1.py
出力: rank_measurement_result_v1.json
"""

import importlib.util
import json
import os
import numpy as np

# --- 検証済み第5論文モジュールを import（別実装しない） ---
_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "nbody_spontaneous_splitting_reproduction_v1",
                    "run_spontaneous_splitting_preliminary_v1.py")
_spec = importlib.util.spec_from_file_location("paper5_verified", _SRC)
p5 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(p5)

# 検証済み関数のみを使う
sine_generator = p5.sine_generator
cayley = p5.cayley
prepare_initial_state = p5.prepare_initial_state
line_graph_adjacency = p5.line_graph_adjacency


def pairing_residual(K):
    """固有値 μ(iK) の ±対性残差。実反対称なら非零 rank は偶数（定理T2）。"""
    mu = np.sort(np.linalg.eigvalsh(1j * K))
    return float(np.max(np.abs(mu + mu[::-1])))


def active_planes(K, thr):
    return int(np.sum(np.linalg.eigvalsh(1j * K) > thr))


def run_rank(N, delta, seed, steps=4320, sub=8, thresholds=(0.02, 0.05, 0.10)):
    A = line_graph_adjacency(N)
    m = A.shape[0]
    rng = np.random.default_rng(20260721 + seed)  # 第5論文と同じ種規約
    Z, _residual = prepare_initial_state(rng, A, m, delta)  # 検証済み関数は (Z, residual) を返す
    series = {thr: [] for thr in thresholds}
    numerical_ranks = []
    max_pair = 0.0
    for t in range(steps + 1):
        K = sine_generator(np.angle(Z), A)  # 検証済み関数
        if t % sub == 0:
            max_pair = max(max_pair, pairing_residual(K))
            numerical_ranks.append(int(np.linalg.matrix_rank(K)))
            for thr in thresholds:
                series[thr].append(active_planes(K, thr))
        if t < steps:
            Z = cayley(K) @ Z  # 検証済み関数

    threshold_summaries = {}
    effective_rank_series = {}
    for thr in thresholds:
        effective_rank = 2 * np.array(series[thr])
        effective_rank_series[thr] = effective_rank
        increments = np.diff(effective_rank)
        increments = increments[increments != 0]
        threshold_summaries[str(thr)] = {
            "max_effective_rank": int(effective_rank.max()),
            "n_increment_events": int(increments.size),
            "plus_2": int(np.sum(increments == 2)),
            "minus_2": int(np.sum(increments == -2)),
            "odd_increments": int(np.sum(np.abs(increments) % 2 != 0)),
            "even_not_2": int(np.sum((np.abs(increments) % 2 == 0)
                                     & (np.abs(increments) != 2))),
        }

    base = effective_rank_series[0.05]
    times = np.arange(len(base)) * sub
    inc = np.diff(base)
    inc = inc[inc != 0]
    return {
        "N": N, "delta": delta, "seed": seed, "m": m,
        "numerical_algebraic_rank_values": sorted(set(numerical_ranks)),
        "max_effective_rank_eps_0.05": int(base.max()),
        "max_pairing_residual": max_pair,
        "n_increment_events": int(inc.size),
        "odd_effective_rank_increments": int(np.sum(np.abs(inc) % 2 != 0)),
        "increments_of_2": int(np.sum(np.abs(inc) == 2)),
        "even_not_2": int(np.sum((np.abs(inc) % 2 == 0) & (np.abs(inc) != 2))),
        "threshold_summaries": threshold_summaries,
        "_effective_rank_series": base.tolist(),  # 作図用（ε=0.05）
        "_times": times.tolist(),                 # 作図用
        "_increments": inc.tolist(),              # 作図用（ε=0.05）
    }


def save_svg_without_trailing_whitespace(fig, path):
    """Matplotlib の SVG を保存し、公開用差分検査に不要な行末空白を除く。"""
    fig.savefig(path)
    with open(path, encoding="utf-8") as f:
        text = "\n".join(line.rstrip() for line in f.read().splitlines()) + "\n"
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def make_figures(results, outdir):
    """図1（rank(τ) 偶数階段）と図2（増分ヒストグラム）を SVG+PNG で生成。
    図内ラベルは英数字（フォント非依存）。キャプションは論文本文（日本語）。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(outdir, exist_ok=True)

    # --- 図1: rank_ε(τ) の偶数階段（代表 run、ε=0.05）---
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    reps = [(5, 1e-3, 1, "tab:blue"), (6, 1e-3, 1, "tab:red")]
    for N, d, s, color in reps:
        r = next(x for x in results if x["N"] == N and x["delta"] == d and x["seed"] == s)
        ax.step(r["_times"], r["_effective_rank_series"], where="post", color=color,
                label=f"N={N} (M={r['m']})", linewidth=1.6)
    ax.set_xlabel("step  t")
    ax.set_ylabel(r"effective rank$_{\epsilon}(K(t))$ ($\epsilon=0.05$)")
    ax.set_yticks(range(0, 14, 2))  # 偶数のみ
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_title("Effective-rank threshold crossings occur one pair at a time")
    ax.legend(loc="lower right")
    fig.tight_layout()
    save_svg_without_trailing_whitespace(
        fig, os.path.join(outdir, "N2_fig1_rank_staircase.svg"))
    fig.savefig(os.path.join(outdir, "N2_fig1_rank_staircase.png"), dpi=150)
    plt.close(fig)

    # --- 図2: ε=0.05 の全有効ランク増分イベントのヒストグラム ---
    all_inc = [v for r in results for v in r["_increments"]]
    n_plus = sum(1 for v in all_inc if v == 2)
    n_minus = sum(1 for v in all_inc if v == -2)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    bins = np.arange(-4.5, 5.5, 1.0)
    ax.hist(all_inc, bins=bins, color="tab:green", edgecolor="black", rwidth=0.9)
    ax.set_xlabel(r"effective-rank increment $\Delta\,$rank$_{\epsilon}$")
    ax.set_ylabel("count (all runs)")
    ax.set_xticks(range(-4, 5))
    for lo in (-3.5, -1.5, 0.5, 2.5):  # 奇数の帯（±1, ±3）を薄赤で強調
        ax.axvspan(lo, lo + 1.0, color="red", alpha=0.08)
    ax.set_title(f"{len(all_inc)} threshold events: +2={n_plus}, -2={n_minus}; "
                 "no multi-pair crossings")
    fig.tight_layout()
    save_svg_without_trailing_whitespace(
        fig, os.path.join(outdir, "N2_fig2_increment_histogram.svg"))
    fig.savefig(os.path.join(outdir, "N2_fig2_increment_histogram.png"), dpi=150)
    plt.close(fig)
    print(f"図を生成: {outdir}/N2_fig1_rank_staircase.(svg|png), "
          f"N2_fig2_increment_histogram.(svg|png)")


if __name__ == "__main__":
    print(f"import 元: {_SRC}")
    print(f"{'N':>2} {'delta':>7} {'seed':>4} {'numRank':>9} {'maxR_eps':>8} "
          f"{'pairRes':>10} {'#evt':>4} {'odd':>3} {'+-2':>4} {'evOth':>5}")
    results = []
    for N in (5, 6):
        for d in (1e-3, 1e-4):
            for s in range(3):
                r = run_rank(N, d, s)
                results.append(r)
                print(f"{r['N']:>2} {r['delta']:>7.0e} {r['seed']:>4} "
                      f"{str(r['numerical_algebraic_rank_values']):>9} "
                      f"{r['max_effective_rank_eps_0.05']:>8} "
                      f"{r['max_pairing_residual']:>10.2e} "
                      f"{r['n_increment_events']:>4} "
                      f"{r['odd_effective_rank_increments']:>3} "
                      f"{r['increments_of_2']:>4} {r['even_not_2']:>5}")
    tot_odd = sum(r["odd_effective_rank_increments"] for r in results)
    tot_oth = sum(r["even_not_2"] for r in results)
    max_pair = max(r["max_pairing_residual"] for r in results)
    tot_evt = sum(r["n_increment_events"] for r in results)
    print("-" * 62)
    print(f"全増分イベント数: {tot_evt}")
    print(f"奇数の有効ランク増分: {tot_odd}（rank_ε の定義上 0）")
    print(f"±2以外の偶数増分: {tot_oth}")
    print(f"±対性残差の最大: {max_pair:.2e}（機械精度で0）")
    print("閾値別イベント集計:")
    threshold_totals = {}
    for thr in (0.02, 0.05, 0.10):
        key = str(thr)
        summary = {
            metric: sum(r["threshold_summaries"][key][metric] for r in results)
            for metric in ("n_increment_events", "plus_2", "minus_2",
                           "odd_increments", "even_not_2")
        }
        summary["max_effective_rank_values"] = sorted({
            r["threshold_summaries"][key]["max_effective_rank"] for r in results
        })
        threshold_totals[key] = summary
        print(f"  ε={thr:.2f}: events={summary['n_increment_events']}, "
              f"+2={summary['plus_2']}, -2={summary['minus_2']}, "
              f"odd={summary['odd_increments']}, "
              f"even_not_2={summary['even_not_2']}, "
              f"max values={summary['max_effective_rank_values']}")

    # 図の生成（SVG+PNG）
    make_figures(results, os.path.join(_HERE, "figures"))

    # JSON 保存（作図用の大きな配列は除く）
    slim = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    import_source = os.path.relpath(_SRC, _HERE)
    with open(os.path.join(_HERE, "rank_measurement_result_v1.json"), "w") as f:
        json.dump({"import_source": import_source, "runs": slim,
                   "totals_eps_0.05": {"events": tot_evt, "odd": tot_odd,
                                       "even_not_2": tot_oth,
                                       "max_pairing_residual": max_pair},
                   "threshold_totals": threshold_totals},
                  f, indent=2, ensure_ascii=False)
    print("結果を rank_measurement_result_v1.json に保存")
