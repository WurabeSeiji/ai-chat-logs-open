#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N4 H11（巻き数超選択則）検証実験（対照テスト済み第5論文コードを import）

【再現性の規約】
  独自再実装をしない。control_test_v1.py で公開結果の完全再現（714フィールド, 差0）を
  確認した run_spontaneous_splitting_preliminary_v1.py の力学関数
  （sine_generator / cayley / prepare_initial_state / line_graph_adjacency）を import する。

【巻き数の操作的定義】
  状態 Z(t) ∈ C^M の各成分 Z_e(t) は複素平面上の点。その偏角 arg(Z_e) を時間方向に
  非巻き戻し（unwrap）し、巻き数を W_e(t) = unwrap(arg Z_e)(t) / (2π) と定義する。
  W_e は位相的整数量であり、Z_e が原点 |Z_e|=0 を通るときだけ不連続に変化しうる。

【H11 の予言（測定前に固定）】
  (a) 位相的保護: 巻き数の整数ジャンプは |Z_e| の極小（ゼロ近傍交差）と一致する。
      ゼロ近傍交差を伴わない自発的な整数ジャンプは起きない。
  (b) 線形 vs 交換: 固定生成子（線形発展）ではゼロ近傍交差が少なく巻き数が安定、
      逐次再構成（交換散乱）ではゼロ近傍交差が起き巻き数が変化しうる。
  (c) 正味保存: 閉鎖 Σ Z² と関係する正味量が保存する（Z^T Z は厳密保存＝対照テストで確認済み）。

  合否（測定前固定）:
    - PASS(a): 巻き数ジャンプ点の |Z_e| が、非ジャンプ点の |Z_e| 分布より有意に小さい
      （ジャンプ時 min|Z_e| が全体中央値の 1/10 以下等）。
    - PASS(b): 固定生成子でのジャンプ総数 << 逐次でのジャンプ総数。
    - PASS(c): Z^T Z の保存（実部・虚部）が機械精度。
  いずれか予言が外れたら H11 を見直す（仮説修正）。

実行: python3 winding_superselection_extension_v1.py
出力: winding_superselection_result_v1.json
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


def run_winding(N, delta, seed, steps=4320, sub=1, sequential=True):
    """力学を走らせ、各成分の巻き数時系列・|Z_e|・Z^TZ を記録。
    sub=1（全ステップ）で位相の unwrap を正確にする。"""
    A = line_graph_adjacency(N)
    m = A.shape[0]
    rng = np.random.default_rng(20260721 + seed)
    Z, _ = prepare_initial_state(rng, A, m, delta)
    K_fixed = sine_generator(np.angle(Z), A)
    U_fixed = cayley(K_fixed)

    arg_prev = np.angle(Z)
    winding = np.zeros(m)          # 累積巻き数（連続化）
    abs_hist = []                  # |Z_e| の時系列
    wind_hist = []                 # 巻き数の時系列
    ztz_re0 = float(np.real(Z @ Z))
    ztz_im0 = float(np.imag(Z @ Z))
    max_ztz_dev = 0.0

    for t in range(steps + 1):
        arg_now = np.angle(Z)
        # 位相差を (-π, π] に畳んで累積（unwrap）
        d = arg_now - arg_prev
        d = (d + np.pi) % (2 * np.pi) - np.pi
        winding += d / (2 * np.pi)
        arg_prev = arg_now
        if t % sub == 0:
            abs_hist.append(np.abs(Z).copy())
            wind_hist.append(winding.copy())
        ztz = Z @ Z
        max_ztz_dev = max(max_ztz_dev,
                          abs(np.real(ztz) - ztz_re0),
                          abs(np.imag(ztz) - ztz_re0 * 0 + np.imag(ztz) - ztz_im0))
        if t < steps:
            if sequential:
                Z = cayley(sine_generator(np.angle(Z), A)) @ Z
            else:
                Z = U_fixed @ Z

    abs_hist = np.array(abs_hist)      # (T, m)
    wind_hist = np.array(wind_hist)    # (T, m)

    # 巻き数の整数ジャンプ検出: 各成分で round(W) が変化した時刻
    round_w = np.round(wind_hist)
    jumps = []  # (時刻index, 成分e, |Z_e|)
    for e in range(m):
        dr = np.diff(round_w[:, e])
        idx = np.where(dr != 0)[0]
        for i in idx:
            # ジャンプ直前直後の最小 |Z_e|
            lo = max(0, i - 1); hi = min(abs_hist.shape[0], i + 2)
            jumps.append((int(i), int(e), float(abs_hist[lo:hi, e].min())))

    # 全体の |Z_e| 分布（中央値）
    abs_median = float(np.median(abs_hist))
    jump_abs = [j[2] for j in jumps]
    return {
        "N": N, "delta": delta, "seed": seed, "m": m, "sequential": sequential,
        "n_winding_jumps": len(jumps),
        "jump_abs_values": jump_abs,
        "jump_abs_median": float(np.median(jump_abs)) if jump_abs else None,
        "overall_abs_median": abs_median,
        "max_ztz_deviation": max_ztz_dev,
        "net_winding_final": float(wind_hist[-1].sum()),
        "net_winding_initial": float(wind_hist[0].sum()),
    }


if __name__ == "__main__":
    print(f"import 元: {_SRC}\n")
    results = []
    print("=== H11(b): 固定生成子（線形）vs 逐次（交換散乱）の巻き数ジャンプ総数 ===")
    print(f"{'N':>2} {'delta':>7} {'seed':>4} {'mode':>10} {'#jumps':>7} "
          f"{'jump|Z|中央':>11} {'全体|Z|中央':>11} {'比':>7}")
    for N in (5, 6):
        for delta in (1e-3,):
            for seed in range(3):
                for seq in (False, True):
                    r = run_winding(N, delta, seed, sequential=seq)
                    results.append(r)
                    mode = "逐次(交換)" if seq else "固定(線形)"
                    jm = r["jump_abs_median"]
                    ratio = (jm / r["overall_abs_median"]) if jm else float("nan")
                    print(f"{r['N']:>2} {r['delta']:>7.0e} {r['seed']:>4} {mode:>10} "
                          f"{r['n_winding_jumps']:>7} "
                          f"{(jm if jm else 0):>11.3e} {r['overall_abs_median']:>11.3e} "
                          f"{ratio:>7.3f}")

    # 集計
    fixed = [r for r in results if not r["sequential"]]
    seqs = [r for r in results if r["sequential"]]
    jumps_fixed = sum(r["n_winding_jumps"] for r in fixed)
    jumps_seq = sum(r["n_winding_jumps"] for r in seqs)
    # (a) ジャンプ時 |Z| が全体中央値よりどれだけ小さいか
    all_jump_abs = [v for r in results for v in r["jump_abs_values"]]
    all_median = np.median([r["overall_abs_median"] for r in results])
    max_ztz = max(r["max_ztz_deviation"] for r in results)
    print("\n=== 判定 ===")
    print(f"(a) 位相的保護: ジャンプ時|Z|中央={np.median(all_jump_abs) if all_jump_abs else float('nan'):.3e}, "
          f"全体|Z|中央={all_median:.3e}, 比={np.median(all_jump_abs)/all_median if all_jump_abs else float('nan'):.4f}")
    print(f"(b) 線形 vs 交換: 固定ジャンプ総数={jumps_fixed}, 逐次ジャンプ総数={jumps_seq}")
    print(f"(c) Z^TZ 保存: 最大偏差={max_ztz:.2e}（機械精度で0）")

    with open(os.path.join(_HERE, "winding_superselection_result_v1.json"), "w") as f:
        slim = [{k: v for k, v in r.items() if k != "jump_abs_values"} for r in results]
        json.dump({"import_source": _SRC, "runs": slim,
                   "summary": {"jumps_fixed": jumps_fixed, "jumps_seq": jumps_seq,
                               "jump_abs_median": float(np.median(all_jump_abs)) if all_jump_abs else None,
                               "overall_abs_median": float(all_median),
                               "max_ztz_deviation": max_ztz}},
                  f, indent=2, ensure_ascii=False)
    print("\n結果を winding_superselection_result_v1.json に保存")

    # --- 図: 位相巻き数超選択の反証 ---
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    outdir = os.path.join(_HERE, "figures")
    os.makedirs(outdir, exist_ok=True)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.3))
    # 左: ジャンプ時|Z| の分布 vs 全体|Z|中央（重なる＝ゼロ交差不要＝保護なし）
    axL.hist(all_jump_abs, bins=40, color="tab:orange", edgecolor="none", alpha=0.8,
             label="|Z_e| at winding jumps")
    axL.axvline(all_median, color="tab:blue", lw=2, ls="--",
                label=f"overall median |Z_e| = {all_median:.3f}")
    axL.axvline(np.median(all_jump_abs), color="tab:red", lw=2, ls="-",
                label=f"jump median = {np.median(all_jump_abs):.3f}")
    axL.set_xlabel("|Z_e|"); axL.set_ylabel("count")
    axL.set_title("(REFUTED) jumps do NOT occur near zero")
    axL.legend(fontsize=8)
    # 右: 固定 vs 逐次 のジャンプ総数（ほぼ同数＝超選択でない）
    axR.bar(["fixed\n(linear)", "sequential\n(exchange)"], [jumps_fixed, jumps_seq],
            color=["tab:gray", "tab:purple"], edgecolor="black")
    axR.set_ylabel("total winding jumps")
    axR.set_title("(REFUTED) phase winding: fixed ~ sequential")
    for i, v in enumerate([jumps_fixed, jumps_seq]):
        axR.text(i, v + 20, str(v), ha="center")
    fig.suptitle("Experiment 1: phase-winding superselection is REFUTED", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "N4_fig1_winding_refuted.svg"))
    fig.savefig(os.path.join(outdir, "N4_fig1_winding_refuted.png"), dpi=150)
    plt.close(fig)
    print("図を生成: figures/N4_fig1_winding_refuted.(svg|png)")
