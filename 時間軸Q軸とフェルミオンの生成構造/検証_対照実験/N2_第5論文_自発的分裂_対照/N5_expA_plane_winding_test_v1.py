#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N5 予備実験A: 個別平面の巻き数は超選択量か（対照テスト済み第5論文コードを import）

N4 で成分位相巻き数 arg(Z_e) は超選択でないと分かった。粒子マッピングが量子数として
使うのは「個別平面の巻き数」（xy巻=スピン, zR巻=質量, tQ巻=電荷 の想定）なので、
これを N4 と同型に直接検定する。

【平面巻き数の操作的定義】
  初期生成子 K0 の固有平面 (p_j, q_j) を固定参照枠とする。状態 Z を射影し
  ζ_j = (p_j·Z) + i(q_j·Z) とおくと、平面が角 φ 回転すると ζ_j → e^{iφ} ζ_j。
  ゆえに W_j(t) = unwrap(arg ζ_j)(t) / 2π が平面 j の巻き数（回転数）。

【予言（N4 の当初 H11 と同型・測定前固定）】
  (a) 位相的保護: 平面巻き数の整数ジャンプは |ζ_j| のゼロ近傍でのみ。
  (b) 線形 vs 交換: 固定生成子（線形）でジャンプ総数 << 逐次（交換）。
  → 成立なら平面巻き数は超選択量（＝量子数の候補として生存）。
  → 反証なら平面巻き数も超選択でなく、「量子数＝巻き数」を捨て、根ラベルへ移る。

実行: python3 N5_expA_plane_winding_test_v1.py
出力: N5_expA_plane_winding_result_v1.json
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


def reference_frames(K0, tol=1e-9):
    """初期生成子の固有平面 (p,q) を固定参照枠として返す（活性のみ）。"""
    mu, V = np.linalg.eigh(1j * K0)
    order = np.argsort(-mu)
    frames = []
    for idx in order:
        if mu[idx] <= tol:
            continue
        v = V[:, idx]
        p = np.sqrt(2.0) * v.real; p /= np.linalg.norm(p)
        q = np.sqrt(2.0) * v.imag; q = q - (q @ p) * p; q /= np.linalg.norm(q)
        frames.append((p, q))
    return frames


def run_plane_winding(N, delta, seed, steps=4320, sequential=True):
    A = line_graph_adjacency(N)
    m = A.shape[0]
    rng = np.random.default_rng(20260721 + seed)
    Z, _ = prepare_initial_state(rng, A, m, delta)
    K0 = sine_generator(np.angle(Z), A)
    U_fixed = cayley(K0)
    frames = reference_frames(K0)
    nfr = len(frames)

    def zetas(Zc):
        out = np.empty(nfr, dtype=complex)
        for j, (p, q) in enumerate(frames):
            out[j] = (p @ Zc) + 1j * (q @ Zc)
        return out

    z_prev = zetas(Z)
    arg_prev = np.angle(z_prev)
    winding = np.zeros(nfr)
    absz_hist = []
    wind_hist = []
    for t in range(steps + 1):
        zc = zetas(Z)
        arg_now = np.angle(zc)
        d = (arg_now - arg_prev + np.pi) % (2 * np.pi) - np.pi
        winding += d / (2 * np.pi)
        arg_prev = arg_now
        absz_hist.append(np.abs(zc))
        wind_hist.append(winding.copy())
        if t < steps:
            Z = (cayley(sine_generator(np.angle(Z), A)) if sequential else U_fixed) @ Z
    absz_hist = np.array(absz_hist)
    wind_hist = np.array(wind_hist)
    round_w = np.round(wind_hist)
    jumps_abs = []
    for j in range(nfr):
        idx = np.where(np.diff(round_w[:, j]) != 0)[0]
        for i in idx:
            lo = max(0, i - 1); hi = min(absz_hist.shape[0], i + 2)
            jumps_abs.append(float(absz_hist[lo:hi, j].min()))
    return {
        "N": N, "delta": delta, "seed": seed, "sequential": sequential, "nfr": nfr,
        "n_jumps": len(jumps_abs),
        "jump_absz_median": float(np.median(jumps_abs)) if jumps_abs else None,
        "overall_absz_median": float(np.median(absz_hist)),
        "_jump_absz": jumps_abs,
    }


if __name__ == "__main__":
    print(f"import 元: {_SRC}\n")
    print("=== 実験A: 個別平面の巻き数の超選択検定 ===")
    print(f"{'N':>2} {'seed':>4} {'mode':>10} {'#jumps':>7} {'jump|ζ|中央':>11} "
          f"{'全体|ζ|中央':>11} {'比':>7}")
    results = []
    for N in (5, 6):
        for seed in range(3):
            for seq in (False, True):
                r = run_plane_winding(N, 1e-3, seed, sequential=seq)
                results.append(r)
                mode = "逐次(交換)" if seq else "固定(線形)"
                jm = r["jump_absz_median"]
                ratio = (jm / r["overall_absz_median"]) if jm else float("nan")
                print(f"{r['N']:>2} {r['seed']:>4} {mode:>10} {r['n_jumps']:>7} "
                      f"{(jm if jm else 0):>11.3e} {r['overall_absz_median']:>11.3e} "
                      f"{ratio:>7.3f}")

    fixed = [r for r in results if not r["sequential"]]
    seqs = [r for r in results if r["sequential"]]
    jf = sum(r["n_jumps"] for r in fixed)
    js = sum(r["n_jumps"] for r in seqs)
    all_jump = [v for r in results for v in r["_jump_absz"]]
    all_med = np.median([r["overall_absz_median"] for r in results])
    print("\n=== 判定 ===")
    print(f"(a) 位相的保護: ジャンプ時|ζ|中央={np.median(all_jump) if all_jump else float('nan'):.3e}, "
          f"全体中央={all_med:.3e}, 比={np.median(all_jump)/all_med if all_jump else float('nan'):.4f}")
    print(f"(b) 線形 vs 交換: 固定ジャンプ総数={jf}, 逐次ジャンプ総数={js}")
    protected = (all_jump and np.median(all_jump) / all_med < 0.2)
    superselected = (jf < js * 0.2)
    print(f"⇒ 平面巻き数は超選択量か: "
          f"{'PASS（保護あり・線形不変）' if (protected and superselected) else 'FAIL（超選択でない＝量子数=巻き数を捨てる）'}")

    with open(os.path.join(_HERE, "N5_expA_plane_winding_result_v1.json"), "w") as f:
        slim = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
        json.dump({"import_source": _SRC, "runs": slim,
                   "jumps_fixed": jf, "jumps_seq": js,
                   "jump_absz_median": float(np.median(all_jump)) if all_jump else None,
                   "overall_absz_median": float(all_med),
                   "plane_winding_superselected": bool(protected and superselected)},
                  f, indent=2, ensure_ascii=False)
    print("\n結果を N5_expA_plane_winding_result_v1.json に保存")
