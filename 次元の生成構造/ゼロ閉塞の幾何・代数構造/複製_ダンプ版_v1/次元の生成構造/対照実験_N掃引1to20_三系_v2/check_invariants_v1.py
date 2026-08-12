#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""主軸の向きの安定性と、保存する二次不変量の検査 v1（主張18 の導出プログラム）

----------------------------------------------------------------------
何を問うているか
----------------------------------------------------------------------
主張12 は上位3方向への集中を述べ、§0B は主軸ベクトルが tau について連続
しないことを述べた。ではこの向きは tau を長く掃引すれば安定するのか。

**重要なのは向きでも符号でもなく、内積の取れる保存量が存在するかどうか
である。** 本プログラムは両方を測る。

----------------------------------------------------------------------
検査A：向きは安定するか（18-a）
----------------------------------------------------------------------
個々の主軸は固有値が縮退すればいくらでも回るので、判定すべきは上位 k 本が
張る **部分空間** である。射影子の重なり

    ov(P1, P2) = tr(P1 P2) / k = || U1^T U2 ||_F^2 / k

で測る。1 なら同一部分空間、0 なら直交。**乱数の k 次元部分空間どうしの
期待値は k/(N-1)** であり、これが「向きを完全に失った」水準である。

後期を三つの窓に分け、ラグ別に重なりを測る。安定化しているなら後の窓ほど
重なりが高くなるはずである。

----------------------------------------------------------------------
検査B：内積の取れる保存量はあるか（18-b, 18-c, 18-d）
----------------------------------------------------------------------
関係ごとの複素振幅を
    x_e(tau) = sum over the trailing axes of C2[tau, e]
とする（距離の読出し lengths_from_C2 と同じ集約）。次を測る。

    sum_e |x_e|^2      エルミート内積（共役を取る）
    sum_e  x_e^2       双線形（共役を取らない。零閉鎖そのもの）。複素数
    tr(B), tr(B^2), tr(B^3)     グラム行列のスペクトル不変量

さらに **集約の順序** を変えた量

    sum_{e,k,j} C2[e,k,j]^2    成分ごとに二乗してから和（モデル内部の
                               closure 診断と一致する層）

も測り、どの層で保存が成り立つかを分ける。

**注意**：sum_e x_e^2 は複素数である。絶対値だけでなく実部・虚部を別々に
見ること。偏角が動けば保存ではない。

----------------------------------------------------------------------
使い方
----------------------------------------------------------------------
    python3 check_invariants_v1.py <stem> [--N 16] [--late 10000]
                                   [--windows 10000,20000,30000,40000]

例:
    python3 check_invariants_v1.py electron_T40000_d0.1_rep-dump40k16_N16 --N 16
    python3 check_invariants_v1.py electron_T40000_d0.1_rep-dump40k_N12  --N 12

出力: figures_tau/invariants_{stem}_v1.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(stem: str, side: str):
    X = np.load(HERE / f"dump_C2_{stem}_{side}_v1.npy", mmap_mode="r")
    taus = np.load(HERE / f"dump_meta_{stem}_{side}_v1.npz")["dump_taus"]
    return X, taus


def amplitudes(C: np.ndarray) -> np.ndarray:
    """関係ごとの複素振幅。倍音を先に和する（距離の読出しと同じ集約）。"""
    return C.reshape(C.shape[0], -1).sum(axis=1)


def gram(d: np.ndarray, N: int) -> np.ndarray:
    """長さ -> 二重中心化グラム行列。"""
    ia, ib = np.triu_indices(N, k=1)
    D2 = np.zeros((N, N))
    D2[ia, ib] = d ** 2
    D2[ib, ia] = d ** 2
    J = np.eye(N) - np.ones((N, N)) / N
    B = -0.5 * (J @ D2 @ J)
    return 0.5 * (B + B.T)


def top_basis(C: np.ndarray, N: int) -> np.ndarray:
    """固有値の大きい順に並べた主軸。"""
    B = gram(np.abs(amplitudes(C)), N)
    lam, U = np.linalg.eigh(B)
    return U[:, np.argsort(-lam)]


def overlap(U1: np.ndarray, U2: np.ndarray, k: int) -> float:
    """tr(P1 P2)/k。1 なら同一部分空間、0 なら直交。"""
    return float(np.sum((U1[:, :k].T @ U2[:, :k]) ** 2) / k)


def check_orientation(X, taus, N, windows, lags, ks):
    """検査A：向きは安定するか。"""
    out = {"random_baseline": {str(k): k / (N - 1) for k in ks}, "windows": []}
    for lo, hi in zip(windows[:-1], windows[1:]):
        idx = np.flatnonzero((taus >= lo) & (taus < hi))[::2]
        if len(idx) < 4:
            continue
        step = int(taus[idx[1]] - taus[idx[0]])
        Us = [top_basis(np.asarray(X[int(f)]), N) for f in idx]
        row = {"tau_lo": int(lo), "tau_hi": int(hi), "n": len(idx),
               "sample_step": step, "k": {}}
        for k in ks:
            row["k"][str(k)] = {
                str(l * step): float(np.mean(
                    [overlap(Us[i], Us[i + l], k) for i in range(len(Us) - l)]))
                for l in lags if l < len(Us)}
        out["windows"].append(row)
    return out


def check_invariants(X, taus, N, late, windows):
    """検査B：保存する二次不変量はどれか。"""
    idx = np.flatnonzero(taus >= late)[::4]
    bil, herm, trB, trB2, trB3, comp = [], [], [], [], [], []
    for f in idx:
        C = np.asarray(X[int(f)])
        x = amplitudes(C)
        bil.append((x ** 2).sum())
        herm.append((np.abs(x) ** 2).sum())
        lam = np.linalg.eigvalsh(gram(np.abs(x), N))
        trB.append(lam.sum()); trB2.append((lam ** 2).sum())
        trB3.append((lam ** 3).sum())
        comp.append((C ** 2).sum())

    def spread(v):
        v = np.asarray(v)
        m = np.abs(v).mean()
        return float((np.abs(v).max() - np.abs(v).min()) / m) if m else float("nan")

    bil = np.asarray(bil); comp = np.asarray(comp)
    res = {
        "late_from": int(late), "n": len(idx),
        "hermitian_sum_abs_x2": {"mean": float(np.mean(herm)), "spread": spread(herm)},
        "bilinear_sum_x2": {"real": float(bil.real.mean()),
                            "imag": float(bil.imag.mean()),
                            "abs": float(np.abs(bil).mean()),
                            "spread": spread(bil),
                            "angle_range": float(np.ptp(np.angle(bil)))},
        "trB": {"mean": float(np.mean(trB)), "spread": spread(trB)},
        "trB2": {"mean": float(np.mean(trB2)), "spread": spread(trB2)},
        "trB3": {"mean": float(np.mean(trB3)), "spread": spread(trB3)},
        "component_layer_sum_C2": {"abs": float(np.abs(comp).mean()),
                                   "spread": spread(comp)},
        "across_transition": [],
    }
    # 転移を貫通するか（区間ごとの実部・虚部）
    for lo, hi in zip(windows[:-1], windows[1:]):
        ii = np.flatnonzero((taus >= lo) & (taus < hi))
        if len(ii) == 0:
            continue
        ii = ii[:: max(1, len(ii) // 12)]
        z = np.array([(amplitudes(np.asarray(X[int(f)])) ** 2).sum() for f in ii])
        res["across_transition"].append(
            {"tau_lo": int(lo), "tau_hi": int(hi), "n": len(z),
             "real": float(z.real.mean()), "imag": float(z.imag.mean()),
             "spread": spread(z)})
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stem")
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--late", type=int, default=10000)
    ap.add_argument("--windows", default="")
    ap.add_argument("--outdir", default="figures_tau")
    ns = ap.parse_args()

    N = ns.N
    lags = (1, 2, 4, 8, 16, 32)
    ks = (1, 3)
    report = {"stem": ns.stem, "N": N}

    for side in ("m", "v"):
        try:
            X, taus = load(ns.stem, side)
        except FileNotFoundError:
            report[side] = {"error": "dump が無い"}
            continue
        wins = ([int(v) for v in ns.windows.split(",")] if ns.windows
                else [ns.late, ns.late + 10000, ns.late + 20000, int(taus[-1]) + 1])
        report[side] = {
            "orientation": check_orientation(X, taus, N, wins, lags, ks),
            "invariants": check_invariants(X, taus, N, ns.late, wins),
        }

    outdir = HERE / ns.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"invariants_{ns.stem}_v1.json"
    path.write_text(json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8")

    for side, lab in (("m", "物質"), ("v", "真空")):
        r = report.get(side, {})
        if "error" in r:
            print(f"=== {lab}: {r['error']}"); continue
        iv = r["invariants"]
        print(f"=== {lab}  N={N}  tau>={iv['late_from']}  {iv['n']} 点 ===")
        print(f"  Σ|x|²(エルミート) {iv['hermitian_sum_abs_x2']['mean']:.8f}"
              f"  変動 {iv['hermitian_sum_abs_x2']['spread']:.2e}")
        b = iv["bilinear_sum_x2"]
        print(f"  Σx² (双線形)  実 {b['real']:+.8e}  虚 {b['imag']:+.8e}"
              f"  変動 {b['spread']:.2e}  偏角幅 {b['angle_range']:.2e}")
        print(f"  tr(B) {iv['trB']['mean']:.8f} 変動 {iv['trB']['spread']:.2e}"
              f"   tr(B²) 変動 {iv['trB2']['spread']:.2e}"
              f"   tr(B³) 変動 {iv['trB3']['spread']:.2e}")
        print(f"  成分層 Σ_{{e,k,j}}C²  |·| {iv['component_layer_sum_C2']['abs']:.4e}"
              f"  変動 {iv['component_layer_sum_C2']['spread']:.2e}")
        for w in iv["across_transition"]:
            print(f"    τ {w['tau_lo']:6d}–{w['tau_hi']:6d}: "
                  f"実 {w['real']:+.8e}  虚 {w['imag']:+.8e}  変動 {w['spread']:.1e}")
        o = r["orientation"]
        print(f"  向き（乱数基準 k=3: {o['random_baseline']['3']:.4f}"
              f" / k=1: {o['random_baseline']['1']:.4f}）")
        for w in o["windows"]:
            for k in ("3", "1"):
                vals = w["k"][k]
                print(f"    τ {w['tau_lo']:6d}–{w['tau_hi']:6d} k={k}: "
                      + "  ".join(f"{lag}:{v:.3f}" for lag, v in vals.items()))
    print(f"-> {path}")


if __name__ == "__main__":
    main()
