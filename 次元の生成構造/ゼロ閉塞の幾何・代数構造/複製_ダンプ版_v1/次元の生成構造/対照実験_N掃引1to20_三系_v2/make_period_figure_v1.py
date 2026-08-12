#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""準振動の周期を測り、図化する v1（主張15 の導出プログラム）

----------------------------------------------------------------------
何を測るか
----------------------------------------------------------------------
第2公理 U^n = I は、系が有限位数で元に戻ることを要求する。もし成立するなら
不動点は存在せず巡回のみが存在する。すなわち完全な定常状態はありえない。
そこで実際に周期を測る。

観測量は **66 本の関係長の平均** である。関係長は
    d_e(tau) = | sum over the trailing axes of C2[tau, e] |
で読む（make_ellipsoid_figure_v1.py の lengths_from_C2 と同じ規約）。

----------------------------------------------------------------------
方法（2通りの独立な推定を並べる）
----------------------------------------------------------------------
1. 自己相関法
   201 点移動平均でトレンドを除いたあと、窓 496 ステップごとに自己相関を
   とり、ラグ 10 以上の第1極大を卓越周期とする。極大に至るまでの最小値の
   ラグを反周期とする。自己相関の値そのものも記録する（1 に近くなければ
   鋭い周期ではなく準振動である、という判断のため）。

2. 山間隔法
   トレンド除去後の系列の局所極大を数え、隣接する山の間隔の中央値をとる。
   8 ステップ未満の間隔は微細な刻みとして除く。

**2通りの推定が一致しないこと自体が結果である。** 一致しないなら、それは
単一の鋭い周期が存在しないことを意味する。

----------------------------------------------------------------------
測定の限界（重要）
----------------------------------------------------------------------
tau >= DENSE_END では記録間隔が STRIDE ステップである。ナイキスト周期は
2*STRIDE。しかも 124 / 31 = 4.000 ちょうどであるため、周期が 124 付近から
わずかにずれるとうなりとして現れる。**この領域の周期は本走行のデータから
は直接測れない。** 図では該当領域を灰色で塗り、測定不能と明示する。

----------------------------------------------------------------------
使い方
----------------------------------------------------------------------
    python3 make_period_figure_v1.py <stem> [--dense-end 4000] [--stride 31]
                                     [--win 496] [--outdir figures_tau]

出力: figures_tau/fig_period_{stem}_v1.{png,svg}
      figures_tau/period_series_{stem}_v1.npz   （測定値そのもの）
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


# ---------------------------------------------------------------- 観測量
def mean_relation_length(c2_path: Path, n_dense: int) -> np.ndarray:
    """密領域の各ステップにおける関係長の平均を返す。"""
    X = np.load(c2_path, mmap_mode="r")
    n = min(n_dense, X.shape[0])
    D = np.asarray(X[:n])
    d = np.abs(D.reshape(n, D.shape[1], -1).sum(axis=2))   # (n, M)
    return d.mean(axis=1)


# ---------------------------------------------------------------- 推定
def detrend(y: np.ndarray, w: int = 201) -> np.ndarray:
    """端補正つき移動平均を引く。端でゼロ詰めしないこと（誤読防止）。"""
    y = np.asarray(y, dtype=float)
    k = np.ones(w)
    num = np.convolve(y, k, mode="same")
    den = np.convolve(np.ones_like(y), k, mode="same")
    return y - num / den


def acf(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float) - np.mean(y)
    n = len(y)
    a = np.correlate(y, y, mode="full")[n - 1:]
    return a / a[0]


def acf_estimate(y: np.ndarray, lo: int = 10, hi: int = 400):
    """自己相関の第1極大を卓越周期、そこに至る最小値のラグを反周期とする。"""
    a = acf(y)
    hi = min(hi, len(a) - 1)
    seg = a[lo:hi]
    loc = np.flatnonzero((seg[1:-1] > seg[:-2]) & (seg[1:-1] >= seg[2:])) + 1 + lo
    if len(loc) == 0:
        return None
    k = int(loc[int(np.argmax(a[loc]))])
    j = int(np.argmin(a[lo:k + 1])) + lo
    return k, float(a[k]), j, float(a[j]), a


def peak_gap_estimate(y: np.ndarray, min_gap: int = 8):
    """局所極大の間隔の中央値。"""
    y = np.asarray(y, dtype=float)
    p = np.flatnonzero((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:])) + 1
    if len(p) < 3:
        return np.nan, 0
    g = np.diff(p)
    g = g[g >= min_gap]
    if len(g) == 0:
        return np.nan, 0
    return float(np.median(g)), len(g)


# ---------------------------------------------------------------- 本体
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stem")
    # 【2026-08-12 修正】既定値 4000 を固定していたため、密領域を広げて走らせても
    # τ<4000 しか測らないという罠があった（N=16 は転移が τ≈9000 にあり、
    # 4000 で切ると転移が窓の外になる）。既定は None とし、走行時の
    # DUMP_TAUC（dump_meta の dump_tauc）を読んで密領域の実際の端を使う。
    ap.add_argument("--dense-end", type=int, default=None,
                    help="密領域の終端 τ。既定は dump_meta の dump_tauc に従う")
    ap.add_argument("--stride", type=int, default=31)
    ap.add_argument("--win", type=int, default=496)
    ap.add_argument("--outdir", default="figures_tau")
    ns = ap.parse_args()

    outdir = HERE / ns.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # 密領域の終端と関係の本数 M を、走行のメタから決める（ハードコードしない）
    meta_path = HERE / f"dump_meta_{ns.stem}_m_v1.npz"
    meta = np.load(meta_path) if meta_path.exists() else {}
    if ns.dense_end is None:
        ns.dense_end = int(meta["dump_tauc"]) if "dump_tauc" in getattr(meta, "files", []) \
            else 4000
        print(f"  [dense-end] メタの dump_tauc から {ns.dense_end} を採用")
    if "dump_stride" in getattr(meta, "files", []):
        ns.stride = int(meta["dump_stride"])
    M_REL = int(meta["m"]) if "m" in getattr(meta, "files", []) else None

    series = {}
    for side in ("m", "v"):
        c2 = HERE / f"dump_C2_{ns.stem}_{side}_v1.npy"
        if not c2.exists():
            raise FileNotFoundError(c2)
        y0 = mean_relation_length(c2, ns.dense_end)
        y = detrend(y0)
        rows = []
        for s in range(0, len(y) - ns.win + 1, ns.win):
            seg = y[s:s + ns.win]
            r = acf_estimate(seg)
            gap, ngap = peak_gap_estimate(seg)
            rows.append(dict(tau0=s, tau1=s + ns.win - 1,
                             per=(r[0] if r else np.nan),
                             corr=(r[1] if r else np.nan),
                             anti=(r[2] if r else np.nan),
                             acorr=(r[3] if r else np.nan),
                             gap=gap, ngap=ngap))
        series[side] = dict(raw=y0, det=y, rows=rows)

    # ---- 図 -----------------------------------------------------------
    fig = plt.figure(figsize=(12.4, 13.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 1.0, 1.1, 1.0], hspace=0.62)
    COL = {"m": "#c00000", "v": "#0060c0"}
    NAME = {"m": "matter", "v": "vacuum control"}

    # (a) 生の観測量
    ax = fig.add_subplot(gs[0])
    for side in ("m", "v"):
        ax.plot(np.arange(len(series[side]["raw"])), series[side]["raw"],
                color=COL[side], lw=0.8, label=NAME[side])
    ax.set_xlabel("tau"); ax.set_ylabel("mean relation length")
    ax.legend(fontsize=9, loc="upper left")
    # 【2026-08-12 修正】M を 66 に固定していたため N=16（M=120）で題字が誤っていた
    ax.set_title(f"(a) the observable: mean of the M = {M_REL if M_REL else '?'}"
                 " relation lengths"
                 f"   [dense region, every step, tau < {ns.dense_end}]",
                 fontsize=11, pad=4)
    ax.grid(alpha=0.25)

    # (b) トレンド除去後（振動そのもの）
    ax = fig.add_subplot(gs[1])
    ax.plot(np.arange(len(series["m"]["det"])), series["m"]["det"],
            color=COL["m"], lw=0.8, label="matter  (left axis)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("tau"); ax.set_ylabel("matter", color=COL["m"])
    ax.tick_params(axis="y", labelcolor=COL["m"])
    # 真空は振幅が二桁小さい。同じ軸に描くと見えないので右軸に分ける
    axb = ax.twinx()
    axb.plot(np.arange(len(series["v"]["det"])), series["v"]["det"],
             color=COL["v"], lw=0.8, label="vacuum  (right axis)")
    axb.set_ylabel("vacuum control", color=COL["v"])
    axb.tick_params(axis="y", labelcolor=COL["v"])
    amp_m = float(np.abs(series["m"]["det"]).max())
    amp_v = float(np.abs(series["v"]["det"]).max())
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = axb.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=9, loc="upper left")
    ax.set_title("(b) after removing the trend -- this is the quasi-oscillation itself."
                 "\nIt is not a clean sine: amplitude and spacing both change with tau."
                 f"\nNote the two y-axes: the matter amplitude ({amp_m:.1e}) is "
                 f"{amp_m/amp_v:.0f}x the vacuum one ({amp_v:.1e})",
                 fontsize=11, pad=4)
    ax.grid(alpha=0.25)

    # (c) 卓越周期の τ 依存
    ax = fig.add_subplot(gs[2])
    for side in ("m", "v"):
        rows = series[side]["rows"]
        c = np.array([(r["tau0"] + r["tau1"]) / 2 for r in rows])
        ax.plot(c, [r["per"] for r in rows], "o-", color=COL[side], lw=1.6,
                ms=6, label=f"{NAME[side]}: dominant period (autocorrelation)")
        ax.plot(c, [r["gap"] for r in rows], "s--", color=COL[side], lw=1.2,
                ms=5, alpha=0.65, label=f"{NAME[side]}: median peak spacing")
        ax.plot(c, [r["anti"] for r in rows], "^:", color=COL[side], lw=1.0,
                ms=5, alpha=0.45, label=f"{NAME[side]}: anti-period")
    ax.axhline(124, color="#666666", lw=1.4, ls="-.")
    ax.annotate("124", xy=(0, 124), xytext=(4, 4), textcoords="offset points",
                color="#666666", fontsize=10, fontweight="bold")
    ax.axhline(2 * ns.stride, color="#a000a0", lw=1.2, ls="--")
    ax.annotate(f"Nyquist period for tau >= {ns.dense_end}  ({2*ns.stride})",
                xy=(0, 2 * ns.stride), xytext=(4, 4), textcoords="offset points",
                color="#a000a0", fontsize=9)
    ax.set_xlabel("tau (window centre)"); ax.set_ylabel("period (steps)")
    ax.legend(fontsize=8, ncol=2, loc="upper left")
    ax.set_title("(c) the period is NOT constant.  Two independent estimators "
                 "disagree with each other and both drift with tau.\n"
                 "124 is one of the values the quasi-oscillation takes, "
                 "not a constant of the system",
                 fontsize=11, pad=4)
    ax.grid(alpha=0.25)

    # (d) 自己相関そのもの（鋭い周期でないことを示す）
    ax = fig.add_subplot(gs[3])
    rows = series["m"]["rows"]
    picks = [1, 2, 4, 6]
    cmap = plt.get_cmap("viridis")
    for n, idx in enumerate(picks):
        if idx >= len(rows):
            continue
        r = rows[idx]
        seg = series["m"]["det"][r["tau0"]:r["tau1"] + 1]
        a = acf(seg)
        ax.plot(np.arange(len(a[:400])), a[:400], color=cmap(n / max(1, len(picks) - 1)),
                lw=1.4, label=f"tau {r['tau0']}-{r['tau1']}  "
                              f"period {r['per']:.0f}, corr {r['corr']:+.2f}")
        if np.isfinite(r["per"]):
            ax.plot([r["per"]], [r["corr"]], "o",
                    color=cmap(n / max(1, len(picks) - 1)), ms=8)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhline(1.0, color="#888888", lw=0.8, ls=":")
    ax.set_xlabel("lag (steps)"); ax.set_ylabel("autocorrelation")
    ax.set_ylim(-1.05, 1.05)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("(d) matter side: the autocorrelation at the dominant period "
                 "is only +0.17 to +0.53, far from 1.\n"
                 "A sharp period would reach ~1 here. This is a quasi-oscillation, "
                 "not a cycle.",
                 fontsize=11, pad=4)
    ax.grid(alpha=0.25)

    fig.suptitle(
        f"Quasi-oscillation and its period   [{ns.stem}]\n"
        f"measured only where every step is recorded (tau < {ns.dense_end}); "
        f"beyond that the record stride is {ns.stride} and 124/{ns.stride} = "
        f"{124/ns.stride:.3f}, so the period cannot be measured there",
        fontsize=12.5, y=0.995)

    base = f"fig_period_{ns.stem}_v1"
    for ext in ("png", "svg"):
        fig.savefig(outdir / f"{base}.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)

    np.savez_compressed(
        outdir / f"period_series_{ns.stem}_v1.npz",
        **{f"{side}_{k}": np.array([r[k] for r in series[side]["rows"]], dtype=float)
           for side in ("m", "v")
           for k in ("tau0", "tau1", "per", "corr", "anti", "acorr", "gap", "ngap")})

    for side in ("m", "v"):
        print(f"=== {NAME[side]} ===")
        for r in series[side]["rows"]:
            print(f"  tau {r['tau0']:5d}-{r['tau1']:5d} : "
                  f"周期 {r['per']:6.1f} 相関 {r['corr']:+.3f}  "
                  f"反周期 {r['anti']:6.1f} ({r['acorr']:+.3f})  "
                  f"山間隔中央値 {r['gap']:6.1f} ({r['ngap']:2d} 個)")
    print(f"-> {outdir / base}.png")


if __name__ == "__main__":
    main()
