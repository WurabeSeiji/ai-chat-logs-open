#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""指定した τ の 3 次元楕円体図を描く v1

全 τ ダンプ `dump_C2_*.npy` から τ を指定して読み、M 本の関係を線分の長さと読み、
N 個の頂点の配置を復元して、頂点が乗る楕円体を描く。半軸は t, R, Q。

出力はファイル名に τ を含めるので上書きされない。
  figures_tau/fig_ellipsoid_tau{τ:05d}_{stem}_{side}_v1.png
  figures_tau/fig_ellipsoid_tau{τ:05d}_{stem}_{side}_v1.svg

使い方:
  python3 make_ellipsoid_figure_v1.py --tau 0 2200 2500 3999
  python3 make_ellipsoid_figure_v1.py --tau 1000 --side v
  python3 make_ellipsoid_figure_v1.py --tau 2500 --elev 25 --azim 40

----------------------------------------------------------------------
構成（本シリーズの定義に従う）
----------------------------------------------------------------------
1. 関係 e の複素振幅       z_e = Σ_{k,η} C2[e,k,η]
   長さ                    d_e = |z_e|
   辺の並びは正本と同じ    (i,j) = np.triu_indices(N, k=1)

2. 距離行列 D から二重中心化でグラム行列を作る（重心を原点に取る）
       B = −(1/2) J D∘D J,     J = I − (1/N) 11ᵀ
   B の固有値が正の方向が実の座標、負の方向は実配置では実現できない方向。
   負固有値の本数を「虚方向の本数」として図に併記する。

3. 上位 3 主軸へ射影して N 頂点の 3 次元座標 V を得る。
   慣性テンソル T = Vᵀ V の固有値の平方根が半軸で、大きい順に t, R, Q。
   楕円体は xᵀ T⁻¹ x = 3/N の等位面（N 頂点の広がりに一致する面）。

4. 描くもの: N 個の頂点、M 本の線分（関係）、楕円体、半軸 t, R, Q。
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

HERE = Path(__file__).resolve().parent
DEFAULT_STEM = "electron_d0.1_rep-dumpall_N12"


def lengths_from_C2(C2_tau: np.ndarray) -> np.ndarray:
    """(M, Nn, Nη) 複素 → (M,) 実。関係ごとの長さ |z_e|。"""
    z = C2_tau.reshape(C2_tau.shape[0], -1).sum(axis=1)
    return np.abs(z)


def configuration_from_lengths(d: np.ndarray, N: int):
    """M 本の長さ → N 頂点の座標（重心が原点）、固有値、虚方向の本数。"""
    ia, ib = np.triu_indices(N, k=1)
    assert len(ia) == len(d), f"辺数が合わない: {len(ia)} vs {len(d)}"
    D2 = np.zeros((N, N))
    D2[ia, ib] = d ** 2
    D2[ib, ia] = d ** 2
    J = np.eye(N) - np.ones((N, N)) / N
    B = -0.5 * J @ D2 @ J
    B = 0.5 * (B + B.T)
    lam, U = np.linalg.eigh(B)
    o = np.argsort(-lam)
    lam, U = lam[o], U[:, o]
    n_imag = int((lam < -1e-10 * max(1.0, abs(lam[0]))).sum())
    pos = np.maximum(lam, 0.0)
    V = U * np.sqrt(pos)                      # (N, N) 各列が主軸方向の座標
    return V, lam, n_imag, (ia, ib)


def draw(tau: int, side: str, stem: str, outdir: Path,
         elev: float, azim: float) -> None:
    c2_path = HERE / f"dump_C2_{stem}_{side}_v1.npy"
    if not c2_path.exists():
        raise FileNotFoundError(c2_path)
    X_all = np.load(c2_path, mmap_mode="r")
    T_total, M = X_all.shape[0], X_all.shape[1]
    if not (0 <= tau < T_total):
        raise IndexError(f"τ={tau} は範囲外（0..{T_total - 1}）")
    N = int(round((1 + np.sqrt(1 + 8 * M)) / 2))
    assert N * (N - 1) // 2 == M, f"M={M} が三角数でない"

    d = lengths_from_C2(np.asarray(X_all[tau]))
    Vfull, lam, n_imag, (ia, ib) = configuration_from_lengths(d, N)
    V = Vfull[:, :3]                                  # 上位3主軸へ射影

    T = V.T @ V
    ev, U3 = np.linalg.eigh(T)
    o = np.argsort(-ev)
    ev, U3 = ev[o], U3[:, o]
    semi = np.sqrt(np.maximum(ev, 0.0) * 3.0 / N)
    t_, R_, Q_ = semi

    # 楕円体メッシュ
    u = np.linspace(0, 2 * np.pi, 90)
    w = np.linspace(0, np.pi, 45)
    ex = np.outer(np.cos(u), np.sin(w)) * t_
    ey = np.outer(np.sin(u), np.sin(w)) * R_
    ez = np.outer(np.ones_like(u), np.cos(w)) * Q_
    P = np.stack([ex, ey, ez], axis=-1) @ U3.T

    fig = plt.figure(figsize=(8.4, 7.6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_wireframe(P[..., 0], P[..., 1], P[..., 2],
                      rstride=6, cstride=6, color="#9a9a9a",
                      linewidth=0.45, alpha=0.5)

    # M 本の関係（線分）。色は長さ。
    segs = [[V[i], V[j]] for i, j in zip(ia, ib)]
    lc = Line3DCollection(segs, cmap="viridis", linewidths=1.0, alpha=0.55)
    lc.set_array(d)
    ax.add_collection3d(lc)

    # N 個の頂点
    ax.scatter(V[:, 0], V[:, 1], V[:, 2], c="k", s=55,
               depthshade=False, zorder=20)

    # 半軸 t, R, Q
    for i, (lab, c) in enumerate(zip(["t", "R", "Q"],
                                     ["#c00000", "#0060c0", "#00a000"])):
        e = U3[:, i] * semi[i]
        ax.plot(*zip(-e, e), color=c, linewidth=3.0, zorder=30)
        ax.text(*(e * 1.18), lab, color=c, fontsize=19,
                fontweight="bold", ha="center")

    r = max(np.abs(V).max(), semi.max()) * 1.2
    ax.set_xlim(-r, r); ax.set_ylim(-r, r); ax.set_zlim(-r, r)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=elev, azim=azim)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.fill = False
        pane.pane.set_edgecolor("#e6e6e6")

    side_en = "matter" if side == "m" else "vacuum control"
    ax.set_title(
        f"tau = {tau}   ({side_en})\n"
        f"N = {N} vertices,  M = {M} relations\n"
        f"t = {t_:.6g}    R = {R_:.6g}    Q = {Q_:.6g}\n"
        f"imaginary directions = {n_imag} / {N - 1}",
        fontsize=12, pad=6)
    cb = fig.colorbar(lc, ax=ax, shrink=0.62, pad=0.02)
    cb.set_label("relation length  |z|", fontsize=10)

    outdir.mkdir(parents=True, exist_ok=True)
    base = f"fig_ellipsoid_tau{tau:05d}_{stem}_{side}_v1"
    for ext in ("png", "svg"):
        fig.savefig(outdir / f"{base}.{ext}", dpi=190, bbox_inches="tight")
    plt.close(fig)
    print(f"  tau={tau:>5}  N={N} M={M}  t={t_:.6g} R={R_:.6g} Q={Q_:.6g}"
          f"  imag={n_imag}  -> {base}.png / .svg")


def main() -> None:
    ap = argparse.ArgumentParser(description="指定 τ の 3 次元楕円体図")
    ap.add_argument("--tau", type=int, nargs="+", required=True,
                    help="描く τ（複数可）")
    ap.add_argument("--side", choices=["m", "v"], default="m",
                    help="m=物質側 / v=真空対照（既定 m）")
    ap.add_argument("--stem", default=DEFAULT_STEM,
                    help=f"ダンプのファイル語幹（既定 {DEFAULT_STEM}）")
    ap.add_argument("--outdir", default=str(HERE / "figures_tau"),
                    help="出力先（既定 ./figures_tau）")
    ap.add_argument("--elev", type=float, default=20.0)
    ap.add_argument("--azim", type=float, default=34.0)
    a = ap.parse_args()
    out = Path(a.outdir)
    for tau in a.tau:
        draw(tau, a.side, a.stem, out, a.elev, a.azim)


if __name__ == "__main__":
    main()
