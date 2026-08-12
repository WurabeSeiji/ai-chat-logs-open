#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""指定した τ の 3 次元楕円体図を描く v1

軸の名前は固有値の大きい順に付けた **序数的な抽象名** である。
  A, B, C : 第1〜第3主軸（3次元射影に現れる＝読み出される3方向）
  D, E, F : 第4〜第6主軸（実だが読み出されない）
  h, i, j, k, l : 第7〜第11主軸（符号は τ により実／虚が入れ替わる）
これらが物理的時空の軸（時間・空間・R・Q など）に対応するかどうかは
**同定の根拠がなく、今後の課題である**。名前に物理的意味は無い。

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
   慣性テンソル T = Vᵀ V の固有値の平方根が半軸で、大きい順に A, B, C。
   楕円体は xᵀ T⁻¹ x = 3/N の等位面（N 頂点の広がりに一致する面）。

4. 描くもの: N 個の頂点、M 本の線分（関係）、楕円体、半軸 A, B, C。
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
    # 二重中心化により B·1 = 0 が厳密に成り立つ。この自明ゼロは主軸ではない
    # ので、実／虚の数え上げから除く。数値誤差で -1e-17 になったものを虚方向
    # として数えてしまうと、(d) の見た目と 1 本ずれる。
    ones = np.ones(N) / np.sqrt(N)
    j0 = int(np.argmax(np.abs(U.T @ ones)))
    keep = np.ones(N, dtype=bool)
    keep[j0] = False
    lam_nt = lam[keep]                        # 自明ゼロを除いた N-1 個
    scale = max(abs(lam_nt[0]), 1e-300)
    n_imag = int((lam_nt < -1e-9 * scale).sum())
    n_bnd = int((np.abs(lam_nt) <= 1e-3 * scale).sum())
    pos = np.maximum(lam, 0.0)
    V = U * np.sqrt(pos)                      # (N, N) 各列が主軸方向の座標
    return V, lam, (n_imag, n_bnd, lam_nt), (ia, ib)


def draw(tau: int, side: str, stem: str, outdir: Path,
         elev: float, azim: float, absmax: float | None = None,
         zoom: float = 2.0) -> None:
    c2_path = HERE / f"dump_C2_{stem}_{side}_v1.npy"
    if not c2_path.exists():
        raise FileNotFoundError(c2_path)
    X_all = np.load(c2_path, mmap_mode="r")
    n_frames, M = X_all.shape[0], X_all.shape[1]
    # 二段サンプリングのとき、フレーム番号 ≠ τ。メタの対応表で引く。
    meta_path = HERE / f"dump_meta_{stem}_{side}_v1.npz"
    if meta_path.exists() and "dump_taus" in np.load(meta_path).files:
        taus = np.load(meta_path)["dump_taus"]
        hit = np.flatnonzero(taus == tau)
        if len(hit) == 0:
            near = taus[np.argmin(np.abs(taus - tau))]
            raise IndexError(
                f"τ={tau} は保存されていない（最も近い保存 τ={near}）")
        frame = int(hit[0])
    else:
        if not (0 <= tau < n_frames):
            raise IndexError(f"τ={tau} は範囲外（0..{n_frames - 1}）")
        frame = tau
    N = int(round((1 + np.sqrt(1 + 8 * M)) / 2))
    assert N * (N - 1) // 2 == M, f"M={M} が三角数でない"

    d = lengths_from_C2(np.asarray(X_all[frame]))
    Vfull, lam, (n_imag, n_bnd, lam_nt), (ia, ib) = \
        configuration_from_lengths(d, N)
    V = Vfull[:, :3]                                  # 上位3主軸へ射影

    T = V.T @ V
    ev, U3 = np.linalg.eigh(T)
    o = np.argsort(-ev)
    ev, U3 = ev[o], U3[:, o]
    semi = np.sqrt(np.maximum(ev, 0.0) * 3.0 / N)
    t_, R_, Q_ = semi

    # ---- 楕円体からの外れ ---------------------------------------------------
    # 楕円体は xᵀT⁻¹x = c（c = 3/N）の等位面。頂点 v に対し
    #     s = sqrt( vᵀT⁻¹v / c )
    # は「同じ方向で楕円体面まで測ったときの倍率」。s=1 なら面上、
    # s>1 なら外、s<1 なら内。全頂点で s=1 なら厳密に楕円体に乗る。
    c_lvl = 3.0 / N
    Tinv = np.linalg.pinv(T)
    s = np.sqrt(np.array([v @ Tinv @ v for v in V]) / c_lvl)
    dev_ratio = float(s.max() / s.min())
    dev_cv = float(s.std() / s.mean())
    dev_med = float(np.median(s))
    foot = V / s[:, None]                       # 同方向で楕円体面上に落とした点

    # 正規化の基準は、閉鎖が固定する二乗平均半径（主張5）
    r_rms = float(np.sqrt(np.trace(T) / N))
    frac3 = float(np.sum(np.maximum(lam, 0)[:3]) / np.sum(np.maximum(lam, 0)))
    lmax = float(np.abs(np.log10(s)).max())

    def panel(ax, k, fixed_half=None, title=""):
        """k で割った座標で1枚描く。fixed_half を与えると軸範囲を固定。"""
        Vk, footk, semik = V / k, foot / k, semi / k
        u = np.linspace(0, 2 * np.pi, 90)
        w = np.linspace(0, np.pi, 45)
        P = np.stack([np.outer(np.cos(u), np.sin(w)) * semik[0],
                      np.outer(np.sin(u), np.sin(w)) * semik[1],
                      np.outer(np.ones_like(u), np.cos(w)) * semik[2]],
                     axis=-1) @ U3.T
        ax.plot_wireframe(P[..., 0], P[..., 1], P[..., 2], rstride=6, cstride=6,
                          color="#9a9a9a", linewidth=0.45, alpha=0.5)
        lcx = Line3DCollection([[Vk[i], Vk[j]] for i, j in zip(ia, ib)],
                               cmap="viridis", linewidths=1.0, alpha=0.40)
        lcx.set_array(d)
        ax.add_collection3d(lcx)
        ax.add_collection3d(Line3DCollection(
            [[footk[i], Vk[i]] for i in range(N)],
            colors="#d62728", linewidths=2.0, alpha=0.9))
        ax.scatter(footk[:, 0], footk[:, 1], footk[:, 2], c="#d62728", s=14,
                   marker="x", depthshade=False, zorder=19)
        scx = ax.scatter(Vk[:, 0], Vk[:, 1], Vk[:, 2], c=np.log10(s),
                         cmap="coolwarm", vmin=-lmax, vmax=lmax, s=70,
                         depthshade=False, zorder=20,
                         edgecolors="k", linewidths=0.5)
        for i, (lab, col) in enumerate(zip(["A", "B", "C"],
                                           ["#c00000", "#0060c0", "#00a000"])):
            e = U3[:, i] * semik[i]
            ax.plot(*zip(-e, e), color=col, linewidth=3.0, zorder=30)
            ax.text(*(e * 1.18), lab, color=col, fontsize=17,
                    fontweight="bold", ha="center")
        h = (fixed_half if fixed_half
             else max(np.abs(Vk).max(), semik.max()) * 1.2) / zoom
        n_clip = int((np.abs(Vk) > h).any(axis=1).sum())
        if n_clip:
            ax.text2D(0.01, 0.945, f"{n_clip} vertices outside the frame",
                      transform=ax.transAxes, fontsize=8, color="#b00000",
                      va="top")
        ax.set_xlim(-h, h); ax.set_ylim(-h, h); ax.set_zlim(-h, h)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=elev, azim=azim)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
            pane.pane.fill = False
            pane.pane.set_edgecolor("#e6e6e6")
        ax.set_title(title, fontsize=11, pad=0)
        return lcx, scx

    fig = plt.figure(figsize=(8.6, 20.4))
    gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 1.0, 0.60, 0.60],
                          hspace=0.30)
    ax1 = fig.add_subplot(gs[0], projection="3d")
    lc, sc2 = panel(
        ax1, 1.0, fixed_half=absmax,
        title=("(a) absolute scale"
               + (f"   [half-range {absmax:g}]" if absmax else "   [auto]")
               + f"\nA = {t_:.5g}   B = {R_:.5g}   C = {Q_:.5g}"
                 f"   r_rms = {r_rms:.5g}"))
    ax2 = fig.add_subplot(gs[1], projection="3d")
    panel(ax2, r_rms, fixed_half=2.2,
          title=("(b) normalised by r_rms   [half-range 2.2]"
                 f"\nA/r_rms = {t_/r_rms:.4f}   B/r_rms = {R_/r_rms:.4f}"
                 f"   C/r_rms = {Q_/r_rms:.4f}"))

    side_en = "matter" if side == "m" else "vacuum control"
    fig.suptitle(
        f"tau = {tau}   ({side_en})   N = {N} vertices,  M = {M} relations\n"
        f"off-ellipsoid:  max/min = {dev_ratio:.2f}   CV = {dev_cv:.3f}"
        f"   median s = {dev_med:.3f}   (1.00 / 0.000 / 1.000 = exactly on)\n"
        f"principal directions: real / imaginary = {N - 1 - n_imag} / {n_imag}"
        f"   (of {N - 1}"
        + (f"; {n_bnd} within 1e-3 of the boundary)" if n_bnd else ")")
        + f"   top-3 share = {frac3:.3f}",
        fontsize=12, y=0.985)
    # ---- (c) スケール履歴（対数）/ (d) 11方向の推移 --------------------------
    def smooth(y, w=41):
        """端補正つき移動平均。

        mode="same" の素の畳み込みは配列外をゼロとみなすため、先頭と末尾の
        w/2 点が 0 に向かって落ちる。これは物理ではなく図化上のエイリアスで、
        「最後に急減した」という誤読を招く。窓内に実際に入ったサンプル数で
        割ることで、端でも平均の意味を保たせる。"""
        y = np.asarray(y, dtype=float)
        if len(y) < w:
            return y
        k = np.ones(w)
        num = np.convolve(y, k, mode="same")
        den = np.convolve(np.ones_like(y), k, mode="same")
        return num / den

    ax3 = fig.add_subplot(gs[2])
    ax4 = fig.add_subplot(gs[3])
    AXNAMES = ["A", "B", "C", "D", "E", "F", "h", "i", "j", "k", "l"]
    ser_path = HERE / f"scale_series_{stem}_{side}_v1.npz"
    if ser_path.exists():
        S = np.load(ser_path)
        tt = S["tau"]
        for key, col, lw, ls, lab in [
                ("t", "#c00000", 1.6, "-", "A  (1st principal axis)"),
                ("R", "#0060c0", 1.6, "-", "B  (2nd)"),
                ("Q", "#00a000", 1.6, "-", "C  (3rd)"),
                ("r_rms_3d", "#7f4fbf", 2.0, "-", "r_rms of the 3D projection"),
                ("r_rms_full", "#7f7f7f", 1.8, ":",
                 "sqrt(sum of POSITIVE eigenvalues / N)  -- NOT conserved"),
                ("r_rms_sgn", "#000000", 2.4, "--",
                 "sqrt(SIGNED tr(T) / N)  -- exactly conserved")]:
            ax3.plot(tt, smooth(S[key]), color=col, lw=lw, ls=ls, label=lab)
        ax3.axvline(tau, color="#ff8800", lw=2.5, alpha=0.9)
        ax3.set_yscale("log"); ax3.set_xlim(0, tt.max())
        ax3.set_xlabel("tau", fontsize=10)
        ax3.set_ylabel("scale (log)", fontsize=10)
        ax3.grid(alpha=0.25, which="both")
        ax3.legend(fontsize=7.6, loc="lower right", framealpha=0.92, ncol=1)
        ax3.set_title("(c) scale history   —   the SIGNED trace (black dashed) is "
                      "exactly conserved,\nwhile the real content (grey dotted) and "
                      "the 3 observable directions both grow.\n"
                      "Real and imaginary content grow together and cancel "
                      "in the signed sum.",
                      fontsize=10.5, pad=4)
        ax3.text(tau, ax3.get_ylim()[1], f" tau = {tau}", color="#cc6600",
                 fontsize=9, va="top", ha="left")

        sp = S["spec"]
        n_dir = sp.shape[1]
        cmap = plt.get_cmap("turbo")
        for i in range(n_dir):
            nm = AXNAMES[i] if i < len(AXNAMES) else str(i + 1)
            ax4.plot(tt, smooth(sp[:, i]), color=cmap(i / max(1, n_dir - 1)),
                     lw=1.3, label=nm)
            ax4.annotate(nm, xy=(tt[-1], smooth(sp[:, i])[-1]),
                         xytext=(4, 0), textcoords="offset points",
                         fontsize=9, color=cmap(i / max(1, n_dir - 1)),
                         va="center", fontweight="bold")
        ax4.axhline(0.0, color="k", lw=1.2)
        ax4.axvline(tau, color="#ff8800", lw=2.5, alpha=0.9)
        ax4.set_xlim(0, tt.max())
        ax4.set_xlabel("tau", fontsize=10)
        ax4.set_ylabel("signed sqrt(lambda)", fontsize=10)
        ax4.grid(alpha=0.25)
        late = tt >= 10000
        pos_frac = (sp[late] > 0).mean(axis=0) if late.sum() else np.zeros(n_dir)
        names = [AXNAMES[i] if i < len(AXNAMES) else str(i + 1)
                 for i in range(n_dir)]
        always_r = " ".join(n for n, f in zip(names, pos_frac) if f >= 0.999)
        cross = " ".join(n for n, f in zip(names, pos_frac) if 0.001 < f < 0.999)
        always_i = " ".join(n for n, f in zip(names, pos_frac) if f <= 0.001)
        ax4.set_title(
            f"(d) all {n_dir} principal directions, named by rank only:  "
            + " ".join(names) + "\n"
            "above 0 = real, below 0 = imaginary.   "
            f"for tau >= 10000:  always real = {always_r or '-'};   "
            f"crossing = {cross or '-'};   always imaginary = {always_i or '-'}\n"
            "Whether these correspond to physical spacetime axes"
            " is an open question.",
            fontsize=10, pad=4)
        ax4.text(tau, ax4.get_ylim()[1], f" tau = {tau}", color="#cc6600",
                 fontsize=9, va="top", ha="left")
    else:
        for a in (ax3, ax4):
            a.text(0.5, 0.5, f"scale series not found: {ser_path.name}\n"
                   "run make_scale_series_v1.py first",
                   ha="center", va="center", transform=a.transAxes, fontsize=10)
            a.set_axis_off()

    cb = fig.colorbar(lc, ax=[ax1, ax2], shrink=0.40, pad=0.02)
    cb.set_label("relation length  |z|", fontsize=10)
    cb2 = fig.colorbar(sc2, ax=[ax1, ax2], shrink=0.40, pad=0.06)
    cb2.set_label("log10  s   (0 = on the ellipsoid)", fontsize=10)
    for _a in (ax1, ax2):                       # 3D 軸の余白を詰めて拡大表示
        _p = _a.get_position()
        _a.set_position([_p.x0 - 0.06 * _p.width, _p.y0 - 0.08 * _p.height,
                         _p.width * 1.20, _p.height * 1.20])

    outdir.mkdir(parents=True, exist_ok=True)
    base = f"fig_ellipsoid_tau{tau:05d}_{stem}_{side}_v1"
    for ext in ("png", "svg"):
        fig.savefig(outdir / f"{base}.{ext}", dpi=190, bbox_inches="tight")
    plt.close(fig)
    print(f"  tau={tau:>6}  t={t_:.5g} R={R_:.5g} Q={Q_:.5g}"
          f"  外れ 最大/最小={dev_ratio:8.2f} CV={dev_cv:.3f} 中央値s={dev_med:.3f}"
          f"  虚方向={n_imag}  上位3={frac3:.3f}  -> {base}")


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
    ap.add_argument("--absmax", type=float, default=None,
                    help="上段(絶対スケール)の軸半幅を固定する値。"
                         "τ 間で大きさを比べるときに指定する")
    ap.add_argument("--zoom", type=float, default=2.0,
                    help="軸半幅をこの値で割って拡大表示する（既定 2.0）")
    ap.add_argument("--elev", type=float, default=20.0)
    ap.add_argument("--azim", type=float, default=34.0)
    a = ap.parse_args()
    out = Path(a.outdir)
    for tau in a.tau:
        draw(tau, a.side, a.stem, out, a.elev, a.azim, a.absmax, a.zoom)


if __name__ == "__main__":
    main()
