#!/usr/bin/env python3
"""make_parent（倍音対応）可視化 v1

(1) 生成した Z ∈ C^{M×H} の全複素成分を複素平面にプロット（段 n で色分け）
(2) 閉塞の検算 δ = Σ(a+ib)² を計算し、δ の実部・虚部を同じ図に表示
    - 左パネル: 複素平面散布図＋δ の位置（赤★、原点近傍）＋数値注記
    - 右パネル: 段ごとの δₙ = v⁽ⁿ⁾ᵀv⁽ⁿ⁾ の実部・虚部と総和 δ の棒グラフ

出力: plot_N5_H8.png, plot_N40_H4.png（本フォルダ内）
テストと同一 seed を使用（40260801 / 40260802）。生成器本体は変更しない。
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from make_parent_harmonic_v1 import make_parent_harmonic

HERE = Path(__file__).resolve().parent

plt.rcParams["font.family"] = ["Hiragino Sans", "Hiragino Kaku Gothic ProN",
                                 "Yu Gothic", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False


def plot_case(n, H, seed, fname):
    # テスト v2 と同一の確定パラメータ（iters=2000, restarts=10, tol=1e-12）
    Z, info = make_parent_harmonic(n, H, seed, iters=2000, restarts=10, tol=1e-12)
    m = Z.shape[0]

    # 検算 δ = Σ(a+ib)²（全成分）と段ごとの δₙ
    delta = complex(np.sum(Z * Z))
    delta_n = [complex(Z[:, k] @ Z[:, k]) for k in range(H)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))
    cmap = plt.get_cmap("viridis", H)

    # ---- 左: 複素平面 ----
    for k in range(H):
        ax1.scatter(Z[:, k].real, Z[:, k].imag, s=28, color=cmap(k),
                    label=f"段 n={k+1} (σ₁={info['levels'][k]['sigma1']:.4f})",
                    alpha=0.85, edgecolors="none")
    r_amp = float(np.mean(np.abs(Z)))
    th = np.linspace(0, 2 * np.pi, 256)
    ax1.plot(r_amp * np.cos(th), r_amp * np.sin(th), ls="--", lw=0.7,
             color="gray", label=f"平均振幅円 r={r_amp:.4f}")
    ax1.plot([delta.real], [delta.imag], marker="*", color="red", ms=16,
             ls="none", label="δ = ΣZ²（赤★）")
    ax1.axhline(0, color="black", lw=0.5)
    ax1.axvline(0, color="black", lw=0.5)
    ax1.set_aspect("equal")
    ax1.set_xlabel("実部 a")
    ax1.set_ylabel("虚部 b")
    ax1.set_title(f"複素平面プロット  N={n}, M={m}, H={H}, seed={seed}\n"
                  f"全成分 {m}×{H}={m*H} 点（段で色分け）")
    ax1.annotate(f"δ = Σ(a+ib)²\nRe δ = {delta.real:+.3e}\nIm δ = {delta.imag:+.3e}\n"
                 f"|δ| = {abs(delta):.3e}",
                 xy=(0.02, 0.02), xycoords="axes fraction",
                 fontsize=10, ha="left", va="bottom",
                 bbox=dict(boxstyle="round", fc="lightyellow", ec="red", alpha=0.9))
    ax1.legend(fontsize=7, loc="upper right")

    # ---- 右: 段ごとの δₙ 実部・虚部 ----
    xs = np.arange(1, H + 1)
    w = 0.38
    ax2.bar(xs - w / 2, [d.real for d in delta_n], w, label="Re δₙ", color="tab:blue")
    ax2.bar(xs + w / 2, [d.imag for d in delta_n], w, label="Im δₙ", color="tab:orange")
    ax2.axhline(0, color="black", lw=0.7)
    ax2.bar([H + 1 - w / 2], [delta.real], w, color="navy", label="Re δ（総和）")
    ax2.bar([H + 1 + w / 2], [delta.imag], w, color="darkorange", label="Im δ（総和）")
    ax2.set_xticks(list(xs) + [H + 1])
    ax2.set_xticklabels([f"n={k}" for k in xs] + ["総和"])
    ax2.set_ylabel("δ の値")
    ax2.set_title("閉塞検算  δₙ = v⁽ⁿ⁾ᵀv⁽ⁿ⁾（段別）と δ = ΣZ²（総和）\n"
                  f"最大 |δₙ| = {max(abs(d) for d in delta_n):.2e}   |δ| = {abs(delta):.2e}")
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax2.legend(fontsize=8)

    fig.tight_layout()
    out = HERE / fname
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"{fname}: δ = {delta.real:+.3e} {delta.imag:+.3e}i  |δ|={abs(delta):.3e}  保存済み")
    return delta


def main():
    plot_case(5, 8, 40260801, "plot_N5_H8.png")
    plot_case(40, 4, 40260802, "plot_N40_H4.png")
    plot_case(300, 4, 40260803, "plot_N300_H4.png")


if __name__ == "__main__":
    main()
