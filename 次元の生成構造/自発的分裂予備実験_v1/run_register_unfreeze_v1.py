#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
実験T: レジスタの凍結解除（第6論文・開幕予備実験）

動機: Σx_n^2 = R^2 を「左辺＝右辺」の関係式として読めば、R^2 は左辺の
再配分に応答するフィードバック量である。従来実装は直交更新で
Z^T Z = R^2 を構成上の定数にし、このループを切断していた（レジスタの凍結）。

最小の修復: レジスタを明示的な成分 Z_R として状態に含め、全体の零閉鎖
    Σ_e Z_e^2 + Z_R^2 = 0
だけを保存する。生成子は位相差正弦結合のまま、レジスタを全関係波と結合する
仮想頂点として追加（K_{e,R} = g·sin(θ_R - θ_e)、零次同次＝公理0.5適合）。

測定:
  r(τ) = |Z_R|^2 / ‖W‖^2  …… レジスタ占有率（交換重みの読出し）
  - 統制 g=0（レジスタ切断）: r は厳密に凍結されるはず
  - g=1: r が動くか。正帰還で暴走するか、特定値に捕捉されるか、空転するか
  - 全過程で |W^T W|（零閉鎖）と ‖W‖² の保存を機械精度で検査
  - 後期窓の r 分布とプラトー（捕捉候補）の検出

出力: register_unfreeze_result_v1/ に JSON と図。
"""

import json
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE_DIR, "register_unfreeze_result_v1")

GAMMA = math.tan(math.pi / 144.0)
N_LIST = [8, 12]
SEEDS = [0, 1, 2, 3]
STEPS = 30000
G_COUPLE = 1.0


def build_adjacency_ext(n, g):
    """辺×辺の線グラフ隣接＋レジスタ仮想頂点（全辺と結合、重み g）。"""
    ea, eb = np.triu_indices(n, k=1)
    m = len(ea)
    A = np.zeros((m + 1, m + 1))
    for i in range(m):
        share = (ea == ea[i]) | (ea == eb[i]) | (eb == ea[i]) | (eb == eb[i])
        A[i, :m][share] = 1.0
    np.fill_diagonal(A, 0.0)
    A[:m, m] = g
    A[m, :m] = g
    return A, m


def step(W, A):
    """凍結解除 Cayley 更新（拡張状態、スペクトルノルム正規化）。"""
    theta = np.angle(W)
    K = A * np.sin(theta[None, :] - theta[:, None])
    sig = np.linalg.eigvalsh(1j * K).max()
    if sig < 1e-15:
        return W
    Kn = K / sig
    I = np.eye(len(W))
    return np.linalg.solve(I - GAMMA * Kn, W + GAMMA * (Kn @ W))


def make_initial(rng, m):
    """一般複素状態＋厳密零閉鎖レジスタ: Z_R = i·sqrt(Z^T Z)。"""
    Z = rng.normal(size=m) + 1j * rng.normal(size=m)
    c = complex(Z @ Z)
    zr = 1j * np.sqrt(c)
    W = np.concatenate([Z, [zr]])
    return W / np.linalg.norm(W)


def plateaus(r, min_len=2000, band=0.01):
    """r(τ) の捕捉候補（band 内に min_len 以上留まる区間）の値を列挙。"""
    out = []
    i = 0
    while i < len(r):
        j = i
        lo = hi = r[i]
        while j < len(r):
            lo = min(lo, r[j])
            hi = max(hi, r[j])
            if hi - lo > band:
                break
            j += 1
        if j - i >= min_len:
            out.append((i, j, float(np.mean(r[i:j]))))
            i = j
        else:
            i += 1
    return out


def run_one(n, seed, g):
    A, m = build_adjacency_ext(n, g)
    rng = np.random.default_rng(90260722 + 1000 * n + seed)
    W = make_initial(rng, m)
    wtw0 = complex(W @ W)
    r_hist = np.empty(STEPS + 1)
    max_dev_closure = 0.0
    for t in range(STEPS + 1):
        r_hist[t] = abs(W[m]) ** 2 / float(np.real(np.conj(W) @ W))
        if t % 500 == 0:
            max_dev_closure = max(max_dev_closure, abs(complex(W @ W) - wtw0))
        if t < STEPS:
            W = step(W, A)
    return {
        "n": n,
        "m": m,
        "seed": seed,
        "g": g,
        "abs_wtw_initial": abs(wtw0),
        "max_dev_closure": max_dev_closure,
        "r_initial": float(r_hist[0]),
        "r_final": float(r_hist[-1]),
        "r_min": float(np.min(r_hist)),
        "r_max": float(np.max(r_hist)),
        "r_late_mean": float(np.mean(r_hist[len(r_hist) // 2 :])),
        "r_late_std": float(np.std(r_hist[len(r_hist) // 2 :])),
        "plateaus": plateaus(r_hist),
    }, r_hist


def main():
    os.makedirs(RESULT_DIR, exist_ok=True)
    results = []
    curves = {}

    for n in N_LIST:
        for seed in SEEDS:
            entry, r = run_one(n, seed, G_COUPLE)
            results.append(entry)
            curves[(n, seed, G_COUPLE)] = r
            print(
                f"[g={G_COUPLE}] N={n} seed={seed}: r {entry['r_initial']:.4f}"
                f" → {entry['r_final']:.4f}"
                f" (min {entry['r_min']:.4f}, max {entry['r_max']:.4f},"
                f" 後期平均 {entry['r_late_mean']:.4f}±{entry['r_late_std']:.4f})"
                f" 捕捉候補={len(entry['plateaus'])}"
                f" 閉鎖偏差={entry['max_dev_closure']:.1e}",
                flush=True,
            )
        entry0, r0 = run_one(n, SEEDS[0], 0.0)
        results.append(entry0)
        curves[(n, SEEDS[0], 0.0)] = r0
        print(
            f"[統制 g=0] N={n}: r 変動幅 = {entry0['r_max'] - entry0['r_min']:.2e}"
            f"（凍結確認）閉鎖偏差={entry0['max_dev_closure']:.1e}",
            flush=True,
        )

    with open(os.path.join(RESULT_DIR, "summary_v1.json"), "w") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)

    tau = np.arange(STEPS + 1)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6), sharey=True)
    for ax, n in zip(axes, N_LIST):
        for seed in SEEDS:
            ax.plot(tau, curves[(n, seed, G_COUPLE)], alpha=0.7, lw=0.8,
                    label=f"seed {seed}")
        ax.plot(tau, curves[(n, SEEDS[0], 0.0)], "k--", alpha=0.8,
                label="control g=0 (frozen)")
        ax.set_xlabel(r"$\tau$")
        ax.set_title(rf"$N={n}$ ($M={n*(n-1)//2}$): register occupancy", fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
    axes[0].set_ylabel(r"$r(\tau)=|Z_R|^2/\|W\|^2$")
    fig.suptitle("Experiment T: unfreezing the register "
                 r"($\sum Z_e^2 + Z_R^2 = 0$ preserved exactly)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "register_occupancy_v1.png"), dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    allr = np.concatenate(
        [curves[(n, s, G_COUPLE)][STEPS // 2 :] for n in N_LIST for s in SEEDS]
    )
    ax.hist(allr, bins=120, density=True, alpha=0.8)
    ax.set_xlabel(r"$r$ (late half, all runs)")
    ax.set_ylabel("density")
    ax.set_title("Late-window distribution of register occupancy", fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "register_histogram_v1.png"), dpi=160)
    plt.close(fig)
    print(f"出力: {RESULT_DIR}")


if __name__ == "__main__":
    main()
