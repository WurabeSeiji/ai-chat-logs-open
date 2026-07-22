#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
実験T2: レジスタ占有率による結合重みの自己変調（α型捕捉の検査）

実験Tは固定結合 g=1 で、占有率が次数重み平衡 r* = N/(9N-8) に鋭く捕捉される
ことを示した。本実験では結合重み自体を力学量にする:

    g(τ) = ( |Z_R| / RMS(|Z_e|) )^β

外部関数を持ち込まない純粋な比率変調（零次同次＝公理0.5適合、パラメータフリー）。
「レジスタの結合の強さは、レジスタの相対振幅そのものである」という
Σx² = R² の左辺＝右辺フィードバックの最小実装。

判定:
  - 暴走: r→1（セクタがレジスタへ全没入＝一本波への再吸収）か r→0 か
  - 捕捉: g(τ) が途中で停留するか。停留値が cos²(πm/n)（U^n=I 型の算術族）
    に乗るか
  - 統制: β=0（g=1 固定＝実験T、既知の収束 r*=N/(9N-8)）
  - 全過程で零閉鎖 ΣZ² + Z_R² = 0 の保存を検査

出力: register_modulation_result_v2/ に JSON と図。
"""

import json
import math
import os
from fractions import Fraction

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE_DIR, "register_modulation_result_v2")

GAMMA = math.tan(math.pi / 144.0)
N_LIST = [8, 12]
SEEDS = [0, 1]
BETAS = [1.0, 0.5]
STEPS = 30000


def build_edge_adjacency(n):
    ea, eb = np.triu_indices(n, k=1)
    m = len(ea)
    A = np.zeros((m, m))
    for i in range(m):
        share = (ea == ea[i]) | (ea == eb[i]) | (eb == ea[i]) | (eb == eb[i])
        A[i, share] = 1.0
    np.fill_diagonal(A, 0.0)
    return A, m


def step(W, A_edge, m, g):
    """拡張状態の凍結解除 Cayley 更新。レジスタ結合重み g は今ステップの値。"""
    theta = np.angle(W)
    K = np.zeros((m + 1, m + 1))
    K[:m, :m] = A_edge * np.sin(theta[None, :m] - theta[:m, None])
    b = g * np.sin(theta[m] - theta[:m])  # K[e,R] = g sin(θ_R - θ_e)
    K[:m, m] = b
    K[m, :m] = -b
    sig = np.linalg.eigvalsh(1j * K).max()
    if sig < 1e-15:
        return W
    Kn = K / sig
    I = np.eye(m + 1)
    return np.linalg.solve(I - GAMMA * Kn, W + GAMMA * (Kn @ W))


def make_initial(rng, m):
    Z = rng.normal(size=m) + 1j * rng.normal(size=m)
    c = complex(Z @ Z)
    W = np.concatenate([Z, [1j * np.sqrt(c)]])
    return W / np.linalg.norm(W)


def coupling_from_state(W, m, beta):
    """g = (|Z_R| / RMS|Z_e|)^β。比率のみ＝スケール無名。"""
    reg = abs(W[m])
    rms = math.sqrt(float(np.mean(np.abs(W[:m]) ** 2)))
    if rms < 1e-300:
        return 1e6
    return (reg / rms) ** beta


def stalls(x, min_len=3000, band_rel=0.005):
    """x(τ) の停留区間（相対幅 band_rel 内に min_len 以上）を列挙。"""
    out = []
    i = 0
    while i < len(x):
        j = i
        lo = hi = x[i]
        while j < len(x):
            lo = min(lo, x[j])
            hi = max(hi, x[j])
            mid = 0.5 * (lo + hi)
            if mid > 0 and (hi - lo) / mid > band_rel:
                break
            j += 1
        if j - i >= min_len:
            out.append((i, j, float(np.mean(x[i:j]))))
            i = j
        else:
            i += 1
    return out


def nearest_cos2(value, n_max=24):
    """value に最も近い cos²(πm/n) と誤差。"""
    best = (None, None, np.inf)
    for n in range(2, n_max + 1):
        for mm in range(1, n):
            c = math.cos(math.pi * mm / n) ** 2
            e = abs(value - c)
            if e < best[2]:
                best = (mm, n, e)
    return {"m": best[0], "n": best[1], "cos2": math.cos(math.pi * best[0] / best[1]) ** 2,
            "err": best[2]}


def run_one(n, seed, beta):
    A_edge, m = build_edge_adjacency(n)
    rng = np.random.default_rng(95260722 + 1000 * n + seed)
    W = make_initial(rng, m)
    wtw0 = complex(W @ W)
    r_hist = np.empty(STEPS + 1)
    g_hist = np.empty(STEPS + 1)
    max_dev = 0.0
    for t in range(STEPS + 1):
        r_hist[t] = abs(W[m]) ** 2 / float(np.real(np.conj(W) @ W))
        g = 1.0 if beta == 0.0 else coupling_from_state(W, m, beta)
        g_hist[t] = g
        if t % 500 == 0:
            max_dev = max(max_dev, abs(complex(W @ W) - wtw0))
        if t < STEPS:
            W = step(W, A_edge, m, g)
    late = slice(len(r_hist) // 2, None)
    st_r = stalls(r_hist)
    st_g = stalls(g_hist)
    entry = {
        "n": n, "m": m, "seed": seed, "beta": beta,
        "max_dev_closure": max_dev,
        "r_initial": float(r_hist[0]), "r_final": float(r_hist[-1]),
        "r_late_mean": float(np.mean(r_hist[late])),
        "r_late_std": float(np.std(r_hist[late])),
        "g_final": float(g_hist[-1]),
        "g_late_mean": float(np.mean(g_hist[late])),
        "g_late_std": float(np.std(g_hist[late])),
        "r_stalls": [(a, b, v, nearest_cos2(v)) for a, b, v in st_r],
        "g_stalls": [(a, b, v) for a, b, v in st_g],
    }
    return entry, r_hist, g_hist


def main():
    os.makedirs(RESULT_DIR, exist_ok=True)
    results = []
    curves = {}
    for n in N_LIST:
        for beta in BETAS + [0.0]:
            for seed in SEEDS if beta != 0.0 else [0]:
                entry, r, g = run_one(n, seed, beta)
                results.append(entry)
                curves[(n, beta, seed)] = (r, g)
                tag = "統制(g=1固定)" if beta == 0.0 else f"β={beta}"
                rs = entry["r_stalls"]
                stall_txt = ""
                if rs:
                    a, b, v, near = rs[-1]
                    stall_txt = (f" 最終停留 r={v:.5f}"
                                 f"（最近接 cos²({near['m']}π/{near['n']})="
                                 f"{near['cos2']:.5f}, 誤差 {near['err']:.1e}）")
                print(f"[{tag}] N={n} seed={seed}: "
                      f"r {entry['r_initial']:.4f}→{entry['r_final']:.5f}"
                      f" (後期 {entry['r_late_mean']:.5f}±{entry['r_late_std']:.5f})"
                      f" g終={entry['g_final']:.4f}"
                      f" 閉鎖偏差={entry['max_dev_closure']:.1e}{stall_txt}",
                      flush=True)

    with open(os.path.join(RESULT_DIR, "summary_v2.json"), "w") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False, default=float)

    tau = np.arange(STEPS + 1)
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 7.6))
    for col, n in enumerate(N_LIST):
        ax = axes[0][col]
        for beta in BETAS:
            for seed in SEEDS:
                r, g = curves[(n, beta, seed)]
                ax.plot(tau, r, alpha=0.7, lw=0.8, label=rf"$\beta={beta}$ s{seed}")
        r0, _ = curves[(n, 0.0, 0)]
        ax.plot(tau, r0, "k--", alpha=0.8, lw=1.0, label="control $g{=}1$")
        ax.set_title(rf"$N={n}$: register occupancy $r(\tau)$", fontsize=10)
        ax.set_xlabel(r"$\tau$")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
        ax = axes[1][col]
        for beta in BETAS:
            for seed in SEEDS:
                r, g = curves[(n, beta, seed)]
                ax.plot(tau, g, alpha=0.7, lw=0.8, label=rf"$\beta={beta}$ s{seed}")
        ax.set_title(rf"$N={n}$: self-modulated coupling $g(\tau)$", fontsize=10)
        ax.set_xlabel(r"$\tau$")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
    axes[0][0].set_ylabel(r"$r$")
    axes[1][0].set_ylabel(r"$g$")
    fig.suptitle("Experiment T2: register-modulated coupling "
                 r"$g=(|Z_R|/\mathrm{RMS}|Z_e|)^\beta$", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(RESULT_DIR, "register_modulation_v2.png"), dpi=160)
    plt.close(fig)
    print(f"出力: {RESULT_DIR}")


if __name__ == "__main__":
    main()
