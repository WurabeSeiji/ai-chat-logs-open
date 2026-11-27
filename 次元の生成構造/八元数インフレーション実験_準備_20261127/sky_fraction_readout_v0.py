#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
結合係数 k (全天空に占める相手の面積の割合) の読みだけ。動力学なし。v0 (粗い最初の一本)。

仕様は同じフォルダの「思考実験.md」§11。思考実験 2 (2026-11-22) §26 の定義を使う。
  距離        : 中心位相の差 Δ (0〜180°)
  相手の大きさ : δ = 360°/M (M = 分解能 = 倍音の上限)
  自分自身     : k_aa = 0
  天空 S^n     : 見かけの半角 α の円が覆う割合 k_n(α) = ∫_0^α sin^(n-1) t dt / ∫_0^π sin^(n-1) t dt
                 n = 1, 2, 4, 6 の四通り

作業者 (Claude) が粗く置いたもの。著者の式ではない。
  [置き 1] 見かけの半角 α の作り方を二通り:
     A (平ら)  : α = (δ/2)/Δ          (角度の比をそのまま半角[rad]とする。π で頭打ち)
     B (球面)  : sin α = sin(δ/2)/sin Δ (sin Δ <= sin(δ/2) のときは α = π。対蹠点で再び全天空を覆う)
  [置き 2] 配置は、円周の上に等間隔に並んだ N 個、全員が同じ M (高対称な初期配置の最も単純な形)。
  [置き 3] 行の和 Σ_b k_ab = 1 になる M を、実数として二分法で求める (M は本来は整数)。

割り算と正規化: k を Σk で割って 1 に合わせる操作はしていない。α と k の定義の中の比だけである。

使い方: python3 sky_fraction_readout_v0.py
"""
import numpy as np

DIMS = (1, 2, 4, 6)


def cap_fraction(alpha, n):
    """S^n の上で、半角 alpha の円が覆う割合。"""
    a = np.clip(alpha, 0.0, np.pi)
    c = np.cos(a)
    if n == 1:
        return a / np.pi
    if n == 2:
        return (1 - c) / 2
    if n == 4:
        return (2 - 3 * c + c ** 3) / 4
    if n == 6:
        return (8 - 15 * c + 10 * c ** 3 - 3 * c ** 5) / 16
    raise ValueError(n)


def half_angle(delta, Delta, kind):
    """見かけの半角 α [rad]。delta, Delta は rad。"""
    if kind == "A":
        return np.minimum((delta / 2) / Delta, np.pi)
    s = np.sin(delta / 2) / np.maximum(np.sin(Delta), 1e-300)
    return np.where(s >= 1.0, np.pi, np.arcsin(np.minimum(s, 1.0)))


def distances(N):
    """等間隔 N 個のうちの一つから見た、他の N-1 個までの距離 [rad]。"""
    j = np.arange(1, N)
    d = 2 * np.pi * j / N
    return np.minimum(d, 2 * np.pi - d)


def row_sum(N, M, n, kind):
    delta = 2 * np.pi / M
    return float(cap_fraction(half_angle(delta, distances(N), kind), n).sum())


def solve_M(N, n, kind, lo=1.0, hi=1e9):
    """Σk = 1 になる M。Σk は M について単調減少。"""
    if row_sum(N, lo, n, kind) < 1.0:
        return None
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        if row_sum(N, mid, n, kind) > 1.0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


def main():
    print("=== (1) 二体。相手の大きさ δ = 360°/M、距離 Δ のときの k ===")
    for M in (6, 24):
        delta = 2 * np.pi / M
        print("  M = %d (δ = %.1f°)" % (M, 360 / M))
        for kind in ("A", "B"):
            print("    半角の作り方 %s" % kind)
            print("      Δ[度]   α[度]   " + "   ".join("S^%d      " % n for n in DIMS))
            for D in (1, 5, 10, 15, 30, 60, 90, 120, 150, 170, 179, 180):
                a = float(half_angle(delta, np.radians(D), kind))
                print("      %5d  %7.2f   " % (D, np.degrees(a))
                      + "   ".join("%.3e" % cap_fraction(a, n) for n in DIMS))

    print("\n=== (2) 等間隔 N 個、全員が同じ M。一つから見た行の和 Σk ===")
    for kind in ("A", "B"):
        print("  半角の作り方 %s" % kind)
        print("      N     M    " + "   ".join("S^%d      " % n for n in DIMS))
        for N in (3, 4, 6, 8, 12, 24, 48, 96):
            for M in (6, 24, 96):
                print("    %4d  %4d   " % (N, M)
                      + "   ".join("%.3e" % row_sum(N, M, n, kind) for n in DIMS))

    print("\n=== (3) 等間隔 N 個で Σk = 1 になる M (実数) と、そのときの δ/間隔 = N/M ===")
    for kind in ("A", "B"):
        print("  半角の作り方 %s" % kind)
        print("      N   " + "   ".join("S^%d: M*      N/M*  " % n for n in DIMS))
        for N in (2, 3, 4, 6, 8, 12, 24, 48, 96, 192, 384):
            cells = []
            for n in DIMS:
                Ms = solve_M(N, n, kind)
                cells.append("      なし          " if Ms is None else "%10.4f  %7.4f" % (Ms, N / Ms))
            print("    %4d   " % N + "   ".join(cells))


if __name__ == "__main__":
    main()
