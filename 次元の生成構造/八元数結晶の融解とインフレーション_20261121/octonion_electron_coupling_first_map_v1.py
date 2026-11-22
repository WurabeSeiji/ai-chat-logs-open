#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
結合係数を先に読んでから演算して描き戻す写像(候補 2)と、その比較。仕様は思考実験2 の §24。

  (A) 結合係数が全部 1 の写像(候補 1) a' = a·(a+b) で描き戻した状態 (§24.1)
  (B) 結合係数の式 √M·k = Σ_{m=1}^{M-1}(1-m/M)e^{imΔ} と弦の長さの逆数 (§24.2)
  (C) 候補 2: 先に k_ab = <a, a·b>/E_a を描き戻す前の状態から読み、a' = a·(k_aa+k_ab) (§24.4)
  (D) 候補 2 を繰り返したときの、ベース波の中心位相。初期の差 0°・180° とその近く (§24.6)

octonion_pulse_map_O2_to_O_v1.py と同じフォルダで実行する。格子点数は 1024 に固定。
状態は倍音 1〜6 を等しい振幅で積んだ C 型((1,e7) の面)と D 型(Q_k の面)。
"""
import sys, io, contextlib
import numpy as np
sys.path.insert(0, '.')
with contextlib.redirect_stdout(io.StringIO()):      # 読み込むときの出力は捨てる
    import octonion_electron_candidates_h6_v1 as h
m = h.m; mul = h.mul; psi = h.psi; NG = h.NG
Cst = h.Cst; Dst = h.Dst; fmt = h.fmt; spec = h.spec


def k(a, b):
    """結合係数。描き戻す前の状態だけから読む(自分の密度で相手の波形を平均した八元数)。"""
    return m.inner(a, mul(a, b)) / m.power(a)


def bc(o):
    return np.broadcast_to(o, (NG, 8))


def centroid(X):
    """円周 ψ の上の密度 |Q|² の重心[度]と山の位置[度]。"""
    d = (X ** 2).sum(1)
    pk = np.degrees(psi[d.argmax()])
    return np.degrees(np.angle((d * np.exp(1j * psi)).sum())), pk - 360 if pk > 180 else pk


def base_phase(X, s):
    """(1,e7) の面のベース波(巻数 s)の位相[度]。"""
    Z = np.fft.fft(X[:, 0] + 1j * X[:, 7]) / NG
    return np.degrees(np.angle(Z[s % NG]))


def step(S):
    """候補 2: Q_i' = Q_i·(Σ_j k_ij)。k_ij は描き戻す前の状態から先に読む。"""
    return [mul(a, bc(sum(k(a, b) for b in S))) for a in S]


PAIRS = (("C·C", lambda D: (Cst(0), Cst(D))), ("C·C逆巻", lambda D: (Cst(0), Cst(D, +1))),
         ("D1·D1", lambda D: (Dst(1, 0), Dst(1, D))), ("D1·D2", lambda D: (Dst(1, 0), Dst(2, D))))

print("===== (A) 結合係数が全部 1 の写像 a'=a·(a+b), b'=b·(a+b) で描き戻した状態 =====")
for nm, fn in PAIRS:
    print("==", nm)
    for Dg in (10, 30, 60, 90, 120, 180):
        a, b = fn(np.radians(Dg)); P = a + b
        ar, br, al = mul(a, P), mul(b, P), mul(P, a)
        ca, cb = centroid(ar), centroid(br)
        print(" Δ=%3d  a': 密度の重心 %+7.2f 山 %+6.1f パワー %6.3f | b': 重心 %+7.2f 山 %+6.1f | 重心の間隔 %6.2f"
              " | 左右の波形の差 %.3f 密度の差 %.1e | a' の巻数"
              % (Dg, ca[0], ca[1], m.power(ar), cb[0], cb[1], (cb[0] - ca[0]) % 360, np.abs(al - ar).max(),
                 np.abs((al ** 2).sum(1) - (ar ** 2).sum(1)).max()), spec(ar))

print("\n===== (B) √M·k = Σ_{m=1}^{M-1}(1-m/M)e^{imΔ} と 1/(2sin(Δ/2)) (式の計算) =====")
for D in (2, 5, 10, 30, 60, 90, 120, 150, 180):
    d = np.radians(D); row = []
    for M in (6, 24, 96, 384):
        mm = np.arange(1, M); s = ((1 - mm / M) * np.exp(1j * mm * d)).sum()
        row.append("M=%d: %7.3f∠%4.0f°" % (M, abs(s), np.degrees(np.angle(s)) % 360))
    print(" Δ=%3d° " % D, "  ".join(row), "  1/弦 %7.3f  90+Δ/2=%5.1f°" % (1 / (2 * np.sin(d / 2)), 90 + D / 2))
a = Cst(0)
print(" 八元数の計算との照合 (M=6): " + "  ".join(
    "Δ=%d° %.3f" % (D, np.linalg.norm(k(a, Cst(np.radians(D)))) * np.sqrt(6)) for D in (10, 30, 60, 90, 120, 180)))

print("\n===== (C) 候補 2: 先に結合係数を読み、a' = a·(k_aa+k_ab) (右)、(k_aa+k_ab)·a (左) =====")
print("  4 面 = (1,e7) Q1 Q2 Q3 の 大きさ∠位相")
for nm, fn in PAIRS:
    print("==", nm)
    for Dg in (10, 30, 60, 90, 120, 180):
        a, b = fn(np.radians(Dg)); kaa, kab = k(a, a), k(a, b); H = kaa + kab
        ar, al = mul(a, bc(H)), mul(bc(H), a)
        print(" Δ=%3d k_aa:" % Dg, fmt(kaa)); print("       k_ab:", fmt(kab))
        print("       右: パワー %.4f 巻数" % m.power(ar), spec(ar), "  <a,a'>:", fmt(m.inner(a, ar)))
        print("       左: パワー %.4f 巻数" % m.power(al), spec(al), "  <a,a'>:", fmt(m.inner(a, al)))

print("\n===== (D) 候補 2 を繰り返す。位置 = ベース波(m=1)の中心位相 =====")
for nm, sb in (("C 型と同じ巻きの C 型", -1), ("C 型と逆巻きの C 型", +1)):
    print("==", nm)
    for Dg in (0, 1, 5, 90, 175, 179, 180):
        S = [Cst(0), Cst(np.radians(Dg), sb)]
        print(" 初期の中心位相の差 %3d°" % Dg)
        for it in range(4):
            pa, pb = base_phase(S[0], -1), base_phase(S[1], sb)
            print("   %d 回: a のベース波の位相 %+8.3f | b のベース波の位相 %+8.3f | 差 %+8.3f | パワー %.4g %.4g"
                  % (it, pa, pb, (pb - pa + 180) % 360 - 180, m.power(S[0]), m.power(S[1])))
            S = step(S)
