#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3 個の八元数だけで、一回の演算で 3^3 の配列が変わらない状態 (初期値の単位) を探す。v0。
仕様は「思考実験.md」§16。

状態     : 八元数 3 個 O_1, O_2, O_3 (倍音なし。実数 8 個ずつ)
配列     : X_ijk = (O_i O_j) O_k  (3^3 = 27 個。値は八元数)
演算     : P = O_1 + O_2 + O_3 として、三つの積を面積 conj(P)P = |P|^2 で一度だけ割る。三つの形:
             右    O_i' = O_i (P P) / |P|^2
             真ん中 O_i' = (P O_i) P / |P|^2
             左    O_i' = (P P) O_i / |P|^2
           相手の重みは全部同じ (結合係数は入れていない)。作業者 (Claude) が置いた最も単純な形である。
割り算   : |P|^2 で一度だけ。P = 0 (線形和が閉じる) では定義できない。そのまま報告する。
使い方   : python3 three_octonion_fixed_point_v0.py
"""
import numpy as np

rng = np.random.default_rng(20261127)
TRI = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
T = np.zeros((8, 8, 8))
T[0, 0, 0] = 1
for i in range(1, 8):
    T[0, i, i] = 1; T[i, 0, i] = 1; T[i, i, 0] = -1
for a, b, c in TRI:
    for (i, j, k), s in [((a, b, c), 1), ((b, c, a), 1), ((c, a, b), 1), ((b, a, c), -1), ((a, c, b), -1), ((c, b, a), -1)]:
        T[i, j, k] = s
mul = lambda x, y: np.einsum('ijk,...i,...j->...k', T, x, y)
nrm = np.linalg.norm


def step(S, form):
    P = S.sum(axis=0); A = P @ P
    if form == "右":
        return np.array([mul(o, mul(P, P)) for o in S]) / A
    if form == "真ん中":
        return np.array([mul(mul(P, o), P) for o in S]) / A
    return np.array([mul(mul(P, P), o) for o in S]) / A


def array3(S):
    return np.array([[[mul(mul(a, b), c) for c in S] for b in S] for a in S])


def gram(S):
    return S @ S.T


def angle_P(S):
    P = S.sum(axis=0)
    return np.degrees(np.arctan2(nrm(P[1:]), P[0]))


def jac(S, form, h=1e-6):
    x0 = S.ravel(); n = x0.size; J = np.zeros((n, n))
    for i in range(n):
        d = np.zeros(n); d[i] = h
        J[:, i] = (step((x0 + d).reshape(S.shape), form).ravel() - step((x0 - d).reshape(S.shape), form).ravel()) / (2 * h)
    return J


def make_closed_imag(N=3):
    """虚部の和が 0、実部の和は 0 でない 3 個。大きさはばらばら。"""
    S = rng.normal(size=(N, 8))
    S[:, 1:] -= S[:, 1:].mean(axis=0)
    return S


print("=== (1) 一般の 3 個: 一回の演算で配列は変わる ===")
S = rng.normal(size=(3, 8))
for form in ("右", "真ん中", "左"):
    S1 = step(S, form)
    print("  %-4s 配列の変化 %.3f  グラム行列の変化 %.1e  総量の比 %.12f  P の角度 %.2f° → %.2f°"
          % (form, nrm(array3(S1) - array3(S)) / nrm(array3(S)), nrm(gram(S1) - gram(S)), (S1 ** 2).sum() / (S ** 2).sum(), angle_P(S), angle_P(S1)))

print("\n=== (2) 虚部の和が 0 の 3 個 (P が実数): 一回の演算で何も変わらない ===")
S = make_closed_imag()
P = S.sum(axis=0)
print("  P = %s" % np.round(P, 12))
for form in ("右", "真ん中", "左"):
    S1 = step(S, form)
    print("  %-4s 状態の変化 %.1e  配列の変化 %.1e" % (form, nrm(S1 - S), nrm(array3(S1) - array3(S))))

print("\n=== (3) その状態のまわりで、ずれが一回ごとに何倍になるか (24 方向) ===")
for form in ("右", "真ん中", "左"):
    ev = np.linalg.eigvals(jac(S, form))
    vals, cnt = np.unique(np.round(np.abs(ev), 4), return_counts=True)
    print("  %-4s 倍率と個数: %s" % (form, ", ".join("%.4f × %d" % (v, c) for v, c in zip(vals, cnt))))

print("\n=== (4) P の角度は、一回ごとに 3 倍になる (右の形) ===")
for th0 in (1.0, 10.0, 40.0, 45.0, 60.0, 90.0, 120.0):
    S = make_closed_imag(); n7 = rng.normal(size=7); n7 /= nrm(n7)
    Pr = S.sum(axis=0)[0]
    S[0, 1:] += n7 * abs(Pr) * np.tan(np.radians(th0)) if th0 != 90 else 0
    if th0 == 90:
        S[0, 0] -= Pr; S[0, 1:] += n7
    if Pr < 0 and th0 != 90:
        S[0, 0] -= 2 * Pr
    seq = [angle_P(S)]
    X0 = array3(S); dist = []
    for _ in range(6):
        S = step(S, "右"); seq.append(angle_P(S)); dist.append(nrm(array3(S) - X0) / nrm(X0))
    print("  初期 %6.2f°: 角度 " % th0 + " ".join("%7.2f" % a for a in seq) + " | 配列の初期からの離れ " + " ".join("%.2f" % d for d in dist))
print("  (角度は 0〜180° の範囲に折り返して表示している。3 倍して 360° を法として読み、180° を超えたら 360° から引く)")

print("\n=== (5) 丸め誤差だけから離れていくか (虚部の和が 0 の状態から、何も足さずに繰り返す) ===")
for form in ("右", "真ん中", "左"):
    S = make_closed_imag(); S0 = S.copy(); X0 = array3(S); out = []
    for it in range(1, 61):
        S = step(S, form)
        if it in (1, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60):
            out.append("%d:%.1e" % (it, nrm(array3(S) - X0) / nrm(X0)))
    print("  %-4s 配列の初期からの離れ (回数:値) " % form + "  ".join(out))
    print("        P の角度 %.3e° (最初) → 総量の比 %.12f" % (angle_P(S0), (S ** 2).sum() / (S0 ** 2).sum()))

print("\n=== (6) P = 0 (線形和が閉じた 3 個) ===")
S = rng.normal(size=(3, 8)); S -= S.mean(axis=0)
print("  |P|^2 = %.1e。割る数が 0 になり、演算が定義できない。" % (S.sum(axis=0) @ S.sum(axis=0)))
