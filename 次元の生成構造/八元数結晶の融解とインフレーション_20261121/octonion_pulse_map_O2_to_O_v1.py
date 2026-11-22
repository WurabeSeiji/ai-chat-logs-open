#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N 個に戻す写像 O^2 |-> O (候補 1) : Q_i' = Q_i · (Σ_j Q_j)

仕様は「思考実験2_ヒッグス場のQ位相ベクトルと頂点の検算_20261122.md」の §21。
この写像が正しいかどうかは決まっていない。書いて計算した最初の候補として保存する。

状態     : 円周 ψ∈[0,360°) の上の八元数の波形 Q_i(ψ)。NG 点で持つ (配列 NG×8)。
           倍音 m=1..M を等しい振幅で積んだパルス
             Q(ψ) = (A/M) Σ_m [cos α_m + n sin α_m],  α_m = φ0 + m(ψ−θ)
           θ=中心位相(位置)、n=7 次元の単位の虚数の軸、M=分解能(幅 360°/M)、φ0=位相のずれ。
           A=√M と置くと、パワー(円周平均の |Q|²)が 1 になる。
写像     : Q_i'(ψ) = Q_i(ψ) · P(ψ),  P(ψ) = Σ_j Q_j(ψ)   (八元数の積。円周の各点で掛ける)
           全部の組 Q_i·Q_j (自分との積 j=i を含む) を作って j について足したものと同じ。
           外から与える数、重みの配列、相手を選ぶ手続き、正規化は、ない。
読み出し : パワー E_i = 円周平均の |Q_i|²
           自分の変化 h_i = <Q_i, Q_i'>/E_i   (<A,B> = 円周平均の conj(A)·B。値は八元数)
           結合係数 k_ij = <Q_i, Q_i·Q_j>/E_i  (Σ_j k_ij = h_i)
           八元数 o の「大きさ」= |o|、「回転角」= atan2(|虚部|, 実部)
格子点数 : 円周平均が正確になるのは、積の次数の合計が NG より小さいときだけ。
           M=24 の三つの積は 72 次、M=48 を含むと 96 次になる。NG=64 では足りない。
           既定は NG=1024。

使い方   : python3 octonion_pulse_map_O2_to_O_v1.py [NG]
"""
import sys
import numpy as np

NG = int(sys.argv[1]) if len(sys.argv) > 1 else 1024
rng = np.random.default_rng(5)          # 軸 n を作る乱数。種は固定
psi = 2 * np.pi * np.arange(NG) / NG

# ---- 八元数の積の表 -------------------------------------------------------
# e_i e_j = -δ_ij + f_ijk e_k。f は完全反対称で、次の 7 組で +1。
TRI = [(1, 2, 3), (1, 4, 5), (1, 7, 6), (2, 4, 6), (2, 5, 7), (3, 4, 7), (3, 6, 5)]
f = np.zeros((8, 8, 8))
for a, b, c in TRI:
    for (i, j, k), s in [((a, b, c), 1), ((b, c, a), 1), ((c, a, b), 1),
                         ((b, a, c), -1), ((a, c, b), -1), ((c, b, a), -1)]:
        f[i, j, k] = s
T = np.zeros((8, 8, 8))                 # (xy)_k = T[i,j,k] x_i y_j
T[0, 0, 0] = 1
for i in range(1, 8):
    T[0, i, i] = 1
    T[i, 0, i] = 1
    T[i, i, 0] = -1
T[1:, 1:, 1:] += f[1:, 1:, 1:]


def mul(x, y):
    """八元数の積。先頭の軸(円周の点)ごとに掛ける。"""
    return np.einsum('ijk,...i,...j->...k', T, x, y)


def conj(x):
    return x * np.r_[1, -np.ones(7)]


# ---- 状態 -----------------------------------------------------------------
def axis():
    v = rng.normal(size=7)
    return v / np.linalg.norm(v)


def pulse(theta, n, M, phi0=0.0, A=1.0):
    """倍音 1..M を等しい振幅で積んだ八元数のパルス。"""
    Q = np.zeros((NG, 8))
    for m in range(1, M + 1):
        al = phi0 + m * (psi - theta)
        Q[:, 0] += np.cos(al)
        Q[:, 1:] += np.outer(np.sin(al), n)
    return A * Q / M


def unit_pulse(theta_deg, n, M, phi0=0.0):
    """パワー 1 のパルス。"""
    return pulse(np.radians(theta_deg), n, M, phi0) * np.sqrt(M)


# ---- 写像 -----------------------------------------------------------------
def F(states):
    """写像 O^2 |-> O : Q_i' = Q_i · (Σ_j Q_j)。右から掛ける。"""
    P = sum(states)
    return [mul(S, P) for S in states]


def F_left(states):
    """比較用 : Q_i' = (Σ_j Q_j) · Q_i。左から掛ける。"""
    P = sum(states)
    return [mul(P, S) for S in states]


# ---- 読み出し ---------------------------------------------------------------
def power(A):
    return (A ** 2).sum() / NG


def inner(A, B):
    """<A,B> = 円周平均の conj(A)·B。値は八元数。"""
    return mul(conj(A), B).sum(axis=0) / NG


def read(o):
    """八元数 → (大きさ, 回転角[度])"""
    return np.linalg.norm(o), np.degrees(np.arctan2(np.linalg.norm(o[1:]), o[0]))


def coupling(Si, Sj):
    """k_ij = <Q_i, Q_i·Q_j>/E_i"""
    return inner(Si, mul(Si, Sj)) / power(Si)


def max_order(X):
    """波形に入っている倍音の次数の上限。"""
    pw = (np.abs(np.fft.fft(X, axis=0) / NG) ** 2).sum(axis=1)
    kk = np.minimum(np.arange(NG), NG - np.arange(NG))
    return kk[pw > 1e-12 * pw.max()].max()


# ---- 計算 -----------------------------------------------------------------
def main():
    print("格子点数 NG = %d" % NG)
    M = 6
    th = [40, 75, 200, 310]
    ax = [axis() for _ in th]
    ph = [0.3, 1.1, 2.0, 0.7]
    S = [unit_pulse(t, a, M, p) for t, a, p in zip(th, ax, ph)]
    S2 = F(S)
    print("=== 写像 Q_i' = Q_i·(Σ_j Q_j) を一回 ===")
    print("パルスの位置(度):", th, " 分解能 M=%d (幅 %d°)" % (M, 360 // M))
    print("位相のずれ φ0:", ph)
    for i, a in enumerate(ax):
        print("   軸 n_%d = [%s]" % (i + 1, ", ".join("%+.4f" % v for v in a)))

    print("\n(1) 各状態のパワー 前→後、 自分の変化の読み <Q_i,Q_i'>/E_i (大きさ, 回転角)")
    for i in range(4):
        r, a = read(inner(S[i], S2[i]) / power(S[i]))
        print("   状態%d: パワー %.4f → %.4f   読み 大きさ %.4f 回転角 %7.2f°"
              % (i + 1, power(S[i]), power(S2[i]), r, a))
    print("   全体のパワー: %.4f → %.4f" % (sum(map(power, S)), sum(map(power, S2))))

    print("\n(2) 結合係数 k_ij = <Q_i, Q_i·Q_j>/E_i  (i が j から受ける分。大きさ / 回転角)")
    K = np.zeros((4, 4, 8))
    for i in range(4):
        row = []
        for j in range(4):
            K[i, j] = coupling(S[i], S[j])
            row.append("%.4f/%6.1f°" % read(K[i, j]))
        print("   i=%d: " % (i + 1) + "   ".join(row))
    print("   Σ_j k_ij と (1) の読みの差: %.1e"
          % max(np.abs(K[i].sum(axis=0) - inner(S[i], S2[i]) / power(S[i])).max() for i in range(4)))
    dens = lambda A: (A ** 2).sum(axis=1)
    print("   k_ij と <|Q_i|², Q_j>/E_i (自分の密度で相手の波形を平均したもの) の差: %.1e"
          % max(np.abs(K[i, j] - (dens(S[i])[:, None] * S[j]).sum(axis=0) / NG / power(S[i])).max()
                for i in range(4) for j in range(4)))

    print("\n(3) 二体で距離と分解能を変える: |k_12| (相手から受ける分) と |k_11| (自分から受ける分)")
    n1 = np.eye(7)[0]
    n2 = np.eye(7)[1]
    for M in (6, 24):
        a = unit_pulse(0, n1, M)
        print("   分解能 M=%d (幅 %d°)  |k_11|=%.4f" % (M, 360 // M, np.linalg.norm(coupling(a, a))))
        for d in (5, 15, 30, 60, 120, 180):
            b = unit_pulse(d, n1, M)
            c = unit_pulse(d, n2, M)
            print("     距離 %3d°: 同じ軸 %.4f   直交する軸 %.4f"
                  % (d, np.linalg.norm(coupling(a, b)), np.linalg.norm(coupling(a, c))))

    print("\n(4) 相手の大きさだけを変える(自分は M=24、距離 30°)")
    a = unit_pulse(0, n1, 24)
    for Mb in (3, 6, 12, 24, 48):
        b = unit_pulse(30, n1, Mb)
        print("     相手の分解能 %2d (幅 %5.1f°): |k_12| = %.4f"
              % (Mb, 360 / Mb, np.linalg.norm(coupling(a, b))))

    print("\n(5) 左から掛ける形 Q_i' = (Σ_j Q_j)·Q_i との違い")
    L = F_left(S)
    print("   Σ_i Q_i' の差: %.1e   個々の状態の差の最大: %.4f"
          % (np.abs(sum(L) - sum(S2)).max(), max(np.abs(l - r).max() for l, r in zip(L, S2))))
    print("   パワーの差の最大: %.1e" % max(abs(power(l) - power(r)) for l, r in zip(L, S2)))

    print("\n(6) 繰り返したとき(戻すたびに何が変わるか)")
    X = [s.copy() for s in S]
    for it in range(4):
        print("   %d 回: 各パワー" % it, ["%.3e" % power(x) for x in X],
              " 全体 %.4e" % sum(map(power, X)), " 次数の上限", max_order(sum(X)))
        X = F(X)


if __name__ == "__main__":
    main()
