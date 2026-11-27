#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
高対称な初期値の汎用初期化。成分は ±1 だけ。N は 4 の倍数。仕様は「思考実験.md」§31。

高対称の条件 (この三つだけ。行 = 元、列 = 成分):
  重心        : どの列も、足すと 0
  モーメント  : 列の長さがすべて同じ
  等方性      : 違う二列がすべて直交
検算: どの二列 (a, b) を選んで z_n = (a 列の成分) + i (b 列の成分) と読んでも
        Σ z_n = 0   かつ   Σ z_n² = 0
      が成り立つこと。Σ z_n² の実部 = 列 a と列 b の長さの二乗の差、虚部 = 2 × (列 a と列 b の内積) なので、
      全部の二列の組で Σ z_n² = 0 であることは、モーメントと等方性の条件そのものである。

出すもの: balanced_init(N) は、N 個の元 (行) × N-1 個の成分 (列) の表を返す。成分はすべて ±1。
  作り方: 成分が ±1 で行どうし・列どうしが直交する N×N の配列 (既知の構成) を作り、行の符号をそろえて
          一列目を全部 +1 にし、その列を取り除く。残った N-1 本の列は、足すと 0 で、互いに直交する。
  N×N の配列の作り方 (どれも式で決まる。乱数なし):
    [二進] N が 2 の冪: 行と列の番号を二進数で書き、同じ桁で両方が 1 になる個数が偶数なら +1、奇数なら -1。
    [平方剰余 1] N = p + 1、p は 4 で割って 3 余る素数: 番号の差が p を法として平方数かどうかで ±1。
    [平方剰余 2] N = 2(p + 1)、p は 4 で割って 1 余る素数。
    [積] N = N1 × N2 で、N1 と N2 のそれぞれが作れるとき: 二つの配列の成分どうしを全部の組で掛けて並べる。
  上のどれにも当たらない 4 の倍数 (例: 92) は、このプログラムでは作れない。その場合はそう報告して止まる。

使い方: python3 balanced_init_pm1_v0.py [N]      (N を省くと 8。表と検算結果を標準出力に出す)
        python3 balanced_init_pm1_v0.py --survey  (4 から 100 までの 4 の倍数について、作れるかと検算結果の一覧)
"""
import sys


def is_prime(x):
    if x < 2:
        return False
    k = 2
    while k * k <= x:
        if x % k == 0:
            return False
        k += 1
    return True


def jacobsthal(p):
    sq = {(x * x) % p for x in range(1, p)}
    chi = lambda x: 0 if x % p == 0 else (1 if x % p in sq else -1)
    return [[chi(j - i) for j in range(p)] for i in range(p)]


def kron(A, B):
    return [[a * b for a in ra for b in rb] for ra in A for rb in B]


def square_array(N):
    """成分 ±1、行どうし・列どうしが直交する N×N の配列と、使った作り方の名前。作れなければ (None, None)。"""
    if N == 1:
        return [[1]], "1"
    if N == 2:
        return [[1, 1], [1, -1]], "二進"
    if N % 4 != 0:
        return None, None
    if N & (N - 1) == 0:
        return [[-1 if bin(i & j).count("1") % 2 else 1 for j in range(N)] for i in range(N)], "二進"
    p = N - 1
    if is_prime(p) and p % 4 == 3:
        Q = jacobsthal(p)
        return [[1] * N] + [[-1] + [Q[i][j] + (1 if i == j else 0) for j in range(p)] for i in range(p)], "平方剰余 1 (p=%d)" % p
    if N % 2 == 0:
        p = N // 2 - 1
        if is_prime(p) and p % 4 == 1:
            Q = jacobsthal(p)
            n = p + 1
            C = [[0] + [1] * p] + [[1] + Q[i] for i in range(p)]
            I = lambda i, j: 1 if i == j else 0
            top = [[C[i][j] + I(i, j) for j in range(n)] + [C[i][j] - I(i, j) for j in range(n)] for i in range(n)]
            bot = [[C[i][j] - I(i, j) for j in range(n)] + [-C[i][j] - I(i, j) for j in range(n)] for i in range(n)]
            return top + bot, "平方剰余 2 (p=%d)" % p
    for a in range(2, N):
        if N % a == 0 and a <= N // a:
            A, na = square_array(a)
            B, nb = square_array(N // a)
            if A is not None and B is not None:
                return kron(A, B), "積 %d×%d [%s] × [%s]" % (a, N // a, na, nb)
    return None, None


def balanced_init(N):
    """N 個の元 × N-1 個の成分の表 (成分は ±1) と、作り方の名前。"""
    if N % 4 != 0:
        raise ValueError("N は 4 の倍数でなければならない: N = %d" % N)
    H, how = square_array(N)
    if H is None:
        raise ValueError("N = %d の作り方は、このプログラムに入っていない" % N)
    assert all(v in (1, -1) for r in H for v in r)
    assert all(sum(H[i][t] * H[j][t] for t in range(N)) == (N if i == j else 0) for i in range(N) for j in range(N))
    H = [[v * r[0] for v in r] for r in H]            # 行の符号をそろえて一列目を全部 +1 に
    return [r[1:] for r in H], how


def verify(T):
    """検算。全部の列について Σz = 0、全部の二列の組について Σz² = 0 を調べ、最大の大きさを返す。"""
    R, C = len(T), len(T[0])
    max_sum = max(abs(sum(T[n][a] for n in range(R))) for a in range(C))
    max_re = max_im = 0
    for a in range(C):
        for b in range(a + 1, C):
            re = sum(T[n][a] ** 2 - T[n][b] ** 2 for n in range(R))
            im = 2 * sum(T[n][a] * T[n][b] for n in range(R))
            max_re, max_im = max(max_re, abs(re)), max(max_im, abs(im))
    row_len = sorted(set(sum(v * v for v in r) for r in T))
    row_dot = sorted(set(sum(x * y for x, y in zip(T[i], T[j])) for i in range(R) for j in range(i + 1, R)))
    return max_sum, max_re, max_im, C * (C - 1) // 2, row_len, row_dot


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--survey":
        print("  N    作り方                                      Σz の最大   Σz² の実部・虚部の最大   判定")
        for N in range(4, 101, 4):
            try:
                T, how = balanced_init(N)
            except ValueError:
                print("  %3d  このプログラムでは作れない" % N)
                continue
            s, re, im, npairs, _, _ = verify(T)
            print("  %3d  %-42s  %d           %d, %d  (%d 組)          %s" % (N, how, s, re, im, npairs, "合格" if s == re == im == 0 else "不合格"))
        return
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    T, how = balanced_init(N)
    print("N = %d : %d 個の元 (行) × %d 個の成分 (列)。成分はすべて ±1。作り方: %s" % (N, N, N - 1, how))
    for n, r in enumerate(T, 1):
        print("  元%-3d " % n + " ".join("+1" if v > 0 else "-1" for v in r))
    s, re, im, npairs, row_len, row_dot = verify(T)
    print("検算")
    print("  Σ z_n   : 全 %d 列で、和の大きさの最大 = %d" % (N - 1, s))
    print("  Σ z_n²  : 全 %d 組の二列で、実部の大きさの最大 = %d、虚部の大きさの最大 = %d" % (npairs, re, im))
    print("  判定     : %s" % ("合格 (Σz = 0、Σz² = 0)" if s == re == im == 0 else "不合格"))
    print("参考 (条件ではない): 元の長さの二乗 = %s、違う二つの元の内積 = %s" % (row_len, row_dot))


if __name__ == "__main__":
    main()
