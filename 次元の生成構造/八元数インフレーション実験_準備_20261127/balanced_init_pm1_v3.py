#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
高対称な初期データの作成。成分は +1 と -1 だけ。N は 4 の倍数。v3 (2026-11-27 の修正版)。

条件 (この三つだけ): 重心、モーメント、等方性。
検算 (これだけ)    : Σ z_n² = 0。
  表は、行 = 元 (n = 1..N)、列 = 成分。どの二つの成分 (a, b) を選んで z_n = (成分 a) + i (成分 b) と読んでも Σ_n z_n² = 0。
成分は区別しない。どの成分も同じ扱いで、「実部」という考えは初期データの中に無い。

返すもの: balanced_init(N) は、N 個の元 × N-1 個の成分の表 (成分はすべて ±1) を返す。これが条件を満たす初期データである。
  このプログラムが外へ渡すのは、この表だけである。

v2 までとの違い (修正した点):
  v1・v2 には、条件を満たす表を作る途中の N×N の配列 (全部の行で同じ符号になる列を含む。重心がずれている) を
  そのまま外へ渡す関数 square_array_shifted / square_array_ordered があり、掛け算のプログラムがそれを初期値に使っていた。
  これは Claude が無許可で足したもので、仕様に反する。v3 では、途中の配列を外へ渡す関数を無くした。
  あわせて、検算を Σ z_n² = 0 だけにした (v0〜v2 は Σ z_n = 0 も「検算」として出していたが、これは仕様に無い)。

作り方 (v0 と同じ。乱数なし):
  成分が ±1 の N×N の配列を既知の構成 (二進 / 平方剰余 1 / 平方剰余 2 / 積) で作り、
  行ごとに符号をそろえると全部の行で同じ符号になる列が一本できる。その列は足しても 0 にならないので取り除く。
  残った N-1 本が成分である。どの列を取り除く列に選んでも、条件を満たす表ができる (下の試験で全部確かめる)。
行の並べ替え: order に、ずらす量 (整数) か、並べ替えの順番 ("1,6,4,3,7,2,5,8" の形) を指定できる。行の集まりは変わらない。

使い方: python3 balanced_init_pm1_v3.py [N] [ずらす量 または 並べ替え]   (省くと N = 8、並べ替えなし)
        python3 balanced_init_pm1_v3.py --survey                        (4 から 100 までの 4 の倍数の一覧)
"""
import sys


def _is_prime(x):
    if x < 2:
        return False
    k = 2
    while k * k <= x:
        if x % k == 0:
            return False
        k += 1
    return True


def _jacobsthal(p):
    sq = {(x * x) % p for x in range(1, p)}
    chi = lambda x: 0 if x % p == 0 else (1 if x % p in sq else -1)
    return [[chi(j - i) for j in range(p)] for i in range(p)]


def _kron(A, B):
    return [[a * b for a in ra for b in rb] for ra in A for rb in B]


def _square(N):
    """作る途中の N×N の配列 (内部用。初期データとしては使わない) と、作り方の名前。"""
    if N == 1:
        return [[1]], "1"
    if N == 2:
        return [[1, 1], [1, -1]], "二進"
    if N % 4 != 0:
        return None, None
    if N & (N - 1) == 0:
        return [[-1 if bin(i & j).count("1") % 2 else 1 for j in range(N)] for i in range(N)], "二進"
    p = N - 1
    if _is_prime(p) and p % 4 == 3:
        Q = _jacobsthal(p)
        return [[1] * N] + [[-1] + [Q[i][j] + (1 if i == j else 0) for j in range(p)] for i in range(p)], "平方剰余 1 (p=%d)" % p
    p = N // 2 - 1
    if _is_prime(p) and p % 4 == 1:
        Q = _jacobsthal(p)
        n = p + 1
        C = [[0] + [1] * p] + [[1] + Q[i] for i in range(p)]
        I = lambda i, j: 1 if i == j else 0
        top = [[C[i][j] + I(i, j) for j in range(n)] + [C[i][j] - I(i, j) for j in range(n)] for i in range(n)]
        bot = [[C[i][j] - I(i, j) for j in range(n)] + [-C[i][j] - I(i, j) for j in range(n)] for i in range(n)]
        return top + bot, "平方剰余 2 (p=%d)" % p
    for a in range(2, N):
        if N % a == 0 and a <= N // a:
            A, na = _square(a)
            B, nb = _square(N // a)
            if A is not None and B is not None:
                return _kron(A, B), "積 %d×%d [%s] × [%s]" % (a, N // a, na, nb)
    return None, None


def _parse_order(N, spec):
    if spec is None:
        return list(range(N))
    if isinstance(spec, int) or "," not in str(spec):
        k = int(spec) % N
        return [(i + k) % N for i in range(N)]
    order = [int(v) - 1 for v in str(spec).split(",")]
    if sorted(order) != list(range(N)):
        raise ValueError("並べ替えは 1 から %d までの数を一つずつ、%d 個書く: %s" % (N, N, spec))
    return order


def balanced_init(N, order=None, drop=0):
    """高対称な初期データ: N 個の元 × N-1 個の成分の表 (成分は ±1) と、作り方の名前。
    order: 行の並べ替え。drop: 作る途中の配列で、取り除く列に選ぶ番号 (0 始まり。どれを選んでも条件を満たす)。"""
    if N % 4 != 0:
        raise ValueError("N は 4 の倍数でなければならない: N = %d" % N)
    H, how = _square(N)
    if H is None:
        raise ValueError("N = %d の作り方は、このプログラムに入っていない" % N)
    H = [[v * r[drop] for v in r] for r in H]            # 行ごとに符号をそろえる。drop 列目が全部の行で同じ符号になる
    T = [[v for j, v in enumerate(r) if j != drop] for r in H]   # その列を取り除く
    od = _parse_order(N, order)
    if od != list(range(N)):
        T = [T[od[i]] for i in range(N)]
        how += "、行の並べ替え " + ",".join(str(v + 1) for v in od)
    if not check_conditions(T)[0]:
        raise AssertionError("作った表が条件を満たしていない (プログラムの誤り)")
    return T, how


def check_conditions(T):
    """条件の確認。重心 (各成分の和が 0)、モーメント (各成分の二乗の和が同じ)、等方性 (Σz² = 0 と同じ内容)。"""
    R, C = len(T), len(T[0])
    centroid = all(sum(T[n][a] for n in range(R)) == 0 for a in range(C))
    moment = len({sum(T[n][a] ** 2 for n in range(R)) for a in range(C)}) == 1
    iso = verify(T)[0] == 0
    return centroid and moment and iso, centroid, moment, iso


def verify(T):
    """検算: 全部の二成分の組 (a, b) について Σ_n z_n² (z_n = 成分 a + i 成分 b)。大きさの最大と、組の数を返す。"""
    R, C = len(T), len(T[0])
    worst = 0
    for a in range(C):
        for b in range(a + 1, C):
            re = sum(T[n][a] ** 2 - T[n][b] ** 2 for n in range(R))
            im = 2 * sum(T[n][a] * T[n][b] for n in range(R))
            worst = max(worst, abs(re), abs(im))
    return worst, C * (C - 1) // 2


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--survey":
        print("  N    作り方                                      検算 Σz² = 0")
        for N in range(4, 101, 4):
            try:
                T, how = balanced_init(N)
            except ValueError:
                print("  %3d  このプログラムでは作れない" % N)
                continue
            w, npairs = verify(T)
            print("  %3d  %-42s  %s (全 %d 組、最大 %d)" % (N, how, "合格" if w == 0 else "不合格", npairs, w))
        return
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    order = sys.argv[2] if len(sys.argv) > 2 else None
    T, how = balanced_init(N, order)
    print("N = %d : %d 個の元 (行) × %d 個の成分 (列)。成分はすべて ±1。作り方: %s" % (N, N, N - 1, how))
    for n, r in enumerate(T, 1):
        print("  元%-3d " % n + " ".join("+1" if v > 0 else "-1" for v in r))
    ok, c, m, iso = check_conditions(T)
    print("条件: 重心 %s / モーメント %s / 等方性 %s" % tuple("満たす" if x else "満たさない" for x in (c, m, iso)))
    print("      各成分の +1 と -1 の個数: " + " ".join("%d:%d" % (sum(1 for r in T if r[a] > 0), sum(1 for r in T if r[a] < 0)) for a in range(N - 1)))
    w, npairs = verify(T)
    print("検算: Σ z_n² = 0  ->  全 %d 組の二成分で、Σ z_n² の大きさの最大 = %d  ->  %s" % (npairs, w, "合格" if w == 0 else "不合格"))
    # どの列を取り除く列に選んでも同じか (成分を区別していないことの確認)
    res = []
    for drop in range(N):
        Td, _ = balanced_init(N, order, drop)
        res.append(check_conditions(Td)[0] and verify(Td)[0] == 0)
    print("成分を区別していないことの確認: 取り除く列を 1〜%d 番目のどれに選んでも条件を満たし検算に合格するか -> %s"
          % (N, "はい (全 %d 通り)" % N if all(res) else "いいえ: " + str(res)))


if __name__ == "__main__":
    main()
