# 偶数定理の全数列挙検証（N1 論文 付録B 再現スクリプト）
# 等振幅アルファベット上の非自明閉鎖 sum(x_n^2) = 0 を全数列挙し、
# 閉形式（定理3.2 / 命題4.3）と突き合わせる。整数演算のみ（ガウス整数・円分整数）。
#
# アルファベット1（前提A・四半回転）: {1, i, -1, -i}     -> 二乗は {+1, -1}
# アルファベット2（対照・六分回転）  : {±1, ±ω, ±ω²}    -> 二乗は {1, ω, ω²}
#
# 実行: python3 N1_even_theorem_enumeration.py

from itertools import product
from math import comb, factorial


def count_z4(N):
    """四半回転アルファベット: x^2 ∈ {+1,-1} を (p,q) 計数で全数列挙。"""
    count = 0
    for combo in product([1, -1], repeat=N):  # 各項の二乗値
        if sum(combo) == 0:
            count += 2 ** N  # 各二乗値の逆像は2個（+1<-{1,-1}, -1<-{i,-i}）
    return count


def count_z4_closed(N):
    """閉形式: N偶数なら C(N, N/2)·2^N、奇数なら0。"""
    return comb(N, N // 2) * 2 ** N if N % 2 == 0 else 0


def count_z6(N):
    """六分回転アルファベット: x^2 ∈ {1, ω, ω²}。
    和が零 ⟺ 三種の個数が等しい（1+ω+ω²=0、{1,ω}はR上独立）。"""
    count = 0
    for combo in product([0, 1, 2], repeat=N):  # 二乗値の ω 指数
        a = combo.count(0)
        b = combo.count(1)
        c = combo.count(2)
        if a == b == c:
            count += 2 ** N  # 各二乗値の逆像は2個（±）
    return count


def count_z6_closed(N):
    """閉形式: N=3k なら N!/(k!)³·2^N、そうでなければ0。"""
    if N % 3 != 0:
        return 0
    k = N // 3
    return factorial(N) // (factorial(k) ** 3) * 2 ** N


def brute_force_z4(N):
    """独立検算: ガウス整数 (re,im) の直積で直接 sum(x^2)=0 を数える。"""
    alphabet = [(1, 0), (0, 1), (-1, 0), (0, -1)]  # 1, i, -1, -i

    def sq(z):
        a, b = z
        return (a * a - b * b, 2 * a * b)

    count = 0
    for tup in product(alphabet, repeat=N):
        s = (0, 0)
        for z in tup:
            q = sq(z)
            s = (s[0] + q[0], s[1] + q[1])
        if s == (0, 0):
            count += 1
    return count


if __name__ == "__main__":
    print("N | Z4列挙 | Z4閉形式 | Z4直接 | Z6列挙 | Z6閉形式")
    for N in range(1, 9):
        z4 = count_z4(N)
        z4c = count_z4_closed(N)
        z4b = brute_force_z4(N) if N <= 8 else None
        z6 = count_z6(N) if N <= 9 else None
        z6c = count_z6_closed(N)
        assert z4 == z4c, f"Z4 mismatch at N={N}"
        assert z4b is None or z4b == z4c, f"Z4 brute mismatch at N={N}"
        assert z6 is None or z6 == z6c, f"Z6 mismatch at N={N}"
        print(f"{N} | {z4} | {z4c} | {z4b} | {z6} | {z6c}")
    print("N=9 Z6:", count_z6(9), "=", count_z6_closed(9))
    print("all checks passed")
