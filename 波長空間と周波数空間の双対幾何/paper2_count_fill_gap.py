"""
論文2 補強: 4次元単位セル完全内接数 N0(R) と 充填率・ギャップ の厳密計算
-------------------------------------------------------------------
完全内接条件:
    sum_{i=1}^4 (|k_i| + 1/2)^2 <= R^2 ,   k in Z^4
両辺 4 倍 (整数演算化):
    sum_{i=1}^4 (2|k_i|+1)^2 <= (2R)^2     （正の奇数平方和）

量の定義:
    V0(R) = 外接 4 次元超球の体積       = (pi^2 / 2) R^4
    V1(R) = 内接単位超立方体の合計体積   = N0(R) * 1^4 = N0(R)
    充填率  fill = V1 / V0   ( < 1.0, R->∞ で 1.0 に漸近 )
    ギャップ gap  = V0 - V1   ( > 0, R とともに拡大; ~R^3 オーダ )
    2ρ(R) = 積み上げセル集合の最大対角線長 = sqrt( max(4-奇数平方和) <= (2R)^2 )

アルゴリズム:
  1 軸の奇数平方分布 d1 を作り、2 軸畳み込み g2 を作る。
  N0 は g2 同士の畳み込み和 (sum<=T) を累積和で O(|g2|) に集計。
  2ρ は g2 のキー集合に対する two-pointer で max(a+b)<=T を求める。
  すべて整数演算なので結果は厳密。R=0.5 刻みなら (2R)^2 は整数。

検算: N0(1,2,3,4,5) = 1, 9, 137, 473, 1545 （論文2 既知値）。
"""
from math import pi


def build_g2(T):
    """2 軸分の (奇数平方和 s -> 重複度) 分布 g2[s], s<=T を返す。"""
    d1 = {}
    m = 0
    while (2 * m + 1) ** 2 <= T:
        v = (2 * m + 1) ** 2
        d1[v] = d1.get(v, 0) + (1 if m == 0 else 2)   # k=0 は重複度1, k=±m は2
        m += 1
    items1 = sorted(d1.items())
    g2 = {}
    for v1, w1 in items1:
        for v2, w2 in items1:
            s = v1 + v2
            if s <= T:
                g2[s] = g2.get(s, 0) + w1 * w2
    return g2


def N0_and_diag(R):
    """完全内接セル数 N0(R) と 積み上げ対角線長 2ρ(R) を厳密に返す。"""
    T = int(round((2 * R) * (2 * R)))          # (2R)^2
    g2 = build_g2(T)
    if not g2:
        return 0, 0.0
    # --- N0: 4 軸畳み込み (sum <= T) を累積和で集計 ---
    maxs = max(g2)
    arr = [0] * (maxs + 1)
    for s, c in g2.items():
        arr[s] = c
    pref = [0] * (maxs + 2)
    for s in range(maxs + 1):
        pref[s + 1] = pref[s] + arr[s]

    def cum_le(x):
        if x < 0:
            return 0
        if x > maxs:
            x = maxs
        return pref[x + 1]

    N0 = 0
    for s1, c1 in g2.items():
        rem = T - s1
        if rem >= 0:
            N0 += c1 * cum_le(rem)
    # --- 2ρ: max(a+b)<=T over achievable 2-sums (two-pointer) ---
    keys = sorted(g2.keys())
    lo, hi = 0, len(keys) - 1
    best = 0
    while lo <= hi:
        s = keys[lo] + keys[hi]
        if s <= T:
            if s > best:
                best = s
            lo += 1
        else:
            hi -= 1
    # best = sum(2|k|+1)^2 = 4 * sum(|k|+1/2)^2 = (2ρ)^2  => 2ρ = sqrt(best)
    return N0, best ** 0.5


def V0(R):
    return (pi * pi / 2.0) * R ** 4


if __name__ == "__main__":
    # 検算
    known = {1: 1, 2: 9, 3: 137, 4: 473, 5: 1545}
    for R in [1, 2, 3, 4, 5]:
        n, _ = N0_and_diag(R)
        assert n == known[R], (R, n)
    print("checksum OK: N0(1..5) = 1, 9, 137, 473, 1545")

    print("\n# R, 2R, N0, 2rho, V0, fill=V1/V0, gap=V0-V1")
    rows = [i * 0.5 for i in range(1, 21)] + [100, 1000, 10000]
    for R in rows:
        n, two_rho = N0_and_diag(R)
        v0 = V0(R)
        print(f"R={R:<8} 2R={2*R:<8} N0={n:<20} 2rho={two_rho:<14.6g} "
              f"V0={v0:.6e} fill={n/v0:.8f} gap={v0-n:.6e}")
