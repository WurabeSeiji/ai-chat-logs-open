#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G9d: 巡回頂点代数の純数学部——到達グラフの分類（査読提案の実施）

対象: Z_n 上の生成規則 倍加 x→2x（海頂点）と交差 (x,y)→2x−y。

定理候補（本スクリプトで全数検証）:
 (M1) 無次数閉包の自明化: {0,a} の (x,y)→2x−y 閉包 = 部分群 ⟨gcd(a,n)⟩。
      ゆえに開閉の物理は「次数つき（graded）到達」にのみ宿る。
 (M2) 単元同変性: u∈Z_n^* は両生成規則と可換（u(2x−y)=2(ux)−(uy)）
      ——E(u·m_A,u·m_B)=E の代数的根拠。
 (M3) 海接続の判定: 倍加軌道 a→2a→… が 0 に到達 ⟺ q|a（q=nの奇数部分）。
      系: n=2^s なら全電荷が有限次数（≤s）で海に接続（n/2 経由）。
      n に奇数因子があれば q∤a の電荷は永続セクター（周期軌道・超選択則）。
 (M4) 倍加合流次数（n=2^s）: 2^{k}a≡2^{k}b (mod 2^s) ⟺ k ≥ s−v₂(a−b)。
      最短合流次数 = s−v₂(a−b)（各側）——2進付値が開度の次数を決める。
 (M5) 次数単調性: G9c 実測開度は graded 到達次数の単調減少関数と整合。
使い方: python3 run_g9d_cyclic_vertex_algebra_math_v1.py
"""
from __future__ import annotations
import json, time
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def closure_pair(a, b, n):
    """{0,a,b} の (x,y)->2x-y 閉包"""
    S = {0, a % n, b % n}
    changed = True
    while changed:
        changed = False
        cur = list(S)
        for x in cur:
            for y in cur:
                z = (2 * x - y) % n
                if z not in S:
                    S.add(z)
                    changed = True
    return S


def subgroup(g, n):
    g = gcd(g, n)
    return {(g * i) % n for i in range(n // g)} if g else {0}


def doubling_orbit(a, n, kmax=None):
    """a, 2a, 4a, ... 最初の再訪まで。0到達なら(到達次数, True)"""
    seen = {}
    x = a % n
    k = 0
    while x not in seen:
        seen[x] = k
        if x == 0:
            return k, True, sorted(seen)
        x = (2 * x) % n
        k += 1
    return None, False, sorted(seen)


def v2(x):
    if x == 0:
        return 10 ** 9
    k = 0
    while x % 2 == 0:
        x //= 2
        k += 1
    return k


def graded_reach(seeds, n, kmax):
    """次数つき到達: R[k] = k回以下の頂点適用で得られる元（海0は常時使用可）"""
    R = [set(seeds) | {0}]
    for _ in range(kmax):
        cur = R[-1]
        nxt = set(cur)
        for x in cur:
            for y in cur:
                nxt.add((2 * x - y) % n)
        R.append(nxt)
    return R


def main():
    t0 = time.time()
    results = {}

    # (M1) 無次数閉包 = 部分群（全 n=3..64, 全 a,b）
    m1_ok = True
    for n in range(3, 65):
        for a in range(1, n):
            S = closure_pair(a, a, n)  # {0,a} の閉包
            if S != subgroup(a, n):
                m1_ok = False
                print(f"M1反例: n={n} a={a}")
    print(f"(M1) 無次数閉包={{0,a}}→⟨gcd(a,n)⟩: {'全数成立 (n≤64)' if m1_ok else '不成立'}")

    # (M2) 単元同変性（恒等式なので構造チェックのみ・n=16で全数）
    n = 16
    m2_ok = all((u * (2 * x - y)) % n == (2 * (u * x) - (u * y)) % n
                for u in range(1, n, 2) for x in range(n) for y in range(n))
    print(f"(M2) 単元同変性 u(2x−y)=2(ux)−(uy): {'全数成立' if m2_ok else '不成立'}")

    # (M3) 海接続 ⟺ q|a（全 n=3..64, 全 a）
    m3_ok = True
    sectors = {}
    for n in range(3, 65):
        q = odd_part(n)
        for a in range(1, n):
            k0, reaches, orb = doubling_orbit(a, n)
            if reaches != (a % q == 0):
                m3_ok = False
                print(f"M3反例: n={n} a={a}")
        sectors[n] = {"odd_part": q,
                      "persistent": sorted(a for a in range(1, n) if a % q != 0)}
    print(f"(M3) 0到達⟺q|a（q=奇数部分）: {'全数成立 (n≤64)' if m3_ok else '不成立'}")
    print(f"     例 n=16: 永続セクター{sectors[16]['persistent']}（空=全電荷が海に溶ける）")
    print(f"     例 n=12: 永続セクター{sectors[12]['persistent']}（3の倍数以外は海に届かない）")

    # (M4) 合流次数（n=2^s）: min k s.t. 2^k a ≡ 2^k b ⟺ k=s−v2(a−b)
    m4_ok = True
    for s in (3, 4, 5, 6):
        n = 2 ** s
        for a in range(1, n):
            for b in range(1, n):
                if a == b:
                    continue
                pred = max(0, s - v2((a - b) % n))
                k = next((kk for kk in range(0, s + 1)
                          if (2 ** kk * a) % n == (2 ** kk * b) % n), None)
                if k != pred:
                    m4_ok = False
    print(f"(M4) 倍加合流次数 = s−v₂(a−b)（n=2^s, s≤6）: {'全数成立' if m4_ok else '不成立'}")

    # (M5) graded到達次数 vs G9c実測開度の単調性（ne=16）
    n = 16
    MEAS = {(1, 2): 0.69, (3, 6): 0.69, (2, 4): 0.22, (1, 4): 0.03,
            (1, 3): 0.006, (1, 5): 0.006, (3, 4): 0.003}  # 開度/|coh|（盲目は符号差の半分程度）
    print("\n対 (a,b): 最小次数 k*（A系とB系のgraded到達が交わる最小合計次数）と実測開度")
    orders = {}
    for (a, b) in MEAS:
        RA = graded_reach([a], n, 6)
        RB = graded_reach([b], n, 6)
        k_star = None
        for kt in range(0, 13):
            for ka in range(0, min(kt, 6) + 1):
                kb = kt - ka
                if kb > 6:
                    continue
                inter = (RA[ka] & RB[kb]) - {0}
                if inter:
                    k_star = kt
                    break
            if k_star is not None:
                break
        orders[f"{a},{b}"] = k_star
        print(f"  ({a},{b}): k*={k_star}  実測開度={MEAS[(a,b)]}")
    # 単調性: k*が小さいほど開度が大きい（同k*内の分裂=約数類・別法則）
    pairs_sorted = sorted(MEAS, key=lambda p: MEAS[p], reverse=True)
    ks = [orders[f"{p[0]},{p[1]}"] for p in pairs_sorted]
    m5_ok = all(ks[i] <= ks[j] for i in range(len(ks)) for j in range(i + 1, len(ks))
                if abs(MEAS[pairs_sorted[i]] - MEAS[pairs_sorted[j]]) > 0.05)
    print(f"(M5) 開度の次数単調性（開度差>0.05の対で）: {'成立' if m5_ok else '不成立'}")

    verdict = ("巡回頂点代数の数学部成立: 開閉則は Z_ne の代数定理（graded到達）・"
               "海接続⟺2べき・奇数因子は超選択セクターを生む"
               if (m1_ok and m2_ok and m3_ok and m4_ok and m5_ok) else "要精査")
    print(verdict)
    out = {"M1": bool(m1_ok), "M2": bool(m2_ok), "M3": bool(m3_ok),
           "M4": bool(m4_ok), "M5": bool(m5_ok),
           "persistent_sectors_n12": sectors[12]["persistent"],
           "persistent_sectors_n16": sectors[16]["persistent"],
           "graded_orders_n16": orders, "measured_openings": {f"{a},{b}": v for (a, b), v in MEAS.items()},
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "result_g9d_cyclic_vertex_algebra_math_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
