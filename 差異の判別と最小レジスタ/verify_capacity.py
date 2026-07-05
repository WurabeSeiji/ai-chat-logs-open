#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
[5] 最終段：容量の数え上げ（定理側の限界までの仕事）

問い: 各 W8 配置型が担い得る巻き数の集合（容量）を、quadrature 実装の候補
スキームごとに組合せ論で列挙し、W8 のスピン表が容量内に無矛盾に収まるかを検定する。
（どの軸がどの粒子か、の割り当ては解釈問題であり本計算の対象外）

quadrature 実装スキーム（組合せの可能性の列挙、D4 計算の帰結より
「1 quadrature 単位 = ±π/2 位相対 = 巻き ±1」）:
  スキームA: 活性な動的軸1本が方位角と対になる —— 各軸が c ∈ {−1, 0, +1} を寄与
  スキームB: 活性な動的軸同士が対を組む —— floor(n/2) 単位まで
フェルミオン（k=2）: 実符号対＋対蹠同一視の自動 1/2（D4 計算で導出済み）が加算

W8 のスピン表（原典 定義2.1・§4-6）:
  ボソン: s = n（n = {t, R, Q遷移} の活性本数、0..3。スピン2は tR/tQ/RQ の3型、スピン3は tRQ）
  フェルミオン: s = n + 1/2（n = {R, Q遷移} の活性本数、0..2）
"""
from itertools import product

def reachable_sums(n, choices=(-1, 0, 1)):
    """n 本の軸が各々 c∈choices を寄与するときの到達可能な巻き数の集合"""
    if n == 0:
        return {0}
    return set(sum(c) for c in product(choices, repeat=n))

def fmt(s):
    return "{" + ", ".join(f"{x:+g}" if x else "0" for x in sorted(s)) + "}"

def main():
    print("=== 容量の数え上げ（スキームA: 各活性軸×方位角 / スキームB: 軸同士の対）===")
    print()
    print("--- ボソン型（k≠2）---")
    boson_types = [
        (0, "（活性なし）", 0),
        (1, "{t} / {R} / {Q遷移}（3型）", 1),
        (2, "{tR} / {tQ} / {RQ}（3型）", 2),
        (3, "{tRQ}（1型）", 3),
    ]
    print(f"{'n':>2} {'型':<24} {'容量A':<22} {'max A':>5} {'容量B':<12} {'max B':>5} {'W8のs':>5} {'A整合':>5} {'B整合':>5}")
    okA = okB = True
    for n, label, s_w8 in boson_types:
        capA = reachable_sums(n)
        capB = reachable_sums(n // 2)
        inA = s_w8 in capA and s_w8 == max(capA)
        inB = s_w8 in capB and s_w8 == max(capB)
        okA &= inA
        okB &= inB
        print(f"{n:>2} {label:<24} {fmt(capA):<22} {max(capA):>+5} {fmt(capB):<12} {max(capB):>+5} {s_w8:>5} {'✓' if inA else '✗':>5} {'✓' if inB else '✗':>5}")
    print()
    print("--- フェルミオン型（k=2、自動 1/2 加算：D4計算で導出済み）---")
    print(f"{'n':>2} {'容量A（+1/2 込み）':<30} {'max A':>7} {'W8のs':>6} {'A整合':>5}")
    for n in [0, 1, 2]:
        cap = {m + 0.5 for m in reachable_sums(n)} | {m - 0.5 for m in reachable_sums(n)}
        s_w8 = n + 0.5
        inA = s_w8 in cap and s_w8 == max(cap)
        okA &= inA
        capstr = "{" + ", ".join(f"{x:+.1f}" for x in sorted(cap)) + "}"
        print(f"{n:>2} {capstr:<30} {max(cap):>+7.1f} {s_w8:>6.1f} {'✓' if inA else '✗':>5}")
    print()
    print("=== 判定 ===")
    print(f"スキームA（各活性軸×方位角）: {'W8 の全スピン表が容量の最大値として整合 ✓' if okA else '不整合 ✗'}")
    print(f"スキームB（軸同士の対）:      {'整合' if okB else '不整合 ✗（スピン2で max 1、スピン3で max 1 となり W8 の表に届かない）'}")
    print()
    print("=== 付随する観察（解釈でなく容量集合の形の記録）===")
    for n in [1, 2]:
        print(f"n={n} の容量集合 {fmt(reachable_sums(n))} は、質量を持つスピン{n}粒子の"
              f"ヘリシティ多重項 {{0, ±1{', ±2' if n==2 else ''}}} と同じ集合")

if __name__ == "__main__":
    main()
