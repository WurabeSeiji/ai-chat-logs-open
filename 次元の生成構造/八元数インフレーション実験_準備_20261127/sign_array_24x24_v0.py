#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
24x24 で、成分がすべて ±1、行どうし・列どうしが直交する配列を一つ作って確かめる。
どの軸が実数軸かは考えない。8x8 の S と同じ条件 (配列 × 転置 = 単位の配列の定数倍) だけを課す。

作り方は平方剰余による標準的な構成 (ペイリーの構成)。乱数は使わない。
  p = 23 (素数、4 で割って 3 余る)。chi(x) = 0 (x = 0)、+1 (x が 23 を法として平方数)、-1 (それ以外)。
  Q[i][j] = chi(j - i)  (23x23)。
  H = [[1, 1...1], [-1.., Q + I]] の形に縁を付けると 24x24 になる。
使い方: python3 sign_array_24x24_v0.py
"""
p = 23
squares = {(x * x) % p for x in range(1, p)}


def chi(x):
    x %= p
    return 0 if x == 0 else (1 if x in squares else -1)


n = p + 1
H = [[0] * n for _ in range(n)]
for j in range(n):
    H[0][j] = 1
for i in range(1, n):
    H[i][0] = -1
    for j in range(1, n):
        H[i][j] = chi((j - 1) - (i - 1)) + (1 if i == j else 0)

assert all(v in (1, -1) for r in H for v in r)
G = [[sum(H[i][t] * H[j][t] for t in range(n)) for j in range(n)] for i in range(n)]
Gc = [[sum(H[t][i] * H[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
ok_r = all(G[i][j] == (n if i == j else 0) for i in range(n) for j in range(n))
ok_c = all(Gc[i][j] == (n if i == j else 0) for i in range(n) for j in range(n))
print("24x24、成分はすべて ±1。")
print("行どうしが直交し、長さが同じ (H × H の転置 = 24 × 単位):", ok_r)
print("列どうしが直交し、長さが同じ (H の転置 × H = 24 × 単位):", ok_c)
print("H =")
for r in H:
    print("  " + " ".join("+" if v > 0 else "-" for v in r))
print("行の和:", [sum(r) for r in H])
print("列の和:", [sum(H[i][j] for i in range(n)) for j in range(n)])
