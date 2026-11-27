#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
8x8 の配列 S (成分はすべて ±1、行どうし・列どうしが直交) を初期値にして、
    O_n = O_1 × O_2 × O_3 × ... × O_{n-1}      (ふつうの配列の積。常に左から順に掛ける)
を n = 20 まで繰り返す。木原さんの指示 (2026-11-27)。仕様は「思考実験.md」§26。

  O_1 = S (初期値)。O_2 = O_1 (式のとおり、掛けるものが O_1 一つだけ)。O_3 = O_1 × O_2。以下同じ。
  割り算・正規化・乱数は入れない。成分は整数のまま厳密に計算する。

出力 (n ごと):
  [配列]  8x8 の全成分。120 桁以下なら整数そのまま。それを超えたら「仮数(17 桁) e 指数」で表示する
          (計算は最後まで整数で厳密。表示だけを短くしている)。
  [読出し] 各行を八元数 1 個として読み、o = A (cosθ + n sinθ) の A・θ・n を出す。今回は第一成分を実部と読む。
          A = sqrt(o0² + ... + o7²)、cosθ = o0 / A (θ は 0〜180 度)、n = (o1..o7) / sqrt(o1² + ... + o7²)。
          読出しは毎回、厳密な整数から一度だけ計算する (mpmath 60 桁)。演算には戻さない。
          虚部が全部 0 の行は n が決まらないので「不定」と出す。

使い方: python3 array_product_cumulative_S8_v0.py > results/array_product_cumulative_S8_v0.txt
"""
import sys

import mpmath as mp

sys.set_int_max_str_digits(0)
mp.mp.dps = 60
N_MAX = 20
FULL_DIGITS = 120

S = [[1, 1, 1, 1, 1, 1, 1, 1],
     [1, -1, 1, -1, 1, -1, -1, 1],
     [1, -1, -1, 1, 1, 1, -1, -1],
     [1, 1, -1, -1, 1, -1, 1, -1],
     [1, -1, -1, -1, -1, 1, 1, 1],
     [1, 1, -1, 1, -1, -1, -1, 1],
     [1, 1, 1, -1, -1, 1, -1, -1],
     [1, -1, 1, 1, -1, -1, 1, -1]]


def matmul(A, B):
    return [[sum(A[i][t] * B[t][j] for t in range(8)) for j in range(8)] for i in range(8)]


def ndigits(v):
    v = abs(v)
    return 1 if v == 0 else (len(str(v)) if v < 10 ** 400 else int(v.bit_length() * 0.30102999566398) + 1)


def show_int(v, full):
    if full:
        return str(v)
    if v == 0:
        return "0"
    return mp.nstr(mp.mpf(v), 17, min_fixed=0, max_fixed=0)


def fnum(x, d=17):
    return mp.nstr(x, d, min_fixed=-5, max_fixed=17, strip_zeros=True)


def readout(row):
    o = [mp.mpf(v) for v in row]
    A = mp.sqrt(mp.fsum(x * x for x in o))
    vn = mp.sqrt(mp.fsum(x * x for x in o[1:]))
    if A == 0:
        return "A=0", None, None, None
    th = mp.degrees(mp.atan2(vn, o[0]))
    n = None if vn == 0 else [x / vn for x in o[1:]]
    return A, th, n, vn


history = []                      # O_1, O_2, ...
rows_summary = []
for n in range(1, N_MAX + 1):
    if n == 1:
        O = [r[:] for r in S]
        how = "初期値 S"
    else:
        O = history[0]
        for k in range(1, n - 1):         # 左から順に O_1 × O_2 × ... × O_{n-1}
            O = matmul(O, history[k])
        how = " × ".join("O_%d" % k for k in range(1, n)) if n <= 6 else "O_1 × O_2 × ... × O_%d" % (n - 1)
    history.append(O)
    n_S = 1 if n <= 2 else 2 ** (n - 2)
    dg = max(ndigits(v) for r in O for v in r)
    full = dg <= FULL_DIGITS
    print("=" * 110)
    print("n = %d    O_%d = %s    (S を %d 個掛けたもの)    成分の最大桁数 %d" % (n, n, how, n_S, dg))
    print("=" * 110)
    print("[配列]" + ("" if full else "  (成分が %d 桁を超えるので、仮数 e 指数 で表示。計算は整数で厳密)" % FULL_DIGITS))
    for i, r in enumerate(O):
        print("  行%d: " % (i + 1) + "  ".join(show_int(v, full) for v in r))
    print("[読出し] o = A (cosθ + n sinθ)、第一成分を実部と読む")
    summ = []
    for i, r in enumerate(O):
        A, th, nv, vn = readout(r)
        if th is None:
            print("  行%d: A=0" % (i + 1))
            summ.append((mp.mpf(0), None))
            continue
        ns = "不定 (虚部が 0)" if nv is None else "(" + ", ".join(fnum(x, 10) for x in nv) + ")"
        print("  行%d: A = %s   θ = %s°   n = %s" % (i + 1, fnum(A), mp.nstr(th, 14, min_fixed=-5, max_fixed=6, strip_zeros=True), ns))
        summ.append((A, th))
    rows_summary.append((n, n_S, summ))
    sys.stdout.flush()

print("=" * 110)
print("まとめ: 各 n の振幅 A と角度 θ (行 1〜8)")
print("=" * 110)
for n, n_S, summ in rows_summary:
    As = sorted(set(fnum(a, 12) for a, _ in summ))
    ths = [("-" if t is None else mp.nstr(t, 8, min_fixed=-5, max_fixed=6, strip_zeros=True)) for _, t in summ]
    print("n=%2d (S %d 個): A の種類 = %s" % (n, n_S, ", ".join(As)))
    print("        θ[度] 行1〜8 = " + "  ".join(ths))
