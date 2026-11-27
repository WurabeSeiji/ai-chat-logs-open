#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
8x8 の配列 S (成分はすべて ±1、行どうし・列どうしが直交) を初期値にして、
    O_n = O_1 × O_2 × O_3 × ... × O_{n-1}      (ふつうの配列の積。常に左から順に掛ける)
を n = 20 まで繰り返す。v1。演算は v0 と同じで、読出しを指数部の形に変えた (木原さんの指示、2026-11-27)。
仕様は「思考実験.md」§26 と §27。

  O_1 = S (初期値)。O_2 = O_1 (式のとおり、掛けるものが O_1 一つだけ)。O_3 = O_1 × O_2。以下同じ。
  3 回目以降は毎回「前回の配列の二乗」になる。割り算・正規化・乱数は入れない。成分は整数のまま厳密に計算する。

出力 (n ごと):
  [配列]      8x8 の全成分。120 桁以下なら整数そのまま。それを超えたら「仮数(17 桁) e 指数」で表示する
              (計算は最後まで整数で厳密。表示だけを短くしている)。
  [行の指数部] 各行を八元数 1 個として o = exp(a + nθ) と読み、指数部の 8 成分を出す。今回は第一成分を実部と読む。
                a  = ln A = ln sqrt(o0² + ... + o7²)                  (実の成分。振幅の対数)
                θn = θ × (o1..o7) / sqrt(o1² + ... + o7²)             (虚の 7 成分。θ はラジアン)
              θ は cosθ = o0 / A で決まる 0〜π の値 (主値) である。一つの行の成分だけからは、
              何周したかは分からないので、行ごとの θ は折り返した値になる。
              虚部が全部 0 の行は向き n が決まらない。θ = 0 なら θn = 0、θ = π なら「向き不定」と出す。
  [配列の指数部] 配列全体については、折り返さない値が出せる。O_n は S を m 個掛けたものなので、
              S の固有値 λ = |λ| exp(iφ) を使って、指数部は m × (ln|λ| + iφ) である。
              S の固有値は +√8 (1 個)、-√8 (1 個)、-1 ± i√7 (3 個ずつ) で、大きさはどれも √8。
              φ は順に 0、π、± arccos(-1/√8)。これを m 倍したものを、折り返さずに出す。
              (固有値は特性多項式 (x² - 8)(x² + 2x + 8)³ から。§26.3)
  読出しは毎回、厳密な整数から一度だけ計算する (mpmath 60 桁)。演算には戻さない。

使い方: python3 array_product_cumulative_S8_v1.py > results/array_product_cumulative_S8_v1.txt
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


def exponent(row):
    """行 (o0..o7) の指数部 (a, θ, θn の 7 成分 または None)。"""
    o = [mp.mpf(v) for v in row]
    A = mp.sqrt(mp.fsum(x * x for x in o))
    vn = mp.sqrt(mp.fsum(x * x for x in o[1:]))
    if A == 0:
        return None, None, None
    a = mp.log(A)
    th = mp.atan2(vn, o[0])
    if vn == 0:
        return a, th, ([mp.mpf(0)] * 7 if th == 0 else None)
    return a, th, [th * x / vn for x in o[1:]]


LN_R8 = mp.log(mp.sqrt(8))
PHI = mp.acos(-1 / mp.sqrt(8))            # 固有値 -1 ± i√7 の位相 (ラジアン)


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
    print("[行の指数部] o = exp(a + nθ)、第一成分を実部と読む。a = ln A、θ は主値 (0〜π)")
    summ = []
    for i, r in enumerate(O):
        a, th, tn = exponent(r)
        if a is None:
            print("  行%d: 全成分 0" % (i + 1))
            summ.append((None, None))
            continue
        tns = "向き不定 (虚部が 0 で θ = π)" if tn is None else "(" + ", ".join(fnum(x, 10) for x in tn) + ")"
        print("  行%d: a = %s   θ = %s rad (%s°)   θn = %s"
              % (i + 1, fnum(a), fnum(th, 14), mp.nstr(mp.degrees(th), 10, min_fixed=-5, max_fixed=6, strip_zeros=True), tns))
        summ.append((a, th))
    m = n_S
    print("[配列の指数部] S を m = %d 個掛けたもの。固有値ごとに m × (ln|λ| + iφ)。折り返さない値" % m)
    print("  実の成分 (全固有値に共通): m × ln√8 = %s" % fnum(m * LN_R8))
    for name, mult, phi in (("+√8", 1, mp.mpf(0)), ("-√8", 1, mp.pi), ("-1 + i√7", 3, PHI), ("-1 - i√7", 3, -PHI)):
        tot = m * phi
        wrapped = mp.fmod(tot + mp.pi, 2 * mp.pi)
        if wrapped < 0:
            wrapped += 2 * mp.pi
        wrapped -= mp.pi
        print("  固有値 %-9s (%d 個): 角度 m × φ = %s rad (%s°)   折り返すと %s°   周回数 %s"
              % (name, mult, fnum(tot, 14), mp.nstr(mp.degrees(tot), 12, min_fixed=-5, max_fixed=12, strip_zeros=True),
                 mp.nstr(mp.degrees(wrapped), 10, min_fixed=-5, max_fixed=6, strip_zeros=True),
                 mp.nstr(mp.floor(abs(tot) / (2 * mp.pi)), 30, min_fixed=0, max_fixed=30).rstrip("0").rstrip(".")))
    rows_summary.append((n, n_S, summ))
    sys.stdout.flush()

print("=" * 110)
print("まとめ: 各 n の行の指数部 (行 1〜8) と、配列の指数部")
print("=" * 110)
print("  n   S の個数 m   log2(m)   行の a = ln A (8 行共通か)        行の θ[度] 行1〜8                      m × ln√8        m × φ [度] (折り返さない)")
for n, n_S, summ in rows_summary:
    As = sorted(set(fnum(a, 12) for a, _ in summ))
    ths = " ".join(("-" if t is None else mp.nstr(mp.degrees(t), 6, min_fixed=-5, max_fixed=6, strip_zeros=True)) for _, t in summ)
    print("  %2d   %-10d   %-7s   %-32s   %-40s   %-14s  %s"
          % (n, n_S, mp.nstr(mp.log(n_S, 2), 6), ", ".join(As), ths, fnum(n_S * LN_R8, 12),
             mp.nstr(mp.degrees(n_S * PHI), 14, min_fixed=-5, max_fixed=14, strip_zeros=True)))
