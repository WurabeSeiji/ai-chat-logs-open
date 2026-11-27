#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
8x8 の配列 S (成分はすべて ±1、行どうし・列どうしが直交) を初期値にして、
    O_n = O_1 × O_2 × O_3 × ... × O_{n-1}      (ふつうの配列の積。常に左から順に掛ける)
を n = 20 まで繰り返す。v3。演算は v0〜v2 と同じ。v2 に、二枚目の図 (全成分の ln) を足した (木原さんの指示、2026-11-27)。
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

図 (results/array_product_cumulative_S8_v2_exponent.png、5 面):
  1. 配列の指数部の大きさ (対数の軸): 実の成分 m ln√8 と、角度 m φ [rad]。折り返さない値。
  2. 指数部の平面: 横が実の成分、縦が角度。両方とも対数の軸。点に回数 n を添える。
  3. 配列の角度を折り返した値 [度]。各点に値を併記。
  4. 行の θ [度] (行を八元数として読んだ主値)。行 1 と行 2〜8。値を併記。
  5. 行 2 の虚の 7 成分 θn。7 成分は三種類の値だけを取るので、その三つを描く。

二枚目の図 (results/array_product_cumulative_S8_v3_ln_components.png):
  各回の 8x8 の全成分 (64 個) について、ln|成分| をそのまま描く。加工はそれだけ。
  符号は色と位置で分ける (正は青で左寄り、負は橙で右寄り)。成分が 0 のものは ln が取れないので描かず、個数を図に書く。
  上の面は縦軸をそのままの目盛り、下の面は同じ値を縦軸だけ対数の目盛りにしたもの (上の面では後半しか見えないため)。

使い方: python3 array_product_cumulative_S8_v3.py > results/array_product_cumulative_S8_v3.txt
        (図は results/ に自動で保存される。図だけ作り直すことはできない。計算と図化は一体)
"""
import sys

import os

import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
plot_rows = []                    # 図化用: (n, m, 行1のθ, 行2のθ, 行2のθn)
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
    plot_rows.append((n, n_S, summ[0][1], summ[1][1], exponent(O[1])[2]))
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


# ---------------------------------------------------------------- 図化
def make_figure(path):
    INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
    C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"          # 検証済みの三色 (青・橙・緑青)。この順で固定
    plt.rcParams["font.family"] = ["Hiragino Sans", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    ns = [r[0] for r in plot_rows]
    ms = [r[1] for r in plot_rows]
    a_arr = [float(m * LN_R8) for m in ms]
    ang_arr = [float(m * PHI) for m in ms]

    def wrap_deg(m):
        w = mp.fmod(m * PHI + mp.pi, 2 * mp.pi)
        if w < 0:
            w += 2 * mp.pi
        return float(mp.degrees(w - mp.pi))
    ang_wrapped = [wrap_deg(m) for m in ms]
    th1 = [float(mp.degrees(r[2])) for r in plot_rows]
    th2 = [float(mp.degrees(r[3])) for r in plot_rows]
    # 行 2 の虚の 7 成分の三つの組 (S を 2 個掛けた配列の行 2 = (0,-4,0,4,0,4,4,0) の位置で分ける)
    g_self, g_a, g_b = [0], [2, 4, 5], [1, 3, 6]
    tn = [r[4] for r in plot_rows]
    v_self = [float(t[g_self[0]]) for t in tn]
    v_a = [float(t[g_a[0]]) for t in tn]
    v_b = [float(t[g_b[0]]) for t in tn]
    for t in tn:                                           # 三種類の値だけ、という前提を確かめる
        assert all(abs(t[i] - t[g_a[0]]) < mp.mpf(10) ** -40 for i in g_a)
        assert all(abs(t[i] - t[g_b[0]]) < mp.mpf(10) ** -40 for i in g_b)

    fig = plt.figure(figsize=(15, 16.5), facecolor=SURF)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.0], hspace=0.42, wspace=0.2,
                          left=0.07, right=0.97, top=0.925, bottom=0.05)

    def style(ax, title, xlabel, ylabel):
        ax.set_facecolor(SURF)
        ax.set_title(title, fontsize=13, color=INK, loc="left", pad=10)
        ax.set_xlabel(xlabel, fontsize=11, color=INK2)
        ax.set_ylabel(ylabel, fontsize=11, color=INK2)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(AXIS)
        ax.tick_params(colors=MUTED, labelsize=10)

    # 1. 指数部の大きさ
    ax = fig.add_subplot(gs[0, 0])
    style(ax, "1. 配列の指数部の大きさ(折り返さない値)", "回数 n", "指数部の成分(対数の軸)")
    ax.semilogy(ns, ang_arr, color=C2, linewidth=2, marker="o", markersize=6, label="角度 m×φ [rad]")
    ax.semilogy(ns, a_arr, color=C1, linewidth=2, marker="o", markersize=6, label="実の成分 m×ln√8")
    ax.set_xticks(ns)
    ax.annotate("角度 m×φ [rad]", (ns[11], ang_arr[11]), textcoords="offset points", xytext=(-12, 16), ha="right", fontsize=10.5, color=INK)
    ax.annotate("実の成分 m×ln√8", (ns[13], a_arr[13]), textcoords="offset points", xytext=(14, -20), ha="left", fontsize=10.5, color=INK)
    ax.legend(loc="upper left", frameon=False, fontsize=10.5, labelcolor=INK2)
    ax.text(0.98, 0.05, "二本は平行。どちらも毎回ちょうど 2 倍\n(n=1 と 2 は同じ値。式のとおり O₂ = O₁)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10, color=INK2)

    # 2. 指数部の平面
    ax = fig.add_subplot(gs[0, 1])
    style(ax, "2. 指数部の平面(横:実の成分、縦:角度。両方とも対数の軸)", "実の成分 m×ln√8", "角度 m×φ [rad]")
    ax.loglog(a_arr, ang_arr, color=C1, linewidth=2, marker="o", markersize=6)
    for n_, x_, y_ in zip(ns, a_arr, ang_arr):
        if n_ == 1:
            continue
        if n_ == 2 or (n_ % 2 == 1 and n_ > 3) or n_ == 20:      # n=3 は n=1, 2 の札と重なるので省く
            ax.annotate("n=1, 2" if n_ == 2 else "n=%d" % n_, (x_, y_), textcoords="offset points",
                        xytext=(10, -2) if n_ == 2 else (8, -12), fontsize=9.5, color=INK2)
    ratio = float(PHI / LN_R8)
    ax.text(0.04, 0.93, "点は一本の直線の上に並ぶ。\n角度 ÷ 実の成分 = φ ÷ ln√8 = %.6f(どの回も同じ)" % ratio,
            transform=ax.transAxes, ha="left", va="top", fontsize=10, color=INK2)

    # 3. 配列の角度を折り返した値
    ax = fig.add_subplot(gs[1, 0])
    style(ax, "3. 配列の角度 m×φ を ±180° に折り返した値", "回数 n", "角度 [度]")
    ax.axhline(0, color=AXIS, linewidth=1)
    ax.plot(ns, ang_wrapped, color=C2, linewidth=1.2, alpha=0.55)
    ax.plot(ns, ang_wrapped, color=C2, linewidth=0, marker="o", markersize=7)
    for n_, y_ in zip(ns, ang_wrapped):
        if n_ == 2:
            continue                                   # n=1 と同じ値なので、n=1 の札にまとめる
        up = y_ >= 0
        ax.annotate(("%.1f°(n=1, 2)" if n_ == 1 else "%.1f°") % y_, (n_, y_), textcoords="offset points",
                    xytext=(-4 if n_ == 1 else 0, 9 if up else -15), ha="left" if n_ == 1 else "center", fontsize=8.8, color=INK)
    ax.set_xticks(ns)
    ax.set_ylim(-250, 215)
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.text(0.98, 0.02, "φ = arccos(−1/√8) ≈ 110.70°。折り返すと、同じ値に戻らない", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=10, color=INK2)

    # 4. 行の θ
    ax = fig.add_subplot(gs[1, 1])
    style(ax, "4. 行を八元数として読んだ θ(主値。第一成分を実部と読む)", "回数 n", "θ [度]")
    ax.plot(ns, th2, color=C1, linewidth=2, marker="o", markersize=6, label="行 2〜8")
    ax.plot(ns, th1, color=C2, linewidth=2, marker="o", markersize=6, label="行 1")
    ax.set_xticks(ns)
    ax.set_ylim(-12, 108)
    ax.set_yticks([0, 30, 60, 90])
    ax.annotate("%.3f°(8 行とも)" % th1[0], (ns[0], th1[0]), textcoords="offset points", xytext=(10, 6), fontsize=10, color=INK)
    ax.annotate("行 2〜8:90.000°", (ns[-1], th2[-1]), textcoords="offset points", xytext=(-4, 8), ha="right", fontsize=10.5, color=INK)
    ax.annotate("行 1:0.000°", (ns[-1], th1[-1]), textcoords="offset points", xytext=(-4, 8), ha="right", fontsize=10.5, color=INK)
    ax.legend(loc="center right", frameon=False, fontsize=10.5, labelcolor=INK2)
    ax.text(0.5, 0.5, "3 回目で 0° と 90° に分かれ、以後動かない", transform=ax.transAxes, ha="center", va="center", fontsize=10, color=INK2)

    # 5. 行 2 の虚の成分
    ax = fig.add_subplot(gs[2, :])
    style(ax, "5. 行 2 の指数部の虚の 7 成分 θn(7 成分は三種類の値だけを取る)", "回数 n", "θn の成分 [rad]")
    ax.axhline(0, color=AXIS, linewidth=1)
    ax.plot(ns, v_self, color=C1, linewidth=2, marker="o", markersize=6, label="e₁ の成分(1 個)")
    ax.plot(ns, v_a, color=C2, linewidth=2, marker="s", markersize=6, label="e₃・e₅・e₆ の成分(3 個、同じ値)")
    ax.plot(ns, v_b, color=C3, linewidth=2, marker="^", markersize=7, label="e₂・e₄・e₇ の成分(3 個、同じ値)")
    ax.set_xticks(ns)
    ax.set_ylim(-1.75, 1.75)
    ax.legend(loc="upper center", ncol=3, frameon=False, fontsize=10.5, labelcolor=INK2)
    ax.text(0.99, 0.04, "θ は 3 回目以降ずっと 90°(π/2 rad)。変わるのは向き n だけ。成分の二乗の和は常に θ²",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10, color=INK2)

    fig.suptitle("8×8 の配列 S の累積の積 Oₙ = O₁×O₂×…×Oₙ₋₁:指数部で読む(n = 1〜20、m = 掛けた S の個数)",
                 fontsize=15, color=INK, x=0.07, ha="left", y=0.975)
    fig.savefig(path, dpi=150, facecolor=SURF)
    plt.close(fig)


OUT_PNG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "array_product_cumulative_S8_v3_exponent.png")
make_figure(OUT_PNG)
print("図を保存: results/" + os.path.basename(OUT_PNG))


# ---------------------------------------------------------------- 二枚目の図: 全成分の ln
def make_figure_ln(path):
    INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
    C1, C2 = "#2a78d6", "#eb6834"
    plt.rcParams["font.family"] = ["Hiragino Sans", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    pos_x, pos_y, neg_x, neg_y, zeros = [], [], [], [], []
    for n, O in enumerate(history, 1):
        z = 0
        for r in O:
            for v in r:
                if v == 0:
                    z += 1
                elif v > 0:
                    pos_x.append(n - 0.13)
                    pos_y.append(float(mp.log(mp.mpf(v))))
                else:
                    neg_x.append(n + 0.13)
                    neg_y.append(float(mp.log(mp.mpf(-v))))
        zeros.append(z)
    ns = list(range(1, len(history) + 1))
    fig, axes = plt.subplots(2, 1, figsize=(15, 13), facecolor=SURF, sharex=True)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.9, bottom=0.07, hspace=0.16)
    for ax, logy in zip(axes, (False, True)):
        ax.set_facecolor(SURF)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(AXIS)
        ax.tick_params(colors=MUTED, labelsize=10)
        ax.scatter(pos_x, pos_y, s=46, color=C1, alpha=0.55, linewidths=0, label="正の成分:ln(成分)")
        ax.scatter(neg_x, neg_y, s=46, color=C2, alpha=0.55, linewidths=0, marker="s", label="負の成分:ln(−成分)")
        ax.set_xticks(ns)
        ax.set_ylabel("ln|成分|" + ("(縦軸は対数の目盛り)" if logy else ""), fontsize=11, color=INK2)
        if logy:
            ax.set_yscale("log")
            ax.set_title("同じ値を、縦軸だけ対数の目盛りにしたもの", fontsize=13, color=INK, loc="left", pad=8)
            ax.set_xlabel("回数 n", fontsize=11, color=INK2)
            ax.text(0.01, 0.95, "n=1, 2 は全成分が ±1 で ln = 0 なので、この目盛りには描けない", transform=ax.transAxes,
                    ha="left", va="top", fontsize=10, color=INK2)
        else:
            ax.set_title("各回の 8×8 の全成分の ln|成分|(加工は ln だけ)", fontsize=13, color=INK, loc="left", pad=8)
            ax.legend(loc="upper left", frameon=False, fontsize=10.5, labelcolor=INK2)
            top = max(pos_y + neg_y)
            for n_, z_ in zip(ns, zeros):
                ax.text(n_, -0.035 * top, "0 が\n%d 個" % z_, ha="center", va="top", fontsize=8, color=MUTED)
            ax.set_ylim(-0.11 * top, 1.05 * top)
            ax.text(0.5, 0.6, "成分が 0 のものは ln が取れないので描いていない(個数は下に記入)", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color=INK2)
    fig.suptitle("8×8 の配列 S の累積の積 Oₙ = O₁×O₂×…×Oₙ₋₁:全成分(各回 64 個)の ln", fontsize=15, color=INK, x=0.07, ha="left", y=0.965)
    fig.savefig(path, dpi=150, facecolor=SURF)
    plt.close(fig)
    return zeros


OUT_PNG2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "array_product_cumulative_S8_v3_ln_components.png")
zs = make_figure_ln(OUT_PNG2)
print("図を保存: results/" + os.path.basename(OUT_PNG2) + "   各回の 0 の成分の個数: " + " ".join(map(str, zs)))
