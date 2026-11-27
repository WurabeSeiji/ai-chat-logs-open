#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N を指定して、N×N の配列 H (成分はすべて ±1、行どうし・列どうしが直交) を初期値に、
    O_n = O_1 × O_2 × ... × O_{n-1}      (ふつうの配列の積。常に左から順に掛ける)
を繰り返す汎用プログラム。v1。8x8 の S に対する array_product_cumulative_S8_v4.py と同じ演算・同じ保存の形。
v0 からの変更: 初期値の行をずらすパラメータ shift を受け取る (balanced_init_pm1_v1.py の square_array_shifted)。
仕様は「思考実験.md」§32。

初期値: balanced_init_pm1_v0.py の square_array(N) が返す N×N の配列 (一列目を取り除く前のもの)。
        掛け算には正方の配列が要るので、N 個 × N-1 成分の表ではなく、こちらを使う。
        この N×N だけでは重心の条件は満たされない (§31.2)。モーメントと等方性は満たす。
O_1 = H。O_2 = O_1 (式のとおり)。3 回目以降は毎回「前回の配列の二乗」になる。
割り算・正規化・乱数は入れない。成分は整数のまま厳密に計算する。

保存 (標準出力。results/ に保存して使う):
  n ごとに [配列] N×N の全成分 (120 桁を超えたら「仮数 e 指数」で表示。計算は厳密)、
           [行の指数部] 各行を o = exp(a + nθ) と読んだ a = ln A、θ (主値 0〜π)、θn の N-1 成分。
           今回は第一成分を実部と読む (読出しの選択)。読出しは演算に戻さない。
図 (results/array_product_cumulative_general_v1_N{N}_shift{shift}.png):
  1. 振幅を ln で読む: 各行の a = ln A を回数に対して (縦軸は対数の目盛り)。
  2. 各行の θ [度] を回数に対して。
  3〜. 向き n を楕円軌道として読む: 一つの行 (行 2) の θn の成分を、回数に対してではなく、成分どうしで描く。
       全回を通じて同じ値を取り続ける成分は一つの組にまとめ、組の代表どうしを描く (組が多いときは最初の三組)。
       点に回数 n を添える。

使い方: python3 array_product_cumulative_general_v1.py N [回数] [shift]     (省くと回数 20、shift 0)
"""
import os
import sys

import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from balanced_init_pm1_v1 import square_array_shifted

sys.set_int_max_str_digits(0)
mp.mp.dps = 60
FULL_DIGITS = 120
N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
N_MAX = int(sys.argv[2]) if len(sys.argv) > 2 else 20
SHIFT = int(sys.argv[3]) if len(sys.argv) > 3 else 0
H, HOW = square_array_shifted(N, SHIFT)
if H is None:
    sys.exit("N = %d の配列は、balanced_init_pm1_v0.py の作り方では作れない" % N)


def matmul(A, B):
    return [[sum(A[i][t] * B[t][j] for t in range(N)) for j in range(N)] for i in range(N)]


def ndigits(v):
    v = abs(v)
    return 1 if v == 0 else (len(str(v)) if v < 10 ** 400 else int(v.bit_length() * 0.30102999566398) + 1)


def show_int(v, full):
    return str(v) if (full or v == 0) else mp.nstr(mp.mpf(v), 17, min_fixed=0, max_fixed=0)


def fnum(x, d=17):
    return mp.nstr(x, d, min_fixed=-5, max_fixed=17, strip_zeros=True)


def exponent(row):
    o = [mp.mpf(v) for v in row]
    A = mp.sqrt(mp.fsum(x * x for x in o))
    vn = mp.sqrt(mp.fsum(x * x for x in o[1:]))
    if A == 0:
        return None, None, None
    a, th = mp.log(A), mp.atan2(vn, o[0])
    if vn == 0:
        return a, th, ([mp.mpf(0)] * (N - 1) if th == 0 else None)
    return a, th, [th * x / vn for x in o[1:]]


print("N = %d。初期値 H の作り方: %s。回数 %d" % (N, HOW, N_MAX))
history, rec = [], []
for n in range(1, N_MAX + 1):
    if n == 1:
        O = [r[:] for r in H]
    else:
        O = history[0]
        for k in range(1, n - 1):
            O = matmul(O, history[k])
    history.append(O)
    m = 1 if n <= 2 else 2 ** (n - 2)
    dg = max(ndigits(v) for r in O for v in r)
    full = dg <= FULL_DIGITS
    print("=" * 110)
    print("n = %d    (H を %d 個掛けたもの)    成分の最大桁数 %d" % (n, m, dg))
    print("=" * 110)
    print("[配列]" + ("" if full else "  (成分が %d 桁を超えるので、仮数 e 指数 で表示。計算は整数で厳密)" % FULL_DIGITS))
    for i, r in enumerate(O):
        print("  行%d: " % (i + 1) + "  ".join(show_int(v, full) for v in r))
    print("[行の指数部] o = exp(a + nθ)、第一成分を実部と読む。a = ln A、θ は主値 (0〜π)")
    ex = []
    for i, r in enumerate(O):
        a, th, tn = exponent(r)
        ex.append((a, th, tn))
        if a is None:
            print("  行%d: 全成分 0" % (i + 1))
            continue
        tns = "向き不定 (虚部が 0 で θ = π)" if tn is None else "(" + ", ".join(fnum(x, 10) for x in tn) + ")"
        print("  行%d: a = %s   θ = %s rad (%s°)   θn = %s"
              % (i + 1, fnum(a), fnum(th, 14), mp.nstr(mp.degrees(th), 10, min_fixed=-5, max_fixed=6, strip_zeros=True), tns))
    rec.append((n, m, ex))
    sys.stdout.flush()

# ---------------------------------------------------------------- 図化
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
C1, C2 = "#2a78d6", "#eb6834"
plt.rcParams["font.family"] = ["Hiragino Sans", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
ns = [r[0] for r in rec]
ROW = 1                                                        # 楕円軌道として読む行 (行 2)

# 全回を通じて同じ値を取り続ける成分を一つの組にまとめる (整数の配列で厳密に比べる)
groups = []
for k in range(1, N):
    for g in groups:
        if all(history[t][ROW][k] == history[t][ROW][g[0]] for t in range(len(history))):
            g.append(k)
            break
    else:
        groups.append([k])
start = 2 if len(rec) > 2 else 0                               # n >= 3 (H を偶数個掛けた回)
valid = [t for t in range(start, len(rec)) if rec[t][2][ROW][2] is not None]
series = [[float(rec[t][2][ROW][2][g[0] - 1]) for t in valid] for g in groups]
vns = [rec[t][0] for t in valid]
n_distinct = len({tuple(round(s[i], 12) for s in series) for i in range(len(vns))})
pairs = [(0, 1), (0, 2), (1, 2)] if len(groups) >= 3 else ([(0, 1)] if len(groups) == 2 else [])
npan = 2 + max(len(pairs), 1)
ncol = 2
nrow = (npan + 1) // 2
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 7.3 * nrow), facecolor=SURF, squeeze=False)
fig.subplots_adjust(left=0.07, right=0.97, top=1 - 1.15 / (7.3 * nrow), bottom=0.5 / (7.3 * nrow), hspace=0.3, wspace=0.2)
flat = [ax for row in axes for ax in row]


def style(ax, title, xl, yl):
    ax.set_facecolor(SURF)
    ax.set_title(title, fontsize=13, color=INK, loc="left", pad=10)
    ax.set_xlabel(xl, fontsize=11, color=INK2)
    ax.set_ylabel(yl, fontsize=11, color=INK2)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=10)


# 1. 振幅を ln で読む
ax = flat[0]
style(ax, "1. 振幅を ln で読む:各行の a = ln A", "回数 n", "a = ln A(縦軸は対数の目盛り)")
a_rows = [[float(rec[t][2][i][0]) for t in range(len(rec))] for i in range(N)]
same_a = all(abs(a_rows[i][t] - a_rows[0][t]) < 1e-12 * max(1.0, abs(a_rows[0][t])) for i in range(N) for t in range(len(rec)))
ax.semilogy(ns, a_rows[0], color=C1, linewidth=2, marker="o", markersize=6)
if not same_a:
    for i in range(1, N):
        ax.semilogy(ns, a_rows[i], color=C1, linewidth=1, alpha=0.5)
ax.set_xticks(ns)
ax.text(0.03, 0.93, ("%d 行とも同じ値。" % N if same_a else "行によって値が違う。") + "1 回目 a = %.6f、%d 回目 a = %.6g" % (a_rows[0][0], ns[-1], a_rows[0][-1]),
        transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=INK2)

# 2. 各行の θ
ax = flat[1]
style(ax, "2. 各行の θ(主値。第一成分を実部と読む)", "回数 n", "θ [度]")
th_rows = [[float(mp.degrees(rec[t][2][i][1])) for t in range(len(rec))] for i in range(N)]
kinds = {}
for i in range(N):
    kinds.setdefault(tuple(round(v, 9) for v in th_rows[i]), []).append(i + 1)
for (key, rows_), col in zip(sorted(kinds.items(), key=lambda kv: kv[1][0]), (C2, C1, "#1baf7a", "#eda100")):
    lab = "行 " + (str(rows_[0]) if len(rows_) == 1 else "%d〜%d" % (rows_[0], rows_[-1]) if rows_ == list(range(rows_[0], rows_[-1] + 1)) else ",".join(map(str, rows_)))
    ax.plot(ns, list(key), color=col, linewidth=2, marker="o", markersize=6, label=lab)
    ax.annotate("%s:%.3f°" % (lab, key[-1]), (ns[-1], key[-1]), textcoords="offset points", xytext=(-4, 8), ha="right", fontsize=10.5, color=INK)
ax.set_xticks(ns)
ax.set_ylim(-12, 192)
ax.set_yticks([0, 45, 90, 135, 180])
ax.legend(loc="upper left", bbox_to_anchor=(0.3, 0.98), frameon=False, fontsize=10.5, labelcolor=INK2)
ax.annotate("1 回目:%.3f°" % th_rows[0][0], (ns[0], th_rows[0][0]), textcoords="offset points", xytext=(8, 8), fontsize=10, color=INK)

# 3〜. 向き n を楕円軌道として読む
gname = lambda g: "成分 " + ",".join(str(k + 1) for k in g) + "(%d 個)" % len(g)
if not pairs:
    ax = flat[2]
    style(ax, "3. 向き n(行 2)", "", "")
    ax.text(0.5, 0.5, "行 2 の θn の成分は、全回を通じて一つの組しか無い\n(成分どうしで描く相手が無い)", transform=ax.transAxes, ha="center", va="center", fontsize=11, color=INK2)
for num, (i, j) in enumerate(pairs, 3):
    ax = flat[num - 1]
    style(ax, "%d. 向き n を軌道として読む(行 2)" % num, "θn の " + gname(groups[i]) + " [rad]", "θn の " + gname(groups[j]) + " [rad]")
    ax.plot(series[i], series[j], linewidth=0, marker="o", markersize=8, color=C1, markeredgecolor=SURF, markeredgewidth=1.2)
    seen = {}
    for n_, px, py in zip(vns, series[i], series[j]):
        seen.setdefault((round(px, 9), round(py, 9)), []).append(n_)
    for (px, py), lst in seen.items():
        lab = str(lst[0]) if len(lst) == 1 else ("n=%d〜%d" % (lst[0], lst[-1]) if lst == list(range(lst[0], lst[-1] + 1)) else ",".join(map(str, lst)))
        ax.annotate(lab, (px, py), textcoords="offset points", xytext=(8, 7), fontsize=9.5, color=INK)
    lim = max(1.8, 1.15 * max(abs(v) for s in series for v in s))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axhline(0, color=AXIS, linewidth=1)
    ax.axvline(0, color=AXIS, linewidth=1)
    if num == 3:
        ax.text(0.03, 0.97, "n=%d〜%d の %d 回分。違う位置は %d 個" % (vns[0], vns[-1], len(vns), n_distinct),
                transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=INK2)
for ax in flat[npan:]:
    ax.axis("off")
fig.suptitle("N = %d(初期値:%s)の累積の積 Oₙ = O₁×O₂×…×Oₙ₋₁:振幅は ln、向きは軌道として読む" % (N, HOW),
             fontsize=14.5, color=INK, x=0.07, ha="left", y=1 - 0.3 / (7.3 * nrow))
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "array_product_cumulative_general_v1_N%d_shift%d.png" % (N, SHIFT))
fig.savefig(out, dpi=150, facecolor=SURF)
print("=" * 110)
print("行 2 の θn の成分の組 (全回を通じて同じ値): " + " / ".join(gname(g) for g in groups))
print("n=%d〜%d の %d 回分で、向きの違う位置は %d 個" % (vns[0], vns[-1], len(vns), n_distinct))
print("図を保存: results/" + os.path.basename(out))
