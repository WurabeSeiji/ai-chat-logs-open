#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v7 (2026-11-27): array_product_cumulative_general_v6.py に、読み方のフラグを一つ足したもの (木原さんの指示)。
  読み方 17 (指定なし): v6 と同じ。各行の N 成分を、実部 1 個 + 虚部 N-1 個と読む。この場合に通る行は v6 のまま。
  読み方 13          : 各行の 8 成分を、成分 1〜4 と成分 5〜8 の二つの四元数と読む (N = 8 のときだけ)。
                       それぞれ q = exp(a + nθ)、実部 1 個 + 虚部 3 個。実部の指定が k なら、成分 k と成分 k+4 を実部と読む (k は 1〜4)。
                       図は v6 と同じ並び: 面 1 が a = ln A、面 2 が θ (どちらも 8 行 × 2 個の 16 本)、
                       面 3〜5 が行 2 の一つ目の四元数の向き θn (3 成分) を成分どうしで描いたもの、面 6〜8 が二つ目の四元数。
  初期値、演算、ダンプの [配列] の部分は v6 と一字も変えていない。変わるのは読出しと図だけ。
  使い方: python3 array_product_cumulative_general_v7.py N [回数] [実部の番号] [plus1: 0 か 1] [shift または 並べ替え] [読み方: 17 か 13]
          (省くと回数 20、実部は 1 番目、plus1 = 0、並べ替えなし (0)、読み方 17)
以下は v6・v5・v4 の説明のまま。

v6 (2026-11-27): array_product_cumulative_general_v5.py から、初期値の取り方だけを変えたもの。
  初期値: rowcol_zero_init_pm1_v1.py の make_init(N, plus1, order)。
    plus1 = 0: 縦横の和が 0 の初期値 (v5 と同じ)。 plus1 = 1: 一行目・一列目が全部 +1 の前回の初期値 (v4 と同じ)。
    order   : 行の並べ替え。初期化のあとに施す。前回 (v4) と同じ書き方・同じ意味。
  演算、保存の形、読出し、図の描き方は v4・v5 と一字も変えていない。
  使い方: python3 array_product_cumulative_general_v6.py N [回数] [実部と読む成分の番号] [plus1: 0 か 1] [shift または 並べ替え]
          (省くと回数 20、実部は 1 番目、plus1 = 0、並べ替えなし)
以下は v5 と v4 の説明のまま (初期値と使い方の記述は v6 には当たらない)。

v5 (2026-11-27): array_product_cumulative_general_v4.py から、初期値の取り方だけを変えたもの。
  初期値: rowcol_zero_init_pm1_v0.py の rowcol_zero_init(N) が返す N×N (±1、縦横の和が全部 0、全部 +1 の行や列が無い)。
  行の並べ替えの指定は無い (初期化プログラムが返す配列をそのまま使う)。
  演算、保存の形、読出し、図の描き方は v4 と一字も変えていない。変えたのは、初期値を取る行、引数、図のファイル名だけ。
  使い方: python3 array_product_cumulative_general_v5.py N [回数] [実部と読む成分の番号]   (省くと回数 20、実部は 1 番目)
以下は v4 の説明のまま (初期値と使い方の記述は v4 のもので、v5 には当たらない)。

N を指定して、N×N の配列 H (成分はすべて ±1、行どうし・列どうしが直交) を初期値に、
    O_n = O_1 × O_2 × ... × O_{n-1}      (ふつうの配列の積。常に左から順に掛ける)
を繰り返す汎用プログラム。v4。8x8 の S に対する array_product_cumulative_S8_v4.py と同じ演算・同じ保存の形。
v0 からの変更: 初期値の行をずらす量 shift (v1)、または行の並べ替えそのもの (v2、例 "1,6,4,3,7,2,5,8") を受け取る
(balanced_init_pm1_v2.py の square_array_ordered)。
v3 の変更: 図の面 2 (各行の θ) で、動き方の違う行が 5 本以上あると 4 本しか描かれなかった不具合を直した。
  全部の行を必ず描く。動き方が 4 種類以下なら色を分けて凡例を付け、5 種類以上なら全部を同じ色の細い線で描き、
  最後の回の値を一覧で書く。あわせて、一回掛けたあと (n = 3) の各行の θ を標準出力の末尾に出す。
v4 の変更: 実部と読む成分を指定できるようにした (木原さんの指示、2026-11-27)。演算は変えない。読出しと図だけが変わる。
  図の面 2 の直し (同日): 「1 回目」の札が行 1 の値だけを全行の値のように書いていた誤りと、線の上の札の重なりを直した。
  行の組ごとに印を変え、1 回目と最後の回の値は凡例に書く。θ の範囲 (0〜180 度) の上に凡例の場所を取った。
  1 なら第一成分を実部と読む (これまでと同じ)。k なら k 番目の成分を実部と読み、残りの N-1 成分を虚部と読む。
仕様は「思考実験.md」§32。

初期値: balanced_init_pm1_v0.py の square_array(N) が返す N×N の配列 (一列目を取り除く前のもの)。
        掛け算には正方の配列が要るので、N 個 × N-1 成分の表ではなく、こちらを使う。
        この N×N だけでは重心の条件は満たされない (§31.2)。モーメントと等方性は満たす。
O_1 = H。O_2 = O_1 (式のとおり)。3 回目以降は毎回「前回の配列の二乗」になる。
割り算・正規化・乱数は入れない。成分は整数のまま厳密に計算する。

保存 (標準出力。results/ に保存して使う):
  n ごとに [配列] N×N の全成分 (120 桁を超えたら「仮数 e 指数」で表示。計算は厳密)、
           [行の指数部] 各行を o = exp(a + nθ) と読んだ a = ln A、θ (主値 0〜π)、θn の N-1 成分。
           実部と読む成分は引数で指定する (読出しの選択)。読出しは演算に戻さない。
図 (results/array_product_cumulative_general_v4_N{N}_{shift または order}.png):
  1. 振幅を ln で読む: 各行の a = ln A を回数に対して (縦軸は対数の目盛り)。
  2. 各行の θ [度] を回数に対して。
  3〜. 向き n を楕円軌道として読む: 一つの行 (行 2) の θn の成分を、回数に対してではなく、成分どうしで描く。
       全回を通じて同じ値を取り続ける成分は一つの組にまとめ、組の代表どうしを描く (組が多いときは最初の三組)。
       点に回数 n を添える。

使い方: python3 array_product_cumulative_general_v4.py N [回数] [shift または 並べ替え] [実部と読む成分の番号]
        (省くと回数 20、並べ替えなし、実部は 1 番目)
"""
import os
import sys

import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from rowcol_zero_init_pm1_v1 import make_init, tag

sys.set_int_max_str_digits(0)
mp.mp.dps = 60
FULL_DIGITS = 120
N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
N_MAX = int(sys.argv[2]) if len(sys.argv) > 2 else 20
PLUS1 = (len(sys.argv) > 4 and sys.argv[4] == "1")
ORDER = sys.argv[5] if len(sys.argv) > 5 else None
H, HOW = make_init(N, PLUS1, ORDER)
RE = (int(sys.argv[3]) - 1) if len(sys.argv) > 3 else 0
assert 0 <= RE < N, "実部と読む成分の番号は 1 から N まで"
IMAG = [k for k in range(N) if k != RE]          # 虚部と読む成分 (番号の小さい順)
TAG = tag(N, PLUS1, ORDER)
READ = sys.argv[6] if len(sys.argv) > 6 else "17"
assert READ in ("17", "13"), "読み方は 17 か 13"
if READ == "13":
    assert N == 8, "読み方 13 (四元数二つ) は N = 8 のときだけ"
    assert 0 <= RE < 4, "読み方 13 では、実部の番号は 1 から 4 まで"
    QBLK = [[0, 1, 2, 3], [4, 5, 6, 7]]                 # 一つ目の四元数 (成分 1〜4)、二つ目の四元数 (成分 5〜8)
    QNAME = ["前", "後"]
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
    im = [o[k] for k in IMAG]
    vn = mp.sqrt(mp.fsum(x * x for x in im))
    if A == 0:
        return None, None, None
    a, th = mp.log(A), mp.atan2(vn, o[RE])
    if vn == 0:
        return a, th, ([mp.mpf(0)] * (N - 1) if th == 0 else None)
    return a, th, [th * x / vn for x in im]


def exponent_q(row, blk):
    """読み方 13: 4 成分 blk を四元数一つと読む。実部は blk[RE]、虚部は残りの 3 成分 (番号順)。式は exponent と同じ。"""
    o = [mp.mpf(row[k]) for k in blk]
    A = mp.sqrt(mp.fsum(x * x for x in o))
    im = [o[k] for k in range(4) if k != RE]
    vn = mp.sqrt(mp.fsum(x * x for x in im))
    if A == 0:
        return None, None, None
    a, th = mp.log(A), mp.atan2(vn, o[RE])
    if vn == 0:
        return a, th, ([mp.mpf(0)] * 3 if th == 0 else None)
    return a, th, [th * x / vn for x in im]


if READ == "13":
    print("N = %d。初期値 H の作り方: %s。回数 %d。読み方 13 (四元数二つ: 成分 1〜4 と成分 5〜8)。実部と読む成分: %d 番目と %d 番目" % (N, HOW, N_MAX, RE + 1, RE + 5))
else:
    print("N = %d。初期値 H の作り方: %s。回数 %d。実部と読む成分: %d 番目" % (N, HOW, N_MAX, RE + 1))
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
    if READ == "17":
        print("[行の指数部] o = exp(a + nθ)、%d 番目の成分を実部と読む。a = ln A、θ は主値 (0〜π)。θn は残りの成分を番号順に" % (RE + 1))
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
    else:
        print("[行の指数部 13] 各行を四元数二つと読む。前 = 成分 1〜4、後 = 成分 5〜8。q = exp(a + nθ)、実部は成分 %d (前) と成分 %d (後)。a = ln A、θ は主値 (0〜π)。θn は残りの 3 成分を番号順に" % (RE + 1, RE + 5))
        ex = []
        for i, r in enumerate(O):
            pair = []
            for q in range(2):
                a, th, tn = exponent_q(r, QBLK[q])
                pair.append((a, th, tn))
                if a is None:
                    print("  行%d %s: 全成分 0" % (i + 1, QNAME[q]))
                    continue
                tns = "向き不定 (虚部が 0 で θ = π)" if tn is None else "(" + ", ".join(fnum(x, 10) for x in tn) + ")"
                print("  行%d %s: a = %s   θ = %s rad (%s°)   θn = %s"
                      % (i + 1, QNAME[q], fnum(a), fnum(th, 14), mp.nstr(mp.degrees(th), 10, min_fixed=-5, max_fixed=6, strip_zeros=True), tns))
            ex.append(pair)
    rec.append((n, m, ex))
    sys.stdout.flush()

# ---------------------------------------------------------------- 図化
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
C1, C2 = "#2a78d6", "#eb6834"
plt.rcParams["font.family"] = ["Hiragino Sans", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
ns = [r[0] for r in rec]
ROW = 1                                                        # 楕円軌道として読む行 (行 2)

if READ == "13":
    # ------------------------------------------------------------ 読み方 13 の図 (並びは読み方 17 と同じ)
    NAN = float("nan")
    npan = 8
    nrow = 4
    fig, axes = plt.subplots(nrow, 2, figsize=(15, 7.3 * nrow), facecolor=SURF, squeeze=False)
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

    ents = [(i, q) for q in range(2) for i in range(N)]            # 16 本: (行, 前か後か)
    ename = lambda e: "%s 行%d" % (QNAME[e[1]], e[0] + 1)

    def group_label(es):
        parts = []
        for q in range(2):
            rows_ = [i + 1 for (i, qq) in es if qq == q]
            if rows_:
                parts.append(QNAME[q] + " 行" + (("%d〜%d" % (rows_[0], rows_[-1])) if (len(rows_) > 1 and rows_ == list(range(rows_[0], rows_[-1] + 1))) else ",".join(map(str, rows_))))
        return " / ".join(parts)

    # 1. 振幅を ln で読む (全成分 0 の四元数は a が定義できないので、その回は点を描かない)
    ax = flat[0]
    style(ax, "1. 振幅を ln で読む:四元数ごとの a = ln A(8 行 × 前・後)", "回数 n", "a = ln A(縦軸は対数の目盛り)")
    a_ser = {e: [(float(rec[t][2][e[0]][e[1]][0]) if rec[t][2][e[0]][e[1]][0] is not None else None) for t in range(len(rec))] for e in ents}
    akinds = {}
    for e in ents:
        akinds.setdefault(tuple(None if v is None else round(v, 9) for v in a_ser[e]), []).append(e)
    for key, es in sorted(akinds.items(), key=lambda kv: ents.index(kv[1][0])):
        ax.semilogy(ns, [NAN if (v is None or v <= 0) else v for v in key], color=C1, linewidth=1.3, alpha=0.7, marker="o", markersize=4)
    n_zero = sum(1 for e in ents for v in a_ser[e] if v is None)
    ax.set_xticks(ns)
    ax.text(0.03, 0.93, "a の動き方は %d 種類(16 本中)。全成分 0 で a が定義できない点は %d 個(描いていない)" % (len(akinds), n_zero),
            transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=INK2)

    # 2. 四元数ごとの θ
    ax = flat[1]
    style(ax, "2. 四元数ごとの θ(主値。成分 %d と成分 %d を実部と読む)" % (RE + 1, RE + 5), "回数 n", "θ [度]")
    th_ser = {e: [(float(mp.degrees(rec[t][2][e[0]][e[1]][1])) if rec[t][2][e[0]][e[1]][1] is not None else None) for t in range(len(rec))] for e in ents}
    kinds = {}
    for e in ents:
        kinds.setdefault(tuple(None if v is None else round(v, 9) for v in th_ser[e]), []).append(e)
    ordered_kinds = sorted(kinds.items(), key=lambda kv: ents.index(kv[1][0]))
    fmt = lambda v: "無し" if v is None else "%.3f°" % v
    if len(ordered_kinds) <= 4:
        for (key, es), col, mk in zip(ordered_kinds, (C2, C1, "#1baf7a", "#eda100"), ("o", "s", "^", "D")):
            ax.plot(ns, [NAN if v is None else v for v in key], color=col, linewidth=2, marker=mk, markersize=6,
                    label="%s:1 回目 %s → %d 回目 %s" % (group_label(es), fmt(key[0]), ns[-1], fmt(key[-1])))
        ax.legend(loc="upper left", frameon=False, fontsize=10, labelcolor=INK2)
    else:
        for key, es in ordered_kinds:
            ax.plot(ns, [NAN if v is None else v for v in key], color=C1, linewidth=1.1, alpha=0.6, marker="o", markersize=3.5)
        items = ["%s %s→%s" % (group_label(es), fmt(key[0]).replace(".000", ""), fmt(key[-1])) for key, es in ordered_kinds]
        txt = "動き方が %d 種類(全 16 本を描画)。1 回目 → %d 回目の θ:\n" % (len(ordered_kinds), ns[-1])
        txt += "\n".join("  ".join(items[q:q + 3]) for q in range(0, min(len(items), 18), 3))
        ax.text(0.02, 0.99, txt, transform=ax.transAxes, ha="left", va="top", fontsize=8.5, color=INK2)
    ax.set_xticks(ns)
    ax.set_ylim(-8, 262)
    ax.set_yticks([0, 45, 90, 135, 180])
    ax.axhline(180, color=AXIS, linewidth=1)

    # 3〜8. 行 2 の向き θn (3 成分) を成分どうしで描く。面 3〜5 が前、面 6〜8 が後
    start = 2 if len(rec) > 2 else 0
    summary = []
    for q in range(2):
        imag = [k for k in QBLK[q] if k != QBLK[q][RE]]
        valid = [t for t in range(start, len(rec)) if rec[t][2][ROW][q][2] is not None]
        ser = [[float(rec[t][2][ROW][q][2][c]) for t in valid] for c in range(3)]
        vns = [rec[t][0] for t in valid]
        n_distinct = len({tuple(round(s[i], 12) for s in ser) for i in range(len(vns))})
        summary.append((q, len(vns), n_distinct, len(rec) - start))
        lim = max(1.8, 1.15 * max([abs(v) for s in ser for v in s] or [0]))
        for pi, (i, j) in enumerate([(0, 1), (0, 2), (1, 2)]):
            num = 3 + 3 * q + pi
            ax = flat[num - 1]
            style(ax, "%d. 向き n を軌道として読む(行 2 の%s:成分 %d〜%d)" % (num, QNAME[q], QBLK[q][0] + 1, QBLK[q][-1] + 1),
                  "θn の 成分 %d [rad]" % (imag[i] + 1), "θn の 成分 %d [rad]" % (imag[j] + 1))
            if not vns:
                ax.text(0.5, 0.5, "描ける点が無い\n(全成分 0、または向きが不定)", transform=ax.transAxes, ha="center", va="center", fontsize=11, color=INK2)
                continue
            ax.plot(ser[i], ser[j], linewidth=0, marker="o", markersize=8, color=C1, markeredgecolor=SURF, markeredgewidth=1.2)
            seen = {}
            for n_, px, py in zip(vns, ser[i], ser[j]):
                seen.setdefault((round(px, 9), round(py, 9)), []).append(n_)
            for (px, py), lst in seen.items():
                lab = str(lst[0]) if len(lst) == 1 else ("n=%d〜%d" % (lst[0], lst[-1]) if lst == list(range(lst[0], lst[-1] + 1)) else ",".join(map(str, lst)))
                ax.annotate(lab, (px, py), textcoords="offset points", xytext=(8, 7), fontsize=9.5, color=INK)
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_aspect("equal")
            ax.axhline(0, color=AXIS, linewidth=1)
            ax.axvline(0, color=AXIS, linewidth=1)
            if pi == 0:
                ax.text(0.03, 0.97, "n=%d〜%d の %d 回分。違う位置は %d 個" % (vns[0], vns[-1], len(vns), n_distinct),
                        transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=INK2)
    fig.suptitle("N = %d(初期値:%s)の累積の積:四元数二つ(成分 1〜4、成分 5〜8)と読む。実部は成分 %d と成分 %d" % (N, HOW, RE + 1, RE + 5),
                 fontsize=14.5, color=INK, x=0.07, ha="left", y=1 - 0.3 / (7.3 * nrow))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "array_product_cumulative_general_v7_N%d_%s_re%d_q13.png" % (N, TAG, RE + 1))
    fig.savefig(out, dpi=150, facecolor=SURF)
    print("=" * 110)
    for q, nv, nd, tot in summary:
        print("行 2 の%s (成分 %d〜%d): n>=3 の %d 回のうち向きが描けるのは %d 回、違う位置は %d 個" % (QNAME[q], QBLK[q][0] + 1, QBLK[q][-1] + 1, tot, nv, nd))
    print("図を保存: results/" + os.path.basename(out))
    if len(rec) >= 3:
        d = lambda v: "無し" if v is None else mp.nstr(mp.degrees(v), 8, min_fixed=-5, max_fixed=6, strip_zeros=True)
        for q in range(2):
            print("一回掛けたあと (n = 3) の θ [度]、%s: " % QNAME[q] + "  ".join("行%d %s" % (i + 1, d(rec[2][2][i][q][1])) for i in range(N)))
        print("面 1 の a の動き方の種類: %d、面 2 の θ の動き方の種類: %d (16 本中)" % (len(akinds), len(ordered_kinds)))
    sys.exit(0)

# 全回を通じて同じ値を取り続ける成分を一つの組にまとめる (整数の配列で厳密に比べる)
groups = []
for k in IMAG:
    for g in groups:
        if all(history[t][ROW][k] == history[t][ROW][g[0]] for t in range(len(history))):
            g.append(k)
            break
    else:
        groups.append([k])
start = 2 if len(rec) > 2 else 0                               # n >= 3 (H を偶数個掛けた回)
valid = [t for t in range(start, len(rec)) if rec[t][2][ROW][2] is not None]
series = [[float(rec[t][2][ROW][2][IMAG.index(g[0])]) for t in valid] for g in groups]
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
style(ax, "2. 各行の θ(主値。%d 番目の成分を実部と読む)" % (RE + 1), "回数 n", "θ [度]")
th_rows = [[float(mp.degrees(rec[t][2][i][1])) for t in range(len(rec))] for i in range(N)]
kinds = {}
for i in range(N):
    kinds.setdefault(tuple(round(v, 9) for v in th_rows[i]), []).append(i + 1)


def rows_label(rows_):
    if len(rows_) == 1:
        return "行 %d" % rows_[0]
    if rows_ == list(range(rows_[0], rows_[-1] + 1)):
        return "行 %d〜%d" % (rows_[0], rows_[-1])
    return "行 " + ",".join(map(str, rows_))


ordered_kinds = sorted(kinds.items(), key=lambda kv: kv[1][0])
if len(ordered_kinds) <= 4:
    # 行の組ごとに色と印を変える。最初の回と最後の回の値は凡例に書く (線の上に札を重ねない)
    for (key, rows_), col, mk in zip(ordered_kinds, (C2, C1, "#1baf7a", "#eda100"), ("o", "s", "^", "D")):
        lab = "%s:1 回目 %.3f° → %d 回目 %.3f°" % (rows_label(rows_), key[0], ns[-1], key[-1])
        ax.plot(ns, list(key), color=col, linewidth=2, marker=mk, markersize=6, label=lab)
    ax.legend(loc="upper left", frameon=False, fontsize=10, labelcolor=INK2)
else:
    for key, rows_ in ordered_kinds:                       # 全部の行を描く (色は一つ。行の区別は下の一覧で)
        ax.plot(ns, list(key), color=C1, linewidth=1.1, alpha=0.6, marker="o", markersize=3.5)
    txt = "動き方の違う行が %d 種類(全 %d 行を描画)。1 回目 → %d 回目の θ:\n" % (len(ordered_kinds), N, ns[-1])
    items = ["%s %.1f°→%.1f°" % (rows_label(rows_), key[0], key[-1]) for key, rows_ in ordered_kinds]
    per = 3
    txt += "\n".join("  ".join(items[q:q + per]) for q in range(0, min(len(items), 15), per)) + ("\n…" if len(items) > 15 else "")
    ax.text(0.02, 0.99, txt, transform=ax.transAxes, ha="left", va="top", fontsize=9, color=INK2)
ax.set_xticks(ns)
ax.set_ylim(-8, 262)
ax.set_yticks([0, 45, 90, 135, 180])
ax.axhline(180, color=AXIS, linewidth=1)

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
fig.suptitle("N = %d(初期値:%s)の累積の積 Oₙ = O₁×O₂×…×Oₙ₋₁:振幅は ln、向きは軌道として読む。実部は %d 番目の成分" % (N, HOW, RE + 1),
             fontsize=14.5, color=INK, x=0.07, ha="left", y=1 - 0.3 / (7.3 * nrow))
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "array_product_cumulative_general_v7_N%d_%s_re%d.png" % (N, TAG, RE + 1))
fig.savefig(out, dpi=150, facecolor=SURF)
print("=" * 110)
print("行 2 の θn の成分の組 (全回を通じて同じ値): " + " / ".join(gname(g) for g in groups))
print("n=%d〜%d の %d 回分で、向きの違う位置は %d 個" % (vns[0], vns[-1], len(vns), n_distinct))
print("図を保存: results/" + os.path.basename(out))
if len(rec) >= 3:
    print("一回掛けたあと (n = 3) の各行の θ [度]: " + "  ".join("行%d %s" % (i + 1, mp.nstr(mp.degrees(rec[2][2][i][1]), 8, min_fixed=-5, max_fixed=6, strip_zeros=True)) for i in range(N)))
    print("面 2 の動き方の種類: %d (行の数 %d)" % (len(ordered_kinds), N))
