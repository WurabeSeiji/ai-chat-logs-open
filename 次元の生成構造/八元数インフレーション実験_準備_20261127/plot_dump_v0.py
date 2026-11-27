#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ダンプのファイルを読んで、図を描くだけのプログラム。計算はしない。2026-11-27。

役目: これまで array_product_cumulative_general_v8.py の中に、演算・読出し・ダンプ・図化が一本で入っていた。
      その図化の部分だけを取り出したもの。読むのは v8 が書き出すダンプ (と、同じ形で書かれたダンプ)。
      図の並びと描き方は v8 と同じ。図化の行は v8 から書き写したもので、変えたのは図の保存先 (引数で受け取る) だけ。

v8 との違い (ダンプを読んで描くことから来るもの):
  - 描く値は、ダンプに書かれた桁の数 (a は 17 桁、θ は 14 桁、θn は 10 桁) である。v8 は計算した値 (60 桁) をそのまま使っていた。
  - 八元数一つと読む場合に「全回を通じて同じ値を取る成分」を一つの組にまとめる所は、
    ダンプの [配列] に書かれた文字が同じかどうかで比べる (120 桁までは整数そのもの、それを超えると仮数 17 桁)。
    v8 は整数で厳密に比べていた。

使い方: python3 plot_dump_v0.py ダンプのファイル 図の保存先(.png)
標準出力に、v8 がダンプの末尾に書くのと同じ集計の行を出す。
"""
import os
import re
import sys

import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

mp.mp.dps = 60
DUMP, OUT = sys.argv[1], sys.argv[2]


def log(msg):
    pass


# ---------------------------------------------------------------- ダンプを読む
lines = open(DUMP, encoding="utf-8").read().split("\n")
m0 = re.match(r"^N = (\d+)。初期値 H の作り方: (.*)。回数 (\d+)。(?:読み方 1(\d+) \(.*\)。)?実部と読む成分: (\d+) 番目", lines[0])
assert m0, "ダンプの 1 行目が読めない: " + lines[0][:80]
N, HOW, N_MAX = int(m0.group(1)), m0.group(2), int(m0.group(3))
BLK = (int(m0.group(4)) + 1) if m0.group(4) else N
RE = (int(m0.group(5)) - 1) % BLK
MULTI = BLK < N
IMAG = [k for k in range(N) if k != RE]
if MULTI:
    NB = N // BLK
    QBLK = [list(range(q * BLK, (q + 1) * BLK)) for q in range(NB)]
    QNAME = ["前", "後"] if NB == 2 else ["第%d" % (q + 1) for q in range(NB)]
    KIND = {4: "四元数", 8: "八元数"}.get(BLK, "%d 成分の元" % BLK)
    CNT = {2: "二つ", 3: "三つ"}.get(NB, " %d 個" % NB)
    RNG = lambda q: "成分 %d〜%d" % (QBLK[q][0] + 1, QBLK[q][-1] + 1)
    NSHOW = min(NB, 3)


def parse_readout(body):
    if body.startswith("全成分 0"):
        return (None, None, None)
    mm = re.match(r"^a = (\S+)   θ = (\S+) rad \((\S+)°\)   θn = (.*)$", body)
    assert mm, "読出しの行が読めない: " + body[:80]
    tn = None if mm.group(4).startswith("向き不定") else [mp.mpf(v) for v in mm.group(4).strip("()").split(", ")]
    return (mp.mpf(mm.group(1)), mp.mpf(mm.group(2)), tn)


history, rec = [], []
k = 1
while k < len(lines):
    mh = re.match(r"^n = (\d+)    \(H を (\d+) 個掛けたもの\)", lines[k])
    if not mh:
        k += 1
        continue
    n, m = int(mh.group(1)), int(mh.group(2))
    k += 2
    assert lines[k].startswith("[配列]"), "[配列] の行が無い (n = %d)" % n
    arr = []
    for i in range(N):
        ml = re.match(r"^  行%d: (.*)$" % (i + 1), lines[k + 1 + i])
        arr.append(ml.group(1).split("  "))
        assert len(arr[-1]) == N
    k += 1 + N
    assert lines[k].startswith("[行の指数部"), "[行の指数部] の行が無い (n = %d)" % n
    k += 1
    ex = []
    for i in range(N):
        if MULTI:
            pair = []
            for q in range(NB):
                ml = re.match(r"^  行%d %s: (.*)$" % (i + 1, QNAME[q]), lines[k])
                assert ml, "行が読めない (n = %d): %s" % (n, lines[k][:60])
                pair.append(parse_readout(ml.group(1)))
                k += 1
            ex.append(pair)
        else:
            ml = re.match(r"^  行%d: (.*)$" % (i + 1), lines[k])
            assert ml, "行が読めない (n = %d): %s" % (n, lines[k][:60])
            ex.append(parse_readout(ml.group(1)))
            k += 1
    history.append(arr)
    rec.append((n, m, ex))
assert len(rec) == N_MAX, "読めた回数 %d が、1 行目の回数 %d と合わない" % (len(rec), N_MAX)

# ---------------------------------------------------------------- 図化
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
C1, C2 = "#2a78d6", "#eb6834"
plt.rcParams["font.family"] = ["Hiragino Sans", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
ns = [r[0] for r in rec]
ROW = 1                                                        # 楕円軌道として読む行 (行 2)

log("演算 終了。図化 開始")
if MULTI:
    # ------------------------------------------------------------ 読み方 13 の図 (並びは読み方 17 と同じ)
    NAN = float("nan")
    npan = 2 + 3 * NSHOW
    nrow = (npan + 1) // 2
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

    ents = [(i, q) for q in range(NB) for i in range(N)]           # N × 塊の数 本: (行, 塊)
    ename = lambda e: "%s 行%d" % (QNAME[e[1]], e[0] + 1)

    def group_label(es):
        parts = []
        for q in range(NB):
            rows_ = [i + 1 for (i, qq) in es if qq == q]
            if rows_:
                parts.append(QNAME[q] + " 行" + (("%d〜%d" % (rows_[0], rows_[-1])) if (len(rows_) > 1 and rows_ == list(range(rows_[0], rows_[-1] + 1))) else ",".join(map(str, rows_))))
        return " / ".join(parts)

    # 1. 振幅を ln で読む (全成分 0 の四元数は a が定義できないので、その回は点を描かない)
    ax = flat[0]
    style(ax, "1. 振幅を ln で読む:%sごとの a = ln A(%d 行 × %s)" % (KIND, N, "・".join(QNAME)), "回数 n", "a = ln A(縦軸は対数の目盛り)")
    a_ser = {e: [(float(rec[t][2][e[0]][e[1]][0]) if rec[t][2][e[0]][e[1]][0] is not None else None) for t in range(len(rec))] for e in ents}
    akinds = {}
    for e in ents:
        akinds.setdefault(tuple(None if v is None else round(v, 9) for v in a_ser[e]), []).append(e)
    for key, es in sorted(akinds.items(), key=lambda kv: ents.index(kv[1][0])):
        ax.semilogy(ns, [NAN if (v is None or v <= 0) else v for v in key], color=C1, linewidth=1.3, alpha=0.7, marker="o", markersize=4)
    n_zero = sum(1 for e in ents for v in a_ser[e] if v is None)
    ax.set_xticks(ns)
    ax.text(0.03, 0.93, "a の動き方は %d 種類(%d 本中)。全成分 0 で a が定義できない点は %d 個(描いていない)" % (len(akinds), len(ents), n_zero),
            transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color=INK2)

    # 2. 四元数ごとの θ
    ax = flat[1]
    style(ax, "2. %sごとの θ(主値。%s を実部と読む)" % (KIND, " と".join("成分 %d" % (QBLK[q][RE] + 1) for q in range(NB))), "回数 n", "θ [度]")
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
        txt = "動き方が %d 種類(全 %d 本を描画)。1 回目 → %d 回目の θ:\n" % (len(ordered_kinds), len(ents), ns[-1])
        txt += "\n".join("  ".join(items[q:q + 3]) for q in range(0, min(len(items), 18), 3))
        ax.text(0.02, 0.99, txt, transform=ax.transAxes, ha="left", va="top", fontsize=8.5, color=INK2)
    ax.set_xticks(ns)
    ax.set_ylim(-8, 262)
    ax.set_yticks([0, 45, 90, 135, 180])
    ax.axhline(180, color=AXIS, linewidth=1)

    # 3〜. 行 2 の向き θn を成分どうしで描く (虚部の最初の三成分)。塊ごとに 3 面。先頭から三つまでの塊
    start = 2 if len(rec) > 2 else 0
    summary = []
    for q in range(NSHOW):
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
    for ax in flat[npan:]:
        ax.axis("off")
    fig.suptitle("N = %d(初期値:%s)の累積の積:%s%s(%s)と読む。実部は%s"
                 % (N, HOW, KIND, CNT, "、".join(RNG(q) for q in range(NB)), " と".join("成分 %d" % (QBLK[q][RE] + 1) for q in range(NB))),
                 fontsize=14.5, color=INK, x=0.07, ha="left", y=1 - 0.3 / (7.3 * nrow))
    out = OUT
    fig.savefig(out, dpi=150, facecolor=SURF)
    print("=" * 110)
    for q, nv, nd, tot in summary:
        print("行 2 の%s (成分 %d〜%d): n>=3 の %d 回のうち向きが描けるのは %d 回、違う位置は %d 個" % (QNAME[q], QBLK[q][0] + 1, QBLK[q][-1] + 1, tot, nv, nd))
    print("図を保存: " + os.path.relpath(out, os.path.dirname(os.path.abspath(__file__))))
    if len(rec) >= 3:
        d = lambda v: "無し" if v is None else mp.nstr(mp.degrees(v), 8, min_fixed=-5, max_fixed=6, strip_zeros=True)
        for q in range(NB):
            print("一回掛けたあと (n = 3) の θ [度]、%s: " % QNAME[q] + "  ".join("行%d %s" % (i + 1, d(rec[2][2][i][q][1])) for i in range(N)))
        print("面 1 の a の動き方の種類: %d、面 2 の θ の動き方の種類: %d (%d 本中)" % (len(akinds), len(ordered_kinds), len(ents)))
    log("終了")
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
out = OUT
fig.savefig(out, dpi=150, facecolor=SURF)
print("=" * 110)
print("行 2 の θn の成分の組 (全回を通じて同じ値): " + " / ".join(gname(g) for g in groups))
print("n=%d〜%d の %d 回分で、向きの違う位置は %d 個" % (vns[0], vns[-1], len(vns), n_distinct))
print("図を保存: " + os.path.relpath(out, os.path.dirname(os.path.abspath(__file__))))
if len(rec) >= 3:
    print("一回掛けたあと (n = 3) の各行の θ [度]: " + "  ".join("行%d %s" % (i + 1, mp.nstr(mp.degrees(rec[2][2][i][1]), 8, min_fixed=-5, max_fixed=6, strip_zeros=True)) for i in range(N)))
    print("面 2 の動き方の種類: %d (行の数 %d)" % (len(ordered_kinds), N))
log("終了")
