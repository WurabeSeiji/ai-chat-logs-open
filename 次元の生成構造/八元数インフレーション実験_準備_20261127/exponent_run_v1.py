#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v1 (2026-11-27): exponent_run_v0.py を、長い走行に使えるようにしたもの。配列の組み直しと読出しの式は一字も変えていない。
  (a) 計算の桁を、指数部のデータの保存の桁 + 20 にする (v0 は 220 桁に固定。データが 200 桁なら v0 と同じ)。
  (b) 桁が足りなくなる回の手前で止まる。掛けた個数 m が大きいと、角度 m × φ の計算に m の桁数ぶんの桁を使う。
      「データの保存の桁 - m の桁数」が 170 を下回る回は計算せず、その旨を標準エラーに書いて止まる (終了コード 3)。
      170 は、Claude が置いた 0 の決まり (ii) (1e-150) が効くだけの余裕を残すための数である。
      これにより、走らせられる回数は おおよそ (保存の桁 - 170) ÷ 0.301 + 2 回になる。
  (c) 見出しに書く「成分の最大桁数」は、2 進の桁数が 10^15 を超えるときは多倍長で数える
      (v0 の式のままだと、桁数が大きくなりすぎたときに計算できなくなるため。10^15 以下では v0 と同じ式・同じ値)。
以下は v0 の説明のまま。

指数部のデータから、積を取らずに、H を m 個掛けた配列と読出しを出す実行系。2026-11-27。

役目 (三段のうちの三段目の前半):
  1. 初期化 (rowcol_zero_init_pm1_v1.py)  2. 指数化 (exponentiate_array_v1.py)
  3. 実行系 (このプログラム): 指数部のデータ (.json) を読み、回ごとに配列を組み直し、読出しをして、ダンプを書く。
     図化は別のプログラム (plot_dump_v0.py) がダンプを読んで行う。

演算: 掛け算はしない。指数部のデータの式
    H を m 個掛けたもの = Σ(実の部分) s^m × exp(m × a) × E + Σ(回る部分) exp(m × a) × (cos(m × φ) × P + sin(m × φ) × J)
  に m を入れるだけ (a と φ を m 倍する足し算)。計算は 220 桁。
回と m の対応は、積による計算 (array_product_cumulative_general_v8.py、O_n = O_1 × … × O_{n-1}) と同じ:
  1 回目と 2 回目は m = 1、3 回目以降は m = 2^(n-2)。
読出しとダンプの形は v8 と同じ。読出しとダンプの行は v8 から書き写したもの (読出しは 60 桁)。
  違うのは [配列] の見出しの但し書きだけ (「計算は整数で厳密」を「指数部のデータから組み直した値」に変えた)。

Claude が置いた決まり (積による計算では要らなかったもの。有限の桁で 0 と整数を扱うため):
  (i)   指数部のデータの配列 E、P、J の成分で、大きさが 1e-150 より小さいものは 0 とする。
  (ii)  組み直した成分で、大きさが「足し合わせた各項の大きさの合計」の 1e-150 倍以下のものは 0 とする
        (項どうしが打ち消し合って 0 になる成分)。
  (iii) その回の成分がすべて 120 桁以内のときは、成分をいちばん近い整数にする。
        整数からのずれが 1e-50 より大きい成分があれば、止まる。
経過のログは標準エラーに 1 行ずつ出す。

使い方: python3 exponent_run_v0.py 指数部のデータ(.json) [回数] [実部の番号] [塊の大きさ]
        (省くと回数 20、実部は 1 番目、塊の大きさ = N)
"""
import json
import re
import sys
import time

import mpmath as mp

sys.set_int_max_str_digits(0)
READ_DPS = 60
FULL_DIGITS = 120
T0 = time.time()


def log(msg):
    sys.stderr.write("[%s] 経過 %8.1f 秒  %s\n" % (time.strftime("%H:%M:%S"), time.time() - T0, msg))
    sys.stderr.flush()


DATA = json.load(open(sys.argv[1], encoding="utf-8"))
CALC_DPS = int(DATA["digits"]) + 20
GUARD = 170                                      # データの保存の桁 - m の桁数 が、これを下回る回は計算しない
mp.mp.dps = CALC_DPS
N = DATA["N"]
N_MAX = int(sys.argv[2]) if len(sys.argv) > 2 else 20
RE = (int(sys.argv[3]) - 1) if len(sys.argv) > 3 else 0
BLK = int(sys.argv[4]) if len(sys.argv) > 4 else N
mh = re.search(r"作り方: (.*)、生成: ", DATA["description"])
assert mh, "指数部のデータの説明から、初期値の作り方が読めない"
HOW = mh.group(1)
assert 1 <= BLK <= N and N % BLK == 0, "塊の大きさは N を割り切る数"
MULTI = BLK < N
assert 0 <= RE < BLK, "実部の番号は 1 から塊の大きさまで"
IMAG = [k for k in range(N) if k != RE]          # 虚部と読む成分 (番号の小さい順)
if MULTI:
    assert BLK >= 4, "塊の大きさは 4 以上"
    NB = N // BLK
    QBLK = [list(range(q * BLK, (q + 1) * BLK)) for q in range(NB)]
    QNAME = ["前", "後"] if NB == 2 else ["第%d" % (q + 1) for q in range(NB)]
    KIND = {4: "四元数", 8: "八元数"}.get(BLK, "%d 成分の元" % BLK)
    CNT = {2: "二つ", 3: "三つ"}.get(NB, " %d 個" % NB)
    RNG = lambda q: "成分 %d〜%d" % (QBLK[q][0] + 1, QBLK[q][-1] + 1)

# ---- 指数部のデータ (決まり (i): 1e-150 より小さい成分は 0)
EPS_C = mp.mpf(10) ** -150
EPS_Z = mp.mpf(10) ** -150
EPS_I = mp.mpf(10) ** -50


def load(M):
    out = [[mp.mpf(v) for v in r] for r in M]
    return [[(mp.mpf(0) if abs(v) < EPS_C else v) for v in r] for r in out]


PARTS = []
for p in DATA["parts"]:
    q = dict(kind=p["kind"], a=mp.mpf(p["a"]))
    if p["kind"] == "real":
        q.update(sign=p["sign"], E=load(p["E"]))
    else:
        q.update(phi=mp.mpf(p["phi"]), P=load(p["P"]), J=load(p["J"]))
    PARTS.append(q)


def rebuild(m):
    """H を m 個掛けた配列 (220 桁)。決まり (ii) で 0 を決める。"""
    S = [[mp.mpf(0)] * N for _ in range(N)]
    T = [[mp.mpf(0)] * N for _ in range(N)]          # 各項の大きさの合計
    for p in PARTS:
        gmag = mp.exp(m * p["a"])
        if p["kind"] == "real":
            c = (p["sign"] ** m) * gmag
            for i in range(N):
                for j in range(N):
                    t = c * p["E"][i][j]
                    S[i][j] += t
                    T[i][j] += abs(t)
        else:
            cc, ss = gmag * mp.cos(m * p["phi"]), gmag * mp.sin(m * p["phi"])
            for i in range(N):
                for j in range(N):
                    t1, t2 = cc * p["P"][i][j], ss * p["J"][i][j]
                    S[i][j] += t1 + t2
                    T[i][j] += abs(t1) + abs(t2)
    return [[(mp.mpf(0) if abs(S[i][j]) <= EPS_Z * T[i][j] else S[i][j]) for j in range(N)] for i in range(N)]


def ndigits(v):
    """成分の桁数。v8 の ndigits と同じ値になるように数える (10^400 以上は 2 進の桁数からの概算。v8 と同じ式)。"""
    if v == 0:
        return 1
    v = abs(v)
    if isinstance(v, int):
        return len(str(v)) if v < 10 ** 400 else int(v.bit_length() * 0.30102999566398) + 1
    with mp.workdps(READ_DPS):
        w = +v                                       # 60 桁に丸めてから数える (2 の冪ちょうどの成分で、桁が 1 つ少なく数えられるのを防ぐ)
    if w < mp.mpf(10) ** 400:
        return len(str(int(mp.nint(w))))
    bits = w.man.bit_length() + w.exp                                   # 2 進の桁数 = 仮数の桁数 + 指数
    if bits < 10 ** 15:
        return int(bits * 0.30102999566398) + 1
    with mp.workdps(len(str(bits)) + 30):
        return int(mp.floor(bits * mp.log10(2))) + 1


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
    im = [o[k] for k in range(BLK) if k != RE]
    vn = mp.sqrt(mp.fsum(x * x for x in im))
    if A == 0:
        return None, None, None
    a, th = mp.log(A), mp.atan2(vn, o[RE])
    if vn == 0:
        return a, th, ([mp.mpf(0)] * (BLK - 1) if th == 0 else None)
    return a, th, [th * x / vn for x in im]


if MULTI:
    print("N = %d。初期値 H の作り方: %s。回数 %d。読み方 1%d (%s%s: %s)。実部と読む成分: %s"
          % (N, HOW, N_MAX, BLK - 1, KIND, CNT, " と".join(RNG(q) for q in range(NB)), "と ".join("%d 番目" % (QBLK[q][RE] + 1) for q in range(NB))))
else:
    print("N = %d。初期値 H の作り方: %s。回数 %d。実部と読む成分: %d 番目" % (N, HOW, N_MAX, RE + 1))
rec = []
for n in range(1, N_MAX + 1):
    m = 1 if n <= 2 else 2 ** (n - 2)
    if int(DATA["digits"]) - len(str(m)) < GUARD:
        log("n = %d (m は %d 桁) は計算しない。指数部のデータの桁 (%d) が足りない。n = %d までで止める" % (n, len(str(m)), int(DATA["digits"]), n - 1))
        sys.stdout.flush()
        sys.exit(3)
    mp.mp.dps = CALC_DPS
    O = rebuild(m)
    dg = max(ndigits(v) for r in O for v in r)
    full = dg <= FULL_DIGITS
    if full:                                                   # 決まり (iii)
        R = [[int(mp.nint(v)) for v in r] for r in O]
        dev = max(abs(O[i][j] - R[i][j]) for i in range(N) for j in range(N))
        assert dev < EPS_I, "n = %d: 整数からのずれが大きい: %s" % (n, mp.nstr(dev, 5))
        O = R
        dg = max(ndigits(v) for r in O for v in r)
    mp.mp.dps = READ_DPS
    if not full:
        O = [[(0 if v == 0 else +v) for v in r] for r in O]   # 読出しの桁 (60 桁) に丸める。0 は整数の 0 にする (ダンプに 0 と書くため)
    print("=" * 110)
    print("n = %d    (H を %d 個掛けたもの)    成分の最大桁数 %d" % (n, m, dg))
    print("=" * 110)
    print("[配列]" + ("" if full else "  (成分が %d 桁を超えるので、仮数 e 指数 で表示。指数部のデータから組み直した値)" % FULL_DIGITS))
    for i, r in enumerate(O):
        print("  行%d: " % (i + 1) + "  ".join(show_int(v, full) for v in r))
    if not MULTI:
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
        print("[行の指数部 1%d] 各行を%s%sと読む。%s。q = exp(a + nθ)、実部は%s。a = ln A、θ は主値 (0〜π)。θn は残りの %d 成分を番号順に"
              % (BLK - 1, KIND, CNT, "、".join("%s = %s" % (QNAME[q], RNG(q)) for q in range(NB)),
                 " と".join("成分 %d (%s)" % (QBLK[q][RE] + 1, QNAME[q]) for q in range(NB)), BLK - 1))
        ex = []
        for i, r in enumerate(O):
            pair = []
            for q in range(NB):
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
    log("n = %d (m = %d) 終了。成分の最大桁数 %d" % (n, m, dg))
    sys.stdout.flush()
log("終了")
