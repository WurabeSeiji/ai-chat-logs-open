#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2x2x2 の立方体 U に 120 度の三つの数 (1, ω, ω²) を置き、積を n = 100 まで繰り返して U^n = I になるかを調べる。v1。
初期値を含め、n 回目ごとに全成分をダンプする。v0 からの変更は出力の追加だけで、演算は同じ。

各 n について三通りで出力する。
  厳密 : アイゼンシュタイン整数 x + yω (x, y は整数、ω = exp(2πi/3)、ω² = -1 - ω)。v0 と同じ。
  直交 : 複素数 Re + Im·i に読み戻した値。Re = (2x - y)/2、Im = (√3/2)·y。
  極   : 振幅 A と位相 θ [度、-180 < θ <= 180]。A = sqrt(x² - xy + y²)、θ = atan2(Im, Re)。
直交と極は、厳密な整数から最後に一度だけ変換する (途中の計算は整数のままなので桁落ちは無い)。
変換は mpmath の 50 桁で行い、倍精度と同じ有効 17 桁で表示する。倍精度の範囲 (約 1e308) を超える大きさも指数つきで表示できる。

升目は 4 種類 (d: 対角、a: 1 番目の軸だけ違う、b: 2 番目、c: 3 番目)、各 2 個。
同じ種類の 2 個が同じ値なら一つにまとめ、違う値なら升目の番号をつけて両方出す。

置き方 (升目 8 個):
  対角 (x,x,x)         -> d
  (y,x,x), y != x      -> a = 1      (1 番目の軸だけが違う)
  (x,y,x), y != x      -> b = ω      (2 番目の軸だけが違う)
  (x,x,y), y != x      -> c = ω²     (3 番目の軸だけが違う)

積:
  三項の積 (Bhattacharya-Mesner):  Prod(A,B,C)[i,j,k] = Σ_t A[i,t,k] B[i,j,t] C[t,j,k]
  二項の積 (t 積, Kilmer-Martin):   (A*B)[i,j,s] = Σ_{x+y=s mod 2} Σ_t A[i,t,x] B[t,j,y]

作業者 (Claude) が決めたこと。指示に無かったので、決め打ちにせず全部走らせる。
  [1] 対角 d の初期値: d = 0 と d = 1 の二通り。
  [2] 「n 回の積」の取り方 (三項の積は結合的でないので取り方で変わる):
        T2  二項の積  P_1 = U, P_{n+1} = P_n * U                         (U を n 個掛けたもの)
        L3  三項の積  P_1 = U, P_{n+1} = Prod(P_n, U, U)                 (左に積む。U は 2n-1 個)
        M3  三項の積  P_1 = U, P_{n+1} = Prod(U, P_n, U)                 (真ん中に積む)
        R3  三項の積  P_1 = U, P_{n+1} = Prod(U, U, P_n)                 (右に積む)
        C3  三項の積  P_1 = U, P_{n+1} = Prod(P_n, P_n, P_n)             (自分自身を三つ掛ける)
  [3] 判定 (読出し。演算には入れない):
        =U      P_n が U と全成分で等しい
        ∝U      P_n が U の定数倍 (定数は 0 でない一つの数)。定数の大きさが 1 なら「位相だけ」
        =I, ∝I  二項の積の単位の立方体 I (0 枚目が単位行列) と等しい / その定数倍
  [4] C3 は成分の桁数が毎回約 3 倍になり、厳密な整数では n = 100 まで持てない。
      桁数が EXACT_LIMIT を超えた回から先は、厳密をやめて mpmath の 50 桁の浮動小数 (指数は無制限) で続ける。
      その区間は「近似」と明記し、判定はしない。他の取り方は n = 100 まで全部厳密。

使い方: python3 cube222_powers_n100_v1.py   (結果は標準出力。results/ に保存して使う)
"""
import itertools
import sys

import mpmath as mp

mp.mp.dps = 50
SQ3H = mp.sqrt(3) / 2

N_MAX = 100
EXACT_LIMIT = 3000        # 厳密な整数で持つ桁数の上限 (C3 だけが超える)
n = 2
IDX = list(itertools.product(range(n), repeat=3))


class E:
    """アイゼンシュタイン整数 x + yω。"""
    __slots__ = ("x", "y")

    def __init__(self, x=0, y=0):
        self.x, self.y = x, y

    def __add__(self, o):
        return E(self.x + o.x, self.y + o.y)

    def __mul__(self, o):
        # (x1 + y1 ω)(x2 + y2 ω) = x1x2 + (x1y2 + y1x2) ω + y1y2 ω²,  ω² = -1 - ω
        return E(self.x * o.x - self.y * o.y, self.x * o.y + self.y * o.x - self.y * o.y)

    def __eq__(self, o):
        return self.x == o.x and self.y == o.y

    def is_zero(self):
        return self.x == 0 and self.y == 0

    def norm(self):
        """|z|² = x² - xy + y² (整数)。"""
        return self.x * self.x - self.x * self.y + self.y * self.y

    def __str__(self):
        if self.y == 0:
            return str(self.x)
        if self.x == 0:
            return "%dω" % self.y if abs(self.y) != 1 else ("ω" if self.y == 1 else "-ω")
        yy = ("+%dω" % self.y) if self.y > 0 else ("%dω" % self.y)
        if abs(self.y) == 1:
            yy = "+ω" if self.y == 1 else "-ω"
        return "%d%s" % (self.x, yy)

    def digits(self):
        """10 進の桁数。大きい数は 2 進の桁数から見積もる (±1 の誤差がありうる)。"""
        def dg(v):
            v = abs(v)
            return len(str(v)) if v < 10 ** 300 else int(v.bit_length() * 0.30102999566398) + 1
        return max(dg(self.x), dg(self.y))


ZERO, ONE, W = E(0, 0), E(1, 0), E(0, 1)
W2 = W * W


def make_U(d):
    U = {}
    for i, j, k in IDX:
        if i == j == k:
            U[i, j, k] = d
        elif j == k:
            U[i, j, k] = ONE
        elif i == k:
            U[i, j, k] = W
        else:
            U[i, j, k] = W2
    return U


def bm(A, B, C):
    R = {}
    for i, j, k in IDX:
        s = ZERO
        for t in range(n):
            s = s + A[i, t, k] * B[i, j, t] * C[t, j, k]
        R[i, j, k] = s
    return R


def tprod(A, B):
    R = {}
    for i, j, s in IDX:
        v = ZERO
        for t in range(n):
            for x in range(n):
                v = v + A[i, t, x] * B[t, j, (s - x) % n]
        R[i, j, s] = v
    return R


IDT = {idx: (ONE if (idx[0] == idx[1] and idx[2] == 0) else ZERO) for idx in IDX}


def proportional(P, Q):
    """P = λ Q となる λ があるか (割り算をせず、たすき掛けで判定)。あれば (True, |λ|² の分子, 分母)。"""
    ref = None
    for idx in IDX:
        if Q[idx].is_zero():
            if not P[idx].is_zero():
                return False, None, None
        elif ref is None:
            ref = idx
    if ref is None or P[ref].is_zero():
        return False, None, None
    for idx in IDX:
        if not (P[idx] * Q[ref] == P[ref] * Q[idx]):
            return False, None, None
    return True, P[ref].norm(), Q[ref].norm()


def judge(P, U):
    out = []
    if all(P[i] == U[i] for i in IDX):
        out.append("=U")
    else:
        ok, pn, qn = proportional(P, U)
        if ok:
            out.append("∝U(位相だけ)" if pn == qn else "∝U(大きさも変わる)")
    if all(P[i] == IDT[i] for i in IDX):
        out.append("=I")
    else:
        ok, pn, qn = proportional(P, IDT)
        if ok:
            out.append("∝I(位相だけ)" if pn == qn else "∝I(大きさも変わる)")
    return " ".join(out) if out else "-"


def dump(P):
    return "  ".join("(%d%d%d)=%s" % (i, j, k, P[i, j, k]) for i, j, k in IDX)


CLASS_CELLS = {"d": [(0, 0, 0), (1, 1, 1)], "a": [(0, 1, 1), (1, 0, 0)],
               "b": [(0, 1, 0), (1, 0, 1)], "c": [(0, 0, 1), (1, 1, 0)]}


def to_mpc(v):
    """厳密な x + yω、または近似の mpc を、mpc (Re + Im i) にする。"""
    if isinstance(v, E):
        return mp.mpc(mp.mpf(2 * v.x - v.y) / 2, SQ3H * v.y)
    return v


def fnum(x):
    return mp.nstr(x, 17, min_fixed=-5, max_fixed=17, strip_zeros=True)


def rect_str(z):
    return "%s %s %s i" % (fnum(z.real), "+" if z.imag >= 0 else "-", fnum(abs(z.imag)))


def polar_str(z):
    A = abs(z)
    if A == 0:
        return "A=0, θ=(不定)"
    th = mp.degrees(mp.arg(z))
    return "A=%s, θ=%s°" % (fnum(A), mp.nstr(th, 12, min_fixed=-5, max_fixed=6, strip_zeros=True))


def same(u, v):
    if isinstance(u, E):
        return u == v
    return u == v


def class_lines(P):
    rect, pol = [], []
    for name in "dabc":
        c1, c2 = CLASS_CELLS[name]
        if same(P[c1], P[c2]):
            items = [(name, P[c1])]
        else:
            items = [("%s(%d%d%d)" % ((name,) + c1), P[c1]), ("%s(%d%d%d)" % ((name,) + c2), P[c2])]
        for lab, v in items:
            z = to_mpc(v)
            rect.append("%s = %s" % (lab, rect_str(z)))
            pol.append("%s: %s" % (lab, polar_str(z)))
    return "       直交: " + " ;  ".join(rect), "       極  : " + " ;  ".join(pol)


def bm_mp(A, B, C):
    R = {}
    for i, j, k in IDX:
        R[i, j, k] = mp.fsum(A[i, t, k] * B[i, j, t] * C[t, j, k] for t in range(n))
    return R


SCHEMES = [
    ("T2", "二項の積  P_{n+1} = P_n * U", lambda P, U: tprod(P, U)),
    ("L3", "三項の積  P_{n+1} = Prod(P_n, U, U)", lambda P, U: bm(P, U, U)),
    ("M3", "三項の積  P_{n+1} = Prod(U, P_n, U)", lambda P, U: bm(U, P, U)),
    ("R3", "三項の積  P_{n+1} = Prod(U, U, P_n)", lambda P, U: bm(U, U, P)),
    ("C3", "三項の積  P_{n+1} = Prod(P_n, P_n, P_n)", lambda P, U: bm(P, P, P)),
]

summary = []
for dname, d in (("d=0", ZERO), ("d=1", ONE)):
    U = make_U(d)
    for key, desc, step in SCHEMES:
        print("=" * 100)
        print("初期値 %s, (a,b,c) = (1, ω, ω²)   取り方 %s: %s" % (dname, key, desc))
        print("=" * 100)
        P = U
        hits = []
        approx_from = None
        ended = None
        for m in range(1, N_MAX + 1):
            if m > 1:
                if approx_from is None:
                    P = step(P, U)
                else:
                    P = bm_mp(P, P, P)          # 近似に切り替わるのは C3 だけ
            if approx_from is None:
                dg = max(P[i].digits() for i in IDX)
                j = judge(P, U) if m > 1 else "(初期値)"
                if m > 1 and j != "-":
                    hits.append((m, j))
                print("n=%3d  判定 %-22s 最大桁数 %d" % (m, j, dg))
                if dg > 400:
                    print("       厳密: (成分が 400 桁を超えるので省略。成分ごとの桁数: %s)"
                          % " ".join(str(P[i].digits()) for i in IDX))
                else:
                    print("       厳密: " + dump(P))
                zero = all(P[i].is_zero() for i in IDX)
            else:
                print("n=%3d  判定 (近似の区間。判定しない)   mpmath 50 桁" % m)
                zero = False
            r, p = class_lines(P)
            print(r)
            print(p)
            if zero:
                ended = m
                print("   -> 全成分が 0 になった (n=%d)。以後も 0 のまま。" % m)
                break
            if approx_from is None and dg > EXACT_LIMIT:
                approx_from = m + 1
                P = {i: to_mpc(P[i]) for i in IDX}
                print("   -> 桁数が %d を超えた。n=%d から先は厳密をやめ、mpmath 50 桁の近似で続ける。" % (EXACT_LIMIT, m + 1))
        summary.append((dname, key, hits, approx_from, ended))
        sys.stdout.flush()

print("=" * 100)
print("まとめ: n = 2..%d で、U や I に等しい、または定数倍になった回" % N_MAX)
print("=" * 100)
for dname, key, hits, approx_from, ended in summary:
    tail = ""
    if approx_from is not None:
        tail += "  (n=%d から近似。判定は n=%d まで)" % (approx_from, approx_from - 1)
    if ended is not None:
        tail += "  (n=%d で全成分 0)" % ended
    if hits:
        ns = [h[0] for h in hits]
        kinds = sorted(set(h[1] for h in hits))
        body = "%d 回 (n = %s%s)  種類: %s" % (len(ns), ", ".join(map(str, ns[:8])), ", ..." if len(ns) > 8 else "", " / ".join(kinds))
    else:
        body = "なし"
    print("  %s  %s : %s%s" % (dname, key, body, tail))
