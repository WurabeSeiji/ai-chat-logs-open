#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v2 (2026-11-27): exponentiate_array_v1.py の桁数を、引数で指定できるようにしたもの。計算の方法は一字も変えていない。
  保存の桁 D を三つ目の引数で受け取る (省くと 200 = v1 と同じ)。計算の桁は D + 20、検算の合格の基準は 10^-(D-20)。
  保存先は results/exponent_data_v2/。D = 200 のときのデータは、プログラム名の欄を除いて v1 と同じになる。
  理由: 2 倍ずつの走行を長く続けると、掛けた個数 m の桁数ぶんだけ角度の桁を使う。長く走らせるには桁が要る。
  使い方: python3 exponentiate_array_v2.py 初期値のファイル 名前 [保存の桁]
以下は v1・v0 の説明のまま。

v1 (2026-11-27): exponentiate_array_v0.py の桁数を上げたもの。計算の方法は一字も変えていない。
  変えたのは次の設定だけ: 計算の桁 80 → 220、保存の桁 70 → 200、検算の合格の基準 1e-60 → 1e-180、
  保存先 results/exponent_data/ → results/exponent_data_v1/、データに書くプログラム名。
  理由: 積による計算のダンプは成分が 120 桁までは整数をそのまま書いている。それを指数部から再現して照合するため。
以下は v0 の説明のまま (桁数の記述は v0 のもの)。

汎用の指数化プログラム。2026-11-27。

役目 (三段のうちの二段目):
  1. 初期化      : ±1 の N×N を作る (rowcol_zero_init_pm1_v1.py。このプログラムは触れない)
  2. 指数化      : 配列 H を読み、指数部のデータに直す  ← このプログラム
  3. 読出しと図化: 指数部のデータと m から、H を m 個掛けた配列と、その読出しを出す (別のプログラム)

することは一つ。正方の整数の配列 H を、次の形に分ける。 m は掛けた個数 (m >= 1)。

    H を m 個掛けたもの  =  Σ(実の部分)  s^m × exp(m × a) × E
                          + Σ(回る部分)  exp(m × a) × ( cos(m × φ) × P + sin(m × φ) × J )

    a = ln|λ| (1 個掛けるごとの伸び)、φ = 角度 (1 個掛けるごとの回転、0 < φ < π)、s = ±1 (実の固有値の符号)。
    E、P、J は H から一度だけ求める N×N の実数の配列。固有値 0 の向きは、1 個掛けた時点で消えるので式に現れない。
  指数部で言えば、掛け算は「a と φ を m 倍する」足し算になる。

方法:
  - 特性多項式を整数係数のまま因数分解する (厳密。sympy)。これが記録として残す厳密な量である。
  - 「重複のない因子の積」に H を代入して 0 になるかを、整数のまま確かめる (厳密)。
    0 にならない配列は、上の形に分けられない (掛けても消えも回りもしない、ずれて残る成分を持つ)。
    その場合は、無理に当てはめず、分けられないと報告して止まる (終了コード 2)。
  - 固有値は、因子の根を mpmath で 80 桁で求める。固有値 λ ごとの配列は Π_{μ≠λ} (H - μ) / (λ - μ)。
    複素共役の組 λ、λ~ は、P = 2 Re(E_λ)、J = -2 Im(E_λ) の実数の配列二つにまとめる。
  - 大きさと角度は割り切れない数なので、有限の桁 (80 桁で計算、70 桁で保存) で持つ。
  割り算・正規化・乱数は、演算としては入っていない (上の分け方の計算の中の割り算だけ)。

検算 (このプログラムの中で必ず行う。通らなければデータを書き出さずに止まる):
  - 部分の配列を全部足すと、「固有値 0 の向きを消す配列」になり、それを H に掛けても H のままであること。
  - 上の式で m = 1、2、3、8 の配列を組み直し、H を実際に整数で掛けたものと比べる。
    差を「その回の成分の最大値に対する比」で出す。1e-60 より大きければ不合格。

入力: 初期値のファイル (1 行目が # で始まる説明、2 行目以降が配列。成分は空白区切りの整数)。
出力: results/exponent_data/{名前}.json  指数部のデータ (数は 70 桁の文字列)
      標準出力                            報告 (results/exponent_data/{名前}.txt に保存して使う)
使い方: python3 exponentiate_array_v0.py 初期値のファイル 名前
"""
import json
import os
import sys

import mpmath as mp
import sympy as sp

SAVE_DIGITS = int(sys.argv[3]) if len(sys.argv) > 3 else 200
assert SAVE_DIGITS >= 100, "保存の桁は 100 以上"
mp.mp.dps = SAVE_DIGITS + 20
TOL = mp.mpf(10) ** -(SAVE_DIGITS - 20)


def read_array(path):
    rows, head = [], ""
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            head = line[1:].strip()
            continue
        rows.append([int(v) for v in line.split()])
    n = len(rows)
    assert n > 0 and all(len(r) == n for r in rows), "正方の配列ではない"
    return rows, head


def matmul_int(A, B):
    n = len(A)
    return [[sum(A[i][t] * B[t][j] for t in range(n)) for j in range(n)] for i in range(n)]


def main():
    path, name = sys.argv[1], sys.argv[2]
    Hi, head = read_array(path)
    N = len(Hi)
    x = sp.symbols("x")
    Hs = sp.Matrix(Hi)
    print("入力: %s" % path)
    print("説明: %s" % head)
    print("N = %d" % N)

    # ---- 特性多項式の因数分解 (厳密)
    chi = Hs.charpoly(x).as_expr()
    const, facs = sp.factor_list(chi, x)
    facs = sorted(facs, key=lambda fe: (sp.degree(fe[0], x), str(fe[0])))
    print("特性多項式 (厳密): " + " × ".join("(%s)^%d" % (f, e) if e > 1 else "(%s)" % f for f, e in facs))

    # ---- 重複のない因子の積が H を消すか (厳密)
    mu0 = sp.Integer(1)
    for f, e in facs:
        mu0 = sp.expand(mu0 * f)
    coeffs = sp.Poly(mu0, x).all_coeffs()
    acc = sp.zeros(N, N)
    for c in coeffs:                                  # ホーナー法。整数のまま
        acc = acc * Hs + c * sp.eye(N)
    semisimple = acc.is_zero_matrix
    print("重複のない因子の積 %s に H を代入すると 0 か (厳密): %s" % (sp.factor(mu0), "はい" if semisimple else "いいえ"))
    if not semisimple:
        print("この配列は、実の部分と回る部分だけには分けられない。指数部のデータは書き出さない。")
        sys.exit(2)

    # ---- 固有値 (80 桁) と、固有値ごとの配列
    lams = []                                         # (λ, 重複度, 因子)
    for f, e in facs:
        cf = [mp.mpf(int(c)) for c in sp.Poly(f, x).all_coeffs()]
        roots = [mp.mpf(0)] if (len(cf) == 2 and cf[1] == 0) else mp.polyroots(cf, maxsteps=500, extraprec=400)
        for r in roots:
            lams.append((mp.mpc(r), int(e), f))
    H = mp.matrix(Hi)
    I = mp.eye(N)

    def projector(lam):
        E = I.copy()
        for mu, _, _ in lams:
            if mu is not lam and mu != lam:
                E = E * (H - mu * I) / (lam - mu)
        return E

    parts = []
    used = [False] * len(lams)
    zero_mult = 0
    for i, (lam, e, f) in enumerate(lams):
        if used[i]:
            continue
        used[i] = True
        if lam == 0:
            zero_mult = e
            continue
        E = projector(lam)
        if abs(mp.im(lam)) < mp.mpf(10) ** -70 * max(1, abs(lam)):          # 実の固有値
            lr = mp.re(lam)
            parts.append(dict(kind="real", sign=(1 if lr > 0 else -1), a=mp.log(abs(lr)), mult=e, factor=str(f),
                              E=E.apply(mp.re)))
        else:                                                              # 複素共役の組
            for j in range(i + 1, len(lams)):
                if not used[j] and abs(lams[j][0] - mp.conj(lam)) < mp.mpf(10) ** -60 * abs(lam):
                    used[j] = True
                    break
            else:
                raise AssertionError("複素共役の相手が見つからない: %s" % lam)
            if mp.im(lam) < 0:
                lam = mp.conj(lam)
                E = E.apply(mp.conj)
            parts.append(dict(kind="pair", a=mp.log(abs(lam)), phi=mp.arg(lam), mult=e, factor=str(f),
                              P=E.apply(lambda z: 2 * mp.re(z)), J=E.apply(lambda z: -2 * mp.im(z))))
    parts.sort(key=lambda p: (-p["a"], p.get("phi", mp.mpf(0)), p.get("sign", 0)))

    def rebuild(m):
        S = mp.zeros(N, N)
        for p in parts:
            g = mp.exp(m * p["a"])
            if p["kind"] == "real":
                S += (p["sign"] ** m) * g * p["E"]
            else:
                S += g * (mp.cos(m * p["phi"]) * p["P"] + mp.sin(m * p["phi"]) * p["J"])
        return S

    # ---- 報告
    f70 = lambda v: mp.nstr(v, 30)
    print("固有値 0 の向き: %d 本 (1 個掛けた時点で消える)" % zero_mult)
    print("部分の数: %d" % len(parts))
    for k, p in enumerate(parts, 1):
        if p["kind"] == "real":
            print("  部分 %d [実]  符号 %+d  大きさ %s  伸び a = ln|λ| = %s  本数 %d  因子 %s"
                  % (k, p["sign"], f70(mp.exp(p["a"])), f70(p["a"]), p["mult"], p["factor"]))
        else:
            print("  部分 %d [回る] 大きさ %s  伸び a = ln|λ| = %s  角度 φ = %s rad (%s 度)  回る平面の数 %d  因子 %s"
                  % (k, f70(mp.exp(p["a"])), f70(p["a"]), f70(p["phi"]), f70(mp.degrees(p["phi"])), p["mult"], p["factor"]))

    # ---- 検算
    ok = True
    Ptot = mp.zeros(N, N)
    for p in parts:
        Ptot += p["E"] if p["kind"] == "real" else p["P"]
    d0 = max(abs(v) for v in (Ptot * H - H))
    print("検算 1: 部分の配列の和を H に掛けても H のまま。差の最大 = %s" % mp.nstr(d0, 5))
    ok &= d0 < TOL * N
    Pw = [r[:] for r in Hi]
    powers = {1: Pw}
    for m in range(2, 9):
        Pw = matmul_int(Pw, Hi)
        powers[m] = Pw
    for m in (1, 2, 3, 8):
        S = rebuild(m)
        big = max(abs(v) for r in powers[m] for v in r)
        d = max(abs(S[i, j] - powers[m][i][j]) for i in range(N) for j in range(N))
        rel = d / big if big else d
        print("検算 2: m = %d の配列を式から組み直し、整数で掛けたものと比べる。差の最大 ÷ 成分の最大 = %s" % (m, mp.nstr(rel, 5)))
        ok &= rel < TOL
    print("判定: %s" % ("合格" if ok else "不合格"))
    if not ok:
        print("検算に通らないので、指数部のデータは書き出さない。")
        sys.exit(1)

    # ---- 書き出し
    s = lambda v: mp.nstr(v, SAVE_DIGITS)
    mat = lambda M: [[s(M[i, j]) for j in range(N)] for i in range(N)]
    out = dict(program="exponentiate_array_v2.py", source=path, description=head, N=N,
               charpoly=" × ".join("(%s)^%d" % (f, e) for f, e in facs), zero_multiplicity=int(zero_mult),
               formula="H^m = sum_real sign^m exp(m a) E + sum_pair exp(m a) (cos(m phi) P + sin(m phi) J), m >= 1",
               digits=SAVE_DIGITS, H=Hi, parts=[])
    for p in parts:
        q = dict(kind=p["kind"], a=s(p["a"]), mult=p["mult"], factor=p["factor"])
        if p["kind"] == "real":
            q.update(sign=p["sign"], E=mat(p["E"]))
        else:
            q.update(phi=s(p["phi"]), P=mat(p["P"]), J=mat(p["J"]))
        out["parts"].append(q)
    odir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "exponent_data_v2")
    os.makedirs(odir, exist_ok=True)
    with open(os.path.join(odir, name + ".json"), "w", encoding="utf-8", newline="\n") as fo:
        json.dump(out, fo, ensure_ascii=False, indent=0)
    print("指数部のデータを保存: results/exponent_data_v2/%s.json (保存の桁 %d)" % (name, SAVE_DIGITS))


if __name__ == "__main__":
    main()
