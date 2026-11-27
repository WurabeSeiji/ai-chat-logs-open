#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
縦横の和が 0 になる、±1 だけの N×N の初期値 (N は 4 の倍数)。2026-11-27。

仕様 (木原さんの指示):
  - 成分は +1 と -1 だけ。N×N。N 個の成分を区別しない (実部・虚部という考えを入れない)。
  - どの行も、どの列も、足すと 0。全部 +1 の行や列を作らない。
  - 検算は Σ z_n² = 0。二つの成分 a, b を選び z_n = (n 行目の成分 a) + i (n 行目の成分 b) と読む。
    ただし「反転の相手どうし」(下記) を実部と虚部に組む読み方は範囲に入れない (木原さんの了承済み)。

作り方: N/2 本の列 A と、その符号を全部反転した -A を横に並べる。 M = [ A | -A ]
  成分 j と成分 j + N/2 が「反転の相手どうし」になる。
  どの行も「A の一行と、その反転」なので、足すと 0。A の列を、足すと 0 で互いの内積が 0 になるように取る。

  (甲) N = 4、または N が 8 の倍数のとき:
       半分の大きさの表 B (N/2 × N/2、±1、列どうしの内積が 0) を作り、  M = [ B  -B ]
                                                                           [ -B   B ]
       B の作り方は square_array(N/2)。行も、反転の相手 (i 行目と i + N/2 行目) を一つずつ持つ。
  (乙) それ以外の 4 の倍数 (12、20、28、…。半分の大きさの B が存在しない):
       N×N の表 H = square_array(N) の各行に自分の一列目の値を掛けて一列目を全部 +1 にし、
       後ろの N/2 本の列を A とする。 M = [ A | -A ]
       **この「後ろの N/2 本を取る」は Claude が決めた選び方である** (ほかの N/2 本でも条件は満たす)。
       乙では、列は反転の相手を持つが、行は反転の相手を持たない (N/2 が 4 の倍数でないと、持たせることができない)。
  N が 2 の冪のときは、甲と乙は同じ配列になる。

  square_array (±1 で行と列の内積が 0 になる正方の表。既知の作り方: 二進、平方剰余 1、平方剰余 2、積) は、
  balanced_init_pm1_v2.py の同名の関数を書き写したもの。この関数が返す表そのものは初期値に使わない
  (一行目や一列目が全部 +1 になるため)。作れない大きさのときは、でたらめな値を返さず、止まる。

乱数、割り算、正規化は入っていない。

使い方:
  python3 rowcol_zero_init_pm1_v0.py 8          N = 8 の初期値と検算を表示
  python3 rowcol_zero_init_pm1_v0.py --export   4 から 100 までの 4 の倍数で作れるもの全部を
                                                results/rowcol_zero_init/ に書き出し、検算の一覧を表示
データの形式: 1 行目は説明 (# で始まる)。2 行目以降が配列の 1 行ずつ。成分は +1 / -1 を空白区切り。
"""
import hashlib
import os
import sys


def is_prime(x):
    if x < 2:
        return False
    k = 2
    while k * k <= x:
        if x % k == 0:
            return False
        k += 1
    return True


def jacobsthal(p):
    sq = {(x * x) % p for x in range(1, p)}
    chi = lambda x: 0 if x % p == 0 else (1 if x % p in sq else -1)
    return [[chi(j - i) for j in range(p)] for i in range(p)]


def kron(A, B):
    return [[a * b for a in ra for b in rb] for ra in A for rb in B]


def square_array(N):
    """成分 ±1、行どうし・列どうしの内積が 0 の N×N の表と、作り方の名前。作れなければ (None, None)。"""
    if N == 1:
        return [[1]], "1"
    if N == 2:
        return [[1, 1], [1, -1]], "二進"
    if N % 4 != 0:
        return None, None
    if N & (N - 1) == 0:
        return [[-1 if bin(i & j).count("1") % 2 else 1 for j in range(N)] for i in range(N)], "二進"
    p = N - 1
    if is_prime(p) and p % 4 == 3:
        Q = jacobsthal(p)
        return [[1] * N] + [[-1] + [Q[i][j] + (1 if i == j else 0) for j in range(p)] for i in range(p)], "平方剰余 1 (p=%d)" % p
    if N % 2 == 0:
        p = N // 2 - 1
        if is_prime(p) and p % 4 == 1:
            Q = jacobsthal(p)
            n = p + 1
            C = [[0] + [1] * p] + [[1] + Q[i] for i in range(p)]
            I = lambda i, j: 1 if i == j else 0
            top = [[C[i][j] + I(i, j) for j in range(n)] + [C[i][j] - I(i, j) for j in range(n)] for i in range(n)]
            bot = [[C[i][j] - I(i, j) for j in range(n)] + [-C[i][j] - I(i, j) for j in range(n)] for i in range(n)]
            return top + bot, "平方剰余 2 (p=%d)" % p
    for a in range(2, N):
        if N % a == 0 and a <= N // a:
            A, na = square_array(a)
            B, nb = square_array(N // a)
            if A is not None and B is not None:
                return kron(A, B), "積 %d×%d [%s] × [%s]" % (a, N // a, na, nb)
    return None, None


def rowcol_zero_init(N):
    """縦横の和が 0 の N×N (±1) と、作り方の名前。条件を満たさないものは返さない。"""
    if N % 4 != 0 or N < 4:
        raise ValueError("N は 4 の倍数でなければならない: N = %d" % N)
    h = N // 2
    M = None
    if N == 4 or N % 8 == 0:                       # (甲)
        B, how = square_array(h)
        if B is not None:
            M = [[(1 if (i < h) == (j < h) else -1) * B[i % h][j % h] for j in range(N)] for i in range(N)]
            name = "甲: B と -B を縦横に並べる。B は %d×%d、%s" % (h, h, how)
    if M is None:                                  # (乙)
        H, how = square_array(N)
        if H is None:
            raise ValueError("N = %d の作り方は、このプログラムに入っていない" % N)
        H = [[v * r[0] for v in r] for r in H]     # 一列目を全部 +1 に。この列は使わない
        A = [r[h:] for r in H]                     # 後ろの N/2 本 (Claude が決めた選び方)
        M = [r + [-v for v in r] for r in A]
        name = "乙: %d×%d の表 (%s) の後ろ %d 列 A と -A を横に並べる" % (N, N, how, h)
    c = check(M)
    if not c["ok"]:
        raise AssertionError("N = %d: 作った配列が条件を満たさない: %s" % (N, c))
    return M, name


def check(M):
    """検算。整数だけで計算する。"""
    N = len(M)
    h = N // 2
    pm1 = all(len(r) == N for r in M) and all(v in (1, -1) for r in M for v in r)
    row_sums = [sum(r) for r in M]
    col_sums = [sum(r[j] for r in M) for j in range(N)]
    const_line = any(len(set(r)) == 1 for r in M) or any(len({r[j] for r in M}) == 1 for j in range(N))
    partner_cols = all(M[i][j + h] == -M[i][j] for i in range(N) for j in range(h))
    # Σ z_n² = (列 a の二乗の和 - 列 b の二乗の和) + 2i (列 a と列 b の内積)
    in_scope_max = 0      # 反転の相手どうしでない組での |実部|、|虚部| の最大
    in_scope = 0
    partner_vals = set()  # 反転の相手どうしの組での値 (範囲の外。参考)
    for a in range(N):
        for b in range(a + 1, N):
            re = sum(r[a] * r[a] - r[b] * r[b] for r in M)
            im = 2 * sum(r[a] * r[b] for r in M)
            if b == a + h:
                partner_vals.add((re, im))
            else:
                in_scope += 1
                in_scope_max = max(in_scope_max, abs(re), abs(im))
    partner_rows = sum(1 for i in range(N) for k in range(i + 1, N) if all(M[i][j] == -M[k][j] for j in range(N)))
    ok = pm1 and not any(row_sums) and not any(col_sums) and not const_line and partner_cols and in_scope_max == 0
    return dict(ok=ok, pm1=pm1, row_max=max(abs(v) for v in row_sums), col_max=max(abs(v) for v in col_sums),
                const_line=const_line, partner_cols=partner_cols, in_scope=in_scope, in_scope_max=in_scope_max,
                partner_vals=sorted(partner_vals), partner_rows=partner_rows)


def show(M):
    return "\n".join(" ".join("+1" if v > 0 else "-1" for v in r) for r in M)


def report(N, M, name):
    c = check(M)
    h = N // 2
    print("N = %d : %d×%d。成分はすべて ±1。作り方: %s" % (N, N, N, name))
    print(show(M))
    print("検算")
    print("  成分が ±1 だけ                         : %s" % ("はい" if c["pm1"] else "いいえ"))
    print("  行の和 (%d 本) の大きさの最大           : %d" % (N, c["row_max"]))
    print("  列の和 (%d 本) の大きさの最大           : %d" % (N, c["col_max"]))
    print("  全部同じ符号の行または列               : %s" % ("ある" if c["const_line"] else "無い"))
    print("  成分 j と j+%d が反転の相手どうし        : %s" % (h, "はい" if c["partner_cols"] else "いいえ"))
    print("  Σ z_n² (反転の相手どうしでない %d 組)   : 実部・虚部の大きさの最大 = %d" % (c["in_scope"], c["in_scope_max"]))
    print("  参考 (範囲の外): 反転の相手どうしの %d 組の Σ z_n² = %s" % (h, ", ".join("%d%+di" % v for v in c["partner_vals"])))
    print("  参考: 反転の相手を持つ行の組の数 = %d (%d 行中)" % (c["partner_rows"], N))
    print("  判定: %s" % ("合格" if c["ok"] else "不合格"))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    if len(sys.argv) > 1 and sys.argv[1] == "--export":
        out = os.path.join(here, "results", "rowcol_zero_init")
        os.makedirs(out, exist_ok=True)
        sums = []
        print("  N    行の和  列の和  同符号の行列  Σz²(範囲内の最大/組の数)  反転の相手を持つ行の組  判定   作り方")
        for N in range(4, 101, 4):
            try:
                M, name = rowcol_zero_init(N)
            except ValueError:
                print("  %3d  このプログラムでは作れない" % N)
                continue
            c = check(M)
            print("  %3d  %d       %d       %s          %d / %d                     %d                     %s  %s"
                  % (N, c["row_max"], c["col_max"], "ある" if c["const_line"] else "無い", c["in_scope_max"], c["in_scope"],
                     c["partner_rows"], "合格" if c["ok"] else "不合格", name))
            fn = "rowcol_zero_N%03d.txt" % N
            with open(os.path.join(out, fn), "w", encoding="utf-8", newline="\n") as f:
                f.write("# N = %d、%d×%d、作り方: %s、生成: rowcol_zero_init_pm1_v0.py\n" % (N, N, N, name))
                f.write(show(M) + "\n")
            sums.append((hashlib.sha256(open(os.path.join(out, fn), "rb").read()).hexdigest(), fn))
        with open(os.path.join(out, "SHA256SUMS.txt"), "w", encoding="utf-8", newline="\n") as f:
            for hsh, fn in sums:
                f.write("%s  %s\n" % (hsh, fn))
        print("%d 個のファイルを書き出した: results/rowcol_zero_init/" % len(sums))
        return
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    M, name = rowcol_zero_init(N)
    report(N, M, name)


if __name__ == "__main__":
    main()
