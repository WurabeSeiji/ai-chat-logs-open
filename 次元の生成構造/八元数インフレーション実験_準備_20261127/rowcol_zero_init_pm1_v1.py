#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v1 (2026-11-27): rowcol_zero_init_pm1_v0.py に、フラグを二つ足したもの (木原さんの指示)。
  (1) plus1  : 一行目・一列目が全部 +1 の、前回の初期値を作る。前任の対話の balanced_init_pm1_v2.py の
               square_array(N) が返す N×N そのもの (偏りがある。縦横の和は 0 にならない)。
               指定しなければ、v0 と同じ、縦横の和が 0 の初期値。
  (2) order  : 行の並べ替え。初期化のあとに施す。書き方と意味は前回のプログラム (balanced_init_pm1_v2.py の
               parse_order と square_array_ordered) と同じ。整数一つならその数だけ行をずらす
               (新しい i 行目 = もとの (i + k) mod N 行目)、"1,6,4,3,7,2,5,8" なら新しい i 行目 = もとの指定した行。
               今回の配列に合わせた並べ替えは作っていない。前回と同じ操作をそのまま掛けるだけである。
  どちらも指定しなければ v0 と同じ配列を返す。v0 の関数 (square_array、rowcol_zero_init、check) は一字も変えていない。
  使い方:
    python3 rowcol_zero_init_pm1_v1.py N [--plus1] [--order 指定]            初期値と確認を表示
    python3 rowcol_zero_init_pm1_v1.py --export [--plus1] [--order 指定]     4〜100 の 4 の倍数で作れるもの全部を
                                         results/rowcol_zero_init_v1/{種類}_{並べ替え}/ に書き出す
        (order に順番そのものを指定した場合は、その長さの N だけを書き出す)
以下は v0 の説明のまま。

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


def parse_order(N, spec):
    """shift (整数一つ) または並べ替え ("1,6,4,...") を、もとの行の番号 (0 始まり) の列に直す。前回のプログラムと同じ。"""
    if spec is None:
        return list(range(N))
    if isinstance(spec, int) or "," not in str(spec):
        k = int(spec) % N
        return [(i + k) % N for i in range(N)]
    order = [int(v) - 1 for v in str(spec).split(",")]
    if sorted(order) != list(range(N)):
        raise ValueError("並べ替えは 1 から %d までの数を一つずつ、%d 個書く: %s" % (N, N, spec))
    return order


def make_init(N, plus1=False, order=None):
    """初期値の N×N と、作り方の名前。
    plus1=False: 縦横の和が 0 の初期値 (v0 の rowcol_zero_init)。 plus1=True: 前回の初期値 (square_array(N) そのもの)。
    そのあとに、行の並べ替え order を施す (新しい i 行目 = もとの order[i] 行目)。"""
    if plus1:
        H, how = square_array(N)
        if H is None:
            raise ValueError("N = %d の作り方は、このプログラムに入っていない" % N)
        M, name = [r[:] for r in H], how
    else:
        M, name = rowcol_zero_init(N)
    od = parse_order(N, order)
    if od != list(range(N)):
        M = [M[od[i]][:] for i in range(N)]
        name += "、行の並べ替え " + ",".join(str(v + 1) for v in od)
    return M, name


def tag(N, plus1, order):
    kind = "plus1" if plus1 else "rowcolzero"
    if order is not None and "," in str(order):
        return kind + "_order" + "-".join(str(v + 1) for v in parse_order(N, order))
    return kind + "_shift%d" % (int(order) % N if order is not None else 0)


def facts(M):
    """配列についての事実 (判定ではない)。整数だけで計算する。"""
    N = len(M)
    row_sums = [sum(r) for r in M]
    col_sums = [sum(r[j] for r in M) for j in range(N)]
    nz = [(a + 1, b + 1) for a in range(N) for b in range(a + 1, N) if sum(r[a] * r[b] for r in M) != 0]
    sym = all(M[i][j] == M[j][i] for i in range(N) for j in range(N))
    return row_sums, col_sums, nz, sym


def report(N, M, name, plus1):
    row_sums, col_sums, nz, sym = facts(M)
    print("N = %d : %d×%d。成分はすべて ±1。作り方: %s" % (N, N, N, name))
    print(show(M))
    print("確認")
    print("  行の和: " + " ".join(str(v) for v in row_sums))
    print("  列の和: " + " ".join(str(v) for v in col_sums))
    print("  Σ z_n² が 0 にならない二成分の組: %d / %d%s" % (len(nz), N * (N - 1) // 2, ("  " + " ".join("(%d,%d)" % p for p in nz)) if 0 < len(nz) <= 60 else ""))
    print("  転置しても同じ配列か: %s" % ("同じ" if sym else "違う"))
    if plus1:
        print("  種類: 前回の初期値 (一行目・一列目が全部 +1 の表から)。縦横の和が 0 の検算の対象ではない")
    else:
        c = check(M)
        print("  種類: 縦横の和が 0 の初期値。検算 (±1 だけ、行と列の和が全部 0、全部同じ符号の行・列が無い、")
        print("        反転の相手どうしでない %d 組で Σ z_n² = 0): %s" % (c["in_scope"], "合格" if c["ok"] else "不合格"))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    args = sys.argv[1:]
    plus1 = "--plus1" in args
    order = args[args.index("--order") + 1] if "--order" in args else None
    rest = [a for k, a in enumerate(args) if a not in ("--plus1", "--order", "--export") and not (k > 0 and args[k - 1] == "--order")]
    if "--export" in args:
        only = len(str(order).split(",")) if (order is not None and "," in str(order)) else None
        sub = None
        sums = []
        for N in range(4, 101, 4):
            if only is not None and N != only:
                continue
            try:
                M, name = make_init(N, plus1, order)
            except ValueError:
                print("  %3d  このプログラムでは作れない" % N)
                continue
            if sub is None:
                sub = tag(N, plus1, order) if only is not None else tag(N, plus1, order).rsplit("_", 1)[0] + "_shift%s" % (order if order is not None else 0)
                out = os.path.join(here, "results", "rowcol_zero_init_v1", sub)
                os.makedirs(out, exist_ok=True)
            row_sums, col_sums, nz, sym = facts(M)
            verdict = "-" if plus1 else ("合格" if check(M)["ok"] else "不合格")
            print("  %3d  行の和の最大 %d  列の和の最大 %d  Σz²≠0 の組 %d/%d  転置で%s  検算 %s  %s"
                  % (N, max(abs(v) for v in row_sums), max(abs(v) for v in col_sums), len(nz), N * (N - 1) // 2,
                     "同じ" if sym else "違う", verdict, name))
            fn = "init_N%03d.txt" % N
            with open(os.path.join(out, fn), "w", encoding="utf-8", newline="\n") as f:
                f.write("# N = %d、%d×%d、作り方: %s、生成: rowcol_zero_init_pm1_v1.py\n" % (N, N, N, name))
                f.write(show(M) + "\n")
            sums.append((hashlib.sha256(open(os.path.join(out, fn), "rb").read()).hexdigest(), fn))
        with open(os.path.join(out, "SHA256SUMS.txt"), "w", encoding="utf-8", newline="\n") as f:
            for hsh, fn in sums:
                f.write("%s  %s\n" % (hsh, fn))
        print("%d 個のファイルを書き出した: results/rowcol_zero_init_v1/%s/" % (len(sums), sub))
        return
    N = int(rest[0]) if rest else 8
    M, name = make_init(N, plus1, order)
    report(N, M, name, plus1)


if __name__ == "__main__":
    main()
