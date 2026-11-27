#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成した初期データが、重心・モーメント・等方性を保っているかを検算する。生成プログラムとは別のプログラム。
生成プログラム (balanced_init_pm1_v3.py) の関数は一つも使わない。書き出されたデータのファイルだけを読む。

データ: results/balanced_data/balanced_N*.txt (export_balanced_data_v0.py が書き出したもの。行の入れ替えなし)。
        行 = 元 (n = 1..N)、列 = 成分。x[n][a] = n 番目の元の a 番目の成分。

検算する項目 (定義どおりに、整数だけで計算する。浮動小数・乱数なし):
  [値]       成分がすべて +1 か -1 であること。
  [重心]     元を全部足した合計が、どの成分についても 0 であること。      Σ_n x[n][a] = 0  (全部の a)
  [モーメント] ある向き u についての二次のモーメントを M(u) = Σ_n (x[n]·u)² とする (x[n]·u は成分ごとの積の和)。
             成分の軸の向き (u = a 番目の成分だけ 1) について、M(u) がどの成分でも同じ値であること。
  [等方性]   どの向きでも二次のモーメントが同じであること。 M(u) / (u の長さの二乗) が、向きに依らず一定。
             調べる向き: 成分の軸 (C 通り)、二つの成分の和の向き u = e_a + e_b と差の向き u = e_a - e_b (全部の組)、
             および全成分が 1 の向き、成分の番号の偶奇で符号を変えた向き。
             (二次のモーメントは向きについての二次式なので、軸と、全部の組の和・差の向きで値が同じなら、あらゆる向きで同じである。)
  [Σz² = 0]  全部の二成分の組 (a, b) について、z_n = x[n][a] + i x[n][b] として Σ_n z_n² = 0 (実部も虚部も 0)。
  [SHA256]   データのファイルが、書き出したときのものと同じであること。

出力: 標準出力。ファイルごとに各項目の合否と数値、最後に全体のまとめ。一つでも不合格があれば終了コード 1。
使い方: python3 verify_balanced_data_v0.py [ファイル ...]     (省くと results/balanced_data/ の全部)
"""
import glob
import hashlib
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "results", "balanced_data")


def read_table(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            rows.append([int(v) for v in line.split()])
    return rows


def moment(X, u):
    """M(u) = Σ_n (x[n]·u)²"""
    return sum(sum(x * w for x, w in zip(r, u)) ** 2 for r in X)


def check(path, recorded):
    name = os.path.basename(path)
    X = read_table(path)
    R, C = len(X), len(X[0])
    res = {}
    res["形"] = (all(len(r) == C for r in X), "元 %d 個 × 成分 %d 個" % (R, C))
    res["値"] = (all(v in (1, -1) for r in X for v in r), "成分は +1 と -1 だけ")
    cs = [sum(X[n][a] for n in range(R)) for a in range(C)]
    res["重心"] = (all(v == 0 for v in cs), "成分ごとの合計の大きさの最大 = %d" % max(abs(v) for v in cs))
    axis = [moment(X, [1 if k == a else 0 for k in range(C)]) for a in range(C)]
    res["モーメント"] = (len(set(axis)) == 1, "成分の軸についての二次のモーメント = %s" % sorted(set(axis)))
    vals = {Fraction(v, 1) for v in axis}
    ndir = C
    for a in range(C):
        for b in range(a + 1, C):
            for s in (1, -1):
                u = [0] * C
                u[a], u[b] = 1, s
                vals.add(Fraction(moment(X, u), 2))
                ndir += 1
    for u in ([1] * C, [1 if k % 2 == 0 else -1 for k in range(C)]):
        vals.add(Fraction(moment(X, u), C))
        ndir += 1
    res["等方性"] = (len(vals) == 1, "%d 通りの向きで、二次のモーメント ÷ 向きの長さの二乗 = %s" % (ndir, sorted(str(v) for v in vals)))
    worst_re = worst_im = 0
    for a in range(C):
        for b in range(a + 1, C):
            worst_re = max(worst_re, abs(sum(X[n][a] ** 2 - X[n][b] ** 2 for n in range(R))))
            worst_im = max(worst_im, abs(2 * sum(X[n][a] * X[n][b] for n in range(R))))
    res["Σz² = 0"] = (worst_re == 0 and worst_im == 0, "全 %d 組で、実部の大きさの最大 = %d、虚部の大きさの最大 = %d" % (C * (C - 1) // 2, worst_re, worst_im))
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    res["SHA256"] = (recorded.get(name) == h, h[:16] + "… " + ("記録と一致" if recorded.get(name) == h else "記録と不一致、または記録なし"))
    return name, res


def main():
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(DATA, "balanced_N*.txt")))
    recorded = {}
    sp = os.path.join(DATA, "SHA256SUMS.txt")
    if os.path.exists(sp):
        for line in open(sp, encoding="utf-8"):
            h, n = line.split()
            recorded[n] = h
    if not files:
        sys.exit("データのファイルが無い。先に export_balanced_data_v0.py を実行する")
    bad = 0
    keys = ["形", "値", "重心", "モーメント", "等方性", "Σz² = 0", "SHA256"]
    summary = []
    for path in files:
        name, res = check(path, recorded)
        print(name)
        for k in keys:
            okk, msg = res[k]
            print("   [%s] %s   %s" % ("合格" if okk else "不合格", k, msg))
            bad += (not okk)
        summary.append((name, all(res[k][0] for k in keys)))
    print("=" * 60)
    print("まとめ: %d 個のファイルのうち、全項目に合格 %d 個、不合格を含む %d 個"
          % (len(summary), sum(1 for _, s in summary if s), sum(1 for _, s in summary if not s)))
    for name, s in summary:
        if not s:
            print("   不合格を含む: " + name)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
