#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
初期値を一つだけファイルに書き出す。2026-11-27。
rowcol_zero_init_pm1_v1.py の --export は N = 4〜100 しか書き出さないので、それより大きい N のために作った。
初期値そのものは rowcol_zero_init_pm1_v1.py の make_init をそのまま呼ぶ (フラグなし = 縦横の和が 0)。
ファイルの形式と置き場所は --export と同じ: results/rowcol_zero_init_v1/rowcolzero_shift0/init_N{N:03d}.txt
書き出す前に、同じプログラムの検算 (check) に通ることを確かめる。
使い方: python3 export_init_single_v0.py N
"""
import hashlib
import os
import sys

from rowcol_zero_init_pm1_v1 import make_init, check, show, tag

N = int(sys.argv[1])
M, name = make_init(N)
c = check(M)
assert c["ok"], "検算に通らない: %s" % c
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "rowcol_zero_init_v1", tag(N, False, None))
os.makedirs(out, exist_ok=True)
fn = os.path.join(out, "init_N%03d.txt" % N)
with open(fn, "w", encoding="utf-8", newline="\n") as f:
    f.write("# N = %d、%d×%d、作り方: %s、生成: rowcol_zero_init_pm1_v1.py\n" % (N, N, N, name))
    f.write(show(M) + "\n")
print("N = %d。作り方: %s" % (N, name))
print("検算: 行の和の最大 %d、列の和の最大 %d、全部同じ符号の行・列 %s、Σz² (反転の相手どうしでない %d 組) の最大 %d、判定 合格"
      % (c["row_max"], c["col_max"], "ある" if c["const_line"] else "無い", c["in_scope"], c["in_scope_max"]))
print("書き出し: results/rowcol_zero_init_v1/%s/init_N%03d.txt  SHA256 %s" % (tag(N, False, None), N, hashlib.sha256(open(fn, "rb").read()).hexdigest()))
