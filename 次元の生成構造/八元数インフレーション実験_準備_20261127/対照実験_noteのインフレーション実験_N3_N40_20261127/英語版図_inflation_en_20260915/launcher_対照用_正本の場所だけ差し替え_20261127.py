#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照実験用の起動役 (2026-11-27)。regenerate_inflation_figs_from_canonical_20260915.py は書き換えない。
原本は正本データを Google Drive の絶対パス (CANON) から読む。対照実験では再生成したデータを読ませたいので、
原本のソースを読み込み、CANON を決める代入の直後に 1 行だけ足して、この対照フォルダ内の
make_parent型初期値_自己無撞着構造_20260907 を指すようにしてから実行する。変更はその 1 行だけで、追加行数 = 1 を確認する。"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'regenerate_inflation_figs_from_canonical_20260915.py')
src = open(SRC, encoding='utf-8').read()
anchor = "                     'make_parent型初期値_自己無撞着構造_20260907')\n"
assert src.count(anchor) == 1
local = os.path.abspath(os.path.join(HERE, '..', 'make_parent型初期値_自己無撞着構造_20260907'))
new = src.replace(anchor, anchor + f"CANON = {local!r}\n")
assert len(new.splitlines()) == len(src.splitlines()) + 1
exec(compile(new, SRC, 'exec'), {'__name__': '__main__', '__file__': SRC})
