#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ラッパー共通部。第六・補遺のプログラムの読み込みと、実行条件の表だけを持つ。

力学、初期状態の作り方、保存形式、図の描き方は、すべて paper6_programs/ の
プログラムが行う。ここでは何も計算しない。
"""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROGRAMS = HERE / "paper6_programs"

# 実行条件。a の電荷は +n*q0、b の電荷は -n*q0。
# aa は a を相互作用の両脚に、bb は b を両脚に入れる。
# 補遺の初期化プログラムの CASES と同じ形式 (ケース名, 電荷倍数 n, 符号 s_A, 符号 s_B)。
SELF_CASES = (
    ("self_aa_n1", 1, +1, +1),
    ("self_bb_n1", 1, -1, -1),
    ("self_aa_n3", 3, +1, +1),
    ("self_bb_n3", 3, -1, -1),
)

_loaded = {}


def load_program(filename: str):
    """paper6_programs/ のプログラムを、書き換えずにモジュールとして読み込む。同じファイルは 1 回だけ読む。"""
    if filename in _loaded:
        return _loaded[filename]
    path = PROGRAMS / filename
    name = path.stem
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    _loaded[filename] = module
    return module


def load_supplement():
    """補遺の初期化プログラムを読み込み、条件の表 CASES だけを今回の 4 条件に替える。"""
    supplement = load_program("paper6_init_G_q0_c1.py")
    supplement.CASES = SELF_CASES
    return supplement


def case_table():
    """補遺の init_case が返す値から、条件の表を作る。

    戻り値の各行は (ケース名, n, s_A, s_B, lambda_A, lambda_B, C, D)。
    """
    supplement = load_supplement()
    rows = []
    for case_id, n, sign_a, sign_b in SELF_CASES:
        _z, meta = supplement.init_case(case_id)
        rows.append((case_id, n, sign_a, sign_b, meta["lambda_A"], meta["lambda_B"], meta["C"], meta["D"]))
    return rows
