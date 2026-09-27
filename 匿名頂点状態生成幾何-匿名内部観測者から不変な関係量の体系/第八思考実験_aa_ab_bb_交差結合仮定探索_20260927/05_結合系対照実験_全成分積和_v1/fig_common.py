#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""図化の共通部。結合系の観測 row データの読み込みと、第六の行の形式への並べ替え。

第六の行の形式（15 列）
  step, P, E, U_re, U_im, H_re, H_im, Q, N, C, D, t_readout, x_readout, y_readout, q_index
結合系の row は状態 10 項目だけを持つ。読み出し 3 項目は、第六の readout と同じ式
  t = TAU0*log(Q),  x = P*H_re,  y = P*H_im
で計算する。q_index は、保存した行がすべて macrostep の境界なので 0。
"""
from __future__ import annotations
import importlib.util
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
from numba import njit

import unified_engine as eng

HERE = Path(__file__).resolve().parent
PROGRAMS = HERE / "paper6_figure_programs"
RELATIONS = ("aa", "ab", "bb")
PAPER6_COLUMNS = ["step", "P", "E", "U_re", "U_im", "H_re", "H_im", "Q", "N", "C", "D",
                  "t_readout", "x_readout", "y_readout", "q_index"]
STATE_ITEMS = ("P", "E", "U_re", "U_im", "H_re", "H_im", "Q", "N", "C", "D")


def case_id(relation: str) -> str:
    return f"unified_{relation}"


def load_program(filename: str):
    path = PROGRAMS / filename
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[path.stem] = module
    spec.loader.exec_module(module)
    return module


def result_dir(pattern_id: str) -> Path:
    return HERE / "results" / pattern_id


def figure_dir(pattern_id: str) -> Path:
    d = HERE / "figures" / pattern_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_definition(pattern_id: str) -> dict:
    return json.loads((result_dir(pattern_id) / "run_definition.json").read_text(encoding="utf-8"))


def relation_conditions(pattern_id: str):
    """(関係, ケース名, lambda_A, lambda_B, C, D) の一覧。値は実行時に補遺の初期化が返したもの。"""
    metas = {m["case_id"]: m for m in run_definition(pattern_id)["supplement_meta"]}
    return [(r, case_id(r), metas[r]["lambda_A"], metas[r]["lambda_B"], metas[r]["C"], metas[r]["D"]) for r in RELATIONS]


def read_unified_rows(pattern_id: str):
    folder = result_dir(pattern_id)
    manifest = json.loads((folder / "part_manifest.json").read_text(encoding="utf-8"))
    blocks, columns = [], None
    for p in manifest["parts"]:
        with h5py.File(folder / p["file"], "r") as hf:
            cols = json.loads(hf.attrs["columns_json"])
            columns = columns or cols
            blocks.append(hf["raw_macro_trajectory"][:])
    return np.concatenate(blocks), columns


@njit(cache=True)
def _readout(P, H_re, H_im, Q, t, x, y):
    for k in range(P.shape[0]):
        t[k] = eng.TAU0 * math.log(Q[k])
        x[k] = P[k] * H_re[k]
        y[k] = P[k] * H_im[k]


def paper6_rows(data, columns, relation: str):
    """結合系の row から、1 つの関係について第六の行の形式（15 列）を作る。"""
    n = data.shape[0]
    g = np.zeros((n, len(PAPER6_COLUMNS)), dtype=np.float64)
    g[:, 0] = data[:, columns.index("step")]
    for name in STATE_ITEMS:
        g[:, PAPER6_COLUMNS.index(name)] = data[:, columns.index(f"{relation}.{name}")]
    t, x, y = np.empty(n), np.empty(n), np.empty(n)
    _readout(np.ascontiguousarray(g[:, 1]), np.ascontiguousarray(g[:, 5]), np.ascontiguousarray(g[:, 6]),
             np.ascontiguousarray(g[:, 7]), t, x, y)
    g[:, 11], g[:, 12], g[:, 13] = t, x, y
    return g
