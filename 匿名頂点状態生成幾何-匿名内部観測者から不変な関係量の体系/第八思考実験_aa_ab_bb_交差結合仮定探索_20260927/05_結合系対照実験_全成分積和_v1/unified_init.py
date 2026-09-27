#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""初期化と、作用 S の定義。

関係 aa, ab, bb と、1 つの関係の成分の並びを知っているのは、このファイルだけである。
演算プログラム unified_engine.py には、ここで作った psi と S の定義だけを渡す。

初期状態
  第六思考実験 補遺の初期化プログラムを変更せずに読み込み、条件の表 CASES だけを
  aa, ab, bb の 3 条件に替えて init_case を呼ぶ。返ってきた 3 つの状態をつないで psi にする。

S の定義
  S = [[S_aa, 0, 0], [0, S_ab, 0], [0, 0, S_bb]]
  S の全成分に、関数の種類 kind と、その関数が読む成分の番号 arg を置く。
  対角ブロックには、第六思考実験の式  z_i' = sum_j q_j F_ij(z)  の F_ij を置く（掛ける相手は位相 q_j）。
  交差成分には、値 0 の成分を置く。
"""
from __future__ import annotations
import importlib.util
import sys
from pathlib import Path
import numpy as np

import unified_engine as eng

HERE = Path(__file__).resolve().parent
SUPPLEMENT_FILE = HERE / "paper6_supplement_program" / "paper6_init_G_q0_c1.py"

RELATION_NAMES = ("aa", "ab", "bb")

# macrostep ごとに保存する項目。軌道の読み出しに必要な永続状態で、並びは第六思考実験の行と同じ。
SAVED_LOCAL_NAMES = ("P", "E", "U_re", "U_im", "H_re", "H_im", "Q", "N", "C", "D")


def relations(pattern):
    """実行パターンから、3 つの関係の条件を作る。

    a の電荷は sign_a*n*q0、b の電荷は sign_b*n*q0。
    戻り値の各行は (関係の名前, 電荷倍数 n, 符号 s_A, 符号 s_B)。
    """
    n = int(pattern["charge_multiple_n"])
    sa, sb = int(pattern["sign_a"]), int(pattern["sign_b"])
    return (("aa", n, sa, sa), ("ab", n, sa, sb), ("bb", n, sb, sb))

LOCAL_NAMES = (["U_re", "U_im", "P", "E", "H_re", "H_im", "Q", "N", "C", "D"]
               + [f"k{s}_{c}" for s in (1, 2, 3, 4) for c in ("P", "E", "t")]
               + ["P_stage", "E_stage", "d_P", "d_E", "d_t", "P_mid", "E_mid", "dphi"]
               + [f"q{j}" for j in range(11)])


def load_supplement(cases=None):
    spec = importlib.util.spec_from_file_location("paper6_init_G_q0_c1", SUPPLEMENT_FILE)
    module = importlib.util.module_from_spec(spec)
    sys.modules["paper6_init_G_q0_c1"] = module
    spec.loader.exec_module(module)
    if cases is not None:
        module.CASES = cases
    return module


def initial_state(pattern):
    """補遺の init_case で 3 つの関係の初期状態を作り、つないで返す。"""
    cases = relations(pattern)
    sup = load_supplement(cases)
    blocks, metas = [], []
    for name, n, sa, sb in cases:
        z, meta = sup.init_case(name)
        blocks.append(z)
        metas.append(meta)
    psi = np.concatenate(blocks).astype(np.float64)
    names = [f"{name}.{local}" for name in RELATION_NAMES for local in LOCAL_NAMES]
    if len(names) != psi.shape[0]:
        raise RuntimeError("component name count mismatch")
    normalization = {"G": sup.G, "q0": sup.Q0, "c": sup.C_LIGHT, "M_over_q0": str(sup.M_OVER_Q0_EXACT)}
    return psi, names, metas, normalization


def relation_candidates(L):
    """1 つの関係について、(行, 段階) ごとの F_ij の種類と、読む成分（関係の中での番号）を返す。

    第六思考実験 transition(z) の cand[i, j] と同じ内容。
    """
    cand = {}
    for i in range(30):
        for j in range(11):
            cand[(i, j)] = (eng.K_IDENT, (i,))

    def rates(krow, src_p, phase):
        cand[(krow, phase)] = (eng.K_RATE_DP, (src_p, L.N, L.C, L.D))
        cand[(krow + 1, phase)] = (eng.K_RATE_DE, ())
        cand[(krow + 2, phase)] = (eng.K_RATE_DT, (src_p, L.C))
        for a in range(3):
            cand[(krow + a, 10)] = (eng.K_ZERO, ())

    rates(L.K1, L.P, 0)
    cand[(L.PS, 1)] = (eng.K_STAGE_HALF, (L.P, L.K1))
    cand[(L.ES, 1)] = (eng.K_STAGE_HALF, (L.E, L.K1 + 1))
    rates(L.K2, L.PS, 2)
    cand[(L.PS, 3)] = (eng.K_STAGE_HALF, (L.P, L.K2))
    cand[(L.ES, 3)] = (eng.K_STAGE_HALF, (L.E, L.K2 + 1))
    rates(L.K3, L.PS, 4)
    cand[(L.PS, 5)] = (eng.K_STAGE_FULL, (L.P, L.K3))
    cand[(L.ES, 5)] = (eng.K_STAGE_FULL, (L.E, L.K3 + 1))
    rates(L.K4, L.PS, 6)
    for a in range(3):
        cand[(L.DV + a, 7)] = (eng.K_COMBINE, (L.K1 + a, L.K2 + a, L.K3 + a, L.K4 + a))
        cand[(L.DV + a, 10)] = (eng.K_ZERO, ())
    cand[(L.PM, 8)] = (eng.K_MID, (L.P, L.DV))
    cand[(L.EM, 8)] = (eng.K_MID, (L.E, L.DV + 1))
    cand[(L.PM, 10)] = (eng.K_ADD, (L.P, L.DV))
    cand[(L.EM, 10)] = (eng.K_ADD, (L.E, L.DV + 1))
    cand[(L.DPHI, 9)] = (eng.K_CONST_H, ())
    cand[(L.DPHI, 10)] = (eng.K_ZERO, ())
    cand[(L.UR, 10)] = (eng.K_ROT_U_RE, (L.UR, L.UI))
    cand[(L.UI, 10)] = (eng.K_ROT_U_IM, (L.UR, L.UI))
    cand[(L.P, 10)] = (eng.K_ADD, (L.P, L.DV))
    cand[(L.E, 10)] = (eng.K_ADD, (L.E, L.DV + 1))
    cand[(L.HR, 10)] = (eng.K_ROT_H_RE, (L.HR, L.HI, L.DPHI))
    cand[(L.HI, 10)] = (eng.K_ROT_H_IM, (L.HR, L.HI, L.DPHI))
    cand[(L.Q, 10)] = (eng.K_CLOCK, (L.Q, L.DV + 2))
    cand[(L.PS, 10)] = (eng.K_ADD, (L.P, L.DV))
    cand[(L.ES, 10)] = (eng.K_ADD, (L.E, L.DV + 1))
    return cand


def action_definition(total):
    """S の全成分の定義（kind と arg）を作る。交差成分は値 0 の成分になる。"""
    L = load_supplement()
    size = L.NST
    kind = np.full((total, total), eng.K_ZERO, dtype=np.int64)
    arg = np.zeros((total, total, 4), dtype=np.int64)
    cand = relation_candidates(L)
    for r in range(len(RELATION_NAMES)):
        base = r * size
        for (i, j), (k, local_args) in cand.items():
            row, col = base + i, base + L.QB + j
            kind[row, col] = k
            for t, a in enumerate(local_args):
                arg[row, col, t] = base + a
        # 位相の巡回  q0' = q10,  qk' = q(k-1)
        kind[base + L.QB, base + L.QB + 10] = eng.K_ONE
        for j in range(1, 11):
            kind[base + L.QB + j, base + L.QB + j - 1] = eng.K_ONE
    return kind, arg


def saved_components(names):
    """macrostep ごとに保存する成分の番号と名前。"""
    wanted = [f"{r}.{local}" for r in RELATION_NAMES for local in SAVED_LOCAL_NAMES]
    index = np.array([names.index(w) for w in wanted], dtype=np.int64)
    return index, wanted


def describe(kind):
    """S の定義の内訳（種類ごとの成分の数）。"""
    values, counts = np.unique(kind, return_counts=True)
    return {eng.KIND_NAMES[int(v)]: int(c) for v, c in zip(values, counts)}
