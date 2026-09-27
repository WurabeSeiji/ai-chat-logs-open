#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""結合系の演算プログラム。

受け取るのは、状態 psi と、作用 S の定義（成分ごとの関数の種類と、その関数が読む成分の番号）だけ。
状態がいくつの関係からできているか、1 つの関係が何成分か、どこが交差成分かは知らない。
成分の数は psi の長さから取る。

1 回の microstep
  1. 更新前の psi で、S の全成分を評価する（値が 0 の成分も評価する）。
  2. 行ごとに、全成分について S[i,j] * psi[j] を計算して足す。
     積が NaN になった項は足さない。
  3. psi の全成分を同時に書き換える。
どの段階の処理をするかは、psi の中の位相成分が決める。このプログラムは段階を数えない。

S の成分の関数は、第六思考実験の transition(z) の候補値の式を、同じ書き方で書き直したもの。
"""
from __future__ import annotations
import math
import numpy as np
from numba import njit

# 第六思考実験と同じ定数
STEPS_PER_ORBIT = 4000
HSTEP = 2.0 * math.pi / STEPS_PER_ORBIT
PH_RE = math.cos(HSTEP)
PH_IM = math.sin(HSTEP)
TAU0 = 1.0e4
MICROSTEPS_PER_MACRO = 11

# S の成分の関数の種類
K_ZERO = 0        # 0
K_ONE = 1         # 1（位相の巡回）
K_IDENT = 2       # psi[a0]
K_RATE_DP = 3     # dP/dphi   (p=a0, n=a1, c=a2, d=a3)
K_RATE_DE = 4     # dE/dphi = 0
K_RATE_DT = 5     # dt/dphi   (p=a0, c=a1)
K_STAGE_HALF = 6  # psi[a0] + 0.5*HSTEP*psi[a1]
K_STAGE_FULL = 7  # psi[a0] + HSTEP*psi[a1]
K_COMBINE = 8     # (HSTEP/6)*(psi[a0] + 2 psi[a1] + 2 psi[a2] + psi[a3])
K_MID = 9         # psi[a0] + 0.5*psi[a1]
K_ADD = 10        # psi[a0] + psi[a1]
K_CONST_H = 11    # HSTEP
K_ROT_U_RE = 12   # psi[a0]*PH_RE - psi[a1]*PH_IM
K_ROT_U_IM = 13   # psi[a0]*PH_IM + psi[a1]*PH_RE
K_ROT_H_RE = 14   # psi[a0]*cos(psi[a2]) - psi[a1]*sin(psi[a2])
K_ROT_H_IM = 15   # psi[a0]*sin(psi[a2]) + psi[a1]*cos(psi[a2])
K_CLOCK = 16      # psi[a0]*exp(psi[a1]/TAU0)

KIND_NAMES = {
    K_ZERO: "zero", K_ONE: "one", K_IDENT: "ident", K_RATE_DP: "rate_dp", K_RATE_DE: "rate_de",
    K_RATE_DT: "rate_dt", K_STAGE_HALF: "stage_half", K_STAGE_FULL: "stage_full", K_COMBINE: "combine",
    K_MID: "mid", K_ADD: "add", K_CONST_H: "const_h", K_ROT_U_RE: "rot_u_re", K_ROT_U_IM: "rot_u_im",
    K_ROT_H_RE: "rot_h_re", K_ROT_H_IM: "rot_h_im", K_CLOCK: "clock",
}


@njit(cache=True)
def cell_value(kind, a0, a1, a2, a3, psi):
    """S の 1 成分の値を、更新前の psi から計算する。"""
    if kind == K_ZERO:
        return 0.0
    if kind == K_ONE:
        return 1.0
    if kind == K_IDENT:
        return psi[a0]
    if kind == K_RATE_DP:
        p = psi[a0]; n = psi[a1]; c = psi[a2]; d = psi[a3]
        rp = math.sqrt(p); rc = math.sqrt(c)
        return -(4.0/3.0)*n*d*rc/rp -(64.0/5.0)*n*c*rc/(p*rp)
    if kind == K_RATE_DE:
        return 0.0
    if kind == K_RATE_DT:
        p = psi[a0]; c = psi[a1]
        rp = math.sqrt(p); rc = math.sqrt(c)
        return p*rp/rc
    if kind == K_STAGE_HALF:
        return psi[a0]+0.5*HSTEP*psi[a1]
    if kind == K_STAGE_FULL:
        return psi[a0]+HSTEP*psi[a1]
    if kind == K_COMBINE:
        return (HSTEP/6.0)*(psi[a0]+2.0*psi[a1]+2.0*psi[a2]+psi[a3])
    if kind == K_MID:
        return psi[a0]+0.5*psi[a1]
    if kind == K_ADD:
        return psi[a0]+psi[a1]
    if kind == K_CONST_H:
        return HSTEP
    if kind == K_ROT_U_RE:
        return psi[a0]*PH_RE-psi[a1]*PH_IM
    if kind == K_ROT_U_IM:
        return psi[a0]*PH_IM+psi[a1]*PH_RE
    if kind == K_ROT_H_RE:
        cd = math.cos(psi[a2]); sd = math.sin(psi[a2])
        return psi[a0]*cd-psi[a1]*sd
    if kind == K_ROT_H_IM:
        cd = math.cos(psi[a2]); sd = math.sin(psi[a2])
        return psi[a0]*sd+psi[a1]*cd
    if kind == K_CLOCK:
        return psi[a0]*math.exp(psi[a1]/TAU0)
    raise ValueError("unknown kind")


@njit(cache=True)
def build_S(psi, kind, arg, S):
    """S の全成分を評価する。値が 0 の成分も、ほかの成分と同じ手順で評価する。"""
    n = psi.shape[0]
    for i in range(n):
        for j in range(n):
            S[i, j] = cell_value(kind[i, j], arg[i, j, 0], arg[i, j, 1], arg[i, j, 2], arg[i, j, 3], psi)


@njit(cache=True)
def product_sum(S, psi, psi_next):
    """全成分の積和。積が NaN になった項は足さない。足さなかった項の数を返す。"""
    n = psi.shape[0]
    skipped = 0
    for i in range(n):
        acc = 0.0
        for j in range(n):
            term = S[i, j] * psi[j]
            if math.isnan(term):
                skipped += 1
            else:
                acc += term
        psi_next[i] = acc
    return skipped


@njit(cache=True)
def microstep(psi, psi_next, S, kind, arg):
    build_S(psi, kind, arg, S)
    skipped = product_sum(S, psi, psi_next)
    for i in range(psi.shape[0]):
        psi[i] = psi_next[i]
    return skipped


@njit(cache=True)
def generate_chunk(psi, psi_next, S, kind, arg, save_index, start_step, max_rows, last_step):
    """macrostep ごとに、指定された成分を 1 行として保存しながら進める。

    保存する成分の番号は save_index で外から受け取る。どの成分が何を表すかは知らない。
    macrostep は 1 つも省かずに保存する。microstep の途中の状態は保存しない。
    行の保存、終了の判定、macrostep の順序は、第六思考実験の generate_chunk と同じ。
    終了は、外から与えた macrostep の番号 last_step で決める。状態の値では止めない。
    戻り値: 行の配列、行数、足さなかった NaN の項の数、最初に NaN の項が出た macrostep（なければ -1）、終了したか
    """
    m_save = save_index.shape[0]
    rows = np.empty((max_rows, m_save + 1), dtype=np.float64)
    count = 0
    skipped_total = 0
    first_skip_step = -1
    finished = False
    for k in range(max_rows):
        step = start_step + k
        rows[count, 0] = step
        for i in range(m_save):
            rows[count, 1 + i] = psi[save_index[i]]
        count += 1
        if step >= last_step:
            finished = True
            break
        for m in range(MICROSTEPS_PER_MACRO):
            s = microstep(psi, psi_next, S, kind, arg)
            if s > 0:
                if first_skip_step < 0:
                    first_skip_step = step
                skipped_total += s
    return rows, count, skipped_total, first_skip_step, finished
