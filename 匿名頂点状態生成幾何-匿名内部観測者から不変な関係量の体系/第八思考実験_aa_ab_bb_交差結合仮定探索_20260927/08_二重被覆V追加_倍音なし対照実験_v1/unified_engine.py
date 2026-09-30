#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 8 experiment 2 engine: original unified engine + minimal double-cover operations.

The original interaction functions are unchanged. Two new stateless cell functions
rotate a persistent complex cover state V by half of the same phase increment used
for U. No harmonic powers V^3,V^5,... are generated in this control experiment.
"""
from __future__ import annotations
import math
import numpy as np
from numba import njit

STEPS_PER_ORBIT = 4000
HSTEP = 2.0 * math.pi / STEPS_PER_ORBIT
PH_RE = math.cos(HSTEP)
PH_IM = math.sin(HSTEP)
PH_HALF_RE = math.cos(0.5 * HSTEP)
PH_HALF_IM = math.sin(0.5 * HSTEP)
TAU0 = 1.0e4
MICROSTEPS_PER_MACRO = 11

K_ZERO = 0
K_ONE = 1
K_IDENT = 2
K_RATE_DP = 3
K_RATE_DE = 4
K_RATE_DT = 5
K_STAGE_HALF = 6
K_STAGE_FULL = 7
K_COMBINE = 8
K_MID = 9
K_ADD = 10
K_CONST_H = 11
K_ROT_U_RE = 12
K_ROT_U_IM = 13
K_ROT_H_RE = 14
K_ROT_H_IM = 15
K_CLOCK = 16
K_ROT_V_RE = 17
K_ROT_V_IM = 18

KIND_NAMES = {
    K_ZERO: "zero", K_ONE: "one", K_IDENT: "ident", K_RATE_DP: "rate_dp", K_RATE_DE: "rate_de",
    K_RATE_DT: "rate_dt", K_STAGE_HALF: "stage_half", K_STAGE_FULL: "stage_full", K_COMBINE: "combine",
    K_MID: "mid", K_ADD: "add", K_CONST_H: "const_h", K_ROT_U_RE: "rot_u_re", K_ROT_U_IM: "rot_u_im",
    K_ROT_H_RE: "rot_h_re", K_ROT_H_IM: "rot_h_im", K_CLOCK: "clock",
    K_ROT_V_RE: "rot_v_half_re", K_ROT_V_IM: "rot_v_half_im",
}

@njit(cache=True)
def cell_value(kind, a0, a1, a2, a3, psi):
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
    if kind == K_ROT_V_RE:
        return psi[a0]*PH_HALF_RE-psi[a1]*PH_HALF_IM
    if kind == K_ROT_V_IM:
        return psi[a0]*PH_HALF_IM+psi[a1]*PH_HALF_RE
    raise ValueError("unknown kind")

@njit(cache=True)
def build_S(psi, kind, arg, S):
    n = psi.shape[0]
    for i in range(n):
        for j in range(n):
            S[i, j] = cell_value(kind[i, j], arg[i, j, 0], arg[i, j, 1], arg[i, j, 2], arg[i, j, 3], psi)

@njit(cache=True)
def product_sum(S, psi, psi_next):
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
