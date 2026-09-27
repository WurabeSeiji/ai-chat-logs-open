#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paper 6 supplement: initializer under G=q0=c=1.

Only the initialization map is changed. The transition map is intentionally not
included or executed here.
"""
from __future__ import annotations
import numpy as np
from fractions import Fraction

# State layout copied exactly from Paper 6.
UR, UI, P, E, HR, HI, Q, N, C, D = range(10)
K1 = 10; K2 = 13; K3 = 16; K4 = 19
PS = 22; ES = 23; DV = 24; PM = 27; EM = 28; DPHI = 29; QB = 30
NST = 41

# New normalization.
G = 1.0
Q0 = 1.0
C_LIGHT = 1.0

# In Paper 6, the n=1 cases used |lambda|=0.30 with equal masses mA=mB=M/2.
# Under q0=1 and qA=+/- n*q0, lambda = q/m = +/- 2*n*q0/M.
# Matching the n=1 Paper 6 value gives 2*q0/M = 0.30 -> M/q0 = 20/3.
M_OVER_Q0_EXACT = Fraction(20, 3)
M_OVER_Q0 = float(M_OVER_Q0_EXACT)
M_TOTAL = M_OVER_Q0 * Q0
M_A = M_TOTAL / 2.0
M_B = M_TOTAL / 2.0
N_SYMMETRIC_MASS_RATIO = 0.25

CASES = (
    ("n1_attractive", 1, +1, -1),
    ("n1_repulsive",  1, +1, +1),
    ("n3_attractive", 3, +1, -1),
    ("n3_repulsive",  3, +1, +1),
    ("n4_attractive", 4, +1, -1),
)


def charge_to_mass_from_multiple(n: int, sign: int, mass: float) -> float:
    q = sign * n * Q0
    return q / mass


def case_relations(n: int, sign_a: int, sign_b: int):
    # Use exact rational arithmetic for initialization, then convert once to float.
    # q0=1, mA=mB=M/2 and M/q0=20/3 give lambda=+/-3n/10 exactly.
    lambda_a_q = Fraction(sign_a * 3 * n, 10)
    lambda_b_q = Fraction(sign_b * 3 * n, 10)
    c0_q = Fraction(1, 1) - lambda_a_q * lambda_b_q
    d0_q = (lambda_a_q - lambda_b_q) ** 2
    return float(lambda_a_q), float(lambda_b_q), float(c0_q), float(d0_q)


def init_state_from_relations(p0: float, n0: float, c0: float, d0: float) -> np.ndarray:
    """Exact Paper 6 initializer, with c0,d0 supplied by the new normalization."""
    z = np.zeros(NST, dtype=np.float64)
    z[UR] = 1.0; z[P] = p0; z[HR] = 1.0; z[Q] = 1.0
    z[N] = n0; z[C] = c0; z[D] = d0
    z[PS] = p0; z[PM] = p0; z[QB] = 1.0
    return z


def init_case(case_id: str, p0: float = 50.0) -> tuple[np.ndarray, dict]:
    rec = next(x for x in CASES if x[0] == case_id)
    _, n, sign_a, sign_b = rec
    lambda_a, lambda_b, c0, d0 = case_relations(n, sign_a, sign_b)
    z = init_state_from_relations(p0, N_SYMMETRIC_MASS_RATIO, c0, d0)
    meta = {
        "case_id": case_id,
        "charge_multiple_n": n,
        "q0": Q0,
        "M_total": M_TOTAL,
        "m_A": M_A,
        "m_B": M_B,
        "lambda_A": lambda_a,
        "lambda_B": lambda_b,
        "N": N_SYMMETRIC_MASS_RATIO,
        "C": c0,
        "D": d0,
    }
    return z, meta
