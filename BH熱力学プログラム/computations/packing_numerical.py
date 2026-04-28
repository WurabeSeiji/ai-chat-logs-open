#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lattice-point packing count N(k) for the four-dimensional ball.

Problem:
    Count integer points (x1,x2,x3,x4) in Z^4 satisfying
        sum_i (|x_i| + 1/2)^2 <= R^2,   R = 2k+1.

Strategy:
    1. Restrict to ordered non-negative tuples (y1>=y2>=y3>=y4>=0).
    2. For each such tuple, multiply by:
       - sign multiplicity: 2^(#{i: y_i > 0})
       - permutation multiplicity: 4! / prod(factorial of group sizes)
    Reduces 4D enumeration from O(R^4) cells to O(R^4 / 384) ordered tuples.
    Trivially fast up through k = 30.
"""

from collections import Counter
from math import factorial
from fractions import Fraction


def N_count(k: int) -> int:
    """Exact count N(k) of integer 4-tuples satisfying the packing condition."""
    R = 2 * k + 1
    # Use integer arithmetic by multiplying through by 4:
    # condition (|x|+1/2)^2 + ... <= R^2  <=>  sum (2|x|+1)^2 <= (2R)^2
    target = (2 * R) ** 2  # = (4k+2)^2
    bound = R  # upper bound on |x_i| (a safe over-estimate)
    total = 0
    for y1 in range(0, bound + 1):
        s1 = (2 * y1 + 1) ** 2
        if s1 > target:
            break
        for y2 in range(0, y1 + 1):
            s2 = s1 + (2 * y2 + 1) ** 2
            if s2 > target:
                break
            for y3 in range(0, y2 + 1):
                s3 = s2 + (2 * y3 + 1) ** 2
                if s3 > target:
                    break
                for y4 in range(0, y3 + 1):
                    s4 = s3 + (2 * y4 + 1) ** 2
                    if s4 > target:
                        break
                    # sign multiplicity
                    nz = sum(1 for y in (y1, y2, y3, y4) if y > 0)
                    sign_mult = 1 << nz  # 2^nz
                    # permutation multiplicity
                    counts = Counter((y1, y2, y3, y4))
                    perm_mult = 24
                    for c in counts.values():
                        perm_mult //= factorial(c)
                    total += sign_mult * perm_mult
    return total


def V4(R) -> Fraction:
    """4-ball volume V_4(R) = (pi^2 / 2) R^4 ; returns symbolic factor (R^4)/2."""
    # We return numerical via float; symbolic factor preserved separately.
    return Fraction(R ** 4, 2)  # multiply by pi^2 externally


def report(k_max: int = 20, file_out: str = None):
    import math
    pi2 = math.pi ** 2
    rows = []
    rows.append(("k", "R", "N(k)", "a_k=N(k)-N(k-1)", "V4(R)", "rho(k)", "Δ(R)", "c=Δ/(2π²R³)"))
    prev = 0
    for k in range(0, k_max + 1):
        n = N_count(k)
        R = 2 * k + 1
        v4 = (pi2 / 2.0) * R ** 4
        delta = v4 - n
        rho = n / v4
        c = delta / (2 * pi2 * R ** 3) if R > 0 else float("nan")
        a = n - prev
        prev = n
        rows.append((k, R, n, a, f"{v4:.4f}", f"{rho:.6f}", f"{delta:.4f}", f"{c:.6f}"))
    out = "\n".join("\t".join(str(c) for c in row) for row in rows)
    print(out)
    if file_out:
        with open(file_out, "w", encoding="utf-8") as f:
            f.write(out + "\n")
    return rows


if __name__ == "__main__":
    import sys
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    out_file = sys.argv[2] if len(sys.argv) > 2 else None
    report(kmax, out_file)
