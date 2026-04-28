#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verify the Lagrange–Jacobi connection between N(k) and the Jacobi
four-square representation function r_4(N).

Goal:
    Express N(k) in closed form using r_4(N) for N ≡ 4 (mod 8).
    Determine the precise normalization factor.
"""

from packing_numerical import N_count


def divisor_sum(n: int) -> int:
    """Sum of positive divisors of n."""
    if n <= 0:
        return 0
    s = 0
    d = 1
    while d * d <= n:
        if n % d == 0:
            s += d
            if d != n // d:
                s += n // d
        d += 1
    return s


def odd_part(n: int) -> int:
    """Largest odd divisor of n."""
    while n % 2 == 0 and n > 0:
        n //= 2
    return n


def r4(N: int) -> int:
    """Jacobi four-square representation count r_4(N)."""
    if N == 0:
        return 1
    if N % 4 != 0:
        return 8 * divisor_sum(N)
    else:
        return 24 * divisor_sum(odd_part(N))


def r4_odd(N: int) -> int:
    """Number of representations of N as a sum of 4 squares of (positive or negative) odd integers."""
    # If N can be written as sum of 4 odd squares, then N ≡ 4 (mod 8).
    if N % 8 != 4:
        return 0
    # By Jacobi: r_4(N) for N=4(2m+1), this counts sum_p p_i^2 = N for p_i in Z.
    # For ALL p_i odd, we need each p_i^2 odd, i.e., p_i odd. The number of such
    # representations is exactly r_4(N) when N ≡ 4 (mod 8), since any decomposition
    # of N into 4 squares with sum ≡ 4 (mod 8) must have all squares odd.
    return r4(N)


def N_from_jacobi(k: int) -> int:
    """Compute N(k) using the Jacobi four-square machinery, by enumerating
    odd positive integer tuples (p1,...,p4) with sum p_i^2 ≤ (4k+2)^2."""
    target = (4 * k + 2) ** 2
    total = 0
    # The original positive-orthant count weighted by sign multiplicity:
    #   N(k) = Σ over (p1,...,p4) ∈ {1,3,5,...}^4 with Σ p_i^2 ≤ target of 2^{#{i : p_i > 1}}.
    # Equivalently, in terms of integer (positive or negative) odd p_i with sum p_i^2 ≤ target,
    # the count in Z^4 of tuples is r4_odd(N) summed over N. We must relate this to N(k).
    # The mapping is:
    #   c_i ∈ Z, p_i := 2|c_i|+1 ∈ {1,3,5,...}.
    #   For each (c_1,...,c_4), we get one (p_1,...,p_4) with p_i ≥ 1 odd.
    #   Conversely, given (p_1,...,p_4) ∈ {1,3,...}^4 (positive odd), the # of c's is
    #     2^{#{i : p_i > 1}}  (sign flexibility of c_i when c_i ≠ 0).
    # So N(k) = Σ over positive odd 4-tuples of 2^{#{i : p_i > 1}}.
    for p1 in range(1, 4*k+3, 2):
        s1 = p1*p1
        if s1 > target:
            break
        for p2 in range(1, 4*k+3, 2):
            s2 = s1 + p2*p2
            if s2 > target:
                break
            for p3 in range(1, 4*k+3, 2):
                s3 = s2 + p3*p3
                if s3 > target:
                    break
                for p4 in range(1, 4*k+3, 2):
                    s4 = s3 + p4*p4
                    if s4 > target:
                        break
                    nz = sum(1 for p in (p1,p2,p3,p4) if p > 1)
                    total += 1 << nz
    return total


def jacobi_sum(k: int) -> int:
    """Sum r_4(N) over all N with N ≤ (4k+2)^2 and N ≡ 4 (mod 8)."""
    target = (4 * k + 2) ** 2
    s = 0
    N = 4
    while N <= target:
        s += r4(N)
        N += 8
    return s


def verify(k_max: int = 5):
    print(f"{'k':>3} {'N(k)':>10} {'N_jacobi':>10} {'jac_sum':>12} {'jac_sum/16':>12} {'match?':>8}")
    for k in range(k_max + 1):
        n_direct = N_count(k)
        n_jacobi = N_from_jacobi(k)
        js = jacobi_sum(k)
        print(f"{k:>3d} {n_direct:>10d} {n_jacobi:>10d} {js:>12d} {js/16:>12.2f} {n_direct == n_jacobi}")


if __name__ == "__main__":
    import sys
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    verify(kmax)
