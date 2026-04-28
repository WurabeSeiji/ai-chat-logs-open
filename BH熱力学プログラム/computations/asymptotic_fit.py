#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Asymptotic analysis of Δ(R) = V_4(R) - N(k) for R = 2k+1.

Hypothesis (boundary correction of the body Ω_R = {x : Σ(|x_i|+1/2)^2 ≤ R^2}):
    Δ(R) = (16π/3) R^3 - 6π R^2 + O(R)         (continuous boundary correction only)
=>  c := Δ(R) / (2π² R^3) → 8/(3π) ≈ 0.84883…
    subleading: -3/(π R)

Empirically the lattice point error E(R) = N(k) - Vol(Ω_R) is also O(R^?). We fit
Δ(R) = a R^3 + b R^2 + d R + e   with least-squares using k from k_min to k_max.
"""

import math
from collections import OrderedDict


def load_data(path):
    data = []
    with open(path, "r", encoding="utf-8") as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            k = int(parts[0])
            R = int(parts[1])
            N = int(parts[2])
            data.append((k, R, N))
    return data


def fit_polynomial_in_R(data, k_min=10, degree=3):
    """Least-squares fit Δ(R) = sum_{j=0..degree} c_j R^j over the given window.

    Returns coefficients c_j for j = 0..degree (constant ... R^degree)."""
    pts = [(R, math.pi**2 / 2 * R**4 - N) for (k, R, N) in data if k >= k_min]
    # Build matrix X with columns R^0, R^1, ..., R^degree; vector y = Δ
    n = degree + 1
    XtX = [[0.0] * n for _ in range(n)]
    Xty = [0.0] * n
    for R, delta in pts:
        for i in range(n):
            for j in range(n):
                XtX[i][j] += (R ** i) * (R ** j)
            Xty[i] += (R ** i) * delta
    # Solve the linear system XtX @ c = Xty
    coeffs = solve(XtX, Xty)
    return coeffs, pts


def solve(A, b):
    n = len(A)
    # Augment
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for i in range(n):
        # find pivot
        piv = i
        for j in range(i + 1, n):
            if abs(M[j][i]) > abs(M[piv][i]):
                piv = j
        M[i], M[piv] = M[piv], M[i]
        # eliminate
        for j in range(n):
            if j == i:
                continue
            ratio = M[j][i] / M[i][i]
            for k in range(i, n + 1):
                M[j][k] -= ratio * M[i][k]
    return [M[i][n] / M[i][i] for i in range(n)]


def main():
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "computations/results/N_table_k0_60.tsv"
    data = load_data(path)
    print(f"Loaded {len(data)} rows; max k = {data[-1][0]}, max R = {data[-1][1]}.")
    print()

    # Direct values of c(k) = Δ/(2π² R^3) at the largest k
    print("Direct estimate c(k) = Δ(R) / (2π² R^3):")
    for (k, R, N) in data[-10:]:
        delta = math.pi**2 / 2 * R**4 - N
        c = delta / (2 * math.pi**2 * R**3)
        print(f"  k={k:3d}  R={R:4d}  c={c:.6f}")
    print()

    # Theoretical leading c
    theo_leading = 16 * math.pi / 3 / (2 * math.pi**2)
    print(f"Theoretical leading c = 16π/3 / (2π²) = 8/(3π) = {theo_leading:.6f}")
    print()

    # Polynomial fits with increasing windows
    for k_min in (10, 20, 30, 40):
        try:
            coeffs, pts = fit_polynomial_in_R(data, k_min=k_min, degree=3)
        except Exception as e:
            print(f"  k_min={k_min}: fit failed ({e})")
            continue
        c0, c1, c2, c3 = coeffs
        c_lead = c3 / (2 * math.pi**2)
        c_sub  = c2 / (2 * math.pi**2)
        print(f"k_min={k_min:3d}, degree-3 fit:  Δ ≈ {c3:.6f}·R³ + {c2:+.6f}·R² + {c1:+.6f}·R + {c0:+.6f}")
        print(f"         leading c (R³ part) = {c_lead:.6f}    diff from 8/(3π): {c_lead - theo_leading:+.6f}")
        print(f"         subleading (R² part) = {c_sub:.6f}")

    # Also compare to predicted boundary correction  -6π R^2:
    predicted_R2_coef = -6 * math.pi
    print(f"\nPredicted R² coefficient (continuous boundary): -6π = {predicted_R2_coef:.6f}")


if __name__ == "__main__":
    main()
