#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Export deterministic orbit trajectories for E3-A fixed-point v1.

This script regenerates only the trajectories needed for the orbit figures from
exactly the same recurrence and fixed parameters used by E3-A0.  It writes a
CSV; plotting is performed separately by plot_E3A_orbits_v1.py.
"""
from pathlib import Path
import csv
import math
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'E3A2_orbit_trajectories_v1.csv'


def rot(t):
    return np.array([[math.cos(t), -math.sin(t)], [math.sin(t), math.cos(t)]], dtype=float)


def module(n, e, u0):
    th = 2 * math.pi / n
    T = np.diag([math.sqrt(1 - e), math.sqrt(1 + e)])
    S = T @ rot(th) @ np.linalg.inv(T)
    x0 = T @ np.array([math.cos(u0), math.sin(u0)])
    xm1 = np.linalg.inv(S) @ x0
    return 2 * math.cos(th), xm1, x0


def simulate(n_in, n_out, eps, u_in, u_out, steps=1400, e_in=.35, e_out=.25):
    ti, im1, i0 = module(n_in, e_in, u_in)
    to, om1, o0 = module(n_out, e_out, u_out)
    xi = np.empty((steps + 1, 2))
    xo = np.empty((steps + 1, 2))
    xi[0] = i0
    xo[0] = o0
    xi[1] = ti * i0 + eps * o0 - im1
    xo[1] = to * o0 + eps * i0 - om1
    for k in range(1, steps):
        xi[k + 1] = ti * xi[k] + eps * xo[k] - xi[k - 1]
        xo[k + 1] = to * xo[k] + eps * xi[k] - xo[k - 1]

    H = np.array([[ti, eps], [eps, to]])
    lam, V = np.linalg.eigh(H)
    jout = int(np.argmax(np.abs(V[1, :])))
    mout = V[0, jout] * xi + V[1, jout] * xo
    return xo, mout, float(lam[jout])


def main():
    n_in, n_out = 31, 127
    u_in, u_out = 0.37, 1.11
    eps_list = [0.0, 1e-3, 4e-3, 8e-3]
    samples = 4 * (n_out // 2 + 1)

    rows = []
    for eps in eps_list:
        xo, mout, lambda_out = simulate(n_in, n_out, eps, u_in, u_out)
        zraw = (xo[:samples, 0] + 1j * xo[:samples, 1]) ** 2
        zmode = (mout[:samples, 0] + 1j * mout[:samples, 1]) ** 2
        for k in range(samples):
            rows.append({
                'n_in': n_in,
                'n_out': n_out,
                'u_in': u_in,
                'u_out': u_out,
                'eps': eps,
                'lambda_out': lambda_out,
                'step': k,
                'raw_outer_x': float(zraw[k].real),
                'raw_outer_y': float(zraw[k].imag),
                'normal_outer_x': float(zmode[k].real),
                'normal_outer_y': float(zmode[k].imag),
            })

    with OUT.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'Wrote {len(rows)} rows to {OUT.name}')


if __name__ == '__main__':
    main()
