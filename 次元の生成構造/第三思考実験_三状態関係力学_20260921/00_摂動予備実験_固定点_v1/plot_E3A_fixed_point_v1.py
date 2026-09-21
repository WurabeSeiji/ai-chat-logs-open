#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create fixed-point figures from saved E3-A CSV data only.

This program does NOT rerun the experiment. It consumes the frozen CSV outputs:
  E3A0_hierarchical_two_kepler_perturbation_aggregates_v1.csv
  E3A1_perturbation_breakdown_scan_v1.csv
and writes five PNG figures in the same directory.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
A0 = HERE / 'E3A0_hierarchical_two_kepler_perturbation_aggregates_v1.csv'
A1 = HERE / 'E3A1_perturbation_breakdown_scan_v1.csv'


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def f(rows, key):
    return np.array([float(r[key]) for r in rows], dtype=float)


def main():
    a0 = read_csv(A0)
    hier = [r for r in a0 if r['suite'] == 'hierarchical_nonres_31_127' and float(r['eps']) > 0]
    eps = f(hier, 'eps')

    # Fig 1: raw local response and conic departure.
    fig = plt.figure(figsize=(7.4, 5.0))
    ax = fig.add_subplot(111)
    ax.loglog(eps, f(hier, 'q_in_mod_med'), label='inner local-area modulation')
    ax.loglog(eps, f(hier, 'q_out_mod_med'), label='outer local-area modulation')
    ax.loglog(eps, f(hier, 'raw_conic_resid_med'), label='raw outer conic residual')
    ax.set_xlabel('coupling epsilon')
    ax.set_ylabel('dimensionless response')
    ax.set_title('E3-A0: raw hierarchical response')
    ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / 'fig_E3A_01_raw_hierarchical_response_v1.png', dpi=200)
    plt.close(fig)

    # Fig 2: exact invariants / normal-mode Kepler readout remain at numerical floor.
    fig = plt.figure(figsize=(7.4, 5.0))
    ax = fig.add_subplot(111)
    ax.loglog(eps, f(hier, 'q_total_drift_max'), label='total-area drift')
    ax.loglog(eps, f(hier, 'mode_q_drift_max'), label='normal-mode area drift')
    ax.loglog(eps, f(hier, 'mode_conic_resid_max'), label='normal-mode conic residual')
    ax.set_xlabel('coupling epsilon')
    ax.set_ylabel('relative residual')
    ax.set_title('E3-A0: conserved structure and exact normal modes')
    ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / 'fig_E3A_02_invariants_and_normal_modes_v1.png', dpi=200)
    plt.close(fig)

    a1 = read_csv(A1)
    nz = [r for r in a1 if float(r['eps']) > 0]
    e1 = f(nz, 'eps')

    # Fig 3: exact outer-law shift vs second-order approximation.
    fig = plt.figure(figsize=(7.4, 5.0))
    ax = fig.add_subplot(111)
    ax.loglog(e1, f(nz, 'tau_out_shift_exact'), label='exact outer-law shift')
    ax.loglog(e1, f(nz, 'tau_out_shift_Oeps2'), label='O(epsilon^2) approximation')
    ax.set_xlabel('coupling epsilon')
    ax.set_ylabel('delta tau_out')
    ax.set_title('E3-A1: perturbative law shift')
    ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / 'fig_E3A_03_outer_law_shift_v1.png', dpi=200)
    plt.close(fig)

    # Fig 4: perturbation error and hierarchy mixing are separate breakdown diagnostics.
    fig = plt.figure(figsize=(7.4, 5.0))
    ax = fig.add_subplot(111)
    ax.loglog(e1, f(nz, 'shift_Oeps2_rel_error'), label='relative error of O(epsilon^2) shift')
    ax.loglog(e1, f(nz, 'inner_weight_in_outer_mode'), label='inner weight in outer mode')
    ax.set_xlabel('coupling epsilon')
    ax.set_ylabel('dimensionless fraction')
    ax.set_title('E3-A1: approximation error and hierarchy mixing')
    ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / 'fig_E3A_04_approximation_error_and_mixing_v1.png', dpi=200)
    plt.close(fig)

    # Fig 5: exact outer eigenvalue crossing the elliptic/parabolic boundary.
    fig = plt.figure(figsize=(7.4, 5.0))
    ax = fig.add_subplot(111)
    ax.semilogx(e1, f(nz, 'lambda_out'), label='lambda_out')
    ax.semilogx(e1, np.full_like(e1, 2.0), label='elliptic/parabolic boundary')
    ax.set_xlabel('coupling epsilon')
    ax.set_ylabel('outer normal-mode eigenvalue')
    ax.set_title('E3-A1: outer normal-mode type boundary')
    ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / 'fig_E3A_05_outer_mode_type_boundary_v1.png', dpi=200)
    plt.close(fig)

    print('Generated 5 figures from frozen CSV data only.')


if __name__ == '__main__':
    main()
