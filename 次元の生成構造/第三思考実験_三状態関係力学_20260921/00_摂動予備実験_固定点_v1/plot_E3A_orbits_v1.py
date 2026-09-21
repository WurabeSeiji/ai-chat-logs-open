#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot the two fixed-point orbit figures from saved E3-A2 trajectory CSV only."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE / 'E3A2_orbit_trajectories_v1.csv'


def read_rows():
    with DATA.open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def series(rows, eps):
    s = [r for r in rows if float(r['eps']) == eps]
    s.sort(key=lambda r: int(r['step']))
    return s


def draw(rows, xkey, ykey, title, filename):
    fig = plt.figure(figsize=(7.2, 7.2))
    ax = fig.add_subplot(111)
    eps_list = sorted({float(r['eps']) for r in rows})
    for eps in eps_list:
        s = series(rows, eps)
        ax.plot([float(r[xkey]) for r in s], [float(r[ykey]) for r in s], label=f'eps={eps:g}')
    ax.scatter([0.0], [0.0], s=18, label='readout origin')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('readout x')
    ax.set_ylabel('readout y')
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / filename, dpi=200)
    plt.close(fig)


def main():
    rows = read_rows()
    draw(rows, 'raw_outer_x', 'raw_outer_y',
         'E3-A orbit: raw outer readout trajectory',
         'fig_E3A_06_orbit_raw_outer_readout_v1.png')
    draw(rows, 'normal_outer_x', 'normal_outer_y',
         'E3-A orbit: outer normal-mode trajectory',
         'fig_E3A_07_orbit_outer_normal_mode_v1.png')
    print('Generated 2 orbit figures from frozen E3A2_orbit_trajectories_v1.csv only.')


if __name__ == '__main__':
    main()
