#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Re-run all E3-A fixed-point v2 programs and compare frozen data and all seven figures."""
from pathlib import Path
import hashlib
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
PROGRAMS = [
    'run_E3A0_hierarchical_two_kepler_perturbation_v1.py',
    'analyze_E3A1_perturbation_breakdown_v1.py',
    'plot_E3A_fixed_point_v1.py',
    'export_E3A2_orbit_trajectories_v1.py',
    'plot_E3A_orbits_v1.py',
]
TARGETS = [
    'E3A0_hierarchical_two_kepler_perturbation_aggregates_v1.csv',
    'E3A0_hierarchical_two_kepler_perturbation_full_v1.json',
    'E3A1_perturbation_breakdown_scan_v1.csv',
    'E3A1_perturbation_breakdown_analysis_v1.md',
    'E3A2_orbit_trajectories_v1.csv',
    'fig_E3A_01_raw_hierarchical_response_v1.png',
    'fig_E3A_02_invariants_and_normal_modes_v1.png',
    'fig_E3A_03_outer_law_shift_v1.png',
    'fig_E3A_04_approximation_error_and_mixing_v1.png',
    'fig_E3A_05_outer_mode_type_boundary_v1.png',
    'fig_E3A_06_orbit_raw_outer_readout_v1.png',
    'fig_E3A_07_orbit_outer_normal_mode_v1.png',
]


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    failed = []
    with tempfile.TemporaryDirectory(prefix='e3a_verify_v2_') as td:
        td = Path(td)
        for name in PROGRAMS:
            shutil.copy2(HERE / name, td / name)
        for name in PROGRAMS:
            subprocess.run(['python3', name], cwd=td, check=True, stdout=subprocess.DEVNULL)
        for name in TARGETS:
            ok = sha(HERE / name) == sha(td / name)
            print(f'{name}: {"MATCH" if ok else "DIFF"}')
            if not ok:
                failed.append(name)
    if failed:
        raise SystemExit('verification failed: ' + ', '.join(failed))
    print('ALL FIXED-POINT V2 OUTPUTS MATCH')


if __name__ == '__main__':
    main()
