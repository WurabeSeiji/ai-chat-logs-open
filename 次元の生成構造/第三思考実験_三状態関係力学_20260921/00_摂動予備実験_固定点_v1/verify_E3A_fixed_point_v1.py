#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Re-run E3-A fixed-point programs in a temporary directory and compare hashes."""
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
]
TARGETS = [
    'E3A0_hierarchical_two_kepler_perturbation_aggregates_v1.csv',
    'E3A0_hierarchical_two_kepler_perturbation_full_v1.json',
    'E3A1_perturbation_breakdown_scan_v1.csv',
    'E3A1_perturbation_breakdown_analysis_v1.md',
    'fig_E3A_01_raw_hierarchical_response_v1.png',
    'fig_E3A_02_invariants_and_normal_modes_v1.png',
    'fig_E3A_03_outer_law_shift_v1.png',
    'fig_E3A_04_approximation_error_and_mixing_v1.png',
    'fig_E3A_05_outer_mode_type_boundary_v1.png',
]


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    failed = []
    with tempfile.TemporaryDirectory(prefix='e3a_verify_') as td:
        td = Path(td)
        for name in PROGRAMS:
            shutil.copy2(HERE / name, td / name)
        subprocess.run(['python3', PROGRAMS[0]], cwd=td, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(['python3', PROGRAMS[1]], cwd=td, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(['python3', PROGRAMS[2]], cwd=td, check=True, stdout=subprocess.DEVNULL)
        for name in TARGETS:
            ok = sha(HERE / name) == sha(td / name)
            print(f'{name}: {"MATCH" if ok else "DIFF"}')
            if not ok:
                failed.append(name)
    if failed:
        raise SystemExit('verification failed: ' + ', '.join(failed))
    print('ALL FIXED-POINT OUTPUTS MATCH')


if __name__ == '__main__':
    main()
