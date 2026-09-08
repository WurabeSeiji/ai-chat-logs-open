#!/usr/bin/env python3
import csv, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def rows(name):
    with (HERE/name).open(encoding='utf-8') as f: return list(csv.DictReader(f))
s=rows('even_N_harmonic_summary.csv')
assert len(s)==19
assert all(int(r['integer_harmonic_relations_2to4'])==0 for r in s)
assert all(int(r['three_wave_sum_resonances'])==0 for r in s)
assert all(int(r['low_den_rational_relations_q_le_16'])==0 for r in s)
assert all(int(r['unique_positive_frequency_count'])==3 for r in s if int(r['N'])>=6)
assert int(s[0]['unique_positive_frequency_count'])==1
assert max(float(r['floor_eigen_residual']) for r in s)<1e-12
assert [int(r['N']) for r in s if int(r['photon_pair_complete_expected_N8k'])==1]==[8,16,24,32,40]
print('PASS: high-symmetry harmonic-spectrum audit verified')
