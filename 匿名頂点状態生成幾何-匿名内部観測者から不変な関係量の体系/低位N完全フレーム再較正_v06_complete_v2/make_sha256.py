#!/usr/bin/env python3
from pathlib import Path
import hashlib

base = Path(__file__).resolve().parent
targets = [
    'run_lowN_v06_recalibration.py',
    'run_all.sh',
    'make_sha256.py',
    'verify_kappa_exact.py',
    'lowN_v06_samples.csv',
    'lowN_v06_summary.csv',
    'N3_identity_checks.csv',
    'N3_identity_summary.json',
    'kappa_exact_values.csv',
    'kappa_exact_comparison.csv',
    'analysis.md',
    'kappa_exact_distribution.md',
    'README.md',
]
lines=[]
for name in targets:
    p=base/name
    if p.exists():
        lines.append(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {name}')
(base/'SHA256SUMS.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print((base/'SHA256SUMS.txt').read_text(),end='')
