#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
import hashlib, shutil, subprocess, tempfile, sys

HERE=Path(__file__).resolve().parent
scripts=[
 'run_E3B0_symmetric_three_body_rigid_v1.py',
 'analyze_E3B1_symmetric_three_body_v1.py',
 'plot_E3B_symmetric_three_body_v1.py',
]
outputs=[
 'E3B0_symmetric_three_body_summary_v1.csv',
 'E3B0_symmetric_three_body_trajectories_v1.csv',
 'E3B0_symmetric_three_body_results_v1.json',
 'E3B1_symmetric_three_body_analysis_v1.md',
 'fig_E3B_01_rigid_triangle_snapshots_v1.png',
 'fig_E3B_02_quadratic_readout_orbits_v1.png',
 'fig_E3B_03_pair_distances_v1.png',
 'fig_E3B_04_symmetry_residuals_v1.png',
 'fig_E3B_05_masslike_ratios_v1.png',
]

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

with tempfile.TemporaryDirectory() as td:
    td=Path(td)
    for s in scripts: shutil.copy2(HERE/s,td/s)
    for s in scripts:
        subprocess.run([sys.executable,str(td/s)],cwd=td,check=True,stdout=subprocess.DEVNULL)
    lines=[]; ok=True
    for name in outputs:
        a=sha(HERE/name); b=sha(td/name); match=a==b
        lines.append(f'{name}: {"MATCH" if match else "MISMATCH"}')
        ok &= match
    lines.append('ALL E3-B FIXED-POINT OUTPUTS MATCH' if ok else 'E3-B VERIFICATION FAILED')
    text='\n'.join(lines)+'\n'
    (HERE/'verification_E3B_v1.txt').write_text(text)
    print(text,end='')
    if not ok: raise SystemExit(1)
