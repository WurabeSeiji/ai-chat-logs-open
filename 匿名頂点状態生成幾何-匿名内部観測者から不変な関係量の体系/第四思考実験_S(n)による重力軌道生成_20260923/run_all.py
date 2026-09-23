#!/usr/bin/env python3
import subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
cmds=[
 [sys.executable,str(HERE/'run_reference_pn_orbit.py')],
 [sys.executable,str(HERE/'run_sn_gravity_orbit.py'),'--order','12'],
 [sys.executable,str(HERE/'compare_orbits.py'),'--order','12'],
]
for c in cmds:
 print('+',' '.join(c),flush=True); subprocess.run(c,check=True,cwd=HERE)
