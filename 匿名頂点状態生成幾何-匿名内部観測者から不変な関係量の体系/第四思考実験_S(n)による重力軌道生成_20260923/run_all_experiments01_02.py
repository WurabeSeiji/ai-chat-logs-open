#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess,sys
from pathlib import Path
here=Path(__file__).resolve().parent
cmds=[
 [sys.executable,str(here/'run_reference_pn_orbit.py')],
 [sys.executable,str(here/'run_experiment01_sn_conservative.py'),'--order','12'],
 [sys.executable,str(here/'run_experiment02_sn_2p5pn_rr.py'),'--order','12'],
 [sys.executable,str(here/'compare_experiments_01_02.py'),'--order','12'],
]
for c in cmds:
    print('+',' '.join(c),flush=True); subprocess.run(c,check=True)
