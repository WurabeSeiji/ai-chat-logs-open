#!/usr/bin/env python3
import subprocess, sys
from pathlib import Path
here=Path(__file__).resolve().parent
base=here.parent
subprocess.run([sys.executable,str(here/'run_mass_ratio_suite.py'),'--out',str(base)],check=True)
subprocess.run([sys.executable,str(here/'plot_mass_ratio_suite.py'),'--base',str(base)],check=True)
