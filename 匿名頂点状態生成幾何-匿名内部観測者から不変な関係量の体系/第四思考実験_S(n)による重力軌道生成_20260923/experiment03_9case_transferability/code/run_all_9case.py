#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
import subprocess, sys

here = Path(__file__).resolve().parent
base = here.parent
subprocess.run([sys.executable, str(here / "run_9case_suite.py"), "--out", str(base)], check=True)
subprocess.run([sys.executable, str(here / "plot_9case_suite.py"), "--base", str(base)], check=True)
