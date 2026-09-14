#!/usr/bin/env python3
from pathlib import Path
import hashlib
base=Path(__file__).resolve().parent
names=["run_highN_generic_null_N2_N40.py","make_plots.py","run_all.sh","make_sha256.py","highN_N2_N40_summary.csv","highN_reconstruction_conditioning.csv","mean_kappa_vs_N.png","mean_kappa_vs_N.svg","reconstruction_conditioning_vs_N.png","reconstruction_conditioning_vs_N.svg","analysis_highN_N2_N40.md","README.md"]
lines=[]
for n in names:
    p=base/n
    if p.exists(): lines.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {n}")
(base/"SHA256SUMS.txt").write_text("\n".join(lines)+"\n",encoding="utf-8")
