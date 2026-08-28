#!/usr/bin/env python3
"""Run all six comparison summaries. Directories are passed as --p34 ... --p1415."""
from pathlib import Path
import argparse, subprocess, sys
PAIRS=[('p34',3,4),('p67',6,7),('p89',8,9),('p1011',10,11),('p1213',12,13),('p1415',14,15)]
def main():
 ap=argparse.ArgumentParser()
 for key,_,_ in PAIRS: ap.add_argument('--'+key,type=Path,required=True)
 ap.add_argument('--outdir',type=Path,required=True); a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
 here=Path(__file__).resolve().parent/'build_pair_comparison_summary.py'
 for key,n1,n2 in PAIRS:
  subprocess.run([sys.executable,str(here),'--pair-dir',str(getattr(a,key)),'--n1',str(n1),'--n2',str(n2),'--out',str(a.outdir/f'N{n1}_N{n2}_comparison_summary.csv')],check=True)
if __name__=='__main__': main()
