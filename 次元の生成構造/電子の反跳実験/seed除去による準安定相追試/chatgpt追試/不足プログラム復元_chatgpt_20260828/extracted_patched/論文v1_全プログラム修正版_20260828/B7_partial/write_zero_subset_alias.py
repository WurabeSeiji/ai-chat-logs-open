#!/usr/bin/env python3
"""Create the historical filename N3_N16_zero_subset_exact_covers.csv.

Audit found this file is byte-identical to N3_N16_zero_triple_exact_covers.csv.
The original analysis script wrote only the latter name; the former was an alias.
"""
from pathlib import Path
import argparse, shutil

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--source',type=Path,required=True,help='N3_N16_zero_triple_exact_covers.csv'); ap.add_argument('--out',type=Path,required=True); a=ap.parse_args(); a.out.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(a.source,a.out); print(a.out)
if __name__=='__main__':main()
