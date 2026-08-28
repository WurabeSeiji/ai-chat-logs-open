#!/usr/bin/env python3
"""Path-explicit deterministic local refinement of a fixed-cardinality subset."""
from __future__ import annotations
import argparse,itertools
from pathlib import Path
import numpy as np,pandas as pd

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--edges',required=True);ap.add_argument('--max-iter',type=int,default=20);a=ap.parse_args();df=pd.read_csv(a.input);q=df.z2_re.to_numpy()+1j*df.z2_im.to_numpy();w=np.abs(q);mp={(int(r.i),int(r.j)):i for i,r in df.iterrows()};S=set(mp[tuple(map(int,e.split('-')))] for e in a.edges.replace(',',';').split(';') if e)
 def rr(T): idx=np.array(sorted(T));return float(abs(q[idx].sum())/w[idx].sum())
 for _ in range(a.max_iter):
  best=rr(S);bestS=None;sel=list(S);uns=[i for i in range(len(df)) if i not in S]
  for o in sel:
   for x in uns:
    T=(S-{o})|{x};r=rr(T)
    if r<best:best=r;bestS=T
  if bestS is None:break
  S=bestS
 print(f'{rr(S):.17g}\t'+','.join(f'{int(df.iloc[i].i)}-{int(df.iloc[i].j)}' for i in sorted(S)))
if __name__=='__main__':main()
