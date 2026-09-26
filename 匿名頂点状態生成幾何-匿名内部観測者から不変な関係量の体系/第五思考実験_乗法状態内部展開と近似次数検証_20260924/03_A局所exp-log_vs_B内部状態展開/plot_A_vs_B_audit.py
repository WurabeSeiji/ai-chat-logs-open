#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input',required=True); ap.add_argument('--out',required=True); a=ap.parse_args()
    d=json.load(open(a.input)); c=d['csv']
    labels=['p','e','t','xy']; A=[c['A_'+k] for k in labels]; B=[c['B_'+k] for k in labels]
    x=range(len(labels)); fig,ax=plt.subplots(figsize=(7,5))
    w=.36
    ax.bar([i-w/2 for i in x],A,w,label='A: local exp/log')
    ax.bar([i+w/2 for i in x],B,w,label='B: explicit microstate')
    ax.set_yscale('log'); ax.set_xticks(list(x),labels); ax.set_ylabel('max absolute error vs saved Paper-4 C1 data')
    ax.set_title('Audit paths A and B'); ax.grid(True,axis='y',which='both',alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(a.out,bbox_inches='tight'); plt.close(fig)
if __name__=='__main__': main()
