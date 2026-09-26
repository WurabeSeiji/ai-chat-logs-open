#!/usr/bin/env python3
import argparse, csv, math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def read_csv(path):
    with open(path,encoding='utf-8') as f:
        return list(csv.DictReader(f))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--angular',required=True)
    ap.add_argument('--integration',required=True)
    ap.add_argument('--outdir',required=True)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    ar=read_csv(a.angular); ir=read_csv(a.integration)

    # 1: angular truncation error, one line per (p0,e0).
    fig,ax=plt.subplots(figsize=(8,5.5))
    cases=sorted({(float(r['p0']),float(r['e0'])) for r in ar})
    for p,e in cases:
        rr=[r for r in ar if float(r['p0'])==p and float(r['e0'])==e]
        rr=sorted(rr,key=lambda r:int(r['N']))
        x=[int(r['N']) for r in rr]
        y=[max(float(r['phi_abs_error']),1e-16) for r in rr]
        ax.plot(x,y,marker='o',label=f'p={p:g}, e={e:g}')
    ax.set_yscale('log'); ax.set_xlabel('Angular truncation order N'); ax.set_ylabel('|Δφ| after 3 radial orbits')
    ax.set_title('Angular-series convergence without increasing physical state dimension')
    ax.grid(True,which='both',alpha=.25); ax.legend(fontsize=7,ncol=3)
    fig.tight_layout(); fig.savefig(out/'angular_order_convergence.svg',bbox_inches='tight'); plt.close(fig)

    # 2: geometric mean combined error over nine cases vs steps/orbit.
    fig,ax=plt.subplots(figsize=(7.5,5.5))
    for order in (4,6,8):
        xs=[]; ys=[]
        for steps in (8,16,32,64):
            vals=[float(r['combined']) for r in ir if int(r['order'])==order and int(r['steps'])==steps]
            if vals:
                xs.append(steps); ys.append(math.exp(sum(math.log(v) for v in vals)/len(vals)))
        ax.loglog(xs,ys,marker='o',label=f'order {order}')
    ax.set_xlabel('Steps per orbit'); ax.set_ylabel('Geometric mean combined error (9 cases)')
    ax.set_title('Integration order raised with the same five persistent states')
    ax.grid(True,which='both',alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(out/'integration_order_convergence.svg',bbox_inches='tight'); plt.close(fig)

    # 3: empirical order from 16 -> 32 for phi, P, E.
    metrics=[('phi_abs','φ'),('P_rel','P'),('E_abs','E')]
    orders=[4,6,8]
    means={o:[] for o in orders}
    for o in orders:
        for m,_ in metrics:
            vals=[]
            for p,e in sorted({(r['p0'],r['e0']) for r in ir}):
                r16=next(r for r in ir if r['p0']==p and r['e0']==e and int(r['order'])==o and int(r['steps'])==16)
                r32=next(r for r in ir if r['p0']==p and r['e0']==e and int(r['order'])==o and int(r['steps'])==32)
                vals.append(math.log(float(r16[m])/float(r32[m]),2))
            means[o].append(sum(vals)/len(vals))
    fig,ax=plt.subplots(figsize=(7.5,5.5))
    x=np.arange(len(orders),dtype=float); width=.23
    for j,(_,label) in enumerate(metrics):
        ax.bar(x+(j-1)*width,[means[o][j] for o in orders],width,label=label)
    ax.plot(x,orders,marker='o',linestyle='--',label='design order')
    ax.set_xticks(x,[str(o) for o in orders]); ax.set_xlabel('Designed Taylor order'); ax.set_ylabel('Mean empirical order (16→32 steps/orbit)')
    ax.set_title('Observed convergence order'); ax.grid(True,axis='y',alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(out/'integration_empirical_order.svg',bbox_inches='tight'); plt.close(fig)

if __name__=='__main__': main()
