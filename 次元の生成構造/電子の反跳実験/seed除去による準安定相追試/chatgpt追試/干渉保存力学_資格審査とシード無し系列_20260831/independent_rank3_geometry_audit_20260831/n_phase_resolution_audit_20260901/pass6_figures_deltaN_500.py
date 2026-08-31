#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Plotting blocks copied from original pass6_figures.py for the hm-only deltaN control.
# Dynamics are not recomputed here; this is readout-only plotting.
import os, csv
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(ROOT,'figures'); os.makedirs(FIG,exist_ok=True)
rows=list(csv.DictReader(open(os.path.join(ROOT,'timeseries_metrics.csv'))))
D={}
for r in rows:
    N=int(r['N']); D.setdefault(N,[]).append(r)
# fig3: original pass6_figures.py layout/style, hm series only
fig,axs=plt.subplots(4,4,figsize=(17,14)); axs=axs.ravel()
for i,N in enumerate(range(3,17)):
    ax=axs[i]; rr=D[N]
    step=np.array([float(r['step']) for r in rr]); frac=np.array([float(r['Hperp_frac']) for r in rr])
    ax.semilogy(step,np.maximum(frac,1e-40),color='#1f77b4',lw=1.1,label='hm')
    ax.set_ylim(1e-34,3); ax.set_title(f'N={N}',fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6,loc='lower right')
    if i%4==0: ax.set_ylabel('Hperp/H')
    if i>=10: ax.set_xlabel('step')
for j in range(14,16): axs[j].axis('off')
fig.suptitle('Hperp/H over 500 steps (delta=2pi/N, seedless, direct readout, interference-preserving frame): handmade equimodular',fontsize=11)
plt.tight_layout(rect=(0,0,1,0.97)); plt.savefig(os.path.join(FIG,'fig3_Hperp_grid_N3_N16_deltaN_500.png'),dpi=130); plt.close()
# fig8: original pass6_figures.py layout/style, hm series only
fig,axs=plt.subplots(4,4,figsize=(17,14)); axs=axs.ravel()
for i,N in enumerate(range(3,17)):
    ax=axs[i]; rr=D[N]
    step=np.array([float(r['step']) for r in rr]); clo=np.array([float(r['global_closure']) for r in rr])
    ax.semilogy(step,np.maximum(clo,1e-18),color='#1f77b4',lw=1.1,label='hm')
    ax.set_ylim(1e-18,3); ax.set_title(f'N={N}',fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6,loc='lower right')
    if i%4==0: ax.set_ylabel('|sum Z^2| / H')
    if i>=10: ax.set_xlabel('step')
for j in range(14,16): axs[j].axis('off')
fig.suptitle('closure |sum Z^2|/H over 500 steps: delta=2pi/N, interference-preserving frame',fontsize=11)
plt.tight_layout(rect=(0,0,1,0.97)); plt.savefig(os.path.join(FIG,'fig8_closure_grid_N3_N16_deltaN_500.png'),dpi=130); plt.close()
print('PASS6 deltaN 500 plotting OK')
