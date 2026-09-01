#!/usr/bin/env python3
import csv
import matplotlib.pyplot as plt
src='timeseries_den124_500.csv'
out='fig3_Hperp_grid_N3_N16_den124_500.png'
data={N:([],[]) for N in range(3,17)}
with open(src) as f:
    for r in csv.DictReader(f):
        N=int(r['N']); data[N][0].append(int(r['step'])); data[N][1].append(float(r['Hperp_frac']))
fig,axs=plt.subplots(4,4,figsize=(14,11)); axs=axs.ravel()
for k,N in enumerate(range(3,17)):
    ax=axs[k]; x,y=data[N]
    ax.semilogy(x,y,linewidth=0.8)
    ax.set_xlim(0,500); ax.set_ylim(1e-34,3)
    ax.set_title(f'N={N}'); ax.grid(alpha=0.25)
    if k//4==3: ax.set_xlabel('step')
    if k%4==0: ax.set_ylabel('Hperp/H')
for k in range(14,16): axs[k].axis('off')
fig.suptitle('hm N=3..16, denominator=124, 500 steps',y=0.995)
fig.tight_layout(); fig.savefig(out,dpi=180)
