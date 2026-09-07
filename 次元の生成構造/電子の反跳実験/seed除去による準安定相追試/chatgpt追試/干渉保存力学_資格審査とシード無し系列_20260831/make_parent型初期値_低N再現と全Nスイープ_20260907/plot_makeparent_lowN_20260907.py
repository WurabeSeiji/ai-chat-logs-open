#!/usr/bin/env python3
import os,csv
import numpy as np
import matplotlib.pyplot as plt
BASE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(BASE,'lowN_validation','results')
N_LIST=[3,4,5,6,7]
rows=list(csv.DictReader(open(os.path.join(OUT,'timeseries.csv'))))
by={}
for r in rows:
    k=(int(r['N']),r['series']); by.setdefault(k,([],[])); by[k][0].append(int(r['step'])); by[k][1].append(float(r['Hperp_frac']))
order=['N-2','N-1','N','N+1','N+2','124']
fig,axs=plt.subplots(2,3,figsize=(14,8)); axs=axs.ravel()
for k,N in enumerate(N_LIST):
    ax=axs[k]
    for label in order:
        if (N,label) in by:
            x,y=by[(N,label)]; ax.semilogy(x,y,label=('2pi/124' if label=='124' else label),linewidth=1.0)
    ax.set_xlim(0,500); ax.set_ylim(1e-34,3); ax.set_title(f'N={N}'); ax.grid(alpha=.25); ax.set_xlabel('step'); ax.set_ylabel('Hperp/H'); ax.legend(fontsize=7,loc='lower right')
axs[-1].axis('off'); fig.suptitle('Actual make_parent low-N control: Hperp/H vs step (N=3..7)',y=.995); fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig_Hperp_makeparent_actual_lowN.png'),dpi=180); plt.close(fig)
fig,axs=plt.subplots(len(N_LIST),2,figsize=(9,3.2*len(N_LIST)))
for i,N in enumerate(N_LIST):
    Z=np.load(os.path.join(OUT,f'hm_N{N}_den_{N}_states_500.npz'))['Z']
    for j,step in enumerate([0,500]):
        z=Z[step]; ax=axs[i,j]; ax.scatter(z.real,z.imag,s=28); ax.axhline(0,lw=.5); ax.axvline(0,lw=.5)
        m=np.max(np.abs(np.r_[z.real,z.imag])); lim=1.15*m; ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim); ax.set_aspect('equal'); ax.grid(alpha=.2); ax.set_title(f'N={N}, step={step}, den={N}')
        if i==len(N_LIST)-1: ax.set_xlabel('Re')
        if j==0: ax.set_ylabel('Im')
fig.tight_layout(); fig.savefig(os.path.join(OUT,'fig_complex_plane_step0_step500_denN_N3_N7.png'),dpi=180); plt.close(fig)
print('PLOTS DONE')
