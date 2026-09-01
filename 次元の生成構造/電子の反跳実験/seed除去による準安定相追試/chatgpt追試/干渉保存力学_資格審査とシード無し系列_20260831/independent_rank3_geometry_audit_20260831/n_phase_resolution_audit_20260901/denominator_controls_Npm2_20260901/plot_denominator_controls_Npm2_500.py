#!/usr/bin/env python3
import os,csv
import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.abspath(__file__))
FIG=os.path.join(ROOT,'figures'); os.makedirs(FIG,exist_ok=True)
rows=list(csv.DictReader(open(os.path.join(ROOT,'timeseries_metrics_denominator_controls_500.csv'))))
for off in (-2,-1,0,1,2):
    for floor in (1e-34,1e-5):
        fig,axs=plt.subplots(4,4,figsize=(14,11)); axs=axs.ravel()
        for k,N in enumerate(range(3,17)):
            ax=axs[k]; rr=[r for r in rows if int(r['N'])==N and int(r['offset'])==off]
            x=[int(r['step']) for r in rr]; y=[max(float(r['Hperp_frac']),1e-40) for r in rr]
            ax.semilogy(x,y,linewidth=.8); ax.set_xlim(0,500); ax.set_ylim(floor,3); ax.grid(alpha=.25); ax.set_title(f'N={N}')
            if k//4==3: ax.set_xlabel('step')
            if k%4==0: ax.set_ylabel('Hperp/H')
        for k in range(14,16): axs[k].axis('off')
        suffix='' if floor<1e-10 else '_floor1e-5'
        label='N' if off==0 else f'N{off:+d}'
        fig.suptitle(f'hm N=3..16, denominator={label}, 500 steps',y=.995); fig.tight_layout()
        fig.savefig(os.path.join(FIG,f'fig3_Hperp_grid_N3_N16_den_{label}{suffix}.png'),dpi=180); plt.close(fig)
