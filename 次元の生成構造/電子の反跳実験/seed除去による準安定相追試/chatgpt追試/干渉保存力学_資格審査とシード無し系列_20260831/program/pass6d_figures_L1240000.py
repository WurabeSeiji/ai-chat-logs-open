#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス6d（追加図・読出し専用）：L=1240000・500 step 診断走行（fig3c の最初の 5 step の内部）の H⊥/H を
fig3 と同一様式（4×4 グリッド・配色・semilogy・線形横軸）で描く。fig3c は残す（木原指示）。横軸は 0→5 step（同日指示：横軸 500 でなく 5 に）。
出力：figures/fig3d_Hperp_grid_N3_N16_L1240000_first5.png"""
import os
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=os.path.join(ROOT,'figures'); os.makedirs(FIG,exist_ok=True)
meth=['mp','hm','ne','rb']; col={'mp':'k','hm':'#1f77b4','ne':'#ff7f0e','rb':'#2ca02c'}
def load(tag):
    f=os.path.join(ROOT,'data',tag,'treatment_linear1240000_500steps_timeseries.csv')
    return np.genfromtxt(f,delimiter=',',names=True) if os.path.exists(f) else None
fig,axs=plt.subplots(4,4,figsize=(17,14)); axs=axs.ravel()
for i,N in enumerate(range(3,17)):
    ax=axs[i]
    for m in meth:
        t=f'{m}_N{N}'; a=load(t)
        if a is None: continue
        ax.semilogy(a['step'],np.maximum(a['H_perp']/a['H_total'],1e-40),color=col[m],lw=1.1,label=m)
    ax.set_xlim(0,5); ax.set_ylim(1e-34,3); ax.set_title(f'N={N}',fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6,loc='lower right')
    if i%4==0: ax.set_ylabel('H⊥/H')
    if i>=10: ax.set_xlabel('step')
for j in range(14,16): axs[j].axis('off')
fig.suptitle('H⊥/H over the first 5 steps at 10000x finer step (L=1240000, seedless, direct readout, interference-preserving frame)',fontsize=11)
plt.tight_layout(rect=(0,0,1,0.97)); plt.savefig(os.path.join(FIG,'fig3d_Hperp_grid_N3_N16_L1240000_first5.png'),dpi=130); plt.close()
print('PASS6d OK (fig3d)')
