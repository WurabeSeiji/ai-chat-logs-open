#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス6b（追加図・読出し専用）：fig3 と同一様式（4×4 グリッド・配色・凡例・縦軸範囲）で、
横軸だけ 0→1000 step に制限した初期カーブ図。木原指示（2026-08-31）：
「横軸 0→1000 までの初期カーブが見えない。全パターンで横軸を 0→1000 にした図を作れ」。
出力：figures/fig3b_Hperp_grid_N3_N16_first1000.png"""
import os, csv
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=os.path.join(ROOT,'figures'); os.makedirs(FIG,exist_ok=True)
XMAX=1000
meth=['mp','hm','ne','rb']; col={'mp':'k','hm':'#1f77b4','ne':'#ff7f0e','rb':'#2ca02c'}
pred={r['tag']:r for r in csv.DictReader(open(os.path.join(ROOT,'results','parents_predictions.csv')))}
def load(tag):
    f=os.path.join(ROOT,'data',tag,'treatment_linear124_amplitude_aware_timeseries.csv'); f=f if os.path.exists(f) else f+'.gz'
    return np.genfromtxt(f,delimiter=',',names=True,max_rows=XMAX+2) if os.path.exists(f) else None
fig,axs=plt.subplots(4,4,figsize=(17,14)); axs=axs.ravel()
for i,N in enumerate(range(3,17)):
    ax=axs[i]
    for m in meth:
        t=f'{m}_N{N}'; a=load(t)
        if a is None: continue
        sel=a['step']<=XMAX
        ax.semilogy(a['step'][sel],np.maximum(a['H_perp'][sel]/a['H_total'][sel],1e-40),color=col[m],lw=1.1,label=f"{m}: pred {float(pred[t]['pred_lambda_f']):.4f}")
    ax.set_xlim(0,XMAX); ax.set_ylim(1e-34,3); ax.set_title(f'N={N}',fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6,loc='lower right')
    if i%4==0: ax.set_ylabel('H⊥/H')
    if i>=10: ax.set_xlabel('step')
for j in range(14,16): axs[j].axis('off')
fig.suptitle('H⊥/H over the first 1000 steps (L=124, seedless, direct readout, interference-preserving frame): four parent-generation methods',fontsize=11)
plt.tight_layout(rect=(0,0,1,0.97)); plt.savefig(os.path.join(FIG,'fig3b_Hperp_grid_N3_N16_first1000.png'),dpi=130); plt.close()
print('PASS6b OK (fig3b)')
