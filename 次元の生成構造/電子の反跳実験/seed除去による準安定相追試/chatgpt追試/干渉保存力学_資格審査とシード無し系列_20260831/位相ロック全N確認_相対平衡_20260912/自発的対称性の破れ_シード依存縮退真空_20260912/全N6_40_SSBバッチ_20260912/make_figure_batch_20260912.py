#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全N6-40 SSBバッチの集約図（CSV読出しのみ）。35/35確認と真空差のN依存。"""
import csv
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('ssb_batch_summary.csv')))
N=[int(r['N']) for r in rows];med=[float(r['vacuum_spread_med_deg']) for r in rows];mx=[float(r['vacuum_spread_max_deg']) for r in rows]
ok=[r['SSB_confirmed']=='True' for r in rows]
fig,ax=plt.subplots(1,2,figsize=(13,5))
ax[0].plot(N,med,'o-',label='vacuum spread median');ax[0].plot(N,mx,'s--',alpha=.6,label='vacuum spread max')
ax[0].set_xlabel('N');ax[0].set_ylabel('seed-to-seed phase-config spread [deg]')
ax[0].set_title(f'(a) degenerate-vacuum spread vs N (SSB={sum(ok)}/{len(ok)} all confirmed)');ax[0].legend(fontsize=8)
# 確認フラグ（全True）を可視化
labels=['inflation','equal amp','Δ=-2π/N','SSB']
cols=['onset_ok','all_equal_amp','all_lock_2piN','SSB_confirmed']
mat=np.array([[1 if r[c]=='True' else 0 for c in cols] for r in rows]).T
im=ax[1].imshow(mat,aspect='auto',cmap='Greens',vmin=0,vmax=1)
ax[1].set_yticks(range(4));ax[1].set_yticklabels(labels);ax[1].set_xticks(range(0,len(N),3));ax[1].set_xticklabels([N[i] for i in range(0,len(N),3)])
ax[1].set_xlabel('N');ax[1].set_title('(b) per-N checks all TRUE (green) — SSB universal N=6..40')
fig.tight_layout();fig.savefig('fig_ssb_batch_allN.png',dpi=120);print('wrote fig_ssb_batch_allN.png')
