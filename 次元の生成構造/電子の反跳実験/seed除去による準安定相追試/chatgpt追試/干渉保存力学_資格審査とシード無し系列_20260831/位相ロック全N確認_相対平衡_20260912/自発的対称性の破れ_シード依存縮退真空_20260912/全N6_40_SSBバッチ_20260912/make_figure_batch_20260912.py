#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全N6-40 SSBバッチの集約図（per_N実測値の読出しのみ）。読める実測トレンドに刷新。
左: 真空差(中央/最大)のN依存。右: 実測 Δ_tail が −2π/N に一致（全シード平均）＋|Δ+2π/N|誤差・振幅比。"""
import math, os, re, glob
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import csv
# per_N から実測値を回収
rows={}
for f in glob.glob('per_N/N*.txt'):
    N=int(re.search(r'N(\d+)\.txt',f).group(1))
    dh=[];rat=[];dlt=[];derr=[];ls=[]
    for line in open(f,encoding='utf-8'):
        m=re.match(r'\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+(-?[\d.]+)\s+([\d.eE+-]+)',line)
        if m:
            ls.append(int(m.group(2)));dh.append(float(m.group(3)));rat.append(float(m.group(4)))
            dlt.append(float(m.group(5)));derr.append(float(m.group(6)))
    if dlt: rows[N]=dict(delta=np.mean(dlt),derr=np.max(derr),ratio=np.max(rat),lockstep=int(np.mean(ls)))
Ns=sorted(rows)
delta=[rows[N]['delta'] for N in Ns];ideal=[-2*math.pi/N for N in Ns]
derr=[rows[N]['derr'] for N in Ns];ratio=[rows[N]['ratio'] for N in Ns];lockstep=[rows[N]['lockstep'] for N in Ns]
# 真空差
cr=list(csv.DictReader(open('ssb_batch_summary.csv')))
med=[float(r['vacuum_spread_med_deg']) for r in cr];mx=[float(r['vacuum_spread_max_deg']) for r in cr];Nc=[int(r['N']) for r in cr]

fig,ax=plt.subplots(1,2,figsize=(14,5.5))
ax[0].plot(Nc,med,'o-',label='vacuum spread median');ax[0].plot(Nc,mx,'s--',alpha=.5,label='vacuum spread max')
ax[0].set_xlabel('N');ax[0].set_ylabel('seed-to-seed phase-config spread [deg]')
ax[0].set_title(f'(a) degenerate-vacuum spread vs N (SSB {sum(r["SSB_confirmed"]=="True" for r in cr)}/{len(cr)})');ax[0].legend(fontsize=8)
# 右: Δ_tail(実測) と -2π/N の一致（曲線が重なる）＋誤差・振幅比を2軸で
ax[1].plot(Ns,delta,'o',ms=6,label='measured Δ_tail (seed-mean)')
xx=np.linspace(min(Ns),max(Ns),200);ax[1].plot(xx,-2*math.pi/xx,'k-',lw=1,label='-2π/N (theory)')
ax[1].set_xlabel('N');ax[1].set_ylabel('1-step rotation Δ_tail [rad]');ax[1].legend(loc='lower right',fontsize=8)
ax2=ax[1].twinx()
ax2.semilogy(Ns,derr,'^',c='tab:red',ms=5,alpha=.7,label='|Δ+2π/N| error')
ax2.semilogy(Ns,np.array(ratio)-1,'v',c='tab:green',ms=5,alpha=.7,label='amp ratio -1')
ax2.set_ylabel('error (log): |Δ+2π/N|, ampratio-1');ax2.legend(loc='upper right',fontsize=8)
ax[1].set_title('(b) measured Δ_tail lands on -2π/N (err~1e-5), equal amp (ratio-1~1e-3)')
fig.suptitle('SSB batch N=6-40 (measured from per_N): degenerate vacua + universal -2π/N lock')
fig.tight_layout();fig.savefig('fig_ssb_batch_allN.png',dpi=120);print('wrote fig_ssb_batch_allN.png  (N数=%d)'%len(Ns))
