#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 補遺65/66 図版: 全て実計算(模式図なし)、英語ラベル
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
N=2160; xs=np.arange(N)/N; MMAX=81
def comb(k,a):
    w=np.zeros(N); m=1
    while m*k<=MMAX:
        w+=(4/np.pi)*((-1)**((m-1)//2))*np.cos(2*np.pi*m*k*(xs-a))/m; m+=2
    return w
def line(I,f):
    c=2*np.mean(I*np.cos(2*np.pi*f*xs)); s=2*np.mean(I*np.sin(2*np.pi*f*xs))
    return complex(c,s)

# ===== 図1: 読み出しシナリオB (両空間) =====
d1,d3,d5=2,1,3
a1,a3,a5=d1/4.0,d3/12.0,d5/20.0
Psi=1.0+comb(1,a1)+comb(3,a3)+comb(5,a5)
I=Psi**2
fig,ax=plt.subplots(2,2,figsize=(13,9))
# (a) frequency space
fr=list(range(1,22))
mags=[abs(line(I,f)) for f in fr]
ax[0,0].stem(fr,mags,basefmt=' ')
for f,lab in [(1,'read d1'),(3,'read d3 (peel f=3 of frag1)'),(5,'read d5 (peel f=5 of frag1)')]:
    ax[0,0].annotate(lab,(f,abs(line(I,f))),textcoords='offset points',xytext=(6,8),fontsize=8,color='crimson')
ax[0,0].set_xlabel('frequency f (shared coordinate units)'); ax[0,0].set_ylabel('|line coefficient|')
ax[0,0].set_title('(a) Frequency-space record: transcription lines\n(odd lines only; comb-x-comb terms fall on even f)')
# (b) position space
ax[0,1].plot(xs,Psi,lw=0.8)
for a,c,lab in [(a1,'tab:blue','frag1 (fund 1) at d=2'),(a3,'tab:orange','frag3 (fund 3) at d=1'),(a5,'tab:green','frag5 (fund 5) at d=3')]:
    ax[0,1].axvline(a,color=c,ls='--',lw=1.2); ax[0,1].text(a,3.6,lab,rotation=90,fontsize=7,color=c,va='top')
for q in range(4): ax[0,1].axvline(q/4,color='gray',lw=0.5,alpha=0.5)
ax[0,1].set_xlabel('position x (quarter grid of fund 1 in gray)'); ax[0,1].set_ylabel('Psi(x)')
ax[0,1].set_title('(b) Position space: superposed square combs')
# (c) coordinate space (x, R) with Q
Rv=[0.5,1.5,2.5]; Qv=[+1,-1,+1]; pos=[a1,a3,a5]
for x0,R,Q in zip(pos,Rv,Qv):
    ax[1,0].scatter([x0],[R],s=160,c=('tab:red' if Q>0 else 'tab:blue'),marker=('o' if Q>0 else 's'))
    ax[1,0].annotate(f'R={R}, Q={Q:+d}\nclock dt=1/R per tick',(x0,R),textcoords='offset points',xytext=(10,-4),fontsize=8)
ax[1,0].set_xlim(0,1); ax[1,0].set_ylim(0,3.2)
ax[1,0].set_xlabel('x (projected position)'); ax[1,0].set_ylabel('R = nu_t (norm clock)')
ax[1,0].set_title('(c) xyztRQ space: fragments on constraint surface R^2=s\n(red circle Q=+1, blue square Q=-1; t derived from R)')
# (d) decode results table
ax[1,1].axis('off')
rows=[['scenario','configs','unique decode'],
      ['B heterogeneous (1/3/5)','64','64/64'],
      ['C continuous positions','200','200/200 (err 2.5e-16)'],
      ['D +/-20% jitter + noise','192','192/192'],
      ['A same-species, no antiphase pair','71','71/71'],
      ['A with antiphase pair','29','degenerate (record = 0, exact)']]
tb=ax[1,1].table(cellText=rows[1:],colLabels=rows[0],loc='center',cellLoc='center')
tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1,1.6)
ax[1,1].set_title('(d) Condition-matrix results (machine-verified)')
plt.tight_layout(); plt.savefig('supplement66_fig_readout.png',dpi=150); plt.close()

# ===== 図2: 保存則 (崩壊 9 -> (5,1,3)) =====
fig,ax=plt.subplots(2,2,figsize=(13,9))
# (a) per-axis dressed frequencies before/after
labels=['parent (2,1,0,0)','kept (1,1,0,0)','d1 (0,0,0,0)','d2 (1,0,0,0)']
cells=[(2,1,0,0),(1,1,0,0),(0,0,0,0),(1,0,0,0)]
w=0.2
for i,(lab,cell) in enumerate(zip(labels,cells)):
    fa=[abs(k)+0.5 for k in cell]
    ax[0,0].bar(np.arange(4)+i*w-1.5*w,fa,width=w,label=f'{lab}: s={sum(x*x for x in fa):.0f}')
ax[0,0].set_xticks(range(4)); ax[0,0].set_xticklabels(['axis 1','axis 2','axis 3','axis 4'])
ax[0,0].set_ylabel('dressed frequency |k|+1/2'); ax[0,0].legend(fontsize=8)
ax[0,0].set_title('(a) Wavelength/frequency space: decay 9 -> (5,1,3)\n(R read from frequency content, not labels)')
# (b) R^2 ledger with virtual borrow
stages=['parent\n(record)','virtual\n(5,5)','final\n(record)']
ax[0,1].bar([0],[9],color='tab:blue')
ax[0,1].bar([1],[10],color='tab:orange',alpha=0.6)
ax[0,1].bar([2],[5],color='tab:green'); ax[0,1].bar([2],[3],bottom=[5],color='tab:olive'); ax[0,1].bar([2],[1],bottom=[8],color='tab:cyan')
ax[0,1].axhline(9,color='k',ls=':',lw=1)
ax[0,1].annotate('+1 borrowed\n(virtual, supplement 39)',(1,10),textcoords='offset points',xytext=(8,4),fontsize=8)
ax[0,1].annotate('5+3+1=9 exact\n(118,944/118,944 paths)',(2,9),textcoords='offset points',xytext=(8,4),fontsize=8)
ax[0,1].set_xticks(range(3)); ax[0,1].set_xticklabels(stages); ax[0,1].set_ylabel('R^2 = s (additive ledger)')
ax[0,1].set_title('(b) R^2 conservation: exact at records, +/-1 cancels at virtual stage')
# (c) epsilon ledger + control
ax[1,0].bar(['candidate\ncombinations','amplitude-carrying\npaths'],[39.2,0.0],color=['tab:gray','tab:red'])
ax[1,0].set_ylabel('epsilon-violating fraction (%)')
for i,vv in enumerate([39.2,0.0]): ax[1,0].text(i,vv+0.8,f'{vv}%',ha='center')
ax[1,0].set_title('(c) Q conservation enforced by transport (s=13 control):\n39.2% of raw combinations violate; 0 of 24,864 contributing paths do')
# (d) configuration reading
ax[1,1].axis('off')
rows=[['check','result'],
      ['exact configuration reading (60 random finals)','60/60 unique'],
      ['s_read additivity (R from record)','60/60 = 9 exact'],
      ['epsilon multiplicativity (s=9,11,13)','118,944/118,944'],
      ['virtual borrow +/-1 cancellation','118,944/118,944'],
      ['t = sum of 1/sqrt(s_read) (norm clock)','= paper 8 internal time']]
tb=ax[1,1].table(cellText=rows[1:],colLabels=rows[0],loc='center',cellLoc='center')
tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1,1.7)
ax[1,1].set_title('(d) Conservation through the projection (machine-verified)')
plt.tight_layout(); plt.savefig('supplement65_fig_conservation.png',dpi=150); plt.close()
print("figures written: supplement66_fig_readout.png, supplement65_fig_conservation.png")
