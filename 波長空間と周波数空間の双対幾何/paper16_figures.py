#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 論文16 図版 (機械検証済み数値の提示、英語ラベル)
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

# fig1: 二結果圧縮
fig,ax=plt.subplots(1,2,figsize=(12.5,4.4))
x9=['{531,333}','{333,333}']
der9=[0.98195,0.01805]; seq9=[0.98263,0.01737]
w=0.35
ax[0].bar(np.arange(2)-w/2,der9,w,label='derived (amplitude)')
ax[0].bar(np.arange(2)+w/2,seq9,w,label='sequential (counting)')
ax[0].set_xticks(range(2)); ax[0].set_xticklabels(x9)
ax[0].set_title('(a) s=9 twin: exclusion forces a 2-outcome space\nand both measures pile ~0.98 on X  →  Δ(9)=0.00068')
ax[0].legend(fontsize=8); ax[0].set_ylabel('probability')
x11=['{137,335}','{155,335}','{335,335}']
der11=[0.31422,0.22624,0.45954]; seq11=[0.40415,0.37767,0.21819]
ax[1].bar(np.arange(3)-w/2,der11,w,label='derived (sector A)')
ax[1].bar(np.arange(3)+w/2,seq11,w,label='sequential')
ax[1].set_xticks(range(3)); ax[1].set_xticklabels(x11,fontsize=8)
ax[1].set_title('(b) s=11 twin: 3 outcomes — compression released\n→  Δ(11)=0.24136 (×350 jump): the 0.0007 was an artifact')
ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('paper16_fig1_compression.png',dpi=150); plt.close()

# fig2: 乖離の地形 (Δ系列、機械検証値)
fig,ax=plt.subplots(figsize=(11,4.4))
labels=['s=9\n(2,1,0,0)A','s=11\nA','s=11\nB','s=13\n(2,1,1,1)A','s=13\n(2,1,1,1)B','s=13\n(2,2,0,0)','s=13\n(3,0,0,0)A']
dbin=[0.00068,0.24136,0.43971,0.17392,0.30860,0.90418,0.80348]
dcfg=[None,0.14452,0.14237,None,None,0.90418,None]
xx=np.arange(len(labels))
ax.bar(xx-0.2,dbin,0.4,label='Δ(s): channel-consistent derived vs sequential')
ax.bar([x+0.2 for x,v in zip(xx,dcfg) if v is not None],[v for v in dcfg if v is not None],0.4,color='tab:orange',label="Δ'(s): configuration granularity")
ax.axhline(0.00098,color='gray',ls=':',lw=1)
ax.text(0.1,0.03,'TV(seq, batch) at s=9 = 0.00098 (counting-side internal width)',fontsize=7,color='gray')
ax.set_xticks(xx); ax.set_xticklabels(labels,fontsize=8)
ax.set_ylabel('total variation distance')
ax.set_title('Structural divergence: derived measures do not classicalize to counting measures\n(grows with s; persists under convention width; sector dependence shrinks at configuration granularity)')
ax.legend(fontsize=8)
plt.tight_layout(); plt.savefig('paper16_fig2_divergence.png',dpi=150); plt.close()

# fig3: 粒度の物理 (鏡像親) + 四値ロック
fig,ax=plt.subplots(1,2,figsize=(12.5,4.4))
groups=['m=2 (s=7)\nA','m=2\nmirror B','m=3 (s=13)\nA','m=3\nmirror B','m=4 (s=21)\nA ch(5,7,9)','m=4\nmirror B ch(5,7,9)']
Wch=[52,32,9216,0,147456,0]; Wcf=[36,16,768,768,3072,3072]
xx=np.arange(6)
b1=ax[0].bar(xx-0.2,Wch,0.4,label='W (channel-consistent)')
b2=ax[0].bar(xx+0.2,Wcf,0.4,label="W' (configuration granularity)")
ax[0].set_yscale('symlog'); ax[0].set_xticks(xx); ax[0].set_xticklabels(groups,fontsize=7)
for b,v in zip(list(b1)+list(b2),Wch+Wcf): ax[0].text(b.get_x()+b.get_width()/2,max(v,1)*1.2,str(v),ha='center',fontsize=7)
ax[0].set_title("(a) Granularity is not a convention: mirror cancellation W=0 with W'>0\nfull at m=3; per-CHANNEL at m=4 (ch (5,7,9)); absent at m=2.\nW' mirror-equality holds for ALL channels at m>=3 (suppl. 84)")
ax[0].legend(fontsize=8)
sizes=[12,12,6,6]
lbl=['ratio 0 (χ=−i)','ratio ∞ (χ=+i)','ratio +1 (A-only)','ratio −1 (B-only)']
ax[1].pie(sizes,labels=[f'{l}\n{s} configs' for l,s in zip(lbl,sizes)],autopct='',startangle=90,textprops={'fontsize':8})
ax[1].set_title('(b) Four-value lock (m=3, all 36 configs, zero exceptions):\n$z^-_K/z^+_K\\in\\{0,\\infty,+1,-1\\}$; exclusive families weigh equally (384=384)')
plt.tight_layout(); plt.savefig('paper16_fig3_granularity.png',dpi=150); plt.close()
print("paper16 figures written")
