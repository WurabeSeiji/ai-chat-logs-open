#!/usr/bin/env python3
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
HERE=Path(__file__).resolve().parent
m=pd.read_csv(HERE/'main_results_by_N.csv'); c=pd.read_csv(HERE/'d_est_N_K_curves.csv'); k=pd.read_csv(HERE/'kappaB_null_contrast.csv')

fig,ax=plt.subplots(figsize=(9,6))
ax.fill_between(m.N,m.boot_lo,m.boot_hi,alpha=.2,label='seed-cluster bootstrap 95% CI (A)')
ax.plot(m.N,m.d_consensus_A,'o-',label='d_consensus (A system)')
ax.plot(m.N,m.d_consensus_B,'s--',label='d_consensus (B system)')
ax.plot(m.N,2*m.N-4,':',color='gray',label='generic 2N-4')
ax.axhline(4,color='crimson',lw=1.2,label='H_B1: d=4')
ax.set_xlabel('N'); ax.set_ylabel('d_est (K=N-1)'); ax.grid(alpha=.3); ax.legend()
ax.set_title('Phase B main estimation (5 seeds): d_terminal(N) vs H_B1')
fig.tight_layout(); fig.savefig(HERE/'dimension_vs_N.png',dpi=170); fig.savefig(HERE/'dimension_vs_N.svg'); plt.close(fig)

fig,axs=plt.subplots(2,2,figsize=(11,8))
for ax,N in zip(axs.ravel(),[10,20,30,40]):
    g=c[(c.N==N)&(c.system=='A')]
    ax.plot(g.K,g.d_consensus,'o-',label='consensus A')
    gB=c[(c.N==N)&(c.system=='B')]
    ax.plot(gB.K,gB.d_consensus,'s--',label='consensus B',alpha=.7)
    ax.axhline(4,color='crimson',lw=1)
    ax.set_title(f'N={N} (n={5*N})'); ax.set_xlabel('K'); ax.set_ylabel('d_est'); ax.grid(alpha=.3); ax.legend(fontsize=8)
fig.suptitle('d_est(N,K) saturation curves'); fig.tight_layout()
fig.savefig(HERE/'dimension_vs_K_N10_20_30_40.png',dpi=170); fig.savefig(HERE/'dimension_vs_K_N10_20_30_40.svg'); plt.close(fig)

fig,ax=plt.subplots(figsize=(9,5.5))
ax.semilogy(k.N,k.kappa_B_median,'o-',label='terminal κ_B median')
ax.semilogy(k.N,k.kappa_B_max,'^--',label='terminal κ_B max',alpha=.7)
ax.semilogy(k.N,k.null_exact_mean,':',color='gray',label='generic null E[κ]')
ax.set_xlabel('N'); ax.set_ylabel('κ_B'); ax.grid(alpha=.3,which='both'); ax.legend()
ax.set_title('Terminal Bob-star closure vs generic null (log10 P ≈ −6…−8 per star)')
fig.tight_layout(); fig.savefig(HERE/'kappa_terminal_vs_null_N.png',dpi=170); fig.savefig(HERE/'kappa_terminal_vs_null_N.svg'); plt.close(fig)
print('figures done')
