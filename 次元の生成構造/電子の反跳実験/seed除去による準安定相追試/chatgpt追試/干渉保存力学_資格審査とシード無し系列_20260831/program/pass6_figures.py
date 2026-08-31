#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス6：図化（読出し専用）。`--static` で走行に依存しない図（fig1, fig2）のみ。
fig1 閉塞 step0（本パッケージ 4 生成法 ＋ 8/29 の 4 系統）／fig2 乱数状態の埋め込み（複素 vs 実）／
fig3 H⊥/H 14 パネル（N）× 4 生成法／fig4 最終状態 PR/M と振幅幅／fig5 ρ−1 帯（予測）と実測分類／fig6 λ 予測 vs 実測／fig7 t50 vs N"""
import os, csv, json, sys
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__))); FIG=os.path.join(ROOT,'figures'); os.makedirs(FIG,exist_ok=True)
STATIC='--static' in sys.argv
meth=['mp','hm','ne','rb']; name={'mp':'make_parent equimodular','hm':'handmade equimodular','ne':'non-equimodular','rb':'random balanced'}
col={'mp':'k','hm':'#1f77b4','ne':'#ff7f0e','rb':'#2ca02c'}; mk={'mp':'s','hm':'o','ne':'^','rb':'D'}
# fig1
C=list(csv.DictReader(open(os.path.join(ROOT,'results','closure_step0_4methods.csv'))))
R=list(csv.DictReader(open(os.path.join(ROOT,'data','reference','closure_step0_4systems_20260829.csv'))))
fig,ax=plt.subplots(figsize=(9,4.2))
for m in meth:
    xs=[int(r['N']) for r in C if r['method']==m]; ys=[float(r['global_closure_abs_sum_z2_over_H']) for r in C if r['method']==m]
    ax.scatter(xs,ys,marker=mk[m],color=col[m],s=45,label=name[m]+' (this package)')
kk={k.strip().lstrip('\ufeff').lower():k for k in R[0].keys()}; Ncol=kk['n']; scol=kk['system']; vcol=kk['abs_sum_z2_over_h']   # 参照 csv（2026-08-29）の列（BOM/大文字小文字を吸収）
for i,s in enumerate(['original','fixed_baseline','fixed','fixed_equimodular']):
    xs=[int(r[Ncol])+0.12*(i+1) for r in R if r[scol]==s]; ys=[max(abs(float(r[vcol])),1e-18) for r in R if r[scol]==s]
    ax.scatter(xs,ys,marker='x',s=30,alpha=.7,label=f'{s} (2026-08-29 reference)')
ax.set_yscale('log'); ax.set_ylim(1e-18,1e-9); ax.set_xlabel('N'); ax.set_ylabel('|Σ z²| / H at step 0'); ax.grid(alpha=.3); ax.legend(fontsize=6,ncol=2)
ax.set_title('zero closure holds for every self-consistent parent, regardless of generation method / normalization / amplitude distribution',fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig1_closure_step0.png'),dpi=150); plt.close()
# fig2
S=list(csv.DictReader(open(os.path.join(ROOT,'results','embed_random_summary.csv'))))
fig,axs=plt.subplots(1,2,figsize=(10,3.8))
axs[0].semilogy([int(r['N']) for r in S],[float(r['complex_embed_err_rel_max']) for r in S],'o-'); axs[0].set_ylim(1e-17,1e-13); axs[0].set_xlabel('N'); axs[0].set_ylabel('max relative embedding error (100 random states)'); axs[0].set_title('complex: every random state embeds exactly (Takagi), rank = N−1',fontsize=8); axs[0].grid(alpha=.3)
axs[1].plot([int(r['N']) for r in S],[int(r['real_PSD_count']) for r in S],'s-',color='#d62728'); axs[1].set_ylim(-3,103); axs[1].set_xlabel('N'); axs[1].set_ylabel('# embeddable (B ⪰ 0) of 100'); axs[1].set_title('real reading d²=|z|²: positivity cone excludes almost all states',fontsize=8); axs[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig2_embed_random.png'),dpi=150); plt.close()
print('static figures done')
if STATIC: raise SystemExit
pred={r['tag']:r for r in csv.DictReader(open(os.path.join(ROOT,'results','parents_predictions.csv')))}
D={r['tag']:r for r in csv.DictReader(open(os.path.join(ROOT,'results','dynamics_summary.csv')))}
def load(tag):
    f=os.path.join(ROOT,'data',tag,'treatment_linear124_amplitude_aware_timeseries.csv'); f=f if os.path.exists(f) else f+'.gz'
    return np.genfromtxt(f,delimiter=',',names=True) if os.path.exists(f) else None
# fig3
fig,axs=plt.subplots(4,4,figsize=(17,14)); axs=axs.ravel()
for i,N in enumerate(range(3,17)):
    ax=axs[i]
    for m in meth:
        t=f'{m}_N{N}'; a=load(t)
        if a is None: continue
        ax.semilogy(a['step'],np.maximum(a['H_perp']/a['H_total'],1e-40),color=col[m],lw=1.1,label=f"{m}: pred {float(pred[t]['pred_lambda_f']):.4f}")
    ax.set_ylim(1e-34,3); ax.set_title(f'N={N}',fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6,loc='lower right')
    if i%4==0: ax.set_ylabel('H⊥/H')
    if i>=10: ax.set_xlabel('step')
for j in range(14,16): axs[j].axis('off')
fig.suptitle('H⊥/H over 40000 steps (L=124, seedless, direct readout, interference-preserving frame): four parent-generation methods',fontsize=11)
plt.tight_layout(rect=(0,0,1,0.97)); plt.savefig(os.path.join(FIG,'fig3_Hperp_grid_N3_N16.png'),dpi=130); plt.close()
# fig4
fig,axs=plt.subplots(1,2,figsize=(11,4))
for m in meth:
    rs=[D[t] for t in D if D[t]['method']==m]
    sat=[r for r in rs if r['measured_kind']=='saturated']
    axs[0].scatter([int(r['N']) for r in sat],[float(r['final_PR_over_M']) for r in sat],marker=mk[m],color=col[m],s=45,label=name[m]+' (saturated)')
    axs[1].scatter([int(r['N']) for r in rs],[float(r['final_amp_std'])/np.sqrt(float(r['norm'])**2/(int(r['N'])*(int(r['N'])-1)/2)) for r in rs],marker=mk[m],color=col[m],s=40,label=name[m])
axs[0].axhline(1,color='gray',ls='--',lw=1); axs[0].text(3.2,1.02,'old dynamics (phase-only K): equipartition PR/M = 1',fontsize=7,color='gray')
axs[0].set_xlabel('N'); axs[0].set_ylabel('final PR/M (step 40000)'); axs[0].set_ylim(0,1.1); axs[0].grid(alpha=.3); axs[0].legend(fontsize=7); axs[0].set_title('saturated runs localize (PR/M ≪ 1)',fontsize=9)
axs[1].set_xlabel('N'); axs[1].set_ylabel('final std|z| / rms|z|'); axs[1].grid(alpha=.3); axs[1].legend(fontsize=7); axs[1].set_title('amplitude spread at step 40000 (0 = equimodular)',fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig4_final_state.png'),dpi=150); plt.close()
# fig5
fig,ax=plt.subplots(figsize=(9,4.2))
for m in meth:
    rs=[D[t] for t in D if D[t]['method']==m]
    for r in rs:
        y=max(float(r['pred_rho_minus_1']),1e-9); filled=r['measured_kind']!='floor'
        ax.scatter(int(r['N'])+{'mp':-0.22,'hm':-0.07,'ne':0.07,'rb':0.22}[m],y,marker=mk[m],s=45,color=col[m] if filled else 'white',edgecolors=col[m],linewidths=1.3)
ax.axhspan(1e-9,1e-3,color='green',alpha=.08); ax.axhline(1e-3,color='green',lw=.8,ls=':'); ax.text(3,3e-4,'neutral band (floor within 40000 steps)',fontsize=8,color='green')
ax.set_yscale('log'); ax.set_ylim(1e-9,5e-2); ax.set_xlabel('N'); ax.set_ylabel('ρ − 1 (co-rotating monodromy, predicted before run)'); ax.grid(alpha=.3)
for m in meth: ax.scatter([],[],marker=mk[m],color=col[m],label=name[m])
ax.scatter([],[],marker='o',color='white',edgecolors='k',label='open = measured floor'); ax.legend(fontsize=7,loc='upper right',ncol=2)
ax.set_title('instability prediction vs measured outcome (filled = saturated/growing, open = floor)',fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig5_rho_bands.png'),dpi=150); plt.close()
# fig6
fig,ax=plt.subplots(figsize=(5.2,5.2))
for m in meth:
    rs=[D[t] for t in D if D[t]['method']==m and D[t]['lambda_measured']]
    ax.scatter([float(r['pred_lambda_f']) for r in rs],[float(r['lambda_measured']) for r in rs],marker=mk[m],color=col[m],s=40,label=name[m])
ax.plot([0,0.025],[0,0.025],'--',color='gray'); ax.set_xlim(0,0.025); ax.set_ylim(0,0.025); ax.set_xlabel('predicted λ_f'); ax.set_ylabel('measured λ (ln H⊥ slope)'); ax.grid(alpha=.3); ax.legend(fontsize=7); ax.set_title('growth rate: prediction vs measurement',fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig6_lambda_pred_vs_meas.png'),dpi=150); plt.close()
# fig7
fig,ax=plt.subplots(figsize=(8,4))
for m in meth:
    rs=[D[t] for t in D if D[t]['method']==m and D[t]['t50_measured'] not in ('','None')]
    ax.scatter([int(r['N']) for r in rs],[int(float(r['t50_measured'])) for r in rs],marker=mk[m],color=col[m],s=45,label=name[m]+' measured')
    rp=[D[t] for t in D if D[t]['method']==m and D[t]['pred_t50']]
    ax.scatter([int(r['N'])+0.15 for r in rp],[float(r['pred_t50']) for r in rp],marker=mk[m],facecolors='none',edgecolors=col[m],s=45)
ax.axhline(40000,color='gray',ls='--',lw=.8); ax.text(3,42000,'run length',fontsize=7,color='gray')
ax.set_yscale('log'); ax.set_xlabel('N'); ax.set_ylabel('t50 (steps to H⊥/H = 0.5)'); ax.grid(alpha=.3); ax.legend(fontsize=7); ax.set_title('saturation time: filled = measured, open = predicted from ρ (f0 = 3e-32)',fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG,'fig7_t50_vs_N.png'),dpi=150); plt.close()
# fig8（本フレーム追加図・様式は fig3 と同一）：閉塞率 |ΣZ²|/H。旧フレームでは保存量（〜1e-13）、本フレームでは力学量
fig,axs=plt.subplots(4,4,figsize=(17,14)); axs=axs.ravel()
for i,N in enumerate(range(3,17)):
    ax=axs[i]
    for m in meth:
        t=f'{m}_N{N}'; a=load(t)
        if a is None: continue
        ax.semilogy(a['step'],np.maximum(a['abs_ZT_Z']/a['H_total'],1e-18),color=col[m],lw=1.1,label=m)
    ax.set_ylim(1e-18,3); ax.set_title(f'N={N}',fontsize=10); ax.grid(alpha=.3); ax.legend(fontsize=6,loc='lower right')
    if i%4==0: ax.set_ylabel('|Σ Z²| / H')
    if i>=10: ax.set_xlabel('step')
for j in range(14,16): axs[j].axis('off')
fig.suptitle('closure |Σ Z²|/H over 40000 steps: conserved (~1e-13) in the previous frame, dynamical in the interference-preserving frame',fontsize=11)
plt.tight_layout(rect=(0,0,1,0.97)); plt.savefig(os.path.join(FIG,'fig8_closure_grid_N3_N16.png'),dpi=130); plt.close()
print('PASS6 OK')
