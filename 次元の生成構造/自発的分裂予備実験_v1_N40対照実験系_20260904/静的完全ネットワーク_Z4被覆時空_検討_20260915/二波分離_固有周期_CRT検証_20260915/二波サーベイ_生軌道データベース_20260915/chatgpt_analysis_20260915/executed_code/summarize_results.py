#!/usr/bin/env python3
"""Common plots and summaries, only after all 99 saved trajectories are analyzed."""
from pathlib import Path
import json,csv,os
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(os.environ.get('SURVEY_OUTPUT_DIR','/mnt/data/survey_analysis_20260915'));F=P/'figures';F.mkdir(exist_ok=True)
R=json.loads((P/'all_results.json').read_text());D=pd.read_csv(P/'all_runs_metrics.csv');
INST={x['run_id']:x for x in json.loads((P/'all_instantaneous.json').read_text())}
FRA={x['run_id']:x for x in json.loads((P/'all_frames.json').read_text())}
A=[];W=[]
for r in R:
 rid=r['run_id']; tt=np.load(P/'timeseries'/f'{rid}.npz'); w=r['windows']['W4'];ins=INST[rid]['samples'][-1]; f=FRA[rid]
 o={k:r[k] for k in ['run_id','series','L','M','ma','mb','parity_class','den','T','chi','analytic_fixed','orbit_abs_max','amp_cv_initial','amp_cv_final','W3_basis_capture_W4','sha256']}
 o.update(r99_W1=r['windows']['W1']['un']['r99'],r99_W4=w['un']['r99'],Reff_W4=w['un']['Reff'],p1_W4=w['un']['p1'],E2_W4=w['un']['E2'],
  centered_r99_W4=w['centered']['r99'],mean_fraction_W4=w['mean_fraction'],eta_initial_W4=w['initial_mode_eta_mean'],
  retained_initial_energy_W4=float(np.mean(1-tt['eta_initial'][3072:4096]**2)),
  instantaneous_dominant_t4096=ins['dominant_weight'],instantaneous_pair_t4096=ins['pair_weight'],
  instantaneous_pair_t0=INST[rid]['samples'][0]['pair_weight'],instantaneous_plane_overlap_3584_4096=ins.get('plane_overlap_previous'),
  fft_peak_period_W4=r['spectra']['W4']['peak_period_steps'],fft_peak_share_W4=r['spectra']['W4']['peak_power_fraction'],
  min_return_abs=r['initial_return_min'],min_return_U1=r['projective_return_min'],step_U1_mean_W4=r['step_projective_W4_mean'])
 if f['frame_valid']:o.update(f['windows']['W4'])
 A.append(o)
 for k,s in r['windows'].items():
  q={z:r[z] for z in ['run_id','series','L','ma','mb','den','parity_class','chi','analytic_fixed']};q['window']=k
  q.update(s['un']);q.update({'centered_'+x:y for x,y in s['centered'].items()});q.update({z:s[z] for z in ['mean_fraction','dynamic_fraction','initial_mode_eta_mean','amp_cv_mean','amp_pr_fraction_mean','generator_norm_mean']});W.append(q)
O=pd.DataFrame(A);O.to_csv(P/'overview_99.csv',index=False);pd.DataFrame(W).to_csv(P/'window_metrics.csv',index=False)
NON=O[~O.analytic_fixed];BASE=NON[NON.series!='den40_control'];MAIN=O[O.series=='main']
summary={'runs':len(O),'fixed':int(O.analytic_fixed.sum()),'nonfixed':len(NON),'all_hashes_verified':True,
 'complex_values':int((O.M*4097).sum()),'snapshots':99*4097,
 'nonfixed_rank99_le2':int((NON.r99_W4<=2).sum()),'nonfixed_rank99_eq1':int((NON.r99_W4==1).sum()),
 'nonfixed_rank99_ge10':int((NON.r99_W4>=10).sum()),
 'all_nonfixed_amplitude_cv_decreased':bool((NON.amp_cv_final<NON.amp_cv_initial).all()),
 'cv_final_initial_ratio_median':float((NON.amp_cv_final/NON.amp_cv_initial).median()),
 'instantaneous_pair_initial_median':float(NON.instantaneous_pair_t0.median()),
 'instantaneous_pair_final_min':float(NON.instantaneous_pair_t4096.min()),
 'instantaneous_pair_final_median':float(NON.instantaneous_pair_t4096.median()),
 'instantaneous_pair_final_max':float(NON.instantaneous_pair_t4096.max()),
 'all_nonfixed_instantaneous_pair_increased':bool((NON.instantaneous_pair_t4096>NON.instantaneous_pair_t0).all()),
 'nonfixed_mean_fraction_W4_max':float(NON.mean_fraction_W4.max()),
 'initial_retained_energy_W4_median':float(NON.retained_initial_energy_W4.median()),
 'nonfixed_min_initial_projective_return':float(NON.min_return_U1.min()),
 'spearman_chi_amp_no_controls':float(BASE.chi.corr(BASE.amp_cv_final,method='spearman')),
 'spearman_chi_rank_no_controls':float(BASE.chi.corr(BASE.r99_W4,method='spearman'))}
summary['within_L']={str(L):{'n':len(g),'chi_amp_spearman':float(g.chi.corr(g.amp_cv_final,method='spearman')),'chi_rank_spearman':float(g.chi.corr(g.r99_W4,method='spearman'))} for L,g in BASE.groupby('L')}
summary['qa_max_H']=float(D.dH_rel_max.max());summary['qa_max_Q2']=float(D.dQ2_rel_max.max())
(P/'summary_statistics.json').write_text(json.dumps(summary,indent=2))
# Compare all nine preselected timestep controls, not a cherry-picked sample.
pairs=[]
for _,c in O[O.series=='den40_control'].iterrows():
 m=O[(O.L==c.L)&(O.ma==c.ma)&(O.mb==c.mb)&(O.den==c.L)].iloc[0]
 pairs.append(dict(L=c.L,ma=c.ma,mb=c.mb,chi=c.chi,r99_main=m.r99_W4,r99_control=c.r99_W4,
  period_main=m.fft_peak_period_W4,period_control=c.fft_peak_period_W4,
  transfer_main=m.W3_basis_capture_W4,transfer_control=c.W3_basis_capture_W4,
  omega_main=m.omega_generator_mean,omega_control=c.omega_generator_mean))
C=pd.DataFrame(pairs);C.to_csv(P/'timestep_comparisons.csv',index=False)
PO=[(1,3),(1,5),(3,5),(2,4),(2,6),(4,6),(1,2),(2,3),(3,6)];LS=[8,10,12,16,20,24,32,40]
def save(name):
 plt.tight_layout();plt.savefig(F/name,dpi=180,bbox_inches='tight');plt.close()
# Every main condition on one complete map.
mat=np.array([[MAIN[(MAIN.ma==a)&(MAIN.mb==b)&(MAIN.L==L)].iloc[0].r99_W4 for L in LS] for a,b in PO],float)
plt.figure(figsize=(10,6));plt.imshow(np.log10(mat),aspect='auto');plt.colorbar(label='log10(r99)')
plt.xticks(range(8),LS);plt.yticks(range(9),[f'{a},{b}  '+('OO' if a%2 and b%2 else 'EE' if a%2==0 and b%2==0 else 'OE') for a,b in PO]);plt.xlabel('L');plt.ylabel('Initial modes (ma, mb)');plt.title('All 72 main runs: r99 on steps 3072-4095')
for i,(a,b) in enumerate(PO):
 for j,L in enumerate(LS):plt.text(j,i,str(int(mat[i,j]))+('*' if (a+b)%L==0 else ''),ha='center',va='center',bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
save('01_main_rank_map.png')
plt.figure(figsize=(8,5))
for tag,mark in [('odd-odd','o'),('even-even','s'),('odd-even','^')]:
 g=BASE[BASE.parity_class==tag];plt.scatter(g.chi,g.amp_cv_final,marker=mark,label=tag,alpha=.8)
plt.xlabel('Conserved chi = |Q2| / H');plt.ylabel('Final edge-amplitude CV');plt.title('85 active main/gcd runs: amplitude spread vs conserved chi');plt.legend();plt.grid(alpha=.2);save('02_amplitude_vs_chi.png')
plt.figure(figsize=(8,5))
for tag,mark in [('odd-odd','o'),('even-even','s'),('odd-even','^')]:
 g=NON[NON.parity_class==tag];plt.scatter(g.r99_W4,g.instantaneous_pair_t4096,marker=mark,label=tag,alpha=.8)
plt.xscale('log');plt.xlabel('Temporal r99 (1024-state late window)');plt.ylabel('Weight of dominant instantaneous +/- eigenpair at t=4096');plt.title('Many temporal directions do not imply many simultaneous eigenmodes');plt.legend();plt.grid(alpha=.2);save('03_instantaneous_vs_temporal.png')
plt.figure(figsize=(8,5))
for tag,mark in [('odd-odd','o'),('even-even','s'),('odd-even','^')]:
 g=NON[NON.parity_class==tag];plt.scatter(g.E2_W4,g.W3_basis_capture_W4,marker=mark,label=tag,alpha=.8)
plt.plot([0,1],[0,1],linestyle='--',label='Equal capture');plt.xlabel('Best rank-2 capture fitted within W4');plt.ylabel('W4 capture by fixed W3 rank-2 basis');plt.title('Local concentration versus persistence of the same subspace');plt.legend();plt.grid(alpha=.2);save('04_subspace_transfer.png')
plt.figure(figsize=(10,5));x=np.arange(len(C));plt.plot(x,C.r99_main,'o-',label='Main den=L');plt.plot(x,C.r99_control,'s-',label='Control den=40');plt.xticks(x,[f'L{int(c.L)}\n({int(c.ma)},{int(c.mb)})' for _,c in C.iterrows()]);plt.ylabel('r99 on steps 3072-4095');plt.title('All nine timestep controls');plt.legend();plt.grid(alpha=.2);save('05_timestep_controls.png')
plt.figure(figsize=(8,5))
for name,g in NON.groupby('series'):plt.scatter(g.generator_cycles_per_step_mean,g.sampled_cycles_per_step,label=name,alpha=.8)
a=np.linspace(0,.5,100);b=np.linspace(.5,1,100);plt.plot(a,a,linestyle='--',label='Principal-branch folding');plt.plot(b,b-1,linestyle='--');plt.xlabel('Instantaneous in-plane generator cycles per update');plt.ylabel('Measured principal-branch frame cycles per update');plt.title('Fast generator rotation and folded sampled phase');plt.legend();plt.grid(alpha=.2);save('06_phase_folding.png')
# Full audit map, all 99 rows.
plt.figure(figsize=(9,24));zz=np.array([[r['windows'][w]['un']['r99'] for w in ['W1','W2','W3','W4']] for r in R]);plt.imshow(np.log10(zz),aspect='auto');plt.colorbar(label='log10(r99)',shrink=.5);plt.xticks(range(4),['0-1023','1024-2047','2048-3071','3072-4095']);plt.yticks(range(99),[r['run_id'] for r in R],fontsize=6);plt.title('All 99 runs, equal-length time windows');plt.xlabel('Stored step window')
for i in range(99):
 for j in range(4):plt.text(j,i,str(zz[i,j]),ha='center',va='center',fontsize=6,bbox=dict(facecolor='white',edgecolor='none',alpha=.65))
save('07_all99_window_ranks.png')
# Late plane evolution for all sizes at fixed pair (2,4), not fitted over a chosen favorable case.
plt.figure(figsize=(9,5))
for L in LS:
 rid=f'L{L}_ma2_mb4_den{L}';s=INST[rid]['samples'];plt.plot([v['t'] for v in s[1:]],[v.get('plane_overlap_previous',np.nan) for v in s[1:]],marker='.',label=f'L={L}')
plt.xlabel('Step');plt.ylabel('Dominant generator-plane overlap with previous 512-step sample');plt.title('The instantaneous dominant plane can move while its weight stays high');plt.legend(ncol=4);plt.grid(alpha=.2);save('08_instantaneous_plane_motion.png')
print(json.dumps(summary,indent=2));print(C.round(6).to_string(index=False))
