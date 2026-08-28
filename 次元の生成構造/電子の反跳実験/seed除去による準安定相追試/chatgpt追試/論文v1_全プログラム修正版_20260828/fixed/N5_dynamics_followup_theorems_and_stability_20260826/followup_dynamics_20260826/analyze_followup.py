#!/usr/bin/env python3
from pathlib import Path
import math,json
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
GAMMA=math.tan(math.pi/144.0)
ANGLE=2.0*math.pi/144.0  # FIX3

# 1) Pump depletion / orthogonal transfer from existing corrected raw-K run
raw=pd.read_csv(ROOT/'N5_sigma_normalization_artifact_test/N5_raw_K_raw_observables.csv')
raw[['step','H_total','H_parallel','H_perp','abs_ztz']].to_csv(HERE/'pump_depletion_timeseries.csv',index=False)
pump_summary={
 'max_abs_Hsum_minus_Htotal':float(np.max(np.abs(raw.H_parallel+raw.H_perp-raw.H_total))),
 'Htotal_min':float(raw.H_total.min()),'Htotal_max':float(raw.H_total.max()),
 'Hparallel_initial':float(raw.H_parallel.iloc[0]),'Hparallel_final':float(raw.H_parallel.iloc[-1]),
 'Hperp_initial':float(raw.H_perp.iloc[0]),'Hperp_final':float(raw.H_perp.iloc[-1]),
 'max_abs_ZtZ':float(raw.abs_ztz.max())}
(HERE/'pump_depletion_summary.json').write_text(json.dumps(pump_summary,indent=2),encoding='utf-8')
plt.figure(figsize=(9,5)); plt.plot(raw.step,raw.H_parallel,label=r'$H_\parallel$'); plt.plot(raw.step,raw.H_perp,label=r'$H_\perp$'); plt.plot(raw.step,raw.H_total,label=r'$H_{total}$',linestyle='--'); plt.xlabel('step'); plt.ylabel('quadratic amplitude'); plt.title('N=5 raw-K: exact redistribution between parent plane and transverse sector'); plt.legend(); plt.tight_layout(); plt.savefig(HERE/'N5_pump_depletion_Hparallel_Hperp.png',dpi=180); plt.close()

# 2) Normalized-vs-raw time reparameterization test
norm=pd.read_csv(ROOT/'N5_sigma_normalization_artifact_test/N5_normalized_K_raw_observables.csv')
for d,kind in [(norm,'normalized'),(raw,'raw')]:
    if kind=='normalized': inc=np.full(len(d),ANGLE)  # FIX3: 線形回転では位相増分は ANGLE（正規化）/ ANGLE·σ（raw）
    else: inc=ANGLE*d.sigma_exact.to_numpy()
    phi=np.r_[0,np.cumsum(inc[:-1])]
    d['cumulative_phase']=phi
# interpolate overlap and compare Hperp in log regime and full regime
lo=max(norm.cumulative_phase.min(),raw.cumulative_phase.min()); hi=min(norm.cumulative_phase.max(),raw.cumulative_phase.max())
grid=np.linspace(lo,hi,20000)
hn=np.interp(grid,norm.cumulative_phase,norm.H_perp); hr=np.interp(grid,raw.cumulative_phase,raw.H_perp)
mask=(hn>1e-18)&(hr>1e-18)&(hn<0.5)&(hr<0.5)
log_rmse=float(np.sqrt(np.mean((np.log10(hn[mask])-np.log10(hr[mask]))**2)))
lin_rmse=float(np.sqrt(np.mean((hn-hr)**2)))
phase_summary={'overlap_phase_min':float(lo),'overlap_phase_max':float(hi),'log10_Hperp_RMSE_growth_regime':log_rmse,'Hperp_RMSE_full_overlap':lin_rmse}
(HERE/'time_reparameterization_summary.json').write_text(json.dumps(phase_summary,indent=2),encoding='utf-8')
pd.concat([norm.assign(branch='normalized'),raw.assign(branch='raw')])[['branch','step','cumulative_phase','H_perp','H_parallel','H_total','sigma_exact']].to_csv(HERE/'time_reparameterization_timeseries.csv',index=False)
plt.figure(figsize=(9,5)); plt.semilogy(norm.cumulative_phase,np.maximum(norm.H_perp,1e-30),label='K/sigma normalization'); plt.semilogy(raw.cumulative_phase,np.maximum(raw.H_perp,1e-30),label='raw K'); plt.xlabel('cumulative dominant-mode Cayley phase'); plt.ylabel(r'$H_\perp$'); plt.title('N=5: normalization comparison after phase-clock reparameterization'); plt.legend(); plt.tight_layout(); plt.savefig(HERE/'N5_normalized_vs_raw_cumulative_phase.png',dpi=180); plt.close()

# 3) tol sweep analysis
ts=pd.read_csv(HERE/'tol_sweep_summary.csv')
x=-np.log(ts.parent_residual.to_numpy()); y=ts['onset_f_ge_1e-8'].to_numpy(float); coef=np.polyfit(x,y,1); pred=np.polyval(coef,x); r2=1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2)
tol_summary={'onset_vs_minus_ln_parent_residual_slope_steps':float(coef[0]),'intercept':float(coef[1]),'R2':float(r2),'growth_rate_mean':float(ts.log_f_growth_rate.mean()),'growth_rate_std':float(ts.log_f_growth_rate.std(ddof=1))}
(HERE/'tol_sweep_fit.json').write_text(json.dumps(tol_summary,indent=2),encoding='utf-8')
plt.figure(figsize=(8,5)); plt.scatter(x,y); xx=np.linspace(x.min(),x.max(),100); plt.plot(xx,np.polyval(coef,xx)); plt.xlabel(r'$-\ln$(make_parent residual)'); plt.ylabel(r'onset step ($H_\perp/H\geq10^{-8}$)'); plt.title(f'N=5 seedless onset vs fixed-point residual (R²={r2:.4f})'); plt.tight_layout(); plt.savefig(HERE/'N5_tol_sweep_onset_log_residual.png',dpi=180); plt.close()
ser=pd.read_csv(HERE/'tol_sweep_timeseries.csv'); plt.figure(figsize=(9,5));
for tol,g in ser.groupby('tol'): plt.semilogy(g.step,np.maximum(g.f,1e-30),label=f'tol={tol:.0e}')
plt.xlabel('step'); plt.ylabel(r'$H_\perp/H$'); plt.title('N=5 seedless: onset shifts with make_parent residual floor'); plt.legend(); plt.tight_layout(); plt.savefig(HERE/'N5_tol_sweep_timeseries.png',dpi=180); plt.close()

# 4) Floquet/Jacobian stability
fl=pd.read_csv(HERE/'floquet_spectrum.csv'); f=fl[np.isclose(fl.fd_eps,1e-7)].copy();
unstable=f[f.modulus>1+1e-5]; floq_summary={'unstable_real_dimensions':int(len(unstable)),'largest_multiplier':float(f.modulus.max()),'largest_log_multiplier':float(f.log_modulus_per_step.max()),'distinct_unstable_moduli':sorted(set(round(x,9) for x in unstable.modulus))}
(HERE/'floquet_summary.json').write_text(json.dumps(floq_summary,indent=2),encoding='utf-8')
plt.figure(figsize=(7,6)); plt.scatter(f.eig_re,f.eig_im); th=np.linspace(0,2*np.pi,400); plt.plot(np.cos(th),np.sin(th),linestyle='--'); plt.axhline(0,linewidth=.5); plt.axvline(0,linewidth=.5); plt.xlabel('Re multiplier'); plt.ylabel('Im multiplier'); plt.title('N=5 rotating-frame one-step Jacobian spectrum'); plt.axis('equal'); plt.tight_layout(); plt.savefig(HERE/'N5_floquet_spectrum.png',dpi=180); plt.close()
# finite-difference stability of leading multipliers
plt.figure(figsize=(8,5));
for rank in [1,2,3,4]:
 g=fl[fl['rank']==rank].sort_values('fd_eps'); plt.semilogx(g.fd_eps,g.modulus,marker='o',label=f'rank {rank}')
plt.xlabel('finite-difference epsilon'); plt.ylabel('|multiplier|'); plt.title('Floquet/Jacobian leading multipliers: finite-difference stability'); plt.legend(); plt.tight_layout(); plt.savefig(HERE/'N5_floquet_fd_stability.png',dpi=180); plt.close()

# 5) N5 moduli seed sweep
ms=pd.read_csv(HERE/'N5_moduli_seed_sweep.csv'); mod_summary={'n_seeds':int(len(ms)),'all_group_counts_3_3_2_2':bool((ms.counts=='3;3;2;2').all()),'family_modulus_mean':float(pd.concat([ms.family1_mod,ms.family2_mod]).mean()),'relative_phase_mod_pi_min':float(ms.relative_phase_mod_pi_rad.min()),'relative_phase_mod_pi_max':float(ms.relative_phase_mod_pi_rad.max()),'relative_phase_mod_pi_std':float(ms.relative_phase_mod_pi_rad.std(ddof=1))}
(HERE/'N5_moduli_seed_sweep_summary.json').write_text(json.dumps(mod_summary,indent=2),encoding='utf-8')
plt.figure(figsize=(8,5)); plt.scatter(ms.seed,ms.relative_phase_mod_pi_rad); plt.axhline(0,linestyle='--'); plt.xlabel('parent random seed index'); plt.ylabel(r'$\arg(v_B/v_A)$ mod $\pi$ [rad]'); plt.title('N=5: relative phase between the two distance families across seeds'); plt.tight_layout(); plt.savefig(HERE/'N5_moduli_relative_phase_seed_sweep.png',dpi=180); plt.close()

# 6) spectral entropy from existing N5 simplex time series
long=pd.read_csv(ROOT/'N5_complex_simplex_complete_analysis_20260826/N5_all_steps_a_b_a2_b2_ab.csv',usecols=['step','r2'])
ent=[]
for step,g in long.groupby('step'):
 p=g.r2.to_numpy(float); p=p/p.sum(); S=-np.sum(np.where(p>0,p*np.log(p),0)); ent.append((step,S,S/math.log(len(p)),p.min(),p.max()))
ent=pd.DataFrame(ent,columns=['step','entropy','entropy_over_lnM','p_min','p_max']); ent.to_csv(HERE/'N5_spectral_entropy_timeseries.csv',index=False)
plt.figure(figsize=(9,5)); plt.plot(ent.step,ent.entropy_over_lnM); plt.xlabel('step'); plt.ylabel(r'$S/\ln M$'); plt.title('N=5 relation-amplitude spectral entropy'); plt.tight_layout(); plt.savefig(HERE/'N5_spectral_entropy.png',dpi=180); plt.close()
summary={'pump_depletion':pump_summary,'time_reparameterization':phase_summary,'tol_sweep':tol_summary,'floquet':floq_summary,'moduli':mod_summary,'entropy_final_over_lnM':float(ent.entropy_over_lnM.iloc[-1])}
(HERE/'FOLLOWUP_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
