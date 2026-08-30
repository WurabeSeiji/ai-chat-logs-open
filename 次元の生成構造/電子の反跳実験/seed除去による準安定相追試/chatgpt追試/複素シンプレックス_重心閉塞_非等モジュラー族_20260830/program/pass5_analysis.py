#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス5：走行結果の集計（読出し専用）。各親について予測（parents_predictions.csv）と実測（data/<tag>/summary.json,
timeseries）を並べる：max f, final f, t50（f≥0.5 に初到達する step）, 実測 λ（ln H⊥ の 1e-10<H⊥<1e-3 区間の傾き）,
f(0) 床, 最終状態の |z| 分布（min/max/std/PR）, 閉塞 |ZᵀZ|/H の最終値, H_total ドリフト。
出力：results/dynamics_summary.csv"""
import os, csv, json, sys, math
import numpy as np
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
pred={r['tag']:r for r in csv.DictReader(open(os.path.join(ROOT,'results','parents_predictions.csv')))}
# 乱数生成の均衡親（パス7）も同じ集計に載せる（予測列名を揃える）
rp=os.path.join(ROOT,'results','balanced_random_parents.csv')
if os.path.exists(rp):
    for r in csv.DictReader(open(rp)):
        pred[r['tag']]=dict(tag=r['tag'],N=r['N'],eps='random',k=r['seed'],pred_kind=r['pred_band'],pred_lambda_f=r['lambda_f'],pred_t50=(f"{(math.log(0.5)+72.6)/float(r['lambda_f']):.0f}" if float(r['lambda_f'])>1e-3 else ''))
rows=[]
for t,p in pred.items():
    dd=os.path.join(ROOT,'data',t); sj=os.path.join(dd,'summary.json')
    if not os.path.exists(sj): print(f"{t}: 未走行"); continue
    S=json.load(open(sj)); T=S['treatment']; g=T['growth_fit']
    fcsv=os.path.join(dd,'treatment_linear124_amplitude_aware_timeseries.csv'); fcsv=fcsv if os.path.exists(fcsv) else fcsv+'.gz'   # 走行後 gzip 済みでも読める
    a=np.genfromtxt(fcsv,delimiter=',',names=True)
    f=a['H_perp']/a['H_total']; i50=np.where(f>=0.5)[0]; t50=int(a['step'][i50[0]]) if len(i50) else None
    lam=None if g is None else g['slope_ln_Hperp_per_step']
    # 成長区間での対数線形性（R²）と、予測との比
    lam_pred=float(p['pred_lambda_f'])
    # 40000 step 後の f の予測（床 f0 から）：ln f(T)=ln f0+λT（飽和前のみ意味を持つ）
    f0=float(f[0]); lnf_pred=math.log(max(f0,1e-300))+lam_pred*40000
    rows.append(dict(tag=t,N=p['N'],eps=p['eps'],k=p['k'],pred_kind=p['pred_kind'],pred_lambda_f=f'{lam_pred:.5f}',pred_t50=p['pred_t50'],
        f0=f'{f0:.2e}',max_f=f'{T["max_Hperp_fraction"]:.3e}',max_f_step=T['max_Hperp_fraction_step'],final_f=f'{f[-1]:.3e}',t50_measured=t50,
        lambda_measured=('' if lam is None else f'{lam:.5f}'),lambda_ratio_meas_over_pred=('' if lam is None else f'{lam/lam_pred:.4f}'),growth_fit_R2=('' if g is None else f'{g["R2"]:.5f}'),
        t50_pred_from_measured_floor=(f'{(math.log(0.5)-math.log(f0))/lam_pred:.0f}' if lam_pred>1e-9 and f0>0 else ''),
        ln_f_40000_pred=f'{lnf_pred:.1f}',ln_f_40000_measured=f'{math.log(max(f[-1],1e-300)):.1f}',
        final_amp_min=f'{T["final_amp_min"]:.4f}',final_amp_max=f'{T["final_amp_max"]:.4f}',final_amp_std=f'{T["final_amp_std"]:.4f}',final_PR_over_M=f'{T["final_PR_over_M"]:.3f}',
        init_amp_min=f'{a["amp_min"][0]:.4f}',init_amp_max=f'{a["amp_max"][0]:.4f}',init_amp_std=f'{a["amp_std"][0]:.4f}',
        final_abs_ZTZ_over_H=f'{a["abs_ZT_Z"][-1]/a["H_total"][-1]:.1e}',Htotal_drift=f'{T["Htotal_max_abs_drift"]:.1e}'))
    print(f"{t}: 予測 {p['pred_kind']} λ={lam_pred:.5f} | 実測 f0={f0:.1e} max f={T['max_Hperp_fraction']:.3e} final f={f[-1]:.3e} t50={t50} λ={lam} 比={'' if lam is None else f'{lam/lam_pred:.4f}'} | 最終|z|∈[{T['final_amp_min']:.4f},{T['final_amp_max']:.4f}] std {T['final_amp_std']:.4f}（初期 std {a['amp_std'][0]:.4f}）PR/M={T['final_PR_over_M']:.3f}")
if rows:
    with open(os.path.join(ROOT,'results','dynamics_summary.csv'),'w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("PASS5 OK")
