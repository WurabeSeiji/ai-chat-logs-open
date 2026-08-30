#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""パス5：走行結果の集計（読出し専用）。親ごとに予測（parents_predictions.csv）と実測（data/<tag>/summary.json, timeseries）を並べる。
実測分類：saturated（max f ≥ 0.5）／growing（指数窓 fit あり・max f > 1e-10・未飽和）／floor（max f < 1e-10）。
予測との照合：inflating ⇔ saturated or growing、neutral ⇔ floor。
出力：results/dynamics_summary.csv、results/matrix_N_by_method.md"""
import os, csv, json, math
import numpy as np
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
pred={r['tag']:r for r in csv.DictReader(open(os.path.join(ROOT,'results','parents_predictions.csv')))}
def load_ts(tag):
    f=os.path.join(ROOT,'data',tag,'treatment_linear124_amplitude_aware_timeseries.csv'); f=f if os.path.exists(f) else f+'.gz'
    return np.genfromtxt(f,delimiter=',',names=True) if os.path.exists(f) else None
rows=[]
for tag,p in pred.items():
    sj=os.path.join(ROOT,'data',tag,'summary.json')
    if not os.path.exists(sj): print(f"{tag}: 未走行"); continue
    S=json.load(open(sj)); T=S['treatment']; g=T['growth_fit']; a=load_ts(tag)
    f=a['H_perp']/a['H_total']; i50=np.where(f>=0.5)[0]; t50=int(a['step'][i50[0]]) if len(i50) else None
    maxf=float(np.nanmax(f)); lam=None if g is None else g['slope_ln_Hperp_per_step']
    if maxf>=0.5: mk='saturated'
    elif lam is not None and lam>1e-3 and maxf>1e-10: mk='growing'
    else: mk='floor'
    lp=float(p['pred_lambda_f']); pk=p['pred_kind']
    agree=(pk=='inflating' and mk in ('saturated','growing')) or (pk=='neutral' and mk=='floor')
    # 後半 1/4 の対数傾き（床の成長形の記録用）
    q=a['step']>=30000; ls=float(np.polyfit(a['step'][q],np.log(np.maximum(f[q],1e-300)),1)[0])
    ll=float(np.polyfit(np.log(a['step'][1000:]),np.log(np.maximum(f[1000:],1e-300)),1)[0])
    rows.append(dict(tag=tag,N=int(p['N']),method=p['method'],design=p['design'],norm=f"{float(p['norm']):.6f}",amp_spread_rel=f"{float(p['amp_spread_rel']):.3e}",local_closed=p['local_closed'],roundness=f"{float(p['roundness']):.3f}",
        pred_rho_minus_1=f"{float(p['pred_rho_minus_1']):.3e}",pred_lambda_f=f"{lp:.5f}",pred_kind=pk,pred_t50=(f"{float(p['pred_t50']):.0f}" if p['pred_t50'] not in ('','None') else ''),
        f0=f"{f[0]:.2e}",max_f=f"{maxf:.3e}",max_f_step=int(a['step'][np.nanargmax(f)]),final_f=f"{f[-1]:.3e}",t50_measured=t50,
        lambda_measured=('' if lam is None else f"{lam:.5f}"),lambda_ratio=('' if lam is None or lp<1e-6 else f"{lam/lp:.4f}"),fit_R2=('' if g is None else f"{g['R2']:.5f}"),
        late_log_slope_per_step=f"{ls:.2e}",loglog_slope_1000_40000=f"{ll:.2f}",measured_kind=mk,prediction_agrees=agree,
        final_PR_over_M=f"{T['final_PR_over_M']:.3f}",final_amp_min=f"{T['final_amp_min']:.4f}",final_amp_max=f"{T['final_amp_max']:.4f}",final_amp_std=f"{T['final_amp_std']:.4f}",
        init_amp_std=f"{a['amp_std'][0]:.4f}",final_abs_ZTZ_over_H=f"{a['abs_ZT_Z'][-1]/a['H_total'][-1]:.1e}",Htotal_drift=f"{T['Htotal_max_abs_drift']:.1e}"))
    print(f"{tag}: 予測 {pk} λ={lp:.5f} | 実測 {mk} max f={maxf:.2e} t50={t50} λ={lam} 比={rows[-1]['lambda_ratio']} PR/M={T['final_PR_over_M']:.3f} {'OK' if agree else '**不一致**'}")
with open(os.path.join(ROOT,'results','dynamics_summary.csv'),'w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
# N × 生成法の表
meth=['mp','hm','ne','rb']; name={'mp':'make_parent 等モジュラー','hm':'手作り等モジュラー','ne':'非等モジュラー','rb':'乱数均衡'}
def cell(r):
    if r is None: return '—'
    s={'saturated':'飽和','growing':'成長中','floor':'床'}[r['measured_kind']]
    if r['measured_kind']=='saturated': s+=f" t50={r['t50_measured']}"
    if r['lambda_measured']: s+=f" λ={float(r['lambda_measured']):.4f}"
    s+=f" (予 {float(r['pred_lambda_f']):.4f} {'○' if r['prediction_agrees'] else '×'})"
    return s
by={(r['method'],r['N']):r for r in rows}
md='# N × 生成法 の結果行列（40000 step、L=124、種なし、直接読出し）\n\n各セル：実測分類、t50、実測 λ（/step）、（予測 λ_f、予測との一致 ○/×）。予測は走行前に固定（`parents_predictions.csv`）。\n\n| N | '+' | '.join(name[m] for m in meth)+' |\n|---|'+'---|'*len(meth)+'\n'
for N in range(3,17): md+=f'| {N} | '+' | '.join(cell(by.get((m,N))) for m in meth)+' |\n'
n_ok=sum(1 for r in rows if r['prediction_agrees']); md+=f'\n予測一致：{n_ok}/{len(rows)}。\n'
open(os.path.join(ROOT,'results','matrix_N_by_method.md'),'w').write(md); print(md); print("PASS5 OK")
