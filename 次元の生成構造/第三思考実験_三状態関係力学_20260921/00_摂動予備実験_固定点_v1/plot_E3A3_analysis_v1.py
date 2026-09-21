#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot E3-A3 analytic figures from fixed/derived CSV data only."""
from pathlib import Path
import csv
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def read_csv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

# Figure 08: exact local-area exchange.
rows = read_csv(HERE/'E3A3_area_exchange_timeseries_v1.csv')
step = np.array([int(r['step']) for r in rows])
dqi = np.array([float(r['dq_in_over_Q']) for r in rows])
dqo = np.array([float(r['dq_out_over_Q']) for r in rows])
dqt = np.array([float(r['dq_total_over_Q']) for r in rows])
plt.figure(figsize=(9,6))
plt.plot(step, dqi, label='inner local area change')
plt.plot(step, dqo, label='outer local area change')
plt.plot(step, dqt, label='total area change')
plt.axhline(0, linewidth=0.8)
plt.xlabel('step k')
plt.ylabel('change / |Q_total(0)|')
plt.title('E3-A3: exact local-area exchange at epsilon=0.004')
plt.legend()
plt.tight_layout()
plt.savefig(HERE/'fig_E3A_08_local_area_exchange_timeseries_v1.png', dpi=180)
plt.close()

# Figure 09: phase robustness from eight fixed phases.
pr = read_csv(HERE/'E3A3_phase_robustness_v1.csv')
plt.figure(figsize=(9,6))
for metric, label in [('q_out_mod','outer local-area modulation'),('raw_outer_conic_resid','raw outer conic residual')]:
    rr = sorted([r for r in pr if r['metric']==metric], key=lambda x: float(x['eps']))
    x = np.array([float(r['eps']) for r in rr])
    lo = np.array([float(r['min']) for r in rr])
    med = np.array([float(r['median']) for r in rr])
    hi = np.array([float(r['max']) for r in rr])
    line, = plt.plot(x, med, marker='o', label=label+' median')
    plt.fill_between(x, lo, hi, alpha=0.16)
# O(eps) guide normalized to first conic median.
conic = sorted([r for r in pr if r['metric']=='raw_outer_conic_resid'], key=lambda x: float(x['eps']))
x0=float(conic[0]['eps']); y0=float(conic[0]['median'])
xguide=np.array([float(r['eps']) for r in conic])
plt.plot(xguide, y0*xguide/x0, linestyle='--', linewidth=1.2, label='O(epsilon) guide')
plt.xscale('log'); plt.yscale('log')
plt.xlabel('coupling epsilon')
plt.ylabel('dimensionless response')
plt.title('E3-A3: robustness across 8 initial relative phases')
plt.legend()
plt.tight_layout()
plt.savefig(HERE/'fig_E3A_09_phase_robustness_v1.png', dpi=180)
plt.close()

# Figure 10: universal control by g = epsilon / |delta_tau|.
scan = read_csv(HERE/'E3A1_perturbation_breakdown_scan_v1.csv')
gpts = np.array([float(r['eps_over_delta']) for r in scan if float(r['eps'])>0])
rpts = np.array([float(r['shift_Oeps2_rel_error']) for r in scan if float(r['eps'])>0])
wpts = np.array([float(r['inner_weight_in_outer_mode']) for r in scan if float(r['eps'])>0])
g = np.logspace(-3, 1, 500)
R = np.sqrt(0.25 + g*g) - 0.5
W = 0.5*(1 - 1/np.sqrt(1+4*g*g))
plt.figure(figsize=(9,6))
plt.plot(g, R, label='relative error of O(epsilon^2) law shift')
plt.plot(g, W, label='inner weight in outer normal mode')
plt.scatter(gpts, rpts, s=18, label='fixed scan: shift error')
plt.scatter(gpts, wpts, s=18, label='fixed scan: inner weight')
for pct in (1,5,10):
    r=pct/100; gt=math.sqrt(r*(1+r))
    plt.axvline(gt, linewidth=0.8, linestyle=':')
    plt.text(gt, r*1.3, f'{pct}% error', rotation=90, va='bottom', ha='right', fontsize=8)
plt.axvline(0.5, linewidth=0.8, linestyle='--', label='2 epsilon / delta_tau = 1')
plt.xscale('log'); plt.yscale('log')
plt.xlabel('g = epsilon / |delta_tau|')
plt.ylabel('dimensionless fraction')
plt.title('E3-A3: perturbative accuracy and mixing are universal functions of g')
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(HERE/'fig_E3A_10_universal_g_control_v1.png', dpi=180)
plt.close()

# Figure 11: degenerate vs nondegenerate hierarchy mixing.
agg = read_csv(HERE/'E3A0_hierarchical_two_kepler_perturbation_aggregates_v1.csv')
plt.figure(figsize=(9,6))
for suite, label in [('hierarchical_nonres_31_127','nondegenerate 31/127'),('degenerate_31_31','degenerate 31/31')]:
    rr=sorted([r for r in agg if r['suite']==suite and float(r['eps'])>0], key=lambda x: float(x['eps']))
    x=np.array([float(r['eps']) for r in rr])
    y=np.degrees(np.array([float(r['mix_angle']) for r in rr]))
    plt.plot(x,y,marker='o',label=label)
plt.axhline(45, linewidth=0.8, linestyle='--', label='maximal 50:50 mixing')
plt.xscale('log')
plt.xlabel('coupling epsilon')
plt.ylabel('mixing angle (deg)')
plt.title('E3-A3: degeneracy destroys the perturbative hierarchy immediately')
plt.legend()
plt.tight_layout()
plt.savefig(HERE/'fig_E3A_11_degenerate_vs_nondegenerate_mixing_v1.png', dpi=180)
plt.close()

print('generated figures 08-11')
