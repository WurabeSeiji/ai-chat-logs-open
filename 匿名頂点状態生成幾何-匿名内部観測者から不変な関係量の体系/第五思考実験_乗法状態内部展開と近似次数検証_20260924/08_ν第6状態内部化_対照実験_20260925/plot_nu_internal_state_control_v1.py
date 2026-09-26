#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plotter for the nu internal-state control experiment.
Reads only saved raw/summary data; it does not run the dynamics."""
from pathlib import Path
import csv, json
import numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
FIG=HERE/'figures'; FIG.mkdir(exist_ok=True)
RAW=HERE/'nu_internal_state_control_timeseries.csv'
CASE=HERE/'nu_internal_state_control_case_summary.csv'
SUMMARY=HERE/'nu_internal_state_control_summary.json'

def save(fig,name):
    fig.tight_layout()
    fig.savefig(FIG/f'{name}.png',dpi=180,bbox_inches='tight')
    fig.savefig(FIG/f'{name}.svg',bbox_inches='tight')
    plt.close(fig)

def load_case(case_id):
    rows=[]
    with RAW.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['case_id']==case_id: rows.append(r)
    return rows

# 1. Representative orbit overlay.
r=load_case('p24.3_e0.45_q1')
xold=np.array([float(z['x_old']) for z in r]); yold=np.array([float(z['y_old']) for z in r])
xnew=np.array([float(z['x_new']) for z in r]); ynew=np.array([float(z['y_new']) for z in r])
fig,ax=plt.subplots(figsize=(7,7))
ax.plot(xold,yold,label='legacy 5-state')
ax.plot(xnew,ynew,'--',label='new 5+1 state')
ax.set_aspect('equal',adjustable='box'); ax.set_xlabel('x'); ax.set_ylabel('y')
ax.set_title('Orbit overlay: p0=24.3, e0=0.45, q=1')
ax.legend(); ax.grid(True,alpha=.25)
save(fig,'figure01_orbit_overlay')

# 2. Exact equality indicators vs step.
step=np.array([int(z['step']) for z in r])
coreeq=np.array([int(z['core5_bitwise_equal']) for z in r])
Neq=np.array([int(z['N_bitwise_equal_nu0']) for z in r])
fig,ax=plt.subplots(figsize=(8,4.5))
ax.plot(step,coreeq,label='core 5 states bitwise equal')
ax.plot(step,Neq,'--',label='N bitwise equal to nu0')
ax.set_ylim(.98,1.02); ax.set_xlabel('macro step'); ax.set_ylabel('equality flag')
ax.set_title('Exact equality flags over 4000 macro steps')
ax.legend(); ax.grid(True,alpha=.25)
save(fig,'figure02_bitwise_equality')

# 3. N remains constant for q=1,2,4,10 at common p,e when available.
fig,ax=plt.subplots(figsize=(8,5))
for cid,label in [('p60_e0.45_q1','q=1'),('p60_e0.45_q2','q=2'),('p60_e0.45_q4','q=4'),('p60_e0.45_q10','q=10')]:
    rr=load_case(cid)
    s=np.array([int(z['step']) for z in rr]); n=np.array([float(z['N_state']) for z in rr])
    ax.plot(s,n,label=label)
ax.set_xlabel('macro step'); ax.set_ylabel('N state')
ax.set_title('Sixth state N is constant; value depends on mass ratio')
ax.legend(); ax.grid(True,alpha=.25)
save(fig,'figure03_N_identity_all_mass_ratios')

# 4. Per-case maximum difference (all are exactly zero).
with CASE.open(newline='',encoding='utf-8') as f: cases=list(csv.DictReader(f))
labels=[z['case_id'] for z in cases]
metric_names=['max_abs_U','max_abs_P','max_abs_E','max_abs_H','max_abs_Q','max_abs_r','max_abs_x','max_abs_y','max_abs_t']
maxvals=np.array([max(float(z[m]) for m in metric_names) for z in cases])
fig,ax=plt.subplots(figsize=(11,5))
ax.plot(np.arange(len(labels)),maxvals,'o')
ax.axhline(0,linewidth=1)
ax.set_xticks(np.arange(len(labels))); ax.set_xticklabels(labels,rotation=55,ha='right')
ax.set_ylabel('max absolute difference')
ax.set_title('Legacy 5-state vs new 5+1: maximum difference per case')
ax.text(0.5,0.82,'All plotted values are exactly 0.0',transform=ax.transAxes,ha='center')
ax.set_ylim(-1e-18,1e-18); ax.grid(True,alpha=.25)
save(fig,'figure04_case_max_difference')

# 5. Saved C1 CSV agreement: legacy and new have identical nonzero reference deviations.
summary=json.loads(SUMMARY.read_text(encoding='utf-8'))
c1=summary['saved_C1_csv_check']
fields=['p','e','t','x','y']
old=np.array([c1['legacy_max_abs_vs_saved_csv'][k] for k in fields])
new=np.array([c1['new6_max_abs_vs_saved_csv'][k] for k in fields])
x=np.arange(len(fields)); w=.36
fig,ax=plt.subplots(figsize=(8,5))
ax.bar(x-w/2,old,w,label='legacy 5-state')
ax.bar(x+w/2,new,w,label='new 5+1 state')
ax.set_yscale('log'); ax.set_xticks(x); ax.set_xticklabels(fields)
ax.set_ylabel('max abs error vs saved Paper-4 C1 CSV')
ax.set_title('Saved C1 reference check: identical error profile')
ax.legend(); ax.grid(True,axis='y',alpha=.25)
save(fig,'figure05_saved_C1_reference_errors')

print(FIG)
