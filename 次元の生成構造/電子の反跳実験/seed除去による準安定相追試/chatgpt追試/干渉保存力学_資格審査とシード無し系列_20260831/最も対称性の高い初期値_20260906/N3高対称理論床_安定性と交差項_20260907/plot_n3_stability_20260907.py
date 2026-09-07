#!/usr/bin/env python3
import csv, os
import numpy as np
import matplotlib.pyplot as plt

ROOT=os.path.dirname(os.path.abspath(__file__))
CSV=os.path.join(ROOT,'n3_stability_timeseries_20260907.csv')
rows=[]
with open(CSV,encoding='utf-8') as f:
    for r in csv.DictReader(f):
        r['epsilon']=float(r['epsilon']); r['step']=int(r['step'])
        for k in ('orbit_distance','norm2','centroid_abs','square_closure_abs','cross_sum_abs'):
            r[k]=float(r[k])
        rows.append(r)

# Figure 1: distance to the N=3 regular-triangle orbit.
fig,ax=plt.subplots(figsize=(9,6))
for kind,label_prefix in [('closure_preserving','closure-preserving'),('nonclosure_control','closure-broken control')]:
    for eps in sorted({r['epsilon'] for r in rows if r['case']==kind}, reverse=True):
        rr=[r for r in rows if r['case']==kind and r['epsilon']==eps]
        ax.plot([r['step'] for r in rr],[max(r['orbit_distance'],1e-16) for r in rr],label=f'{label_prefix}, eps={eps:g}')
ax.set_yscale('log')
ax.set_xlabel('step')
ax.set_ylabel('distance to regular-triangle orbit (global phase removed)')
ax.set_title('N=3 stability test: closure-preserving perturbations return, closure-broken controls do not')
ax.grid(True,which='both',alpha=0.25)
ax.legend(fontsize=8,ncol=2)
fig.tight_layout()
fig.savefig(os.path.join(ROOT,'fig_n3_stability_distance_20260907.png'),dpi=180)
plt.close(fig)

# Figure 2: cross term and square closure for eps=1e-2.
eps=1e-2
fig,ax=plt.subplots(figsize=(9,6))
for kind,label_prefix in [('closure_preserving','closure-preserving'),('nonclosure_control','closure-broken control')]:
    rr=[r for r in rows if r['case']==kind and r['epsilon']==eps]
    ax.plot([r['step'] for r in rr],[max(r['cross_sum_abs'],1e-18) for r in rr],label=f'|sum i<j z_i z_j|: {label_prefix}')
    ax.plot([r['step'] for r in rr],[max(r['square_closure_abs'],1e-18) for r in rr],linestyle='--',label=f'|sum z_i^2|: {label_prefix}')
ax.set_yscale('log')
ax.set_xlabel('step')
ax.set_ylabel('magnitude')
ax.set_title('N=3 cross term and square-closure response (eps=0.01)')
ax.grid(True,which='both',alpha=0.25)
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(ROOT,'fig_n3_cross_terms_20260907.png'),dpi=180)
plt.close(fig)
