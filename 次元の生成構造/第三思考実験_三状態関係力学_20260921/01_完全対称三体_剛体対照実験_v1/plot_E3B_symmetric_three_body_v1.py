#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot E3-B symmetric three-body results from saved CSV only."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
TRAJ=HERE/'E3B0_symmetric_three_body_trajectories_v1.csv'
SUMMARY=HERE/'E3B0_symmetric_three_body_summary_v1.csv'

def load_csv(p):
    with open(p,encoding='utf-8') as f:return list(csv.DictReader(f))
tr=load_csv(TRAJ); su=load_csv(SUMMARY)
rep=0.01
rr=[r for r in tr if abs(float(r['eps'])-rep)<1e-15]

# 1. State-space rigid triangle snapshots
plt.figure(figsize=(8,8))
for k in [0,4,8,12]:
    r=rr[k]
    xs=[float(r['a_x']),float(r['b_x']),float(r['c_x']),float(r['a_x'])]
    ys=[float(r['a_y']),float(r['b_y']),float(r['c_y']),float(r['a_y'])]
    plt.plot(xs,ys,marker='o',label=f'k={k}')
plt.gca().set_aspect('equal',adjustable='box')
plt.xlabel('state coordinate 1'); plt.ylabel('state coordinate 2')
plt.title('E3-B symmetric three-body rigid rotation (eps=0.01)')
plt.legend(); plt.tight_layout(); plt.savefig(HERE/'fig_E3B_01_rigid_triangle_snapshots_v1.png',dpi=180); plt.close()

# 2. Quadratic readout trajectories over one full X period
plt.figure(figsize=(8,8))
period=[r for r in rr if int(r['step'])<=31]
for body in ['a','b','c']:
    x=[float(r[f'{body}_readout_x']) for r in period]
    y=[float(r[f'{body}_readout_y']) for r in period]
    plt.plot(x,y,marker='.',markersize=3,label=body)
plt.gca().set_aspect('equal',adjustable='box')
plt.xlabel('readout x'); plt.ylabel('readout y')
plt.title('E3-B quadratic readout trajectories: same orbit family, phase-shifted')
plt.legend(); plt.tight_layout(); plt.savefig(HERE/'fig_E3B_02_quadratic_readout_orbits_v1.png',dpi=180); plt.close()

# 3. Rigidity: relative pair-distance drift
plt.figure(figsize=(9,5.5))
short=[r for r in rr if int(r['step'])<=124]
steps=[int(r['step']) for r in short]
for key,label in [('d_ab','ab'),('d_bc','bc'),('d_ca','ca')]:
    d=np.array([float(r[key]) for r in short])
    rel=d/d[0]-1.0
    plt.plot(steps,rel,label=label)
plt.xlabel('step'); plt.ylabel('relative change from initial pair distance')
plt.title('E3-B rigid-body check: all three pair distances stay fixed')
plt.legend(); plt.tight_layout(); plt.savefig(HERE/'fig_E3B_03_pair_distances_v1.png',dpi=180); plt.close()

# 4. Symmetry error vs coupling
plt.figure(figsize=(8,5.5))
eps=[float(r['eps']) for r in su]
series=[
 ('centroid_rel_max','centroid closure'),
 ('pair_distance_equality_error','pair equality'),
 ('triangle_area_rel_drift','area drift'),
 ('body_radius_split','body radius split'),
 ('body_q_split','body q split'),
]
for key,label in series:
    vals=[max(float(r[key]),1e-18) for r in su]
    plt.semilogy(eps,vals,marker='o',label=label)
plt.xlabel('epsilon'); plt.ylabel('relative error / split')
plt.title('E3-B symmetry residuals remain at numerical precision')
plt.legend(); plt.tight_layout(); plt.savefig(HERE/'fig_E3B_04_symmetry_residuals_v1.png',dpi=180); plt.close()

# 5. Mass-like ratios: absolute deviation from unity
plt.figure(figsize=(8,5.5))
nz=[r for r in su if float(r['eps'])>0]
eps2=[float(r['eps']) for r in nz]
for key,label in [('mass_ratio_ab','ab / ba'),('mass_ratio_bc','bc / cb'),('mass_ratio_ca','ca / ac')]:
    dev=[max(abs(float(r[key])-1.0),1e-18) for r in nz]
    plt.loglog(eps2,dev,marker='o',label=label)
plt.axhline(1e-10,linestyle='--',label='predeclared tolerance')
plt.xlabel('epsilon'); plt.ylabel('|directed coupling ratio - 1|')
plt.title('E3-B mass-like readout: no nontrivial ratio appears')
plt.legend(); plt.tight_layout(); plt.savefig(HERE/'fig_E3B_05_masslike_ratios_v1.png',dpi=180); plt.close()

print('plots generated')
