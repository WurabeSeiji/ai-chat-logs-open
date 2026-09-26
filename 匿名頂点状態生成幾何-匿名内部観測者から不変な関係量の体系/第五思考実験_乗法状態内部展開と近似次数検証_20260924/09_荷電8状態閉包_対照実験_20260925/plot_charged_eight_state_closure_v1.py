#!/usr/bin/env python3
from pathlib import Path
import json, math, csv
import numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
OUT=HERE/'figures'
OUT.mkdir(exist_ok=True)

d=np.genfromtxt(HERE/'charged_eight_state_benchmark_sampled.csv',delimiter=',',names=True)
ref=np.genfromtxt(HERE/'reference_raw.csv',delimiter=',',names=True)
cv=np.genfromtxt(HERE/'convergence.csv',delimiter=',',names=True)
summary=json.load(open(HERE/'charged_eight_state_summary.json'))

# 01 orbit overlay
fig,ax=plt.subplots(figsize=(7.2,7.2))
ax.plot(ref['x_rel'][::20],ref['y_rel'][::20],linewidth=1,label='saved analytic reference')
ax.plot(d['x_state'],d['y_state'],linewidth=1,label='8-state generator')
ax.set_aspect('equal',adjustable='box'); ax.set_xlabel('x/M'); ax.set_ylabel('y/M')
ax.set_title('Charged binary orbit: 8-state generator vs saved reference')
ax.legend(); ax.grid(True,alpha=.25); fig.tight_layout()
fig.savefig(OUT/'figure01_orbit_overlay.png',dpi=180); fig.savefig(OUT/'figure01_orbit_overlay.svg'); plt.close(fig)

# 02 radius vs time
fig,ax=plt.subplots(figsize=(8.4,5.2))
ax.plot(ref['t'],ref['r'],linewidth=1,label='saved analytic reference')
ax.plot(d['t_state'],d['r_state'],linewidth=1,label='8-state generator')
ax.set_xlabel('t/M'); ax.set_ylabel('r/M'); ax.set_title('Radius evolution')
ax.legend(); ax.grid(True,alpha=.25); fig.tight_layout()
fig.savefig(OUT/'figure02_radius_vs_time.png',dpi=180); fig.savefig(OUT/'figure02_radius_vs_time.svg'); plt.close(fig)

# 03 residuals against independent/saved reference and closed form
fig,ax=plt.subplots(figsize=(8.4,5.2))
ax.semilogy(d['t_state'],np.maximum(d['abs_r_saved_residual'],1e-18),label='|r - r_saved(t)|')
ax.semilogy(d['t_state'],np.maximum(d['abs_H_closed_residual'],1e-18),label='|H - exp(i phi_closed)|')
ax.semilogy(d['t_state'],np.maximum(d['abs_t_closed_residual'],1e-18),label='|t - t_closed(r)|')
ax.set_xlabel('t/M'); ax.set_ylabel('absolute residual'); ax.set_title('Generator residuals')
ax.legend(); ax.grid(True,alpha=.25); fig.tight_layout()
fig.savefig(OUT/'figure03_reference_residuals.png',dpi=180); fig.savefig(OUT/'figure03_reference_residuals.svg'); plt.close(fig)

# 04 conserved internal states
fig,ax=plt.subplots(figsize=(8.4,5.2))
ax.plot(d['step'],d['N_state'],label='N = eta')
ax.plot(d['step'],d['C_state'],label='C = Z')
ax.plot(d['step'],d['D_state'],label='D = (Delta lambda)^2')
ax.set_xlabel('macro step'); ax.set_ylabel('state value'); ax.set_title('Identity-propagated internal states')
ax.legend(); ax.grid(True,alpha=.25); fig.tight_layout()
fig.savefig(OUT/'figure04_identity_states.png',dpi=180); fig.savefig(OUT/'figure04_identity_states.svg'); plt.close(fig)

# 05 closure counterexamples
labels=[]; vals=[]
with open(HERE/'closure_counterexamples.csv',encoding='utf-8',newline='') as f:
    for row in csv.DictReader(f):
        labels.append(row['scenario'].replace('_','\n'))
        vals.append(abs(float(row['dP_dphi'])))
fig,ax=plt.subplots(figsize=(8.6,5.2))
ax.bar(np.arange(len(vals)),vals)
ax.set_xticks(np.arange(len(vals)),labels)
ax.set_ylabel('|dP/dphi| at identical first-six state'); ax.set_title('Six-/seven-state non-closure counterexamples')
ax.grid(True,axis='y',alpha=.25); fig.tight_layout()
fig.savefig(OUT/'figure05_closure_counterexamples.png',dpi=180); fig.savefig(OUT/'figure05_closure_counterexamples.svg'); plt.close(fig)

# 06 convergence at fixed 100 orbits
fig,ax=plt.subplots(figsize=(8.4,5.2))
ax.loglog(cv['steps_per_orbit'],np.maximum(cv['abs_r_error'],1e-18),marker='o',label='|Delta r| at 100 orbits')
ax.loglog(cv['steps_per_orbit'],np.maximum(cv['abs_t_error'],1e-18),marker='o',label='|Delta t| at 100 orbits')
ax.set_xlabel('steps per orbit'); ax.set_ylabel('absolute error'); ax.set_title('Step-size audit (roundoff-limited regime)')
ax.legend(); ax.grid(True,which='both',alpha=.25); fig.tight_layout()
fig.savefig(OUT/'figure06_step_audit.png',dpi=180); fig.savefig(OUT/'figure06_step_audit.svg'); plt.close(fig)

print(OUT)
