#!/usr/bin/env python3
"""Independent algebraic/numerical checks for the analytic charged-binary benchmark."""
from __future__ import annotations
import csv, json, math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
CSV=ROOT/'charged_binary_analytic_reference_v1_raw.csv'
META=ROOT/'charged_binary_analytic_reference_v1_metadata.json'

with META.open(encoding='utf-8') as f:
    meta=json.load(f)
p=meta['initial_conditions']; c=meta['analytic_coefficients']
A=float(c['A_em_in_dr_dt']); B=float(c['B_gw_in_dr_dt']); Z=float(p['Z'])
mu=float(p['mu']); dl=float(p['delta_lambda'])

cols={}
with CSV.open(newline='',encoding='utf-8') as f:
    rd=csv.DictReader(f)
    rows=list(rd)
for k in rows[0]:
    if k=='index': continue
    cols[k]=np.array([float(r[k]) for r in rows])
r=cols['r']; t=cols['t']; phi=cols['phi']; om=cols['omega_lo']; dr=cols['dr_dt']

# Check differential identities away from boundaries.
dt_dr=np.gradient(t,r,edge_order=2)
dphi_dr=np.gradient(phi,r,edge_order=2)
expected_dt_dr=-r**3/(A*r+B)
expected_dphi_dr=-np.sqrt(Z)*r**1.5/(A*r+B)
rel_dt=np.max(np.abs((dt_dr[2:-2]-expected_dt_dr[2:-2])/expected_dt_dr[2:-2]))
rel_dphi=np.max(np.abs((dphi_dr[2:-2]-expected_dphi_dr[2:-2])/expected_dphi_dr[2:-2]))

# Check energy-balance constructed radial speed from powers.
dEdr=mu*Z/(2*r**2)
dr_from_power=-(cols['P_total']/dEdr)
rel_dr=np.max(np.abs((dr_from_power-dr)/dr))

# Check Kepler relation and frequency harmonics.
kepler=np.max(np.abs(om**2*r**3/Z-1))
frel=np.max(np.abs(cols['f_gw_quadrupole']/(2*cols['f_orb'])-1))
ferel=np.max(np.abs(cols['f_em_dipole']/cols['f_orb']-1))

out={
 'max_relative_dt_dr_error':float(rel_dt),
 'max_relative_dphi_dr_error':float(rel_dphi),
 'max_relative_dr_dt_from_energy_balance_error':float(rel_dr),
 'max_abs_kepler_identity_error':float(kepler),
 'max_abs_fgw_over_2forb_minus1':float(frel),
 'max_abs_fem_over_forb_minus1':float(ferel),
 'max_abs_energy_balance_residual':float(np.max(np.abs(cols['E_balance_residual']))),
}
path=ROOT/'validation_results.json'
with path.open('w',encoding='utf-8') as f:
    json.dump(out,f,indent=2); f.write('\n')
print(json.dumps(out,indent=2))
