#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot the saved N=40 diagnostic trajectory without rerunning dynamics."""
import csv
import os
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, 'diagnostic_0_200.csv')
OUT = os.path.join(HERE, 'fig_diagnostic_N40.png')

with open(CSV, newline='') as f:
    rows = list(csv.DictReader(f))

step = [int(r['step']) for r in rows]
hperp = [float(r['Hperp_frac']) for r in rows]
onset = next((int(r['step']) for r in rows if float(r['Hperp_frac']) > 0.05), None)

fig = plt.figure(figsize=(10, 5))
plt.semilogy(step, hperp, label='Hperp/H')
plt.axhline(0.05, linestyle='--', label='threshold = 0.05')
if onset is not None:
    plt.axvline(onset, linestyle=':', label=f'onset = {onset}')
plt.xlabel('step')
plt.ylabel('Hperp/H')
plt.title('N=40, Delta tau=2pi/40 diagnostic')
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
fig.savefig(OUT, dpi=180)
plt.close(fig)
print(OUT)
