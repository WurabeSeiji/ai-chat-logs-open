# -*- coding: utf-8 -*-
"""第十四思考実験 E2 図化：te14_E2_radius_scan.csv から te14_E2_fig.png / .svg"""
import csv, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows = [r for r in csv.DictReader(open('te14_E2_radius_scan.csv')) if r['P'] != 'NA']
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
targets = {'1/6': 1/6, '3/16': 3/16, '5/27': 5/27, '13/70': 13/70, '18/97': 18/97}
for axn, col, mk in (('Z4', 'C0', 'o'), ('Z6', 'C3', 's')):
    s = [r for r in rows if r['axis'] == axn]
    ax[0].scatter([int(r['R'])*(1.03 if axn == 'Z6' else 1) for r in s], [float(r['rho']) for r in s],
                  s=14, color=col, marker=mk, alpha=.6, label=f'axis {axn}')
    ai = [(int(r['R'])*(1.03 if axn == 'Z6' else 1), float(r['alpha_read_inv'])) for r in s]
    ax[1].scatter([a for a, b in ai], [b for a, b in ai], s=14, color=col, marker=mk, alpha=.6, label=f'axis {axn}')
ax[0].axhline(23/124, color='k', lw=1, label='23/124')
for k, v in targets.items():
    ax[0].axhline(v, color='gray', lw=.6, ls=':'); ax[0].text(66, v, k, fontsize=8, va='center')
ax[0].set_xscale('log'); ax[0].set_ylim(0.16, 0.19)
ax[0].set_xlabel('initial radius'); ax[0].set_ylabel('effective rotation rho (turns/step)')
ax[0].set_title('(a) capture by convergents of 23/124'); ax[0].legend(fontsize=8, loc='lower right')
ax[1].axhspan(137.0, 137.036042914599, color='y', alpha=.25, label='window [137, x_23/124]')
ax[1].axhline(137.035999177, color='k', lw=.8, label='CODATA 2022 (reference only)')
ax[1].axhline(137.036042914599, color='gray', lw=.8, ls='--', label='x_23/124')
ax[1].set_xscale('log'); ax[1].set_ylim(136.995, 137.045)
ax[1].set_xlabel('initial radius'); ax[1].set_ylabel('alpha_read^-1 = 4pi/sin^4(pi rho)')
ax[1].set_title('(b) K4 read-out near the window'); ax[1].legend(fontsize=7, loc='lower right')
fig.suptitle('TE14-E2: Hurwitz lattice (4D) + rotation by lambda_a of R_{124,23}')
fig.tight_layout()
for ext in ('png', 'svg'): fig.savefig(f'te14_E2_fig.{ext}', dpi=150)
