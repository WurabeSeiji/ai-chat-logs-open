# -*- coding: utf-8 -*-
"""第十四思考実験 E1 図化：te14_E1_orbits.csv から te14_E1_fig.png / te14_E1_fig.svg を作る"""
import csv, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from te14_E1_lepton_phase_lattice import to_xy, round_a2, norm
rows = list(csv.DictReader(open('te14_E1_orbits.csv')))
MS = (5, 8, 10); COL = {5: 'C0', 8: 'C1', 10: 'C2'}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
# (a) rho - 1/m vs initial radius
for m in MS:
    r = [x for x in rows if int(x['m']) == m]
    ax[0].scatter([math.sqrt(int(x['N0'])) for x in r], [float(x['rho_minus_1_over_m']) for x in r],
                  s=9, color=COL[m], label=f'm={m}', alpha=.7)
ax[0].axhline(0, color='k', lw=.6)
ax[0].set_xlabel('initial radius |z0|'); ax[0].set_ylabel('rho - 1/m  (turns/step)')
ax[0].set_title('(a) effective rotation number vs 1/m'); ax[0].legend()
# (b) fraction with rho != 1/m by radius bin
bins = list(range(0, 45, 5))
for m in MS:
    r = [x for x in rows if int(x['m']) == m]; fr = []; cen = []
    for lo in bins[:-1]:
        s = [x for x in r if lo < math.sqrt(int(x['N0'])) <= lo+5]
        if s:
            fr.append(sum(1 for x in s if float(x['rho_minus_1_over_m']) != 0)/len(s)); cen.append(lo+2.5)
    ax[1].plot(cen, fr, 'o-', color=COL[m], label=f'm={m}')
ax[1].set_xlabel('initial radius bin'); ax[1].set_ylabel('fraction rho != 1/m')
ax[1].set_title('(b) non-locked fraction'); ax[1].legend()
# (c) accumulated phase slip for representative cycles (m=5)
def cyc(a, b, m):
    c, s = math.cos(2*math.pi/m), math.sin(2*math.pi/m); seen = {}; path = []
    while (a, b) not in seen:
        seen[(a, b)] = len(path); path.append((a, b)); x, y = to_xy(a, b); a, b, _ = round_a2(c*x-s*y, s*x+c*y)
    return path[seen[(a, b)]:]
for (a, b, lab) in [(9, 6, 'cycle 11/54'), (10, 0, 'cycle 13/66'), (10, 2, 'cycle 1/5 (locked)')]:
    C = cyc(a, b, 5); P = len(C); ph = 0; ys = [0]
    for i in range(3*P):
        x0, y0 = to_xy(*C[i % P]); x1, y1 = to_xy(*C[(i+1) % P])
        d = (math.atan2(y1, x1)-math.atan2(y0, x0)+math.pi) % (2*math.pi)-math.pi
        ph += d-2*math.pi/5; ys.append(ph/(2*math.pi))
    ax[2].plot(range(len(ys)), ys, lw=1.2, label=lab)
ax[2].set_xlabel('update index k'); ax[2].set_ylabel('accumulated slip (turns)')
ax[2].set_title('(c) phase slip, m=5'); ax[2].legend()
fig.suptitle('TE14-E1: rotation (C_A: 4D-allowed order m) + rounding to A2 lattice (C_B: Z6)')
fig.tight_layout()
for ext in ('png', 'svg'): fig.savefig(f'te14_E1_fig.{ext}', dpi=150)
