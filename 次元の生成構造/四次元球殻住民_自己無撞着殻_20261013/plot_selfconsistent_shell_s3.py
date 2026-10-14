# -*- coding: utf-8 -*-
"""selfconsistent_shell_s3.py の結果（shell3_*.npz, result3_*.json）を図化する。子午線断面 (x, ±rho) を描く"""
import numpy as np, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in font_manager.findSystemFonts():
    if 'NotoSansCJK' in f: font_manager.fontManager.addfont(f)
plt.rcParams['font.family'] = 'Noto Sans CJK JP'; plt.rcParams['svg.fonttype'] = 'none'
NAMES = {'odd': '奇数のみ', 'even': '偶数のみ', 'all': '両方'}
FORMS = {'a': '(a) cos(nψ)', 'b': '(b) U_n(cos ψ)/(n+1)'}

def load(key):
    z = np.load(f'shell3_{key}.npz'); r = json.load(open(f'result3_{key}.json'))
    return z['R']*np.cos(z['chi']), z['R']*np.sin(z['chi']), r, z

rows = [('odd', 17), ('even', 18), ('all', 18), ('all', 4)]
fig, axs = plt.subplots(len(rows), 2, figsize=(9, 4.0*len(rows)))
for i, (s, N) in enumerate(rows):
    for j, f in enumerate('ab'):
        ax = axs[i, j]; x, rho, r, z = load(f'{f}_{s}_N{N}')
        ax.plot(x, rho, color='#1D9E75', lw=1.4); ax.plot(x, -rho, color='#1D9E75', lw=1.4)
        ax.plot(0, 0, '+', color='#444441', ms=11, mew=1.5)
        ax.plot(r['d'], 0, 'o', color='#D85A30', ms=5)
        ax.plot(x[0], 0, 'v', color='#7F77DD', ms=7)
        ax.set_aspect('equal'); ax.grid(alpha=0.25)
        st = '' if r['converged'] else '（未収束）'
        ax.set_title(f'{FORMS[f]}　{NAMES[s]} N={N}{st}\nd = {r["d"]:+.4f}', fontsize=10)
fig.suptitle('S³ 自己無撞着殻の子午線断面（右端▼がピーク、＋が R の中心、●が重心）', fontsize=11)
fig.tight_layout(); fig.savefig('fig1_s3_meridian_shells.png', dpi=150); fig.savefig('fig1_s3_meridian_shells.svg')

fig, axs = plt.subplots(1, 2, figsize=(10, 4))
cols = {'odd': '#7F77DD', 'even': '#D85A30', 'all': '#1D9E75'}
for j, f in enumerate('ab'):
    for s, N in (('odd', 17), ('even', 18), ('all', 18)):
        x, rho, r, z = load(f'{f}_{s}_N{N}'); o = np.argsort(z['psi'])
        axs[j].plot(np.degrees(z['psi'][o]), z['R'][o], color=cols[s], lw=1.3, label=f'{NAMES[s]} N={N}')
    axs[j].set_xlim(0, 180); axs[j].set_xlabel('住民位相 ψ（度、ピーク 0 → 対蹠点 180）'); axs[j].set_ylabel('R')
    axs[j].set_title(FORMS[f], fontsize=11); axs[j].grid(alpha=0.25); axs[j].legend(fontsize=8)
fig.tight_layout(); fig.savefig('fig2_s3_R_vs_psi.png', dpi=150); fig.savefig('fig2_s3_R_vs_psi.svg')
print('ok')
