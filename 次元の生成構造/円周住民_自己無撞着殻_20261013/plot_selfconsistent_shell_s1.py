# -*- coding: utf-8 -*-
"""selfconsistent_shell_s1.py の結果（shell_*.npz, result_*.json）を図化する"""
import numpy as np, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in font_manager.findSystemFonts():
    if 'NotoSansCJK' in f: font_manager.fontManager.addfont(f)
plt.rcParams['font.family'] = 'Noto Sans CJK JP'; plt.rcParams['svg.fonttype'] = 'none'

def load(key):
    z = np.load(f'shell_{key}.npz'); r = json.load(open(f'result_{key}.json'))
    P = z['R'][:, None]*np.c_[np.cos(z['theta']), np.sin(z['theta'])]
    return P, np.array(r['d']), r

# 図1: 収束した殻の形（構成中心と重心）
keys = [('odd_N1', '奇数のみ N=1'), ('all_N4', '両方 N=4'), ('all_N18', '両方 N=18'), ('even_N18', '偶数のみ N=18')]
fig, axs = plt.subplots(2, 2, figsize=(9, 9))
for ax, (k, title) in zip(axs.flat, keys):
    P, d, r = load(k)
    ax.plot(*np.vstack([P, P[:1]]).T, color='#1D9E75', lw=1.6)
    ax.plot(0, 0, '+', color='#444441', ms=12, mew=1.5, label='R を測った中心')
    ax.plot(*d, 'o', color='#D85A30', ms=6, label='殻の重心')
    ax.annotate('', xy=d, xytext=(0, 0), arrowprops=dict(arrowstyle='->', color='#D85A30', lw=1.2))
    ax.plot(P[0, 0], P[0, 1], 'v', color='#7F77DD', ms=8, label='ピーク（φ0）')
    ax.set_aspect('equal'); ax.grid(alpha=0.25)
    ax.set_title(f'{title}\nd = {d[0]:+.4f}　d/L = {r["d_over_L"]:+.4f}', fontsize=11)
axs[0, 0].legend(fontsize=8, loc='lower left')
fig.suptitle('自己無撞着殻の収束形（A = 1、R² = Σcos²(n(φ−φ0)/2)、位相は弧長）', fontsize=12)
fig.tight_layout(); fig.savefig('fig1_converged_shells.png', dpi=160); fig.savefig('fig1_converged_shells.svg')

# 図2: 重心から測り直すたびの並進（両方 N=18）
P, d, r = load('all_N18')
fig, ax = plt.subplots(figsize=(9, 5))
for k in range(5):
    Q = P + k*d
    ax.plot(*np.vstack([Q, Q[:1]]).T, color='#1D9E75', lw=1.4, alpha=1-0.17*k)
    ax.plot(*(k*d), '+', color='#444441', ms=9)
ax.annotate('', xy=4*d+np.array([0, -5.2]), xytext=(0, -5.2), arrowprops=dict(arrowstyle='->', color='#D85A30', lw=1.4))
ax.text(2*d[0], -5.6, '重心から測り直すたびに d だけ並進（形は不変）', ha='center', va='top', fontsize=10, color='#D85A30')
ax.set_aspect('equal'); ax.grid(alpha=0.25); ax.set_ylim(-6.3, 5)
ax.set_title(f'両方 N=18：d = {d[0]:.4f}（ピーク方向）', fontsize=11)
fig.tight_layout(); fig.savefig('fig2_translation_all_N18.png', dpi=160); fig.savefig('fig2_translation_all_N18.svg')

# 図3: 住民位相 φ に対する R
fig, ax = plt.subplots(figsize=(9, 4.2))
cols = {'odd_N1': '#7F77DD', 'all_N4': '#BA7517', 'all_N18': '#1D9E75', 'even_N18': '#D85A30'}
for k, title in keys:
    z = np.load(f'shell_{k}.npz'); o = np.argsort(z['phi'])
    ax.plot(np.degrees(z['phi'][o]), z['R'][o], color=cols[k], lw=1.4, label=title)
ax.set_xlabel('住民位相 φ − φ0（度、弧長で測る）'); ax.set_ylabel('R'); ax.set_xlim(0, 360); ax.grid(alpha=0.25); ax.legend(fontsize=9)
ax.set_title('殻の半径 R と住民位相', fontsize=11)
fig.tight_layout(); fig.savefig('fig3_R_vs_phase.png', dpi=160); fig.savefig('fig3_R_vs_phase.svg')
print('ok')
