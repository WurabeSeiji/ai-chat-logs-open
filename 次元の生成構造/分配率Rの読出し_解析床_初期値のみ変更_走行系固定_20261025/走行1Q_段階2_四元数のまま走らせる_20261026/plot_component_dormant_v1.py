#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""走行1Q 段階2（2026-10-26）：成分ごとの Π の外の割合 f_1, f_i, f_j, f_k を timeseries CSV から描く（読み出しのみ）。
グリッドは走行1 の図と同じ 8×5（N＝3..40、空パネルあり）、semilogy、ylim 1e-34..3。den ごとに 1 枚。
使い方: python plot_component_dormant_v1.py RESULTS_DIR TITLE"""
import csv, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rdir, title = sys.argv[1], sys.argv[2]
by = {}
with open(os.path.join(rdir, 'timeseries_64bit_with124_N3_N40.csv')) as fh:
    for r in csv.DictReader(fh):
        key = (int(r['N']), r['series'])
        by.setdefault(key, {'t': [], '1': [], 'i': [], 'j': [], 'k': []})
        by[key]['t'].append(int(r['step']))
        for c in '1ijk':
            by[key][c].append(max(float(r['f_' + c]), 1e-40))
for series in ('N', '124'):
    fig, axs = plt.subplots(8, 5, figsize=(20, 24)); axs = axs.ravel()
    for k, N in enumerate(range(3, 41)):
        ax = axs[k]
        if (N, series) in by:
            dd = by[(N, series)]
            for c, col in zip('1ijk', ('tab:gray', 'tab:blue', 'tab:orange', 'tab:green')):
                ax.semilogy(dd['t'], dd[c], label=f'f_{c}', color=col, linewidth=1.0)
            ax.legend(fontsize=7, loc='lower right')
        ax.set_xlim(0, 500); ax.set_ylim(1e-34, 3); ax.set_title(f'N={N}'); ax.grid(alpha=.25)
        if k // 5 == 7: ax.set_xlabel('step')
        if k % 5 == 0: ax.set_ylabel('component outside Pi / H')
    for k in range(38, 40): axs[k].axis('off')
    fig.suptitle(f'{title} — den={"N" if series == "N" else "124"}', y=.995)
    fig.tight_layout()
    out = os.path.join(rdir, f'fig_component_dormant_den{series}.png')
    fig.savefig(out, dpi=150); plt.close(fig)
    print('wrote', out)
