#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""H⊥/H のグリッド図を timeseries CSV から描く（読み出しのみ、2026-10-25）。
走行本体 run_sweep_from_parents_v1.py のコピーは、行単位の編集で原本の suptitle 行にあった fig.savefig(...) を落としていたため、
H⊥/H グリッド図（fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png）が出力されていなかった。走行系（one_step・A・保存データ）には影響なし。
本スクリプトは原本 run_N3_N40_stage123_v1.py の図化部（8×5 グリッド、semilogy、ylim 1e-34..3、den=N と 124）を CSV から再現する。
使い方: python plot_hperp_grid_from_timeseries_v1.py RESULTS_DIR TITLE"""
import csv, os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rdir, title = sys.argv[1], sys.argv[2]
by = {}
with open(os.path.join(rdir, 'timeseries_64bit_with124_N3_N40.csv')) as fh:
    for r in csv.DictReader(fh):
        key = (int(r['N']), r['series'])
        by.setdefault(key, ([], []))
        by[key][0].append(int(r['step'])); by[key][1].append(float(r['Hperp_frac']))
fig, axs = plt.subplots(8, 5, figsize=(20, 24)); axs = axs.ravel(); order = ['N-2', 'N-1', 'N', 'N+1', 'N+2', '124']
for k, N in enumerate(range(3, 41)):
    ax = axs[k]
    for label in order:
        if (N, label) in by:
            x, y = by[(N, label)]; ax.semilogy(x, y, label=('2pi/124' if label == '124' else label), linewidth=1.0)
    ax.set_xlim(0, 500); ax.set_ylim(1e-34, 3); ax.set_title(f'N={N}'); ax.grid(alpha=.25)
    if k // 5 == 7: ax.set_xlabel('step')
    if k % 5 == 0: ax.set_ylabel('Hperp/H')
    if any((N, l) in by for l in order): ax.legend(fontsize=7, loc='lower right')
for k in range(38, 40): axs[k].axis('off')
fig.suptitle(title, y=.995); fig.tight_layout()
out = os.path.join(rdir, 'fig_Hperp_denominator_controls_with_124_N3_N40_stage123.png')
fig.savefig(out, dpi=180); plt.close(fig); print('wrote', out)
