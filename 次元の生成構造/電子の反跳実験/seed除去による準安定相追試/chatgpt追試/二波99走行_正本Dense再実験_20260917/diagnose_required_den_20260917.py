#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""必要 den の推移診断（2026-09-17 木原指示・read-only 後処理）。

目的: den の設計指標を「実績」から得る。各走行の全 4097 状態について
生成子 K(Z_t) の全固有スペクトルを毎 step 計算（部分固有分解は使わない・§9）、
最大固有値 ω_max(t) から

  必要 den(t) = 2π·ω_max(t) / c     （c = 1step あたりの許容位相回転）

の時系列を求める。c は事前固定で 0.3 rad（分解基準）と 0.1 rad（十分基準）の
2 本を記録する（グラフは 0.3 基準、0.1 は log 軸上で ×3 平行移動）。

出力（diagnostics/required_den/）:
  - omega_spectrum_<rid>.npz … 全 step の全固有スペクトル（float64、保存規律 §9）
  - required_den_<rid>.csv   … t, omega_max, dtau_omega, den_req_c03, den_req_c01
  - fig_required_den_all_runs_20260917.png … 全走行重ね描き（色=L、対数縦軸）
  - progress log は diagnose_required_den_progress.log

一次データ（states.npz）は read-only。力学・カーネルは一切変更しない。
既に npz が states.npz より新しい run はスキップ（増分実行）。
"""
import os

THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
for _v in THREAD_VARS:
    os.environ[_v] = '1'

import csv                                                     # noqa: E402
import math                                                    # noqa: E402
import sys                                                     # noqa: E402
import time                                                    # noqa: E402

import numpy as np                                             # noqa: E402
import matplotlib                                              # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                # noqa: E402

plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from interaction_kernel_theory_v1 import adjacency, K_of       # noqa: E402
from gen_manifest import build_runs, T, MAIN_L                 # noqa: E402

C_RESOLVE = 0.3    # 事前固定（分解基準）
C_FINE = 0.1       # 事前固定（十分基準）
OUTDIR = os.path.join(HERE, 'diagnostics', 'required_den')
LOG = os.path.join(HERE, 'diagnose_required_den_progress.log')


def log_line(msg):
    line = f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())} UTC] {msg}'
    print(line, flush=True)
    with open(LOG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def process_run(row):
    rid = row['run_id']
    src = os.path.join(HERE, 'runs', rid, 'states.npz')
    if not os.path.exists(src):
        return None
    dst = os.path.join(OUTDIR, f'omega_spectrum_{rid}.npz')
    if os.path.exists(dst) and os.path.getmtime(dst) > os.path.getmtime(src):
        with np.load(dst) as d:
            return np.abs(d['spec']).max(axis=1)
    L, den = row['L'], row['den']
    A = adjacency(L)
    with np.load(src) as d:
        Z = np.array(d['Z'])
    assert Z.shape[0] == T + 1
    M = Z.shape[1]
    spec = np.empty((T + 1, M), np.float64)
    t0 = time.time()
    for t in range(T + 1):
        H = (1j * K_of(Z[t], A)).astype(np.complex128, copy=False)
        spec[t] = np.linalg.eigvalsh(H)      # 全固有値（部分固有分解禁止・§9）
    wall = time.time() - t0
    np.savez(dst, spec=spec, L=np.int64(L), den=np.int64(den), T=np.int64(T))
    wmax = np.abs(spec).max(axis=1)
    dtau = 2.0 * math.pi / den
    with open(os.path.join(OUTDIR, f'required_den_{rid}.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['t', 'omega_max', 'dtau_omega',
                    f'den_req_c{C_RESOLVE}', f'den_req_c{C_FINE}'])
        for t in range(T + 1):
            w.writerow([t, repr(wmax[t]), repr(dtau * wmax[t]),
                        repr(2 * math.pi * wmax[t] / C_RESOLVE),
                        repr(2 * math.pi * wmax[t] / C_FINE)])
    log_line(f'{rid}: M={M} wmax(t=0)={wmax[0]:.1f} wmax最大={wmax.max():.1f} '
             f'den必要(0.3rad)最大={2 * math.pi * wmax.max() / C_RESOLVE:.0f} {wall:.1f}s')
    return wmax


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    rows = build_runs()
    series = {}
    log_line(f'diagnose start: {sum(1 for r in rows if os.path.exists(os.path.join(HERE, "runs", r["run_id"], "states.npz")))} run(s) 対象')
    for row in rows:
        wmax = process_run(row)
        if wmax is not None:
            series[row['run_id']] = (row['L'], row['den'], wmax)

    # 図: 必要 den(t) = 2π·ω_max(t)/0.3、全走行重ね描き（色=L、対数縦軸）
    cmap = plt.cm.turbo(np.linspace(0.05, 0.95, len(MAIN_L)))
    lcol = {L: cmap[i] for i, L in enumerate(MAIN_L)}
    fig, ax = plt.subplots(figsize=(13.5, 8.0))
    tt = np.arange(T + 1)
    for rid, (L, den, wmax) in series.items():
        ax.plot(tt, 2 * math.pi * wmax / C_RESOLVE, color=lcol[L], lw=0.6, alpha=0.55)
    for L in MAIN_L:
        if any(v[0] == L for v in series.values()):
            ax.axhline(L, color=lcol[L], lw=0.8, ls=':', alpha=0.9)
    ax.set_yscale('log')
    ax.set_xlim(0, T)
    ax.set_xlabel('step', fontsize=11)
    ax.set_ylabel(f'必要 den(t) = 2π·ω_max(K(Z_t)) / {C_RESOLVE} rad', fontsize=11)
    ax.set_title(f'必要 den の推移（実測 K スペクトルから、c={C_RESOLVE} rad 基準、'
                 f'{len(series)} run）\n実線=各走行（色=L）、点線=現行 den=L の水準。'
                 f'c={C_FINE} rad 基準は縦軸を×{C_RESOLVE / C_FINE:.0f}（log軸で平行移動）',
                 fontsize=12)
    handles = [plt.Line2D([], [], color=lcol[L], lw=2, label=f'L={L}')
               for L in MAIN_L if any(v[0] == L for v in series.values())]
    ax.legend(handles=handles, fontsize=9, ncol=4, loc='lower right', frameon=False)
    out = os.path.join(OUTDIR, 'fig_required_den_all_runs_20260917.png')
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    log_line(f'figure -> {out} ({len(series)} run)')


if __name__ == '__main__':
    main()
