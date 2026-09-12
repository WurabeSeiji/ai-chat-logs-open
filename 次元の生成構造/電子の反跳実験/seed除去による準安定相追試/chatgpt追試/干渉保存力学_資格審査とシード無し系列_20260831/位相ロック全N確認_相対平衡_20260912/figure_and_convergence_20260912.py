#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全N位相ロックの図化＋「False=未収束か否か」の判別（読出しのみ）。

判別の要点: ロックの終端値ではなく、末尾でのロック指標がまだ減少中（=収束途中）か、
下げ止まり（=真にロックしない）かを、末尾100stepの log(rel_rate) の傾きで見る。
傾き<0 が続けば、onsetが遅く整定時間 T-onset が短いだけで、同じロックへ向かっている。
"""
import csv
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', 'make_parent型初期値_自己無撞着構造_20260907'))
_spec = importlib.util.spec_from_file_location('orig', os.path.join(GEN, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(orig)
NS = list(range(3, 41))


def series_state(N):
    X, Y, fz = orig.series(N)
    S = np.asarray(np.load(os.path.join(orig.STATES, f'hm_N{N}_den_{N}_states_500.npz'))['Z'], np.complex128)
    return S, fz


def rel_rate_series(S):
    phi = np.unwrap(np.angle(S), axis=0)
    rel = phi - phi[:, [0]]
    return np.median(np.abs(np.diff(rel, axis=0)), axis=1)  # (T-1,)


def main():
    rows = list(csv.DictReader(open(os.path.join(HERE, 'phaselock_allN.csv'))))
    conv = []
    for N in NS:
        S, fz = series_state(N)
        rr = rel_rate_series(S)
        onset = int(np.argmax(fz > 1e-3))
        tail = rr[-100:]
        tail = np.maximum(tail, 1e-18)
        x = np.arange(len(tail))
        slope = np.polyfit(x, np.log10(tail), 1)[0]   # decade/step。負=まだ減少=収束中
        settle = 500 - onset
        conv.append((N, onset, settle, float(rr[-1]), float(slope)))

    # 図
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        Ns = [c[0] for c in conv]
        rig_tail = [float(r['rig_tail']) for r in rows]
        rel_tail = [float(r['rel_rate_tail']) for r in rows]
        settle = [c[2] for c in conv]
        slopes = [c[4] for c in conv]
        amp_ratio = [float(r['amp_ratio_tail']) for r in rows]
        delta = [float(r['delta_tail']) for r in rows]
        ideal = [float(r['delta_ideal']) for r in rows]
        fig, ax = plt.subplots(2, 2, figsize=(13, 9))
        ax[0, 0].semilogy(Ns, np.maximum(rig_tail, 1e-18), 'o-', label='rigid-rotation residual (tail)')
        ax[0, 0].semilogy(Ns, np.maximum(rel_tail, 1e-18), 's-', label='rel-phase rate (tail, gauge-inv)')
        ax[0, 0].axhline(1e-6, ls=':', c='gray'); ax[0, 0].set_xlabel('N'); ax[0, 0].set_title('lock metrics at tail vs N'); ax[0, 0].legend(fontsize=8)
        ax[0, 1].plot(Ns, np.array(delta), 'o-', label='Delta_tail (measured)')
        ax[0, 1].plot(Ns, np.array(ideal), 'k--', label='-2pi/N')
        ax[0, 1].set_xlabel('N'); ax[0, 1].set_title('tail 1-step rotation = -2 pi / N'); ax[0, 1].legend(fontsize=8)
        ax[1, 0].plot(settle, np.log10(np.maximum(rel_tail, 1e-18)), 'o')
        for c in conv:
            ax[1, 0].annotate(str(c[0]), (c[2], np.log10(max(c[3], 1e-18))), fontsize=6)
        ax[1, 0].set_xlabel('settling time T-onset [step]'); ax[1, 0].set_ylabel('log10 rel-phase rate (tail)')
        ax[1, 0].set_title('unlocked = short settling time (left)')
        ax[1, 1].plot(Ns, slopes, 'o-'); ax[1, 1].axhline(0, ls=':', c='r')
        ax[1, 1].set_xlabel('N'); ax[1, 1].set_ylabel('tail-100 slope of log(rel_rate) [dec/step]')
        ax[1, 1].set_title('all negative = still decreasing = converging to lock')
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_phaselock_allN.png'), dpi=110)
        print('wrote fig_phaselock_allN.png')
    except Exception as e:
        print('figure skipped:', e)

    n_still_decreasing = sum(1 for c in conv if c[4] < 0)
    print(f'\n末尾でrel_rateがまだ減少中(傾き<0)のN: {n_still_decreasing}/{len(conv)}')
    print(f"{'N':>3}{'onset':>6}{'settle=T-onset':>15}{'rel_rate_tail':>14}{'tail_slope[dec/step]':>22}")
    for N, onset, settle, rrt, slope in conv:
        print(f'{N:>3}{onset:>6}{settle:>15}{rrt:>14.2e}{slope:>22.3e}')
    return conv


if __name__ == '__main__':
    main()
