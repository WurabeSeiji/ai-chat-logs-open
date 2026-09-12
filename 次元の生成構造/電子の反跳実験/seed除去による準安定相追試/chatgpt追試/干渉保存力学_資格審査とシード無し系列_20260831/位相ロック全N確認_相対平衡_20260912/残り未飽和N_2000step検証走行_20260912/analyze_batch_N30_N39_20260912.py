#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=30..39・2000step バッチ走行の位相ロック一括解析＋予測検証（読出しのみ）。

入力: results/hm_N{N}_den_{N}_states_2000.npz（本フォルダ batch wrapper で生成）。
指標は N22/N40 解析と同一（onset・rig残差・Δ_tail=-2π/N・N step周期再帰・ゲージ不変rel_rate lockstep・等半径）。
予測（../未ロックN_2000step到達予測_20260912, モデルB lockstep≈35·N）と実測lockstepを対比。
対照: 各N step0 が旧500step走行と bit一致（step1 は機械ε、以後インフレ指数増幅）。
出力: REPORT用のstdout + phaselock_batch_N30_N39.csv + fig_batch_lockstep_vs_pred.png
"""
import csv
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907'))
_spec = importlib.util.spec_from_file_location('orig', os.path.join(GEN, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(orig)
NS = list(range(30, 40))
THR = 1e-6


def rel_rate(S):
    phi = np.unwrap(np.angle(S), axis=0); rel = phi - phi[:, [0]]
    return np.median(np.abs(np.diff(rel, axis=0)), axis=1)


def hperp(S):
    z0 = S[0]; p = z0.real / np.linalg.norm(z0.real)
    q = z0.imag - np.dot(z0.imag, p) * p; q /= np.linalg.norm(q)
    return np.maximum(np.array([np.vdot(z-p*np.dot(p, z)-q*np.dot(q, z),
                                        z-p*np.dot(p, z)-q*np.dot(q, z)).real/np.vdot(z, z).real
                                for z in S]), orig.FLOOR)


def rig(S):
    z0, z1 = S[:-1], S[1:]
    ov = np.sum(np.conj(z0)*z1, axis=1)
    return 1-np.abs(ov)**2/(np.sum(np.abs(z0)**2, 1)*np.sum(np.abs(z1)**2, 1)), np.angle(ov)


def main():
    rows = []
    for N in NS:
        f = os.path.join(HERE, 'results', f'hm_N{N}_den_{N}_states_2000.npz')
        if not os.path.exists(f):
            print(f'N={N}: npz なし（未走行?） skip'); continue
        S = np.asarray(np.load(f)['Z'], np.complex128); T, M = S.shape
        fz = hperp(S); onset = int(np.argmax(fz > 1e-3))
        r, dang = rig(S); rr = rel_rate(S); pk = int(np.argmax(rr))
        reach = np.where(rr[pk:] < THR)[0]
        lockstep = int(pk+reach[0]) if len(reach) else None
        delta_tail = float(np.median(dang[-200:])); ideal = -2*np.pi/N
        amp = np.abs(S[-1])
        # 対照: 旧500stepとstep0 bit一致
        old = os.path.join(orig.STATES, f'hm_N{N}_den_{N}_states_500.npz')
        s0 = float(np.max(np.abs(S[0]-np.asarray(np.load(old)['Z'][0], np.complex128)))) if os.path.exists(old) else -1
        pred = 35.0*N  # モデルB近似
        rows.append(dict(N=N, M=M, onset=onset, rig_floor=float(np.median(r[:max(onset-20,2)])),
                         rig_tail=float(np.median(r[-200:])), delta_tail=delta_tail, delta_ideal=ideal,
                         delta_err=abs(delta_tail-ideal), recur=delta_tail*N/(2*np.pi),
                         lockstep=lockstep, pred_lockstep=pred, rel_rate_tail=float(np.median(rr[-200:])),
                         amp_ratio=float(amp.max()/amp.min()), step0_ctrl=s0, locked_by_2000=(lockstep is not None)))
    if not rows:
        print('解析対象なし'); return
    with open(os.path.join(HERE, 'phaselock_batch_N30_N39.csv'), 'w', newline='') as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys())); w.writeheader()
        for r_ in rows:
            w.writerow(r_)
    print(f"{'N':>3}{'M':>5}{'onset':>6}{'lockstep':>9}{'pred≈35N':>9}{'2000内':>7}"
          f"{'Δ_tail':>9}{'-2π/N':>9}{'再帰':>8}{'amp比':>8}{'step0対照':>10}")
    for r_ in rows:
        print(f"{r_['N']:>3}{r_['M']:>5}{r_['onset']:>6}{str(r_['lockstep']):>9}{r_['pred_lockstep']:>9.0f}"
              f"{str(r_['locked_by_2000']):>7}{r_['delta_tail']:>9.4f}{r_['delta_ideal']:>9.4f}"
              f"{r_['recur']:>8.4f}{r_['amp_ratio']:>8.4f}{r_['step0_ctrl']:>10.1e}")
    nlock = sum(r_['locked_by_2000'] for r_ in rows)
    print(f'\n2000内ロック: {nlock}/{len(rows)}  Δ=-2π/N一致(全N最大誤差)={max(r_["delta_err"] for r_ in rows):.1e}')

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        Ns = [r_['N'] for r_ in rows]; ls = [r_['lockstep'] or np.nan for r_ in rows]; pr = [r_['pred_lockstep'] for r_ in rows]
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        ax[0].plot(Ns, ls, 'o-', label='actual lockstep')
        ax[0].plot(Ns, pr, 's--', label='predicted ~35N (model B)')
        ax[0].axhline(2000, ls=':', c='r', label='2000 step'); ax[0].set_xlabel('N'); ax[0].set_ylabel('lockstep')
        ax[0].set_title('N=30..39: actual vs predicted lockstep'); ax[0].legend(fontsize=8)
        ax[1].plot(Ns, [r_['delta_tail'] for r_ in rows], 'o-', label='Δ_tail')
        ax[1].plot(Ns, [r_['delta_ideal'] for r_ in rows], 'k--', label='-2π/N')
        ax[1].set_xlabel('N'); ax[1].set_title('rotation locks to -2π/N'); ax[1].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_batch_lockstep_vs_pred.png'), dpi=110)
        print('wrote fig_batch_lockstep_vs_pred.png')
    except Exception as e:
        print('figure skipped:', e)


if __name__ == '__main__':
    main()
