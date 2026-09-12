#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=40・den=40・2000step 走行の位相ロック解析＋予測検証（読出しのみ）。

入力: results/hm_N40_den_40_states_2000.npz（本フォルダ wrapper で生成）。
指標は全N解析・N64解析と同一（rig残差・Δ_tail・ゲージ不変rel_rate・振幅・onset）。
予測（../未ロックN_2000step到達予測_20260912）: N=40 lockstep(1e-6)≈1476、2000でロック。本走行で検証。

★対照テスト（忠実再現の確認）: 既存500step走行 full_N3_N40_sweep/states/hm_N40_den_40_states_500.npz と
  step0 が bit一致、step1 差=O(1e-16)＝機械イプシロン。以後はインフレーションが丸め差を指数増幅する
  （step300 で O(0.1)）。親・物理は同一で、アトラクタ量（onset/Δ=-2π/N/ロック）は軌道非依存で検証可能。
"""
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907'))
_spec = importlib.util.spec_from_file_location('orig', os.path.join(GEN, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(orig)
NPZ = os.path.join(HERE, 'results', 'hm_N40_den_40_states_2000.npz')
OLD500 = os.path.join(orig.STATES, 'hm_N40_den_40_states_500.npz')
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
    ov = np.sum(np.conj(z0) * z1, axis=1)
    return 1 - np.abs(ov)**2/(np.sum(np.abs(z0)**2, 1)*np.sum(np.abs(z1)**2, 1)), np.angle(ov)


def main():
    S = np.asarray(np.load(NPZ)['Z'], np.complex128)
    N = 40; T, M = S.shape
    fz = hperp(S); onset = int(np.argmax(fz > 1e-3))
    r, dang = rig(S)
    rr = rel_rate(S); pk = int(np.argmax(rr))
    reach = np.where(rr[pk:] < THR)[0]
    lockstep = int(pk + reach[0]) if len(reach) else None
    ideal = -2*np.pi/N; delta_tail = float(np.median(dang[-200:]))
    amp = np.abs(S)

    # 対照テスト（bit一致→指数増幅）
    old = np.asarray(np.load(OLD500)['Z'], np.complex128)
    ctrl = {s: float(np.max(np.abs(S[s]-old[s]))) for s in [0, 1, 10, 100, 200, 300, 500]}

    print(f'N={N} den={N} T={T}(step0..{T-1}) M={M}')
    print(f'[対照] 旧500stepとのstep差: ' + ' '.join(f'{s}:{v:.1e}' for s, v in ctrl.items())
          + '  (step0=0 bit一致, step1~1e-16 機械ε, 以後インフレ指数増幅)')
    print(f'onset(H⊥/H>1e-3)=step{onset}   H⊥/H末尾={fz[-1]:.3f}')
    print(f'剛体回転残差 rig: 床={np.median(r[:max(onset-20,2)]):.2e} ピーク={np.max(r[max(onset-5,0):]):.2e} '
          f'末尾200={np.median(r[-200:]):.2e}')
    print(f'Δ_tail={delta_tail:+.5f} rad  -2π/{N}={ideal:+.5f}  |差|={abs(delta_tail-ideal):.2e}')
    print(f'{N}step周期再帰: Δ_tail×{N}/(2π)={delta_tail*N/(2*np.pi):+.4f}')
    print(f'振幅ロック: 末尾 min/max={amp[-1].min():.4f}/{amp[-1].max():.4f} 比{amp[-1].max()/amp[-1].min():.4f}')
    print(f'相対位相変化率: onset={rr[onset]:.2e} 末尾200中央={np.median(rr[-200:]):.2e} rad/step')
    print(f'\n★実測 lockstep(rel_rate<{THR:.0e})={lockstep}   予測≈1476   '
          f'判定: {"予測と整合・2000内ロック確認" if lockstep and lockstep<2000 else "要検討"}')

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        step = np.arange(T)
        fig, ax = plt.subplots(2, 2, figsize=(13, 9))
        ax[0, 0].semilogy(step, fz); ax[0, 0].axvline(onset, ls=':', c='r')
        ax[0, 0].set_title(f'H_perp/H (onset step{onset})'); ax[0, 0].set_xlabel('step')
        ax[0, 1].plot(step[:-1], dang); ax[0, 1].axhline(ideal, ls='--', c='k', label=f'-2pi/{N}')
        ax[0, 1].set_title('1-step rotation -> -2pi/N'); ax[0, 1].legend(); ax[0, 1].set_xlabel('step')
        ax[1, 0].semilogy(step[:-1], np.maximum(rr, 1e-18)); ax[1, 0].axhline(THR, ls='--', c='r')
        if lockstep: ax[1, 0].axvline(lockstep, ls='-', c='g', label=f'lockstep {lockstep}')
        ax[1, 0].axvline(1476, ls=':', c='orange', label='predicted 1476'); ax[1, 0].legend(fontsize=8)
        ax[1, 0].set_title('gauge-inv rel-phase rate & lockstep'); ax[1, 0].set_xlabel('step')
        ax[1, 1].semilogy(step[:-1], np.maximum(r, 1e-18)); ax[1, 1].axvline(onset, ls=':', c='r')
        ax[1, 1].set_title('rigid-rotation residual'); ax[1, 1].set_xlabel('step')
        fig.suptitle(f'N=40 den=40 2000step phase-lock verification (read-only run of faithful copy)')
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_phaselock_N40_2000.png'), dpi=110)
        print('wrote fig_phaselock_N40_2000.png')
    except Exception as e:
        print('figure skipped:', e)


if __name__ == '__main__':
    main()
