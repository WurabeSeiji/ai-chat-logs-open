#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=64・den=64・2000step 走行の位相ロック検証（読出しのみ・物理無変更）。

木原指示（2026-09-12）: N64_2000/N100_2000 走行系が同条件か確認せよ。
  - N100 は未実行（状態npzなし）→解析不能。
  - N64 は既存 2000step 状態 npz を読み、全N(3-40)解析と同一の指標で位相ロックを検証。

入力: ../../make_parent型初期値_自己無撞着構造_20260907/N64_2000step走行_20260907/results/
      hm_N64_den_64_states_2000.npz  （key: Z(2001,2016), N=64, denominator=64, steps=2000）
H⊥/H の式は正本 plot3d_makeparent_selectedN_v1.series と同一（本ファイルで再実装せず同式を適用）。

測定は全N解析 analyze_phaselock_allN と同一:
  剛体回転残差 rig(τ)=1-|⟨z,z'⟩|²/(‖z‖²‖z'‖²)、Δ_tail=arg⟨z,z'⟩、
  ゲージ不変相対位相 φ_m-φ_0 の変化率、振幅の等半径度、onset(H⊥/H>1e-3)。
出力: REPORT 兼 stdout（run_all.sh から run.log / SHA へ）、fig_phaselock_N64_2000.png
"""
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907'))
NPZ = os.path.join(GEN, 'N64_2000step走行_20260907', 'results', 'hm_N64_den_64_states_2000.npz')
_spec = importlib.util.spec_from_file_location('orig', os.path.join(GEN, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(orig)
FLOOR = orig.FLOOR
TAIL = 200


def hperp_over_h(S):
    """orig.series と同一の H⊥/H（初期平面直交成分の割合）。"""
    z0 = S[0]
    p = z0.real.astype(np.float64).copy(); p /= np.linalg.norm(p)
    q = z0.imag.astype(np.float64).copy(); q -= np.dot(q, p) * p; q /= np.linalg.norm(q)
    f = np.empty(S.shape[0])
    for s in range(S.shape[0]):
        z = S[s]; a = np.dot(p, z); b = np.dot(q, z); zp = z - p * a - q * b
        f[s] = np.vdot(zp, zp).real / np.vdot(z, z).real
    return np.maximum(f, FLOOR)


def rig(S):
    z0, z1 = S[:-1], S[1:]
    ov = np.sum(np.conj(z0) * z1, axis=1)
    n0 = np.sum(np.abs(z0) ** 2, axis=1); n1 = np.sum(np.abs(z1) ** 2, axis=1)
    return 1 - np.abs(ov) ** 2 / (n0 * n1), np.angle(ov)


def main():
    d = np.load(NPZ)
    S = np.asarray(d['Z'], np.complex128)
    N = int(d['N']); assert int(d['denominator']) == N
    T, M = S.shape
    fz = hperp_over_h(S)
    onset = int(np.argmax(fz > 1e-3))
    r, dang = rig(S)
    phi = np.unwrap(np.angle(S), axis=0); rel = phi - phi[:, [0]]
    rate = np.median(np.abs(np.diff(rel, axis=0)), axis=1)
    rr = np.abs(S)
    delta_tail = float(np.median(dang[-TAIL:])); ideal = -2 * np.pi / N
    tl = np.maximum(rate[-300:], 1e-18); slope = float(np.polyfit(np.arange(len(tl)), np.log10(tl), 1)[0])

    print(f'N={N} den={N} T={T}(step0..{T-1}) M={M}')
    print(f'onset(H⊥/H>1e-3)=step{onset}  (事前登録9.5×N≈{9.5*N:.0f}, 実測比 {onset/N:.2f}×N)')
    print(f'H⊥/H: step0={fz[0]:.2e} 末尾={fz[-1]:.3f}')
    print(f'剛体回転残差 rig: 床={np.median(r[:max(onset-20,2)]):.2e} ピーク={np.max(r[max(onset-5,0):]):.2e} '
          f'末尾{TAIL}={np.median(r[-TAIL:]):.2e}')
    print(f'Δ_tail={delta_tail:+.5f} rad  -2π/{N}={ideal:+.5f} rad  |差|={abs(delta_tail-ideal):.2e}  '
          f'={np.degrees(abs(delta_tail)):.3f}°')
    print(f'{N}step周期再帰: Δ_tail×{N}/(2π)={delta_tail*N/(2*np.pi):+.4f} (=-1で1周/{N}step)')
    print(f'相対位相変化率(θ0非依存): onset={rate[onset]:.2e} 中間={rate[T//2]:.2e} 末尾{TAIL}={np.median(rate[-TAIL:]):.2e} '
          f'/ 末尾300 log傾き={slope:+.2e}dec/step(負=収束中)')
    print(f'振幅ロック: 末尾 min/max={rr[-1].min():.4f}/{rr[-1].max():.4f} 比{rr[-1].max()/rr[-1].min():.4f} '
          f'変動係数{np.median(rr[-TAIL:].std(0)/rr[-TAIL:].mean(0)):.2e}')

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        step = np.arange(T)
        fig, ax = plt.subplots(2, 2, figsize=(13, 9))
        ax[0, 0].semilogy(step, fz); ax[0, 0].axvline(onset, ls=':', c='r')
        ax[0, 0].set_title(f'H_perp/H (onset step{onset})'); ax[0, 0].set_xlabel('step')
        ax[0, 1].semilogy(step[:-1], np.maximum(r, 1e-18)); ax[0, 1].axvline(onset, ls=':', c='r')
        ax[0, 1].set_title('rigid-rotation residual 1-|<z,z\'>|^2'); ax[0, 1].set_xlabel('step')
        ax[1, 0].plot(step[:-1], dang); ax[1, 0].axhline(ideal, ls='--', c='k', label=f'-2pi/{N}')
        ax[1, 0].set_title('1-step rotation angle -> -2pi/N'); ax[1, 0].set_xlabel('step'); ax[1, 0].legend()
        ax[1, 1].semilogy(step[:-1], np.maximum(rate, 1e-18))
        ax[1, 1].set_title('gauge-inv rel-phase rate [rad/step]'); ax[1, 1].set_xlabel('step')
        fig.suptitle(f'N=64 den=64 2000step phase-lock (read-only, no physics change)')
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_phaselock_N64_2000.png'), dpi=110)
        print('wrote fig_phaselock_N64_2000.png')
    except Exception as e:
        print('figure skipped:', e)


if __name__ == '__main__':
    main()
