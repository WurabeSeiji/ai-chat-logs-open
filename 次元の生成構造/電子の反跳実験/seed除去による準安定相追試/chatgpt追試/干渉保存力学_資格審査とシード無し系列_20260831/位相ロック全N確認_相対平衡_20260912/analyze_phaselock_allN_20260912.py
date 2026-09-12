#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全N（3..40）で「相対平衡→インフレーション→相対平衡ロック」を定量化する（読出しのみ・物理無変更）。

木原指示（2026-09-12）: N=6 で確認した大τでの位相ロックを全Nで確認する。

測定量（すべて既存 den=N・500step 状態 npz の読出しのみ。走行・親生成・物理変更なし）:
  1. 剛体回転残差 rig(τ) = 1 − |⟨z(τ),z(τ+1)⟩|² / (‖z(τ)‖²‖z(τ+1)‖²)
       相対平衡（固定形が一定角速度で回る）なら 0。onset前の床と末尾で小、inflation中に増大。
  2. 末尾50step の rig 中央値 rig_tail、1step回転角 Δ_tail=arg⟨z,z'⟩ の収束値。
  3. ゲージ不変な相対位相 φ_m−φ_0 の末尾変化率（θ0非依存＝物理のロック指標）。
  4. 末尾振幅の等半径度: min/max 比、変動係数中央値。
  5. onset = 最初に H⊥/H>1e-3 となる step。inflation中の rig ピーク。
  θ0・H⊥/H・PC 等の物理式は正本 plot3d_makeparent_selectedN_v1 から取得（再実装しない）。

出力: REPORT.md, phaselock_allN.csv, fig_phaselock_allN.png, run.log, SHA256SUMS.txt（run_all.sh から）
"""
import csv
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', 'make_parent型初期値_自己無撞着構造_20260907'))
_spec = importlib.util.spec_from_file_location(
    'orig', os.path.join(GEN, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(orig)
NS = list(range(3, 41))
TAIL = 50


def rigid_resid(S):
    """各stepの剛体回転残差 1-|overlap|²。"""
    z0, z1 = S[:-1], S[1:]
    ov = np.sum(np.conj(z0) * z1, axis=1)
    n0 = np.sum(np.abs(z0) ** 2, axis=1)
    n1 = np.sum(np.abs(z1) ** 2, axis=1)
    return 1 - np.abs(ov) ** 2 / (n0 * n1), np.angle(ov)


def analyze(N):
    X, Y, fz = orig.series(N)
    d = np.load(os.path.join(orig.STATES, f'hm_N{N}_den_{N}_states_500.npz'))
    S = np.asarray(d['Z'], np.complex128)
    T, M = S.shape
    onset = int(np.argmax(fz > 1e-3)) if np.any(fz > 1e-3) else -1
    rig, dang = rigid_resid(S)
    rig_tail = float(np.median(rig[-TAIL:]))
    delta_tail = float(np.median(dang[-TAIL:]))
    rig_peak = float(np.max(rig[max(onset - 5, 0):] if onset > 0 else rig))
    # 床（onset前）の剛体回転残差
    rig_floor = float(np.median(rig[:max(onset - 10, 2)])) if onset > 12 else float('nan')

    # ゲージ不変相対位相 φ_m-φ_0（θ0非依存）の末尾変化率
    phi = np.unwrap(np.angle(S), axis=0)
    rel = phi - phi[:, [0]]
    rel_rate_onset = float(np.median(np.abs(np.diff(rel, axis=0)), axis=1)[onset]) if onset > 0 else float('nan')
    rel_rate_tail = float(np.median(np.abs(np.diff(rel, axis=0))[-TAIL:]))
    rel_std_tail = float(np.median(rel[-TAIL:].std(axis=0)))

    # 振幅の等半径度（末尾）
    r = np.abs(S)
    ratio_tail = float(r[-1].max() / r[-1].min())
    cv_tail = float(np.median(r[-TAIL:].std(0) / r[-TAIL:].mean(0)))

    # den=N の理想1step回転 -2π/N と Δ_tail の比較
    delta_ideal = -2 * np.pi / N
    # 主閉包内で近い整数倍 k: Δ_tail ≈ -2πk/N
    kfrac = -delta_tail * N / (2 * np.pi)

    return dict(N=N, M=M, onset=onset, rig_floor=rig_floor, rig_peak=rig_peak,
                rig_tail=rig_tail, delta_tail=delta_tail, delta_ideal=delta_ideal,
                k_over_1=kfrac, rel_rate_onset=rel_rate_onset, rel_rate_tail=rel_rate_tail,
                rel_std_tail=rel_std_tail, amp_ratio_tail=ratio_tail, amp_cv_tail=cv_tail,
                locked=(rel_rate_tail < 1e-4 and rig_tail < 1e-6))


def main():
    rows = [analyze(N) for N in NS]
    cols = list(rows[0].keys())
    with open(os.path.join(HERE, 'phaselock_allN.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow(r)
    n_locked = sum(r['locked'] for r in rows)
    print(f'全{len(rows)}本中 ロック成立(rel_rate_tail<1e-4 & rig_tail<1e-6): {n_locked}')
    print(f"{'N':>3}{'onset':>6}{'rig_floor':>11}{'rig_peak':>10}{'rig_tail':>11}"
          f"{'Δ_tail':>9}{'-2π/N':>9}{'rel_rate_tail':>14}{'amp_ratio':>10}{'lock':>6}")
    for r in rows:
        print(f"{r['N']:>3}{r['onset']:>6}{r['rig_floor']:>11.2e}{r['rig_peak']:>10.2e}"
              f"{r['rig_tail']:>11.2e}{r['delta_tail']:>9.4f}{r['delta_ideal']:>9.4f}"
              f"{r['rel_rate_tail']:>14.2e}{r['amp_ratio_tail']:>10.4f}{str(r['locked']):>6}")
    return rows


if __name__ == '__main__':
    main()
