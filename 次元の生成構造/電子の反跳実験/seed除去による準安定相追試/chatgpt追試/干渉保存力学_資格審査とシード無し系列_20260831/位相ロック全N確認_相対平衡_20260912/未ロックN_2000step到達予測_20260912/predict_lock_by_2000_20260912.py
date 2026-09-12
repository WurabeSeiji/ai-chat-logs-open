#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""500step系(N=3..40)で未ロックのNを抽出し、各Nが2000stepでロックするかを予測する（読出しのみ・走行なし）。

木原指示（2026-09-12）: 未ロックN（N=22等）を2000stepまで無停止走行してはどうか。
  大Nが2000でロックするか推測可能か、N=40の推測を基準に。

ロック指標: ゲージ不変な相対位相変化率 rel_rate=median|Δ(φ_m−φ_0)| [rad/step]。THR=1e-6 で"ロック"。
物理: ロック=相対平衡への収束。収束レートは線形化のスペクトルギャップ ∝ 1/N（line-graph 固有値
差 ∝ 1/N）なので、ピーク後の減衰レート r(N) ∝ 1/N を仮定。

★仮説の反証と修正（この解析で確定）:
  仮説A: 「N=40 のピーク直後レート r40=-6.0e-3 dec/step を r(N)=r40·40/N で外挿」。
    → N=64 予測 lockstep≈1428。しかし N=64 実データ(2000step)は step2000 で rel_rate=3.0e-6
      と 1e-6 未到達（実測 lockstep>2000）。**仮説Aは N=64 の実測に反証された**（早すぎる）。
    原因: 減衰はピーク直後が速く固定点接近で遅くなる。N=40 の500stepは速い初期相しか含まない。
  仮説B(採用): 「唯一の長時間走行 N=64 の tail レート r64t=-1.87e-3 (rate×N≈-0.12) を
    r(N)=r64t·64/N で外挿」。N=64 を自己再現し、tail の遅い相を正しく反映。
    lockstep(N) = peak(N) + (log10 peakval − log10 THR)/|r(N)|、peak(N)=peak係数·N（実測≈8.5N）。

入力:
  ../../make_parent型初期値_自己無撞着構造_20260907/full_N3_N40_sweep/states/hm_N{N}_den_{N}_states_500.npz
  ../../make_parent型初期値_自己無撞着構造_20260907/N64_2000step走行_20260907/results/hm_N64_den_64_states_2000.npz
"""
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEN = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907'))
_spec = importlib.util.spec_from_file_location('orig', os.path.join(GEN, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(orig)
N64PATH = os.path.join(GEN, 'N64_2000step走行_20260907', 'results', 'hm_N64_den_64_states_2000.npz')
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


def load(N):
    if N == 64:
        return np.asarray(np.load(N64PATH)['Z'], np.complex128)
    return np.asarray(np.load(os.path.join(orig.STATES, f'hm_N{N}_den_{N}_states_500.npz')['Z']
                              if False else os.path.join(orig.STATES, f'hm_N{N}_den_{N}_states_500.npz'))['Z'],
                      np.complex128)


def main():
    # 1) 未ロックN
    notlocked, onsets, peaks = [], {}, {}
    for N in range(3, 41):
        S = load(N); fz = hperp(S); rr = rel_rate(S)
        onsets[N] = int(np.argmax(fz > 1e-3)); peaks[N] = int(np.argmax(rr))
        if np.median(rr[-50:]) > 1e-4:
            notlocked.append(N)
    print('500stepで未ロック(rel_rate_tail>1e-4)のN:', notlocked)
    onset_a, onset_b = np.polyfit(list(onsets), [onsets[N] for N in onsets], 1)
    peak_coef = np.mean([peaks[N] / N for N in (40,)])  # peak≈係数·N。N=40基準
    print(f'onset(N)≈{onset_a:.2f}N+{onset_b:.1f}   peak(N)≈{peak_coef:.1f}·N (N=40基準)')

    # 2) レート測定
    rr40 = rel_rate(load(40)); pk40 = int(np.argmax(rr40)); pv40 = float(rr40[pk40])
    seg40 = np.log10(np.maximum(rr40[pk40:], 1e-18)); r40 = float(np.polyfit(np.arange(len(seg40)), seg40, 1)[0])
    S64 = load(64); rr64 = rel_rate(S64); pk64 = int(np.argmax(rr64))
    tail = np.log10(np.maximum(rr64[-500:], 1e-18)); r64t = float(np.polyfit(np.arange(500), tail, 1)[0])
    decades = np.log10(pv40) - np.log10(THR)
    print(f'\nN=40: peak@{pk40} peakval={pv40:.2e} 初期レート r40={r40:+.2e} dec/step')
    print(f'N=64: peak@{pk64} tailレート(最後500) r64t={r64t:+.2e} dec/step (rate×N={r64t*64:+.2f})')
    print(f'落とす桁 decades(peakval→1e-6)={decades:.2f}')

    def lockstepA(N):  # 反証済み: N=40初期レート
        rN = r40 * 40.0 / N
        return peak_coef * N + decades / abs(rN)

    def lockstepB(N):  # 採用: N=64 tailレート
        rN = r64t * 64.0 / N
        return peak_coef * N + decades / abs(rN)

    # 3) N=64 での検証（実測 lockstep）
    reach = np.where(rr64[pk64:] < THR)[0]
    actual64 = (pk64 + int(reach[0])) if len(reach) else None
    print(f'\n【検証】N=64 実測lockstep(1e-6)={actual64} (Noneなら2000内未到達, 末尾{rr64[-1]:.2e})')
    print(f'  仮説A予測={lockstepA(64):.0f} → 実測>2000 と矛盾（反証）')
    print(f'  仮説B予測={lockstepB(64):.0f} → 実測(tail外挿≈2250)と整合（採用）')

    # 4) 採用モデルBで予測
    print('\n【予測 採用モデルB】lockstep(1e-6到達):')
    cross = None
    for N in [22, 30, 40, 50, 54, 57, 64, 80, 100]:
        ls = lockstepB(N); ok = ls < 2000
        if cross is None and not ok:
            cross = N
        print(f'  N={N:>3}: lockstep≈{ls:>5.0f}  ({"< 2000 ロック" if ok else ">=2000 未達"})')
    # 交差N
    Ns = np.arange(6, 121); lsB = np.array([lockstepB(N) for N in Ns])
    idx = np.where(lsB >= 2000)[0]
    print(f'\n2000stepでロックできなくなる交差N ≈ {Ns[idx[0]] if len(idx) else ">120"}')
    print(f'→ 未ロック12本(N=22,30..40)は全て2000で解消（ロック）と予測。')
    print(f'  N=40基準の headline: lockstep≈{lockstepB(40):.0f}（余裕≈{2000-lockstepB(40):.0f}step）')

    # しきい依存（N=40, モデルB）
    print('\nしきい依存(N=40, モデルB):')
    for thr in [1e-4, 1e-5, 1e-6]:
        dec = np.log10(pv40) - np.log10(thr)
        print(f'  THR={thr:.0e}: lockstep≈{peak_coef*40 + dec/abs(r64t*64/40):.0f}')

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        ax[0].plot(Ns, [lockstepA(N) for N in Ns], '--', c='gray', label='model A (N=40 init-rate, FALSIFIED)')
        ax[0].plot(Ns, lsB, '-', label='model B (N=64 tail-rate, adopted)')
        ax[0].axhline(2000, ls='--', c='r', label='2000 step'); ax[0].axhline(500, ls=':', c='gray', label='500 step')
        ax[0].plot([64], [actual64 or 2250], 'k*', ms=13, label='N=64 actual (~2250, >2000)')
        ax[0].plot([40], [lockstepB(40)], 'go', ms=9, label=f'N=40 basis (~{lockstepB(40):.0f})')
        ax[0].set_xlabel('N'); ax[0].set_ylabel('predicted lockstep (rel_rate<1e-6)')
        ax[0].set_title('lock-by-2000 prediction'); ax[0].legend(fontsize=7)
        ax[1].semilogy(np.arange(len(rr64)), np.maximum(rr64, 1e-18), label='N=64 actual')
        ax[1].semilogy(np.arange(len(rr40)), np.maximum(rr40, 1e-18), alpha=0.7, label='N=40 (500step)')
        ax[1].axhline(THR, ls='--', c='r'); ax[1].axvline(2000, ls=':', c='gray')
        ax[1].set_xlabel('step'); ax[1].set_ylabel('rel-phase rate'); ax[1].legend(fontsize=8)
        ax[1].set_title('anchor: N=64 reaches 3e-6 at step2000 (not <1e-6)')
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_predict_lock_by_2000.png'), dpi=110)
        print('\nwrote fig_predict_lock_by_2000.png')
    except Exception as e:
        print('figure skipped:', e)


if __name__ == '__main__':
    main()
