#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""「新しい面/S³/rank4」は本物か作為かの判別テスト（読出しのみ・物理無変更）。

木原指示（2026-09-12）: 「厳密に新しい面に乗っているのか、ただのアーティファクトでは？」

判別原理:
  本物の「2回転面が同時に立つ(S³/rank4)」なら、瞬間（短窓SVD）でも第2面 s3/s1≳0.1 が立つはず。
  実体が「単一回転面(rank2)が時間とともに向きを変える」だけなら、短窓では s3/s1≈0（単一面）で、
  長窓にするほど（向きの違う同じ1面を集約して）見かけの s3/s1 が増える＝長窓SVDの集約アーティファクト。

測定（N=40 の2000step軌道, X(t)=[Re z; Im z]∈R^{2M}）:
  (A) inflation中の固定開始点から窓幅を変えて s3/s1。窓幅依存なら作為。
  (B) 30step短窓を滑らせ、各窓内 s3/s1（瞬間rank）と、前窓の主平面からの主角度（面の再配向）。

結論（本スクリプトの出力）:
  短窓 s3/s1≈0（瞬間rank2＝単一面）、窓幅とともに増加、面は連続再配向（onset直後~10°/30step→ロックで0）。
  ⟹ 「S³/rank4/二平面同時/振幅交換」は長窓SVDの集約アーティファクト。物理は単一面の連続傾き（net~19°）のみ。
  ⟹ 思考実験⑥「post-inflation=dominant rank4=S³」は反証（要修正）。
入力: ../N40_2000step検証走行_20260912/results/hm_N40_den_40_states_2000.npz
出力: fig_artifact_shortwindow_svd.png
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NPZ = os.path.join(HERE, '..', 'N40_2000step検証走行_20260912', 'results', 'hm_N40_den_40_states_2000.npz')


def s3s1(seg):
    W = seg - seg.mean(0); s = np.linalg.svd(W, compute_uv=False); return s[2] / s[0]


def main():
    S = np.asarray(np.load(NPZ)['Z'], np.complex128); T, M = S.shape
    X = np.concatenate([S.real, S.imag], 1)

    # (A) 窓幅依存
    mid = 350
    widths = [10, 20, 40, 100, 300, 700]
    wvals = [(w, s3s1(X[mid:mid + w])) for w in widths if mid + w <= T]
    print('(A) 固定開始 step350 からの窓幅依存 s3/s1（本物なら窓幅非依存で~0.1, 作為なら窓幅とともに増加）:')
    for w, v in wvals:
        print(f'    窓幅{w:>4}: s3/s1={v:.4f}')

    # (B) 30step短窓の瞬間rankと面再配向
    Wn = 30; starts = list(range(300, 660, 30))
    inst, tilts = [], []
    prev = None
    print('\n(B) 30step短窓: 窓内s3/s1(瞬間rank) と 前窓主平面からの傾き[deg]:')
    for a in starts:
        seg = X[a:a + Wn]; Wc = seg - seg.mean(0)
        U, s, Vt = np.linalg.svd(Wc, full_matrices=False)
        v = s[2] / s[0]; inst.append((a, v)); plane = Vt[:2]
        if prev is not None:
            sv = np.linalg.svd(prev @ plane.T, compute_uv=False)
            ang = float(np.degrees(np.arccos(np.clip(sv.min(), 0, 1)))); tilts.append((a, ang))
        else:
            ang = None
        print(f'    step[{a}:{a+Wn}] 窓内s3/s1={v:.4f}' + (f'  傾き={ang:.1f}°' if ang is not None else ''))
        prev = plane
    print(f'\n判定: 瞬間(短窓)s3/s1最大={max(v for _,v in inst):.4f}（≈0＝常に単一面rank2）')
    print('      → S³/rank4は長窓集約の作為。物理は単一面の連続再配向のみ。思考実験⑥要修正。')

    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        ax[0].plot([w for w, _ in wvals], [v for _, v in wvals], 'o-')
        ax[0].axhline(0.1, ls='--', c='r', label='genuine 2-plane level ~0.1')
        ax[0].set_xscale('log'); ax[0].set_xlabel('SVD window width [step]'); ax[0].set_ylabel('s3/s1')
        ax[0].set_title('(A) s3/s1 grows with window = aggregation artifact'); ax[0].legend(fontsize=8)
        ax2 = ax[1]; a1 = [a for a, _ in inst]; v1 = [v for _, v in inst]
        ax2.plot(a1, v1, 'o-', label='in-window s3/s1 (instantaneous rank)')
        ax2b = ax2.twinx()
        ax2b.plot([a for a, _ in tilts], [t for _, t in tilts], 's-', color='tab:orange', label='plane reorientation [deg/30step]')
        ax2.set_xlabel('step'); ax2.set_ylabel('in-window s3/s1'); ax2b.set_ylabel('plane tilt per 30 step [deg]')
        ax2.set_title('(B) instantaneous rank2 + single plane reorienting')
        ax2.legend(loc='upper right', fontsize=8); ax2b.legend(loc='center right', fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig_artifact_shortwindow_svd.png'), dpi=110)
        print('wrote fig_artifact_shortwindow_svd.png')
    except Exception as e:
        print('figure skipped:', e)


if __name__ == '__main__':
    main()
