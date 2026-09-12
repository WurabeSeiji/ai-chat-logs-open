#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""追試B: ロック終端状態が最小ゼロ閉塞ブロック（2波±i対／3波）へ分解されるか、N=3..40で監査
（読出しのみ・物理無変更・新規走行なし。既存 10000step 終端状態を使用）。

木原指示（2026-09-12）: 論文9（検討）が「終状態のブロック分解の N 依存＝種の目録が S_N 不変量に立つか」を
既存データで検定可能（未実施）とした。既存 10000step 終端状態でこれを実施する。

方法: 各 N の終端状態 z=S[-1]（ロック済み）について
  - 総閉塞 Σz² を確認（保存されているはず）。
  - 全 unordered pair の二波自己閉塞残差 ε_ij=|z_i²+z_j²|/(|z_i|²+|z_j|²) を総当たり。
    exact（ε<1e-9）・near（ε<1e-6）の本数と最小 ε を記録。
  - 対照として同 N の床（step0）でも同じ監査（論文4：床は N=8k で ±i 対を持つ）。
判定: 終端状態が最小ブロックへ分解されるか、床の分解構造が終端で保たれるか／消えるか。
正本 run_N3_N40_stage123_v1.py（SHA照合）から adjacency を import（辺構造の確認用）。
"""
import csv
import hashlib
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.abspath(os.path.join(HERE, '..', '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py'))
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
LONG = os.path.abspath(os.path.join(HERE, '..', '..', 'N3_N40_long10000_20260905', 'results'))
FLOOR = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907', 'full_N3_N40_sweep', 'states'))

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'


def pair_audit(z):
    M = len(z); z2 = z * z
    num = np.abs(z2[:, None] + z2[None, :])
    den = (np.abs(z) ** 2)[:, None] + (np.abs(z) ** 2)[None, :]
    eps = num / den
    iu = np.triu_indices(M, 1)
    e = eps[iu]
    return int(np.sum(e < 1e-9)), int(np.sum(e < 1e-6)), float(e.min()), M


def main():
    rows = []
    for N in range(3, 41):
        M = N * (N - 1) // 2
        st = os.path.join(LONG, f'hm_N{N}_den_{N}_states_10000.npz')
        if not os.path.exists(st):
            continue
        S = np.asarray(np.load(st)['Z'], np.complex128)
        zT = S[-1]; z0 = S[0]
        clT = float(abs(np.sum(zT * zT)) / np.vdot(zT, zT).real)
        cl0 = float(abs(np.sum(z0 * z0)) / np.vdot(z0, z0).real)
        exT, neT, minT, _ = pair_audit(zT)
        ex0, ne0, min0, _ = pair_audit(z0)
        # 等モジュラー度（終端）
        r = np.abs(zT); eqdev = float((r.max() - r.min()) / r.mean())
        rows.append(dict(N=N, M=M, closure_term=clT, closure_floor=cl0,
                         exact_pair_term=exT, near_pair_term=neT, minEps_term=minT,
                         exact_pair_floor=ex0, near_pair_floor=ne0, minEps_floor=min0,
                         term_eqmod_dev=eqdev))
        print(f'N={N:>2} M={M:>4}: 終端 閉塞={clT:.1e} exact±i対={exT} near={neT} minε={minT:.2e} 等モジュラー偏差={eqdev:.1e}'
              f'  | 床 exact±i対={ex0} minε={min0:.2e}')
    with open(os.path.join(HERE, 'results', 'terminal_block_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    n8k = [r for r in rows if r['N'] % 8 == 0]
    print(f'\n総括: 終端 exact±i対を持つN = {[r["N"] for r in rows if r["exact_pair_term"]>0] or "なし"}')
    print(f'      床 exact±i対を持つN = {[r["N"] for r in rows if r["exact_pair_floor"]>0] or "なし"}（論文4：N=8k のはず）')
    print('判定: 終端に exact ブロックが無ければ、最小ブロック分解は床（高対称）の性質であり、'
          'ロック終端の等モジュラー化で消える＝種の目録は床側に宿る。')


if __name__ == '__main__':
    main()
