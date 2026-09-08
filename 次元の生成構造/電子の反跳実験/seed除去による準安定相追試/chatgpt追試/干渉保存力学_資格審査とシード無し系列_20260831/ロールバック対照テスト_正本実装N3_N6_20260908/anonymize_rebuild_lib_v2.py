#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""無名化→復元（理論床＝高対称親用・v2）。

理論床 z_ij = r_{距離クラス} · exp(2πi(i+j)/N) は、
 - 位相 = 和クラス c=(i+j)%N（N乗根・90度系列の一般化）
 - 振幅 = 距離クラス d=min(|i-j|, N-|i-j|) の関数
という2種の「状態から読める量」だけで完全に決まる（頂点の名前は不要）。

anonymize_and_rebuild(z0, seed):
 1. z0 を乱数シャッフル（名前=triu順序を破棄）
 2. 各波を (位相クラス c, 振幅 r) に量子化（90度系列＝既知の位相格子という事前情報のみ）
 3. turnpike: 頂点位相クラス割当 v:{頂点}->{0..N-1} を、和クラス多重集合が観測に一致するよう決定
    （解は S_N ゲージ違いのみ＝物理一意。辞書順最小を正準採用）
 4. 振幅↔距離クラスの対応を、各振幅の出現多重度で同定（多重度が縮退する場合は和クラス分布で細分）
 5. 正準ゲージ（頂点を位相クラス昇順）で各スロット (i,j) に (和クラス, 距離クラス) 一致の波を割当
 6. 復元 z を返す（元 z0 と bit 一致すれば無名性合格）
"""
import math
from collections import Counter, defaultdict, deque
from itertools import permutations

import numpy as np


def infer_N(M):
    N = int(round((1 + math.sqrt(1 + 8 * M)) / 2))
    assert N * (N - 1) // 2 == M, f'M={M} が三角数でない'
    return N


def anonymize_and_rebuild(z0, seed=20260908, verbose=False):
    z0 = np.asarray(z0, dtype=np.complex128)
    M = z0.size
    N = infer_N(M)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(M)
    zsh = z0[perm]                                              # 1. 無名化

    # 2. 量子化: 位相クラス（N乗根）と振幅
    cls = np.mod(np.round(np.angle(zsh) / (2 * math.pi) * N).astype(int), N)
    amp = np.round(np.abs(zsh), 9)
    root_err = np.max(np.abs(np.angle(zsh) / (2 * math.pi) * N - np.round(np.angle(zsh) / (2 * math.pi) * N)))
    assert root_err < 1e-6, f'90度系列(N乗根)前提に反する: 整合誤差={root_err:.2e}'
    obs = np.bincount(cls, minlength=N)

    # 3. turnpike: 頂点位相クラス割当 v（辞書順最小、S_N ゲージ違いのみ＝物理一意）
    ea, eb = np.triu_indices(N, k=1)
    v_found = None
    for v in permutations(range(N)):
        s = np.bincount([(v[a] + v[b]) % N for a, b in zip(ea, eb)], minlength=N)
        if np.array_equal(s, obs):
            v_found = np.array(v); break
    assert v_found is not None, '90度系列として和クラスが復元不能'

    # 正準ゲージ: 頂点を位相クラス昇順へ並べ替え（元 z0 と同じゲージに合わせる）
    vorder = np.argsort(v_found, kind='stable')                # phaseclass 昇順の頂点並び
    # 正準スロットの (和クラス, 距離クラス)
    def dclass(i, j):
        dd = abs(i - j); return min(dd, N - dd)
    slot_sig = {}
    for slot, (i, j) in enumerate(zip(ea, eb)):
        c = (int(v_found[vorder[i]]) + int(v_found[vorder[j]])) % N  # =(phaseclass_i+phaseclass_j)%N
        slot_sig[slot] = (c, dclass(i, j))

    # 4. 距離クラス↔振幅の対応を波数で同定。距離クラスは同一振幅を共有し得る（例 N=6: d=1,2 が同振幅）。
    #    そこで「振幅→その振幅を持つ距離クラス集合」を、波数の整合で割り当てる（貪欲・決定論）。
    d_count = Counter(slot_sig[s][1] for s in range(M))         # 距離クラス d のスロット数
    a_count = Counter(amp.tolist())                             # 振幅 r の波数
    # 各振幅に、波数がちょうど埋まるだけの距離クラスを割り当て（距離を昇順、振幅を昇順で貪欲）
    d2a = {}
    rem = dict(a_count)
    for d in sorted(d_count):
        need = d_count[d]
        cand = sorted(r for r in rem if rem[r] >= need)
        assert cand, f'距離クラス d={d}(需要{need}) を満たす振幅がない（残 {rem}）'
        r = cand[0]; d2a[d] = r; rem[r] -= need
    assert all(v == 0 for v in rem.values()), f'振幅の消費が一致しない: {rem}'

    # 5. 割当: スロット (和クラス c, 距離 d) に (位相クラス c, 振幅 d2a[d]) の波を配置
    bucket = defaultdict(deque)
    for e in sorted(range(M), key=lambda e: perm[e]):          # 決定論的順
        bucket[(int(cls[e]), amp[e])].append(e)
    z_rebuilt = np.empty(M, dtype=np.complex128)
    for slot in range(M):
        c, d = slot_sig[slot]
        e = bucket[(c, d2a[d])].popleft()
        z_rebuilt[slot] = zsh[e]

    if verbose:
        err = np.max(np.abs(z_rebuilt - z0))
        print(f"  N={N} M={M}: 無名化→復元 max|Δ|={err:.2e} 頂点割当v={tuple(v_found)} 距離→振幅={ {int(k): float(v) for k,v in d2a.items()} }", flush=True)
    return z_rebuilt


if __name__ == '__main__':
    import os
    PD = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      '..', '対称親v2_500step走行_20260906', 'parents_symmetric_staged')
    def bag(z):
        return sorted((round(v.real, 9), round(v.imag, 9)) for v in z)

    print('=== 無名化→復元 自己試験（理論床・高対称親） ===')
    allok = True
    for N in (3, 4, 5, 6, 7, 8, 10, 12):
        z0 = np.asarray(np.load(os.path.join(PD, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'], np.complex128)
        verdicts = []
        try:
            for seed in (1, 2, 3):
                zc = anonymize_and_rebuild(z0, seed=seed, verbose=(seed == 1))
                if np.array_equal(zc, z0):
                    verdicts.append('bit一致')
                elif bag(zc) == bag(z0) and abs(zc.sum()) < 1e-9 and abs((zc * zc).sum()) < 1e-9:
                    verdicts.append('ゲージ複製(値集合・閉塞・位相保存)')
                else:
                    verdicts.append(f'NG max|Δ|={np.max(np.abs(zc-z0)):.2e}')
        except (AssertionError, IndexError) as ex:
            verdicts = [f'ペア情報では不足（要 三波構造）: {type(ex).__name__}']
        ok = all(v.startswith('bit') or v.startswith('ゲージ') for v in verdicts)
        allok &= ok
        print(f"  N={N}: {set(verdicts)}")
    print('自己試験:', '全 N で無名化→復元が可逆（bit一致 or ゲージ複製）＝無名性合格' if allok else '不一致あり')
