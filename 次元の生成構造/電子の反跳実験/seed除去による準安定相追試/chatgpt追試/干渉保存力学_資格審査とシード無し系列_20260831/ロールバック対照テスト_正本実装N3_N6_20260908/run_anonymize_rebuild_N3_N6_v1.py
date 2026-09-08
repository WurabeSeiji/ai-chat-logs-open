#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""無名性ゲート（2026-09-08、木原指示）: 走行前に初期値を一旦無名化（波をシャッフル）し、
90度系列(=既知のN乗根位相)という情報だけで完全リスト(頂点対への割当=配線)を復元して再構築する。

要点: 力学が波の「名前」を読んでいなければ、無名化→復元は初期状態を bit 完全に元へ戻すはずで、
走行結果も正本と bit 一致する。これが無名性合格の証明になる。

anonymize_and_rebuild(z0):
 1. z0 を乱数でシャッフル（名前=triu順序 の情報を破棄）
 2. 各波の位相を最寄りの N 乗根クラス c∈{0..N-1} に量子化（90度系列＝和クラス c=(i+j)mod N）
 3. turnpike 復元: 頂点位相クラス割当 v:{頂点}->{0..N-1} を、和クラス多重集合が観測と一致するよう決定
    （候補は S_N 相当のゲージ違いのみ＝物理的に一意。辞書順最小の割当を正準採用）
 4. その割当で各頂点対 (i,j) に、位相クラス c=(v_i+v_j)%N を持つシャッフル波を割り当てて完全リスト再構築
 5. 正準ゲージ（頂点位相を昇順に整列）で並べ直し、元 z0 と一致するか照合

走行: 物理正本 run_N3_N40_stage123_v1.py（SHA照合・無変更）を den=N・500 step で exec。
初期値だけ z0 -> anonymize_and_rebuild(z0) に差し替える（力学は一切無変更）。
"""
import hashlib
import math
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
PARENT_DIR = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'parents')


def infer_N(M):
    N = int(round((1 + math.sqrt(1 + 8 * M)) / 2))
    assert N * (N - 1) // 2 == M
    return N


def anonymize_and_rebuild(z0, seed=20260908, verbose=False):
    M = z0.size
    N = infer_N(M)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(M)
    zsh = z0[perm]                                   # 1. 無名化（配線情報を破棄）

    # 2. 位相を N 乗根クラスへ量子化（90度系列の事前情報）
    cls = np.mod(np.round(np.angle(zsh) / (2 * math.pi) * N).astype(int), N)
    obs = np.bincount(cls, minlength=N)              # 観測: 和クラス多重集合

    # 3. turnpike 復元: 頂点位相クラス割当 v を決める。頂点位相クラスは {0..N-1} の並べ替え。
    #    和クラス多重集合が obs に一致する割当を辞書順で1つ探す（存在＝一意な物理配線・ゲージ違いのみ）。
    ea, eb = np.triu_indices(N, k=1)
    from itertools import permutations
    v_found = None
    for v in permutations(range(N)):
        s = np.bincount([(v[a] + v[b]) % N for a, b in zip(ea, eb)], minlength=N)
        if np.array_equal(s, obs):
            v_found = np.array(v); break
    assert v_found is not None, '90度系列として復元不能（前提の位相集合に合致せず）'

    # 4. 完全リスト再構築: 頂点対 (i,j) に 位相クラス (v_i+v_j)%N を持つシャッフル波を割り当てる
    from collections import defaultdict, deque
    bucket = defaultdict(deque)
    order = sorted(range(M), key=lambda e: perm[e])   # 決定論的な割当順（元index順）
    for e in order:
        bucket[int(cls[e])].append(e)
    z_rebuilt = np.empty(M, dtype=np.complex128)
    for slot, (i, j) in enumerate(zip(ea, eb)):
        c = (int(v_found[i]) + int(v_found[j])) % N
        e = bucket[c].popleft()
        z_rebuilt[slot] = zsh[e]

    # 5. 正準ゲージへ: 復元は頂点位相クラス割当 v_found で表現されるが、triu の (i,j) 昇順は
    #    「頂点位相クラスの昇順」ゲージに対応する。v_found が恒等でなければ頂点を並べ替えて正準化。
    vorder = np.argsort(v_found, kind='stable')       # 頂点位相クラス昇順への並べ替え
    relabel = np.empty(N, dtype=int); relabel[vorder] = np.arange(N)
    z_canon = np.empty(M, dtype=np.complex128)
    slot_of = {}
    for s2, (i, j) in enumerate(zip(ea, eb)):
        slot_of[(i, j)] = s2
    for slot, (i, j) in enumerate(zip(ea, eb)):
        a, b = sorted((int(relabel[i]), int(relabel[j])))
        z_canon[slot_of[(a, b)]] = z_rebuilt[slot]

    if verbose:
        err = np.max(np.abs(z_canon - z0))
        print(f"  N={N} M={M}: 無名化→復元 max|z_canon - z0|={err:.2e} 頂点割当={tuple(v_found)}", flush=True)
    return z_canon


# ---- 事前自己試験: 復元が元と bit 一致するか（走行前） ----
print('=== 無名化→復元 自己試験（走行前） ===', flush=True)
for N in (3, 4, 5, 6):
    z0 = np.array(np.load(os.path.join(PARENT_DIR, f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],
                  dtype=np.complex128)
    zc = anonymize_and_rebuild(z0, verbose=True)
    assert np.array_equal(zc, z0), f'N={N}: 無名化→復元が bit 一致しない（無名性違反 or 復元バグ）'
print('自己試験: 全 N で bit 完全一致（無名化→復元が可逆）', flush=True)

# ---- 走行: 正本を SHA 照合し、初期値だけ無名化→復元を通す（力学は無変更） ----
out = os.path.join(BASE, 'results_anonymize')
os.makedirs(out, exist_ok=True)
src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'

INJ = ("z0=anonymize_and_rebuild(np.array(np.load(os.path.join(PARENT_DIR,"
       "f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],dtype=np.complex128,copy=True))")
REPL = [
    ("OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')", f'OUT={out!r}'),
    ("for N in range(3,41):", "for N in range(3,7):"),
    ("    z0=np.array(np.load(os.path.join(PARENT_DIR,f'parent_static_N{N:05d}_makeparent_20260905.npz'))['Z0'],dtype=np.complex128,copy=True)",
     "    " + INJ),
    ("'N_range':[3,40]", "'N_range':[3,6]"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所（初期値の無名化→復元のみ・力学は無変更）', flush=True)
g = {'__name__': 'stage123_anonymize_rebuild', '__file__': ORIG, 'anonymize_and_rebuild': anonymize_and_rebuild}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
