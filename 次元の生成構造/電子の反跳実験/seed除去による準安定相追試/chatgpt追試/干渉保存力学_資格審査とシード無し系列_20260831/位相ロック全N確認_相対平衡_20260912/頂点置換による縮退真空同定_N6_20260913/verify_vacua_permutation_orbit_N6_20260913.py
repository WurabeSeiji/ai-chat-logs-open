#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""縮退真空多様体（谷）が有限離散な写像対称性軌道か連続モジュライかの決定的テスト（N=6、2026-09-13）。

論文5 の SSB（メキシカンハット型：不安定な対称極大→縮退した谷→微小シードの自発選択）は実測事実。
本テストは縮退真空多様体（谷）の構造を特定する：谷は (A) 有限離散な写像対称性軌道（S_N＋大域位相）か、
(B) より大きな（連続）モジュライ（傾いた相対平衡の連続族＝メキシカンハットの谷の一般形）か。
どちらでも SSB は成り立つ（谷は本来 連続縮退多様体で有限離散軌道である必要はない）。(A) に尽きるかを判定する。

方法（N=6・6!=720 頂点置換を総当たり・読出しのみ・物理無変更）:
  1. make_parent 床 Z0（full_N3_N40_sweep step0）に向きの異なる微小シード5通り（batch と同一 seed
     rng=default_rng(1000*N+seed), amp=1e-8）を与え、正本 one_step で相対平衡へ強収束させる。
  2. 各終状態 z^(a) を、頂点置換 P∈S_6 が誘導する辺置換で写し、大域位相 θ を最適化して
     res(a,b)=min_P min_θ ‖z^(b) - e^{iθ} P z^(a)‖/‖z^(b)‖ を全ペアで測る。
  3. res~機械精度なら「異なる真空は写像対称性軌道上（厳密縮退＝SSB）」、
     res~O(1e-3以上)で改善しないなら「seed 依存多重アトラクタ（粗視化不変量のみ縮退）」。
  補助: 複素共役（時間反転的対称）も許した場合の残差も併記。

正本 run_N3_N40_stage123_v1.py（SHA照合）の adjacency/one_step を import（物理無変更）。
"""
import hashlib
import itertools
import math
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
STATES = os.path.join(BASE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907',
                      'full_N3_N40_sweep', 'states')
RESULTS = os.path.join(BASE, 'results')
os.makedirs(RESULTS, exist_ok=True)

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
g = {'__name__': 'canonical_funcs', '__file__': ORIG}
exec(compile(src.split('rows=[]; summaries=[]')[0], ORIG, 'exec'), g)
adjacency, one_step = g['adjacency'], g['one_step']

N = 6
NSEED = 5
SEED_AMP = 1e-8
STEPS = 20000          # 相対平衡へ強収束（batch の 3e-4 lock よりはるかに深く）
log = []


def w(*a):
    s = ' '.join(str(x) for x in a); print(s); log.append(s)


def rel_config_rate(z, zprev, k):
    rel = np.angle(z * np.conj(z[0])); relp = np.angle(zprev * np.conj(zprev[0]))
    return np.median(np.abs((rel - relp + np.pi) % (2 * np.pi) - np.pi)) / k


A = adjacency(N)
ea, eb = np.triu_indices(N, 1)
M = len(ea)
edge_index = {}
for e in range(M):
    edge_index[(int(ea[e]), int(eb[e]))] = e

Z0 = np.asarray(np.load(os.path.join(STATES, f'hm_N{N}_den_{N}_states_500.npz'))['Z'][0], np.complex128)

# 1. 5 終状態を強収束
vac = []
w(f'# 縮退真空の頂点置換軌道テスト  N={N}  M={M}  日付 2026-09-13')
w(f'# 正本SHA一致 {PROG_SHA}')
w(f'# Z0=make_parent床 step0, seed=default_rng(1000N+s), amp={SEED_AMP}, one_step×{STEPS}')
w('')
for seed in range(NSEED):
    rng = np.random.default_rng(1000 * N + seed)
    z = Z0 + SEED_AMP * (rng.standard_normal(M) + 1j * rng.standard_normal(M))
    z = z / np.linalg.norm(z) * np.linalg.norm(Z0)
    for t in range(1, STEPS + 1):
        z = one_step(z, A, N)
    # 収束度: 直近1stepの相対配置変化率（相対平衡なら機械ゼロ）
    z1 = one_step(z, A, N)
    rate = rel_config_rate(z1, z, 1)
    r = np.abs(z); amp_ratio = r.max() / r.min()
    vac.append(z.copy())
    w(f'seed {seed}: amp_ratio={amp_ratio:.8f}  相対配置変化率(1step)={rate:.2e} rad/step')

w('')
# 2. 720 頂点置換の誘導辺置換を前計算
perms = list(itertools.permutations(range(N)))
edge_perm = []       # 各頂点置換 P に対する辺置換 dst[e]=Pが辺eを送る先の辺index
for P in perms:
    dst = np.empty(M, dtype=int)
    for e in range(M):
        i, j = int(ea[e]), int(eb[e])
        a2, b2 = P[i], P[j]
        if a2 > b2:
            a2, b2 = b2, a2
        dst[e] = edge_index[(a2, b2)]
    edge_perm.append(dst)


def best_match(za, zb, allow_conj=False):
    """min over 720 vertex-perms (and optional conjugation) and global phase of ‖zb - e^{iθ} P za‖/‖zb‖."""
    best = (np.inf, None, False)
    nb = np.linalg.norm(zb)
    cands = [(za, False)] + ([(np.conj(za), True)] if allow_conj else [])
    for zc, conj in cands:
        for pi, dst in enumerate(edge_perm):
            zp = np.empty(M, dtype=complex)
            zp[dst] = zc            # (P z)[dst[e]] = z[e]
            ip = np.vdot(zp, zb)     # 最適大域位相 θ=arg(ip)
            res = np.linalg.norm(zb - (ip / abs(ip)) * zp) / nb if abs(ip) > 0 else np.inf
            if res < best[0]:
                best = (res, perms[pi], conj)
    return best


# 3. 全ペアで最小残差
w('## 頂点置換＋大域位相のみ（写像対称性）')
w('pair | min残差 | 達成置換 P')
pair_res = {}
for a in range(NSEED):
    for b in range(a + 1, NSEED):
        res, P, _ = best_match(vac[a], vac[b], allow_conj=False)
        pair_res[(a, b)] = res
        w(f'  {a}-{b} | {res:.3e} | P={P}')
w('')
w('## 参考: 頂点置換＋大域位相＋複素共役も許した場合')
for a in range(NSEED):
    for b in range(a + 1, NSEED):
        res, P, conj = best_match(vac[a], vac[b], allow_conj=True)
        w(f'  {a}-{b} | {res:.3e} | P={P} conj={conj}')

worst = max(pair_res.values())
best_ = min(pair_res.values())
w('')
w(f'置換のみ: 全ペア最小残差の最大={worst:.3e}  最小={best_:.3e}')
THRESH = 1e-6
if worst < THRESH:
    verdict = f'全ペアが機械精度({worst:.1e}<{THRESH})で頂点置換＋大域位相により一致 → (A) 谷は有限離散な写像対称性軌道'
else:
    verdict = (f'頂点置換＋大域位相では一致しない（最大残差{worst:.1e}≥{THRESH}）→ (B) 谷は有限離散S_N軌道でなく連続モジュライ'
               '（傾いた相対平衡の連続族）。SSBを否定せず整合（メキシカンハットの谷は本来 連続縮退多様体）。有限離散軌道仮説(A)のみ排除')
w('=== 判定:', verdict, '===')

with open(os.path.join(RESULTS, '実行ログ_20260913.log'), 'w', encoding='utf-8') as f:
    f.write('\n'.join(log) + '\n')
