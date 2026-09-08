#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""合成親実験（2026-09-08、木原指示）: 同一の高対称理論床 2 つの合成 k+k（k=3,4,5,6 → N=6,8,10,12）。

step0 処理（親構成時に 1 回だけ実施 — 以後は波の増減がないため毎 step 因数分解は不要）:
 1. 各理論床 k を頂点位相へ因数分解 φ_i（正準分岐 φ0∈[0,π)、決定論・全単射・走行中固定の正準ゲージ）
 2. 頂点位相閉塞 Σe^{iφ}=Σe^{2iφ}=0 を確認（クロス閉塞の因数分解条件）
 3. 合成 N=2k: ブロック内は元の床の値そのまま、クロス波は系譜位相で ε·exp(i(φ_i+φ_j))
    （ε=1e-8、接触ゲージ=0 デフォルト宣言）。全パワー 1 に正規化。
 4. 監査: 全体・クロスセクターの |Σz|・|Σz²| が機械ゼロであること
走行: 物理正本 run_N3_N40_stage123_v1.py（SHA照合・無変更）を den=N・500 step で exec。
ε>0 なので全スロットが生きており、正本経路のみ（系譜供給分岐は不要）。
"""
import hashlib
import math
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
PAR_SRC = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'parents_symmetric_staged')

EPS = 1e-8      # クロス波の誕生振幅（正規化前）
GAUGE = 0.0     # 凝縮体間の接触相対ゲージ（正準デフォルト 0・宣言済みパラメータ）

pd = os.path.join(BASE, 'parents_composite')
out = os.path.join(BASE, 'results_composite')
os.makedirs(pd, exist_ok=True)
os.makedirs(out, exist_ok=True)

# ---- step0 処理: 合成親の生成（因数分解は親構成時の 1 回のみ）----
for k in (3, 4, 5, 6):
    zk = np.array(np.load(os.path.join(PAR_SRC, f'parent_static_N{k:05d}_makeparent_20260905.npz'))['Z0'],
                  dtype=np.complex128, copy=True)
    ea, eb = np.triu_indices(k, k=1)
    idx = {(int(ea[e]), int(eb[e])): e for e in range(len(ea))}
    th = np.angle(zk)
    phi = np.zeros(k)
    phi0 = 0.5 * (th[idx[(0, 1)]] + th[idx[(0, 2)]] - th[idx[(1, 2)]]) % math.pi   # 正準分岐 [0,π)
    phi[0] = phi0
    for j in range(1, k):
        phi[j] = th[idx[(0, j)]] - phi0
    err = max(abs(np.exp(1j * (phi[int(ea[e])] + phi[int(eb[e])])) - np.exp(1j * th[e])) for e in range(len(ea)))
    assert err < 1e-9, f'k={k}: 床が和型位相に因数分解できない err={err:.3e}'
    s1 = abs(np.sum(np.exp(1j * phi)))
    s2 = abs(np.sum(np.exp(2j * phi)))
    assert s1 < 1e-9 and s2 < 1e-9, f'k={k}: 頂点位相閉塞が不成立 |Σe^iφ|={s1:.3e} |Σe^2iφ|={s2:.3e}'

    N = 2 * k
    EA, EB = np.triu_indices(N, k=1)
    Phi = np.concatenate([phi, phi + GAUGE])
    z = np.empty(len(EA), dtype=np.complex128)
    for e in range(len(EA)):
        i, j = int(EA[e]), int(EB[e])
        if i < k and j < k:
            z[e] = zk[idx[(i, j)]]
        elif i >= k and j >= k:
            z[e] = zk[idx[(i - k, j - k)]]
        else:
            z[e] = EPS * np.exp(1j * (Phi[i] + Phi[j]))      # 系譜位相の誕生
    z /= np.linalg.norm(z)                                    # 全パワー 1
    cross = (EA < k) != (EB < k)
    print(f'k={k}→N={N} M={len(EA)}: 位相分解誤差={err:.2e} 頂点|Σe^iφ|={s1:.2e} |Σe^2iφ|={s2:.2e} | '
          f'全体|Σz|={abs(z.sum()):.2e} |Σz²|={abs((z * z).sum()):.2e} '
          f'クロス|Σz|={abs(z[cross].sum()):.2e} |Σz²|={abs((z[cross] ** 2).sum()):.2e}', flush=True)
    np.savez_compressed(os.path.join(pd, f'parent_composite_N{N:05d}_kk_20260908.npz'),
                        Z0=z, k=np.int64(k), eps=np.float64(EPS), gauge=np.float64(GAUGE))
print('合成親 4 本生成完了', flush=True)

# ---- 走行: 正本を SHA 照合し宣言行のみ置換して exec（物理無変更・全スロット生存＝正本経路のみ）----
src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'

REPL = [
    ("OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')", f'OUT={out!r}'),
    ("PARENT_DIR=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','parents')", f'PARENT_DIR={pd!r}'),
    ("for N in range(3,41):", "for N in [6,8,10,12]:"),
    ("parent_static_N{N:05d}_makeparent_20260905.npz", "parent_composite_N{N:05d}_kk_20260908.npz"),
    ("pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0] + [(124,'124')]",
     "pairs=[(N,'N')]"),
    ("fig,axs=plt.subplots(8,5,figsize=(20,24)); axs=axs.ravel(); order=['N-2','N-1','N','N+1','N+2','124']",
     "fig,axs=plt.subplots(1,4,figsize=(18,4.5)); axs=np.atleast_1d(axs).ravel(); order=['N']"),
    ("for k,N in enumerate(range(3,41)):", "for k,N in enumerate([6,8,10,12]):"),
    ("ax.set_xlim(0,500); ax.set_ylim(1e-34,3); ax.set_title(f'N={N}'); ax.grid(alpha=.25)",
     "ax.set_xlim(0,500); ax.set_ylim(1e-34,3); ax.set_title(f'N={N} = {N//2}+{N//2} composite'); ax.grid(alpha=.25)"),
    ("for k in range(38,40): axs[k].axis('off')", "pass"),
    ("Hperp/H denominator control (float64/complex128): 2pi/(N-2), 2pi/(N-1), 2pi/N, 2pi/(N+1), 2pi/(N+2), 2pi/124; N=3..40",
     "Hperp/H composite parents k+k (k=3,4,5,6 -> N=6,8,10,12): dtau=2pi/N, eps=1e-8, contact gauge 0"),
    ("'N_range':[3,40]", "'N_range':[6,12]"),
    ("'denominators':'N-2,N-1,N,N+1,N+2,124'", "'denominators':'N'"),
    ("'input':'per-N static make_parent parents/parent_static_N*_makeparent_20260905.npz Z0 (N=40 bit-identical to canonical)'",
     "'input':'composite parents parent_composite_N*_kk_20260908.npz Z0 (two identical theoretical floors k+k, genealogical cross phases, eps=1e-8, contact gauge 0)'"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所（出力先・親・N集合・図体裁のみ、物理無変更）', flush=True)
g = {'__name__': 'stage123_composite_kk', '__file__': ORIG}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
