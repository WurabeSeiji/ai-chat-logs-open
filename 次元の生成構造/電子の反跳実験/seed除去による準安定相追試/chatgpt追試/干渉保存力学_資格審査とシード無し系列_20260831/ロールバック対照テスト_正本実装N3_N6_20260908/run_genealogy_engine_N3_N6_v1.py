#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""系譜位相供給エンジンの回帰テスト（2026-09-08、木原指示）: N=3..6（M=3,6,10,15）・500 step。

変更は 1 点のみ: 生成子構成時の位相を genealogy_phases(z) 経由にする。
 - 振幅ゼロの波が無い場合: np.exp(1j*np.angle(z)) をそのまま返す（正本と bit 同一経路）
 - 振幅ゼロの波がある場合: 生きた波を N×N 対称行列に並べ、Takagi 型ランク1因数分解で
   頂点振幅 a_i を読み出し、死スロットの位相を系譜合成 φ_i+φ_j で埋める（規約・seed 混入なし）。
飽和親（既存 N=3..6 親は全波が非ゼロ）では後者は一度も発火しないので、正本と bit 一致が合格条件。
走行後に、合成型の模擬状態（積状態のクロス9波をゼロ化）で位相供給の単体自己試験も行う。

方式: 物理正本 run_N3_N40_stage123_v1.py（SHA照合）を読み取り、宣言行置換＋上記注入のみで exec。
"""
import hashlib
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'

out = os.path.join(BASE, 'results_genealogy')
os.makedirs(out, exist_ok=True)

INJ = '''def genealogy_phases(z):
    # 系譜位相供給: 死スロット（振幅ゼロ）の位相を、生きた波の頂点層因数分解から φ_i+φ_j で埋める。
    dead=(np.abs(z)==0.0)
    if not dead.any():
        return np.exp(1j*np.angle(z))          # 正本と bit 同一経路
    M=z.size; N=int(round((1.0+math.sqrt(1.0+8.0*M))/2.0))
    ea,eb=np.triu_indices(N,k=1)
    Zm=np.zeros((N,N),dtype=np.complex128); Zm[ea,eb]=z; Zm[eb,ea]=z
    live=~dead
    # ランク1補完の反復: 対角（自己波なし）と死スロットを現推定 v_i v_j g2 で埋めて再分解。
    # 生きグラフが連結（非二部）なら系譜位相 φ_i+φ_j に収束。非連結なら成分間の相対ゲージは
    # 原理的に非決定（物理: 接触前の凝縮体間ゲージ）で、呼び出し側がゲージ供給する必要がある。
    g2=1.0+0j; v=np.zeros(N,dtype=np.complex128)
    for _ in range(20):
        Wm=Zm@np.conj(Zm)
        wv,Vv=np.linalg.eigh(Wm); v=Vv[:,-1]*(wv[-1] if wv[-1]>0 else 0.0)**0.25
        num=np.sum(np.conj(v[ea[live]]*v[eb[live]])*z[live])
        g2=num/abs(num) if abs(num)>0 else 1.0+0j
        np.fill_diagonal(Zm,v*v*g2)
        Zm[ea[dead],eb[dead]]=v[ea[dead]]*v[eb[dead]]*g2
        Zm[eb[dead],ea[dead]]=Zm[ea[dead],eb[dead]]
    prod=v[ea]*v[eb]*g2
    outp=np.exp(1j*np.angle(z))
    outp[dead]=np.exp(1j*np.angle(prod[dead]))
    return outp
def one_step(z,A,den):'''

REPL = [
    ("OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')", f'OUT={out!r}'),
    ("for N in range(3,41):", "for N in range(3,7):"),
    ("'N_range':[3,40]", "'N_range':[3,6]"),
    ("def one_step(z,A,den):", INJ),
    ("    H=H_of(np.exp(1j*np.angle(z)),A); H=(1j*np.imag(H)).astype(np.complex128,copy=False)",
     "    H=H_of(genealogy_phases(z),A); H=(1j*np.imag(H)).astype(np.complex128,copy=False)"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所（系譜位相供給の注入のみ・力学は無変更）')
g = {'__name__': 'stage123_genealogy_engine', '__file__': ORIG}
exec(compile(new, ORIG, 'exec'), g)

# 単体自己試験: 積状態 z_ij=a_i a_j の一部をゼロ化し genealogy_phases の系譜位相供給を確認
rng = np.random.default_rng(20260908)
N = 6
av = np.exp(1j * rng.uniform(0, 2 * np.pi, N)) * rng.uniform(0.8, 1.2, N)
ea, eb = np.triu_indices(N, k=1)
z_full = av[ea] * av[eb]
cross = (ea < 3) != (eb < 3)

# ケースA: クロス9波のうち4本だけ死（生きグラフ連結）→ 系譜位相の厳密回復を期待
deadA = np.zeros(len(z_full), dtype=bool); deadA[np.flatnonzero(cross)[:4]] = True
zA = z_full.copy(); zA[deadA] = 0.0
phA = g['genealogy_phases'](zA)
errA = np.max(np.abs(np.exp(1j * np.angle(z_full[deadA])) - phA[deadA]))
errA_live = np.max(np.abs(np.exp(1j * np.angle(zA[~deadA])) - phA[~deadA]))
print(f'自己試験A（連結・死4/15）: 系譜位相誤差 max|Δ|={errA:.3e} / 生き波不変 max|Δ|={errA_live:.3e}')

# ケースB: クロス9波を全て死（生きグラフ非連結=接触前の2凝縮体）→ 相対ゲージ非決定のデモ
zB = z_full.copy(); zB[cross] = 0.0
phB = g['genealogy_phases'](zB)
errB = np.max(np.abs(np.exp(1j * np.angle(z_full[cross])) - phB[cross]))
print(f'自己試験B（非連結・死9/15）: 誤差 max|Δ|={errB:.3e} — 非決定は物理（凝縮体間相対ゲージは状態から読めない。合成親側でゲージ供給が必要）')
print('WRAPPER ALL DONE', flush=True)
