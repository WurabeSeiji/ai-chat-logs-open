#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""残り未飽和N（N=30,31,…,39）単独・den=N のみ・2000 step バッチ走行ラッパー。

目的: 500step で未ロックだった 12 本のうち N=22・N=40 は検証済み。残る N=30〜39 を 2000step まで
      無停止走行し、予測（../未ロックN_2000step到達予測_20260912）を一括検証する。

規約（追加実験は元プログラムの忠実コピーのみ）:
  N=40 版 ../N40_2000step検証走行_20260912/wrapper_run_N40_2000_v1.py を忠実コピーし、
  走らせる N の集合だけを [40] → [30,31,…,39] に変えた（物理正本の元ループ for N in range(3,41) を
  そのまま部分集合に絞る最小変更）。物理正本
  ../../N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py（SHA256=1abf2353… を走行前に厳密照合）
  の数式・dtype・保存形式は一切変更しない。各Nの親は元スイープ parents/ に既存のため PARENT_DIR は元のまま。
  出力は本フォルダ results/ に hm_N{N}_den_{N}_states_2000.npz として N ごとに保存される。

置換（各行 count==1 をアサート）:
 1. OUT → 本フォルダ results/
 2. STEPS=500; OFFSETS=(-2,-1,0,1,2) → STEPS=2000; OFFSETS=(0,)
 3. pairs 行の「+ [(124,'124')]」を除去（den=124 を走らせない）
 4. for N in range(3,41): → for N in [30,31,32,33,34,35,36,37,38,39]:
 5. 保存名 states_500.npz → states_2000.npz
 6. RUN_METADATA の 'N_range':[3,40] → [30,39]
 7. 10歩ごとの進捗print追加（物理式は無変更）
"""
import hashlib
import os

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'

out = os.path.join(BASE, 'results')
os.makedirs(out, exist_ok=True)

REPL = [
    ("OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')", f'OUT={out!r}'),
    ("STEPS=500; OFFSETS=(-2,-1,0,1,2)", "STEPS=2000; OFFSETS=(0,)"),
    ("pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0] + [(124,'124')]",
     "pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0]"),
    ("for N in range(3,41):", "for N in [30,31,32,33,34,35,36,37,38,39]:"),
    ("hm_N{N}_den_{den}_states_500.npz", "hm_N{N}_den_{den}_states_2000.npz"),
    ("'N_range':[3,40]", "'N_range':[30,39]"),
    ("            if t<STEPS: z=one_step(z,A,den)",
     "            if t%10==0: print(f'[progress] N={N} den={den} step {t}/{STEPS}',flush=True)\n"
     "            if t<STEPS: z=one_step(z,A,den)"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所（PARENT_DIR は元のまま。N=30..39 バッチ）', flush=True)
g = {'__name__': 'stage123_N30_N39_2000', '__file__': ORIG}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
