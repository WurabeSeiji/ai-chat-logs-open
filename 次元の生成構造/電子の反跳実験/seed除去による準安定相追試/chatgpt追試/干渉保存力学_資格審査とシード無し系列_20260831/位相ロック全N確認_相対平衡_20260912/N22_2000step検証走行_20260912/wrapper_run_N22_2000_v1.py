#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=22 単独・den=N(=22, Δτ=2π/22)のみ・2000 step 走行ラッパー。

目的: 500step で未ロックだった N=22 を 2000step まで無停止走行し、
      予測（../未ロックN_2000step到達予測_20260912、N=22 lockstep≈812）を検証する。

規約（追加実験は元プログラムの忠実コピーのみ）:
  N=64 版 ../../make_parent型初期値_自己無撞着構造_20260907/N64_2000step走行_20260907/wrapper_run_N64_2000_v1.py
  を忠実コピーし、置換対象だけを N=22 用に変えた。物理正本
  ../../N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py（SHA256=1abf2353… を走行前に厳密照合）
  の one_step / H_of / adjacency / plane / metrics の数式・dtype・保存形式は一切変更しない。

★N=64 版との唯一の設計差:
  N=64 は元スイープに親が無く自前生成したが、N=22 の親は元スイープ parents/ に既存。
  500step 走行と同一の初期条件で忠実延長するため、PARENT_DIR は元のまま（再生成しない・置換しない）。

置換（各行 count==1 をアサート）:
 1. OUT → 本フォルダ results/
 2. STEPS=500; OFFSETS=(-2,-1,0,1,2) → STEPS=2000; OFFSETS=(0,)   ← 歩数と分母系列（den=N のみ）
 3. pairs 行の「+ [(124,'124')]」を除去（den=124 を走らせない）
 4. for N in range(3,41): → for N in [22]:
 5. 保存名 states_500.npz → states_2000.npz
 6. RUN_METADATA の 'N_range':[3,40] → [22,22]
 7. 10歩ごとの進捗print追加（物理式は無変更）
注: 正本内蔵の 8x5 俯瞰図は N=3..40 ループ前提のため空図になる（図はこのフォルダの分析で作る）。
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
    ("for N in range(3,41):", "for N in [22]:"),
    ("hm_N{N}_den_{den}_states_500.npz", "hm_N{N}_den_{den}_states_2000.npz"),
    ("'N_range':[3,40]", "'N_range':[22,22]"),
    ("            if t<STEPS: z=one_step(z,A,den)",
     "            if t%10==0: print(f'[progress] N={N} den={den} step {t}/{STEPS}',flush=True)\n"
     "            if t<STEPS: z=one_step(z,A,den)"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所（PARENT_DIR は元のまま＝500step走行と同一親）')
g = {'__name__': 'stage123_N22_2000', '__file__': ORIG}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
