#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=64 単独・den=N(=64, Δτ=2π/64)のみ・2000 step 走行ラッパー（未実行で保存、実行は指示後）。

物理正本 ../../N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py
（SHA256=1abf2353… を走行前に厳密照合）を読み取り、以下の宣言行だけ을 文字列置換して exec する。
one_step / H_of / adjacency / plane / metrics の数式・dtype・保存形式は一切変更しない。

置換（各行 count==1 をアサート、置換行数=6 を検証）:
 1. PARENT_DIR → 本フォルダ parents_N64（N=64 親、資格審査済み）
 2. OUT        → 本フォルダ results/
 3. STEPS=500; OFFSETS=(-2,-1,0,1,2) → STEPS=2000; OFFSETS=(0,)   ← 歩数と分母系列（den=N のみ）
 4. pairs 行の「+ [(124,'124')]」を除去（den=124 を走らせない）
 5. for N in range(3,41): → for N in [64]:
 6. 保存名 states_500.npz → states_2000.npz（歩数を正しく表す）
 補: RUN_METADATA の 'N_range':[3,40] → [100,100]（記録の整合）
注: 正本内蔵の 8x5 俯瞰図は N=3..40 ループのため空図になる（図はこのフォルダの図化プログラムで作る）。
"""
import hashlib
import os

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'

pd = os.path.join(BASE, 'parents_N64')
out = os.path.join(BASE, 'results')
os.makedirs(out, exist_ok=True)

REPL = [
    ("PARENT_DIR=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','parents')", f'PARENT_DIR={pd!r}'),
    ("OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')", f'OUT={out!r}'),
    ("STEPS=500; OFFSETS=(-2,-1,0,1,2)", "STEPS=2000; OFFSETS=(0,)"),
    ("pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0] + [(124,'124')]",
     "pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0]"),
    ("for N in range(3,41):", "for N in [64]:"),
    ("hm_N{N}_den_{den}_states_500.npz", "hm_N{N}_den_{den}_states_2000.npz"),
    ("'N_range':[3,40]", "'N_range':[64,64]"),
    # 歩数ログ（木原指示 2026-09-07: 10歩ごとに進捗print。物理式は無変更）
    ("            if t<STEPS: z=one_step(z,A,den)",
     "            if t%10==0: print(f'[progress] N={N} den={den} step {t}\u002f{STEPS}',flush=True)\n"
     "            if t<STEPS: z=one_step(z,A,den)"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所（うち1箇所は進捗print追加で+1行）')
g = {'__name__': 'stage123_N64_2000', '__file__': ORIG}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
