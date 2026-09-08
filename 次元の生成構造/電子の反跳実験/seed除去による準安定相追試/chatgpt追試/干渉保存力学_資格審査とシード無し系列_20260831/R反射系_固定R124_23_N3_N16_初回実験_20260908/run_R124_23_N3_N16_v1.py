#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""固定 R_{124,23} 反射散乱 × 現行 Stage1+2+3（対称分割合成）N=3..16・den=N・500 step。

指図書 ClaudeCode_R反射系_固定R124_23_N3_N16_初回実験_指図書_20260908.md と
研究ノート 固定R_23-124交換散乱と現行段123回転の合成方針_研究ノート_20260908.md（正本仕様）準拠。
木原簡素化指示（2026-09-08）: 資格試験群は省略し、走行＋標準図化＋分析に絞る。

方式: 物理正本 run_N3_N40_stage123_v1.py（SHA照合）を読み取り、以下だけを置換して exec。
 - PARENT_DIR → ../対称親v2_500step走行_20260906/parents_symmetric_staged（高対称理論床）
 - OUT → ./results、N範囲 → 3..16、分母 → den=N のみ（STEPS=500 は不変）
 - one_step の呼び出し行を「半散乱 → 正本1歩 → 半散乱」へ（one_step 本体は無変更）
散乱: G = Σ_{A_ef=1} P－^(ef) = ½(D−A)、U_half = exp(iφG/2)、φ=2π·101/124（R=cos²(23π/124) 固定・探索なし）。
"""
import hashlib
import os

BASE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py')
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'

pd = os.path.join(BASE, '..', '対称親v2_500step走行_20260906', 'parents_symmetric_staged')
out = os.path.join(BASE, 'results')
os.makedirs(out, exist_ok=True)

INJ = '''PHI=2.0*math.pi*101.0/124.0          # q=exp(2πi·101/124)、R=cos^2(23π/124) 固定
def build_U_half(A):
    # G_scatter = Σ_{(e,f):A_ef=1} P-（P-=½(|e>-|f>)(<e|-<f|)）を一括合成 = ½(D−A)。順序非依存・状態非依存。
    G=0.5*(np.diag(A.sum(axis=1))-A)
    wg,Vg=np.linalg.eigh(G)
    return ((Vg*np.exp(0.5j*PHI*wg))@Vg.T).astype(np.complex128,copy=False)
def one_step_scatter(z,A,den,U_half):
    z=U_half@z                      # 半散乱
    z=one_step(z,A,den)             # 正本の1歩（無変更）
    return (U_half@z).astype(np.complex128,copy=False)   # 半散乱
rows=[]; summaries=[]'''

REPL = [
    ("PARENT_DIR=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','parents')", f'PARENT_DIR={pd!r}'),
    ("OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')", f'OUT={out!r}'),
    ("STEPS=500; OFFSETS=(-2,-1,0,1,2)", "STEPS=500; OFFSETS=(0,)"),
    ("pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0] + [(124,'124')]",
     "pairs=[(N+o, f'N{o:+d}' if o else 'N') for o in OFFSETS if N+o>0]"),
    ("for N in range(3,41):", "for N in range(3,17):"),
    ("rows=[]; summaries=[]", INJ),
    ("    A=adjacency(N); p,q=plane(z0)", "    A=adjacency(N); p,q=plane(z0); U_half=build_U_half(A)"),
    ("            if t<STEPS: z=one_step(z,A,den)", "            if t<STEPS: z=one_step_scatter(z,A,den,U_half)"),
    ("'N_range':[3,40]", "'N_range':[3,16]"),
]
new = src
for a, b in REPL:
    assert new.count(a) == 1, f'置換対象が一意でない: {a[:60]}'
    new = new.replace(a, b)
print(f'置換適用: {len(REPL)} 箇所')
g = {'__name__': 'stage123_R124_23', '__file__': ORIG}
exec(compile(new, ORIG, 'exec'), g)
print('WRAPPER ALL DONE', flush=True)
