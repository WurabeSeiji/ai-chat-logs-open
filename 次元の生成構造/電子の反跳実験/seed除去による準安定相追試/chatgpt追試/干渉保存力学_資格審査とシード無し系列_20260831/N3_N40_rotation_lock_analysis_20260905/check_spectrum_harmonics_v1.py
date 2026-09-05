#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""初期生成子スペクトルの倍音・組み合わせ共鳴の全系列走査（読み出しのみ）。
問い（木原 2026-09-05）: N=6 あたりから、位相（固有振動数）が倍音関係になる
組み合わせが存在しないか。
方法: 各 N の初期生成子 H=i·K（Z[0] から構成、エンジンは正本の逐語コピー）の
固有振動数 {w_k} について、正の枝 {w>0} で
 (1) 倍音対: w_j = m·w_i (m=2,3,4)
 (2) 組み合わせ共鳴（三波共鳴）: w_i + w_j = w_k (i≤j)
を相対許容 1e-9（厳密級）と 1e-3（近共鳴）の2段で数える。
併せて w² の整数近接（√整数パターンの破れ点）と主固有値の系列を記録。
入力は正本スイープ npz（SHA256 台帳照合）。
出力: check_spectrum_harmonics_v1.json"""
import hashlib
import json
import math
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')

ledger = {}
with open(os.path.join(PKG, 'SHA256SUMS.txt')) as f:
    for line in f:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

# ---- エンジン（run_N3_N40_stage123_v1.py 14-27行の逐語コピー） ----
def edges(N):
    a,b=np.triu_indices(N,k=1); return a.astype(np.int64),b.astype(np.int64)
def adjacency(N):
    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    return A
def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)
# ---- コピーここまで ----

TOL_EXACT = 1e-9
TOL_NEAR = 1e-3
out = {}
for N in range(3, 41):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    assert ledger[rel] == sha256(os.path.join(PKG, rel)), f'INPUT GATE FAIL: {rel}'
    Z0 = np.asarray(np.load(os.path.join(PKG, rel))['Z'][0], dtype=np.complex128)
    A = adjacency(N)
    H = H_of(np.exp(1j * np.angle(Z0)), A); H = (1j * np.imag(H)).astype(np.complex128, copy=False)
    w = np.linalg.eigvalsh(H)
    pos = np.sort(w[w > 1e-12])
    scale = float(pos.max()) if pos.size else 1.0

    def scan(tol):
        mult = []
        for m in (2, 3, 4):
            for i in range(len(pos)):
                for j in range(len(pos)):
                    if abs(pos[j] - m * pos[i]) < tol * scale:
                        mult.append((m, float(pos[i]), float(pos[j])))
        tri = []
        for i in range(len(pos)):
            for j in range(i, len(pos)):
                s = pos[i] + pos[j]
                k = np.searchsorted(pos, s)
                for kk in (k - 1, k):
                    if 0 <= kk < len(pos) and abs(pos[kk] - s) < tol * scale:
                        tri.append((float(pos[i]), float(pos[j]), float(pos[kk])))
        return mult, tri

    mult_e, tri_e = scan(TOL_EXACT)
    mult_n, tri_n = scan(TOL_NEAR)
    wsq = pos ** 2
    wsq_int = [(float(x), int(round(x)), abs(x - round(x))) for x in wsq]
    n_int = sum(1 for _, _, d in wsq_int if d < 1e-9)
    out[str(N)] = {
        'n_pos_modes': int(pos.size),
        'w_max': float(pos.max()),
        'positive_spectrum_first12': [float(x) for x in pos[-12:]],
        'wsq_integer_count_1e-9': n_int,
        'wsq_all_integer': bool(n_int == pos.size),
        'harmonic_pairs_exact': mult_e[:10], 'n_harmonic_pairs_exact': len(mult_e),
        'triad_resonances_exact': tri_e[:10], 'n_triad_exact': len(tri_e),
        'n_harmonic_pairs_near_1e-3': len(mult_n),
        'n_triad_near_1e-3': len(tri_n),
    }
    print(f"N={N}: modes+={pos.size} wsq_int={n_int}/{pos.size} "
          f"harm_exact={len(mult_e)} triad_exact={len(tri_e)} "
          f"harm_1e-3={len(mult_n)} triad_1e-3={len(tri_n)}", flush=True)

with open(os.path.join(BASE, 'check_spectrum_harmonics_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
print('ALL DONE')
