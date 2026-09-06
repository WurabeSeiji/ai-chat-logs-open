#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""高対称「理論床」親の生成・検証・配置（指示書 2026-09-06 準拠）。

理論の正本: ../最も対称性の高い初期値_20260906/高対称理論床初期値の導出_N3_N40_20260906.md（§7 手順）
  1. 位相を θ_ij = 2π(i+j)/N に固定（巡回 Fourier Ansatz、辺順は np.triu_indices(N,k=1)）
  2. その位相から現行生成子 H = i·A∘sin(θ_f−θ_e) を構成（物理正本 one_step と同一の K）
  3. 巡回距離クラス基底 b_d（n_d^{-1/2}·e^{iθ}）へ射影し Q = B†HB（D×D Hermitian）
  4. 全成分同符号に取れる最下位枝 c を選び r_d = c_d/√n_d
  5. z_ij = r_d·e^{iθ_ij} を構成し、独立検算（‖z‖², Σz, Σz², ‖Hz−λz‖/‖z‖）
検証（指示書 §3）: complex128 / M / ノルム<1e-13 / |Σz|<1e-12 / |Σz²|<1e-12 / 残差<1e-12 目標 /
  λ・r_d が導出 md の表（§5・§6）と一致 / 残差 O(1) なら即停止（走行禁止）。
配置: parents_symmetric_staged/parent_static_N{N:05d}_makeparent_20260905.npz（Z0 ほか従来キー）
台帳: theoretical_floor_parent_manifest.json
使い方: python3 make_parents_theoretical_floor_v1.py N_FIRST N_LAST"""
import hashlib
import json
import math
import os
import re
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
DERIV = os.path.join(BASE, '..', '最も対称性の高い初期値_20260906', '高対称理論床初期値の導出_N3_N40_20260906.md')
DST = os.path.join(BASE, 'parents_symmetric_staged')
MANIFEST = os.path.join(BASE, 'theoretical_floor_parent_manifest.json')

# 導出 md から λ（§5 表）と r_d（§6 一覧）の期待値を読む
md = open(DERIV, encoding='utf-8').read()
LAM = {}
for line in md.splitlines():
    m = re.match(r'\| (\d+) \| \d+ \| \d+ \| (-[\d.]+) \|', line)
    if m:
        LAM[int(m.group(1))] = float(m.group(2))
RD = {}
for m in re.finditer(r'\*\*N=(\d+)\*\*: `r_d = \[([^\]]+)\]`', md):
    RD[int(m.group(1))] = [float(x) for x in m.group(2).split(',')]
assert set(LAM) == set(range(3, 41)) and set(RD) == set(range(3, 41)), (len(LAM), len(RD))


def build_floor(N):
    ea, eb = np.triu_indices(N, k=1)           # 正本と同一の辺順
    M = len(ea)
    d = np.minimum((eb - ea) % N, (ea - eb) % N).astype(int)
    D = N // 2
    theta = 2.0 * math.pi * (ea + eb) / N
    A = np.zeros((M, M))
    for e in range(M):
        share = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e])
        share[e] = False
        A[e, share] = 1.0
    K = A * np.sin(theta[None, :] - theta[:, None])
    H = 1j * K                                  # Hermitian
    n_d = np.array([np.sum(d == dd + 1) for dd in range(D)], dtype=float)
    B = np.zeros((M, D), dtype=np.complex128)
    for dd in range(D):
        mask = (d == dd + 1)
        B[mask, dd] = np.exp(1j * theta[mask]) / math.sqrt(n_d[dd])
    Q = B.conj().T @ H @ B
    w, V = np.linalg.eigh(Q)
    chosen = None
    for k in range(len(w)):                     # 昇順: 最下位から、全成分同符号の枝を探す
        v = V[:, k]
        j = int(np.argmax(np.abs(v)))
        v = v * np.exp(-1j * np.angle(v[j]))    # 大域位相を実正に
        if np.max(np.abs(v.imag)) < 1e-10:
            c = v.real
            if np.all(c > 1e-12) or np.all(c < -1e-12):
                chosen = (float(w[k]), np.abs(c))
                break
    assert chosen is not None, f'N={N}: 同符号最下位枝が見つからない'
    lam, c = chosen
    r_d = c / np.sqrt(n_d)
    Z = (r_d[d - 1] * np.exp(1j * theta)).astype(np.complex128)
    # 独立検算
    norm2 = float(np.vdot(Z, Z).real)
    cz = float(abs(Z.sum()))
    s2 = float(abs(np.sum(Z * Z)))
    resid = float(np.linalg.norm(H @ Z - lam * Z) / np.linalg.norm(Z))
    return Z, lam, r_d, dict(N=N, M=M, D=D, lam=lam, norm2=norm2, centroid_abs=cz,
                             square_closure_abs=s2, fixed_point_residual=resid)


def main(n_first, n_last):
    manifest = {}
    if os.path.exists(MANIFEST):
        manifest = json.load(open(MANIFEST))
    for N in range(n_first, n_last + 1):
        Z, lam, r_d, chk = build_floor(N)
        # 資格検査（指示書 §3）——不合格なら即停止・配置しない
        assert Z.dtype == np.complex128
        assert abs(chk['norm2'] - 1) < 1e-13, chk
        assert chk['centroid_abs'] < 1e-12, chk
        assert chk['square_closure_abs'] < 1e-12, chk
        if chk['fixed_point_residual'] >= 1e-12:
            print(f"STOP: N={N} residual={chk['fixed_point_residual']:.2e}（O(1)級またはe-12超）——走行禁止")
            sys.exit(2)
        assert abs(lam - LAM[N]) < 5e-9, (lam, LAM[N])
        assert np.max(np.abs(r_d - np.array(RD[N]))) < 5e-9, f'N={N}: r_d が導出mdと不一致'
        name = f'parent_static_N{N:05d}_makeparent_20260905.npz'
        path = os.path.join(DST, name)
        M = len(Z)
        np.savez_compressed(path, v=Z, g=np.zeros(M, dtype=np.complex128), Z0=Z,
                            sigma=np.array([abs(lam)]), residual=np.array(chk['fixed_point_residual']),
                            n=np.array(N), seed=np.array(-1), delta=np.array(0.0),
                            tol=np.array(0.0), iters=np.array(0),
                            family=np.zeros(M, dtype=np.int64), theta=np.angle(Z),
                            rule=np.array('theoretical_floor_cyclic_fourier'))
        h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
        manifest[name] = dict(N=N, M=M, **{k: chk[k] for k in
                              ('lam', 'fixed_point_residual', 'norm2', 'centroid_abs', 'square_closure_abs')},
                              sha256=h, construction='cyclic_fourier_standard_branch(§7)')
        print(f"N={N:2d} λ={lam:+.12f} 残差={chk['fixed_point_residual']:.2e} |Σz|={chk['centroid_abs']:.1e} "
              f"|Σz²|={chk['square_closure_abs']:.1e} ‖Z‖²={chk['norm2']:.15f} r_d一致OK SHA={h[:12]}…")
    with open(MANIFEST, 'w') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    print('manifest updated:', MANIFEST)


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
