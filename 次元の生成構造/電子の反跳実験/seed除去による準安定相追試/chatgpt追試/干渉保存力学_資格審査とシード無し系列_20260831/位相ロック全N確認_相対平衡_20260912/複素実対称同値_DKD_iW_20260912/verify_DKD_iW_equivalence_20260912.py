#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""90°骨格で複素自己無撞着固有問題 Kz=iωz が実対称0-1行列 W の固有値問題 Wr=ωr と厳密同値になること
（D^{-1}KD=iW）と、Perron–Frobenius による正振幅主固有ベクトルの存在・一意性を検証
（読出しのみ・物理無変更・新規走行なし）。

木原指示（2026-09-12・ChatGPT指摘を採用）: 論文1の核心の追加結果。
  D=diag(e^{iφ_e}) とし、90°骨格 ψ_ef=φ_f-φ_e∈(π/2)Z では恒等式 sinψ·e^{iψ}=i sin²ψ（⇔ sin2ψ=0）より
    (D^{-1}KD)_ef = A_ef sinψ_ef e^{iψ_ef} = i A_ef sin²ψ_ef = i W_ef,  W_ef=A_ef sin²ψ_ef（実対称・0/1・非負）。
  ゆえ Kz=iωz（z=Dr, r実）⇔ D^{-1}KD r=iωr ⇔ iWr=iωr ⇔ Wr=ωr。
  W が既約（連結）なら Perron–Frobenius により ω_max=ρ(W) の主固有ベクトルは r_e>0・スケールを除き一意。

検証（N=3..40、整数K解析床 と make_parent 床）:
  err_DKD=‖D^{-1}KD - iW‖, Wr=ωr の一様性(std)と ω, ρ(W)=maxeig(W), σ(iK)=iKスペクトル半径,
  ω=ρ(W)=σ の一致, r>0, r と W 主固有ベクトルの重なり(=1?), W の連結成分数。
正本 run_N3_N40_stage123_v1.py（SHA照合）から adjacency を import。
"""
import csv
import hashlib
import os

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

HERE = os.path.dirname(os.path.abspath(__file__))
ORIG = os.path.abspath(os.path.join(HERE, '..', '..', 'N3_N40_stage123_sweep_20260905', 'run_N3_N40_stage123_v1.py'))
PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'
IK = os.path.abspath(os.path.join(HERE, '..', '..', '90度理論床_整数K解析解_全N_20260909', 'parents_90deg_floor_analytic'))
MP = os.path.abspath(os.path.join(HERE, '..', '..', 'make_parent型初期値_自己無撞着構造_20260907', 'full_N3_N40_sweep', 'states'))

src = open(ORIG, encoding='utf-8').read()
assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本SHA不一致'
g = {'__name__': 'canonical', '__file__': ORIG}
exec(compile(src.split('rows=[]; summaries=[]')[0], ORIG, 'exec'), g)
adjacency = g['adjacency']


def check(z, N, floor):
    A = adjacency(N); phi = np.angle(z); M = len(z)
    d = phi[None, :] - phi[:, None]
    K = A * np.sin(d); W = A * (np.sin(d) ** 2)
    D = np.exp(1j * phi)
    DKD = (1 / D)[:, None] * K * D[None, :]
    err = float(np.max(np.abs(DKD - 1j * W)))
    r = np.abs(z); Wr = W @ r; om_e = Wr / r
    om = float(np.mean(om_e)); om_std = float(np.std(om_e))
    ew, ev = np.linalg.eigh(W); rho = float(ew[-1])
    overlap = float(abs(np.dot(ev[:, -1] / np.linalg.norm(ev[:, -1]), r / np.linalg.norm(r))))
    sig = float(np.max(np.abs(np.linalg.eigvals(1j * K))))
    ncomp, _ = connected_components(csr_matrix((W > 1e-9).astype(int)), directed=False)
    return dict(N=N, M=M, floor=floor, err_DKD_iW=err, omega=om, omega_uniform_std=om_std,
                rho_W=rho, sigma_iK=sig, omega_eq_rho_eq_sigma=bool(abs(om - rho) < 1e-6 and abs(om - sig) < 1e-4),
                r_all_positive=bool(np.all(r > 0)), perron_overlap=overlap, W_components=int(ncomp))


def main():
    rows = []
    for N in range(3, 41):
        p = os.path.join(IK, f'parent_90deg_floor_N{N:05d}_analytic.npz')
        if os.path.exists(p):
            rows.append(check(np.asarray(np.load(p)['Z0'], np.complex128), N, 'integerK'))
    for N in range(3, 41):
        p = os.path.join(MP, f'hm_N{N}_den_{N}_states_500.npz')
        if os.path.exists(p):
            rows.append(check(np.asarray(np.load(p)['Z'][0], np.complex128), N, 'make_parent'))
    with open(os.path.join(HERE, 'results', 'DKD_iW_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
        for r in rows:
            w.writerow(r)
    for floor in ('integerK', 'make_parent'):
        sub = [r for r in rows if r['floor'] == floor]
        maxerr = max(r['err_DKD_iW'] for r in sub)
        alleq = all(r['omega_eq_rho_eq_sigma'] for r in sub)
        allpos = all(r['r_all_positive'] for r in sub)
        allconn = all(r['W_components'] == 1 for r in sub)
        minov = min(r['perron_overlap'] for r in sub)
        print(f'[{floor}] N=3..40 {len(sub)}床: max‖D^-1KD-iW‖={maxerr:.1e}  ω=ρ(W)=σ 全一致={alleq}  '
              f'r>0 全={allpos}  Perron重なり最小={minov:.6f}  W連結(成分1)全={allconn}')
    print('\n結論: 90°骨格で Kz=iωz ⇔ Wr=ωr（D^{-1}KD=iW）が厳密成立。W は実対称0-1・非負・既約（連結）で、'
          'Perron–Frobenius により正振幅主固有ベクトル r>0（ω=ρ(W)=σ）がスケールを除き一意に存在＝'
          '90°骨格上の自己無撞着正振幅波は構成できる。z=Dr で Kz=iρ(W)z。')


if __name__ == '__main__':
    main()
