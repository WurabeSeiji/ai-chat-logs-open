#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""99-run サーベイ一括分析 v1（収集完了後の最初の解析。全走行に同一パイプライン）。

per-run 指標:
  - Q2(0) の閉形式予測との照合（整合性検査を兼ねる数学的事実）
  - 生成子ノルムの初期/最小/最大（固定点近傍の検出）
  - U(1) 相対・絶対の回帰残差列 eps(q) と最小値・位置（q=1..4096）
  - SVD: full / W1=[0,50] / W3=[200,500] / W_late=[3796,4096] の r99, R_eff, E2
  - 初期2モード部分空間からの逸脱 eta_perp(t)（t=500, 4096, max）と
    late 窓主要2方向との principal angles
集計軸（事前登録 provenance のみ使用）:
  parity_class / gcd(ma,mb) / CRT独立性 gcd(ord_a,ord_b)=1 / L / den 対照 9 対
"""
import csv
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.abspath(os.path.join(HERE, '..'))
PHASEA = os.path.abspath(os.path.join(SURVEY, '..', 'PhaseA_実装一式_v1_20260915'))
sys.path.insert(0, PHASEA)
sys.path.insert(0, SURVEY)
from gen_manifest import build_runs  # noqa: E402

WINDOWS = {'full': (0, 4096), 'W1': (0, 50), 'W3': (200, 500), 'Wlate': (3796, 4096)}


def q2_theory(L, ma, mb, M):
    U = np.exp(2j * np.pi / L)
    def S(r):
        x = U ** r
        return M if abs(x - 1) < 1e-12 else L * x / (1 - x)
    ca2 = np.exp(-1j * np.pi / L)
    cb2 = np.exp(+1j * np.pi / L)
    return ca2 * S(2 * ma) + 2 * S(ma + mb) + cb2 * S(2 * mb)


def rank_indices(s):
    e = s ** 2
    tot = float(np.sum(e))
    cum = np.cumsum(e) / tot
    r99 = int(np.searchsorted(cum, 0.99) + 1)
    return dict(r99=r99, R_eff=float(tot ** 2 / np.sum(e ** 2)),
                E2=float(cum[1]) if len(cum) > 1 else 1.0)


def analyze_run(row):
    rid = row['run_id']
    L, ma, mb, den, M = row['L'], row['ma'], row['mb'], row['den'], row['M']
    Z = np.load(os.path.join(SURVEY, 'runs', rid, 'states.npz'))['Z']  # (4097, M)
    z0 = Z[0]
    H0 = float(np.vdot(z0, z0).real)

    # Q2 理論照合
    Q2m = complex(np.sum(z0 ** 2))
    Q2t = q2_theory(L, ma, mb, M)
    q2_dev = abs(Q2m - Q2t)

    # 生成子ノルム（invariants から。既知の書式バグ: run_survey.py が repr() で
    # 'np.float64(x)' 形式を書いたため剥がして読む。数値情報は無損失。DB は変更しない）
    def _f(s):
        s = s.strip()
        if s.startswith('np.float64('):
            s = s[11:-1]
        return float(s)
    with open(os.path.join(SURVEY, 'runs', rid, 'invariants.csv')) as f:
        next(f)
        gen_vals = [_f(line.rsplit(',', 1)[1]) for line in f]
    gen0, gen_min, gen_max = gen_vals[0], min(gen_vals), max(gen_vals)

    # 回帰残差列
    ip = Z @ z0.conj()                      # <z0, Z(t)> (t=0..4096)
    nrm2 = np.einsum('ij,ij->i', Z.real, Z.real) + np.einsum('ij,ij->i', Z.imag, Z.imag)
    eps_abs = np.sqrt(np.maximum(nrm2 + H0 - 2 * ip.real, 0.0)) / math.sqrt(H0)
    eps_u1 = np.sqrt(np.maximum(nrm2 + H0 - 2 * np.abs(ip), 0.0)) / math.sqrt(H0)
    qa_ = int(np.argmin(eps_abs[2:]) + 2)
    qu_ = int(np.argmin(eps_u1[2:]) + 2)

    # SVD 窓解析
    sv = {}
    for w, (t0, t1) in WINDOWS.items():
        s = np.linalg.svd(Z[t0:t1 + 1].T, compute_uv=False)
        sv[w] = rank_indices(s)
    # late 窓の主要2方向と初期2モード span の principal angles
    ea = np.arange(M)  # placeholder
    from kernel_canonical_dense import edges as _edges
    a_, b_ = _edges(L)
    delta = (b_ - a_).astype(np.float64)
    U = np.exp(2j * np.pi / L)
    B = np.stack([U ** (ma * delta), U ** (mb * delta)], axis=1)
    Qb, _ = np.linalg.qr(B)
    Ul, sl, _ = np.linalg.svd(Z[3796:4097].T, full_matrices=False)
    cosv = np.clip(np.linalg.svd(Qb.conj().T @ Ul[:, :2], compute_uv=False), 0, 1)
    pa = [float(np.arccos(c)) for c in cosv]
    # eta_perp
    G_inv_Bh = np.linalg.solve(B.conj().T @ B, B.conj().T)
    def eta(t):
        c = G_inv_Bh @ Z[t]
        return float(np.linalg.norm(Z[t] - B @ c) / np.linalg.norm(Z[t]))
    eta500, eta_end = eta(500), eta(4096)
    eta_max = max(eta(t) for t in range(0, 4097, 64))

    g = math.gcd(L, math.gcd(ma, mb))
    ord_a, ord_b = row['ord_a'], row['ord_b']
    return dict(run_id=rid, series=row['series'], L=L, M=M, ma=ma, mb=mb, den=den,
                parity=row['parity_class'], gcd_mamb=math.gcd(ma, mb),
                ord_a=ord_a, ord_b=ord_b, crt_indep=int(math.gcd(ord_a, ord_b) == 1),
                pair_order=L // g, q2_dev=q2_dev, absQ2=abs(Q2m), H0=H0,
                gen0=gen0, gen_min=gen_min, gen_max=gen_max,
                eps_abs_min=float(eps_abs[2:].min()), q_abs=qa_,
                eps_u1_min=float(eps_u1[2:].min()), q_u1=qu_,
                D_orbit=float(eps_abs[1:].max()),
                r99_full=sv['full']['r99'], r99_W1=sv['W1']['r99'],
                r99_W3=sv['W3']['r99'], r99_late=sv['Wlate']['r99'],
                Reff_late=sv['Wlate']['R_eff'], E2_late=sv['Wlate']['E2'],
                Reff_full=sv['full']['R_eff'],
                pa1_late=pa[0], pa2_late=pa[1],
                eta500=eta500, eta_end=eta_end, eta_max=eta_max)


def main():
    rows = build_runs()
    out = []
    for i, row in enumerate(rows):
        m = analyze_run(row)
        out.append(m)
        print(f"[{i+1:2d}/99] {m['run_id']:>22} q2dev={m['q2_dev']:.1e} "
              f"eps_u1_min={m['eps_u1_min']:.2e}@{m['q_u1']:>4} "
              f"r99(W1/W3/late)={m['r99_W1']}/{m['r99_W3']}/{m['r99_late']} "
              f"Reff_late={m['Reff_late']:.2f} eta_end={m['eta_end']:.2f}", flush=True)
    with open(os.path.join(HERE, 'survey_metrics.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print('ANALYSIS DONE:', len(out), 'runs')


if __name__ == '__main__':
    main()
