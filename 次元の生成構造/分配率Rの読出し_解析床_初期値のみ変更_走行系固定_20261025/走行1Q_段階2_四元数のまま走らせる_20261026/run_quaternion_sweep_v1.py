#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""走行1Q 段階2（2026-10-26）：四元数のまま走らせる走行系。走行1 の run_sweep_from_parents_v1.py（SHA d7749ae2…）を土台に、変更は one_step と読出しだけ。

状態：Q（M×4 実数、列＝1, i, j, k 成分）。A：走行1 と同じ規則で ea, eb から一度だけ作る（関数 adjacency は走行1 と同一）。
一歩（選択として固定したもの：K は一歩の初めの状態から一度だけ読む、回転は i → j → k の順）：
  u_e = q_e/|q_e|、 R_ef = ū_e u_f（四元数の積）
  K^(c)_ef = A_ef · R_ef の c 成分（c ＝ i, j, k）。実反対称（ū_f u_e ＝ conj(ū_e u_f)）、対角 0。
  Q ← exp(Δτ K^(k)) exp(Δτ K^(j)) exp(Δτ K^(i)) Q 、各 exp は 4 列すべてに同じ実直交行列として掛ける。
  exp の計算は走行1 の one_step と同じ式（H＝iK を eigh、V·diag(e^(−iΔτw))·V^H）を、複素列 a＋ib と c＋id に掛ける。
  K^(c) がすべて 0 の回転は飛ばす（恒等写像。状態が 1–i 面にあれば K^(j)＝K^(k)＝0 で、走行1 と同じ一歩になる）。
状態が 1–i 面（c＝d＝0）にあるとき K^(i)_ef ＝ A_ef sin(θ_f − θ_e) ＝ 走行1 の K。
読出し：Π＝初期状態の 4 列が張る部分空間（SVD、相対 1e-10 で階数を切る）。
  Hperp_frac ＝ ‖Q − P_Π Q‖²/‖Q‖²（走行1 の H⊥/H の拡張。control では Π は Re Z0, Im Z0 の張る面で走行1 と同じ）
  f_1, f_i, f_j, f_k ＝ 各列の Π の外の部分 /‖Q‖²（和が Hperp_frac）、H_total ＝ ‖Q‖²、global_closure ＝ |Σ q_e²|/‖Q‖²
使い方: python run_quaternion_sweep_v1.py STATE_DIR PATTERN OUT TITLE N1 N2 ..."""
import os, math, csv, json, platform, time, sys
import numpy as np

from make_initial_quaternion_v1 import qmul

STATE_DIR = sys.argv[1]; PATTERN = sys.argv[2]; OUT = sys.argv[3]; TITLE = sys.argv[4]
os.makedirs(OUT, exist_ok=True)
STEPS = 500
NS = tuple(int(x) for x in sys.argv[5:])


def adjacency(ea, eb):   # 走行1 と同一
    M = len(ea); A = np.zeros((M, M), dtype=np.float64)
    for e in range(M):
        share = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e]); share[e] = False; A[e, share] = 1.0
    return A


def rotate(K, zs, den):   # 走行1 の one_step と同じ式で exp(Δτ K) を掛ける
    H = (1j * K).astype(np.complex128)
    w, V = np.linalg.eigh(H); phase = np.exp(-1j * np.float64(2.0 * math.pi / den) * w)
    return [(V @ (phase * (V.conj().T @ z))).astype(np.complex128, copy=False) for z in zs]


def one_step_q(Q, A, den):
    U = Q / np.linalg.norm(Q, axis=1)[:, None]
    Ubar = U * np.array([1.0, -1.0, -1.0, -1.0])
    R = qmul(Ubar[:, None, :], U[None, :, :])          # R[e,f] = ū_e u_f
    Ks = []
    for c in (1, 2, 3):
        K = A * R[:, :, c]; np.fill_diagonal(K, 0.0); Ks.append(K)
    z1 = Q[:, 0] + 1j * Q[:, 1]; z2 = Q[:, 2] + 1j * Q[:, 3]
    for K in Ks:                                         # i → j → k
        if np.any(K != 0.0):
            z1, z2 = rotate(K, (z1, z2), den)
    return np.stack([z1.real, z1.imag, z2.real, z2.imag], axis=1)


def basis(Q0):
    U, s, _ = np.linalg.svd(Q0, full_matrices=False)
    return U[:, s > 1e-10 * s[0]]


def metrics(Q, B):
    h = float(np.sum(Q * Q)); P = Q - B @ (B.T @ Q)
    fc = np.sum(P * P, axis=0) / h
    s = np.zeros(4)
    for q in Q:
        s += qmul(q, q)
    return float(fc.sum()), h, float(np.linalg.norm(s) / h), fc


rows = []; summaries = []
for N in NS:
    d = np.load(os.path.join(STATE_DIR, PATTERN.format(N=N)))
    Q0 = np.array(d['Q0'], dtype=np.float64, copy=True); ea, eb = d['ea'], d['eb']; assert len(ea) == Q0.shape[0]
    A = adjacency(ea, eb); B = basis(Q0)
    for den, label in ((N, 'N'), (124, '124')):
        t0 = time.time(); Q = Q0.copy(); vals = np.empty(STEPS + 1); states = np.empty((STEPS + 1,) + Q.shape)
        for t in range(STEPS + 1):
            states[t] = Q; f, h, cl, fc = metrics(Q, B); vals[t] = f
            rows.append((N, label, den, t, f, h, cl, *fc))
            if t < STEPS: Q = one_step_q(Q, A, den)
        np.savez_compressed(os.path.join(OUT, f'hm_N{N}_den_{den}_states_500.npz'), Q=states, N=np.int64(N), denominator=np.int64(den), steps=np.int64(STEPS), rank_Pi=np.int64(B.shape[1]))
        ix = np.flatnonzero(vals > 0.05)
        summaries.append((N, label, den, B.shape[1], int(ix[0]) if ix.size else -1, float(vals[0]), float(vals[1]), float(vals[-1]), float(vals.max())))
        print(f'  N={N} den={den}: rank Π={B.shape[1]} {time.time() - t0:.1f} s', flush=True)
    print('done N', N, flush=True)
with open(os.path.join(OUT, 'timeseries_64bit_with124_N3_N40.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['N', 'series', 'denominator', 'step', 'Hperp_frac', 'H_total', 'global_closure', 'f_1', 'f_i', 'f_j', 'f_k']); w.writerows(rows)
with open(os.path.join(OUT, 'summary_64bit_with124_N3_N40.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['N', 'series', 'denominator', 'rank_Pi', 'onset_gt_0.05', 'initial', 'step1', 'final', 'max']); w.writerows(summaries)
with open(os.path.join(OUT, 'RUN_METADATA.json'), 'w') as f:
    json.dump({'steps': STEPS, 'N_list': list(NS), 'denominators': 'N,124', 'state_dir': STATE_DIR, 'title': TITLE, 'numpy': np.__version__, 'python': platform.python_version(),
               'one_step': 'K^(c)=A*(conj(u_e)u_f)_c, c=i,j,k read once at step start; Q <- exp(dt K^k) exp(dt K^j) exp(dt K^i) Q, same real orthogonal on all 4 columns; zero K skipped'}, f, indent=2)
print('ALL DONE')
