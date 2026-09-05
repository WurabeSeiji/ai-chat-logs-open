#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""回転数ロック全系列調査（N=3..40、den=N 系列、読み出しのみ・新規走行なし）。

入力（コピーせず相対パス参照、実行時に SHA256 入力ゲートで正本台帳と照合）:
  ../N3_N40_stage123_sweep_20260905/results/hm_N{N}_den_{N}_states_500.npz  （38本）
  ../N3_N40_stage123_sweep_20260905/results/summary_64bit_with124_N3_N40.csv （onset 引用）
  ../N3_N40_stage123_sweep_20260905/results/timeseries_64bit_with124_N3_N40.csv（H⊥/H 引用）
  ../N3_N40_stage123_sweep_20260905/SHA256SUMS.txt （入力ゲートの照合先）

エンジン: edges/adjacency/H_of/one_step は run_N3_N40_stage123_v1.py（正本）からの逐語コピー。
対照テスト: 全38系で Z[0] から10歩再ステップし、保存 Z[1..10] と bit 一致を要求
（不一致が1つでもあれば即 abort）。

測定（すべて機械判定）:
 (1) 剛体回転テスト r(s)=min_φ‖Z[s+N]−e^{iφ}Z[s]‖/‖Z[s]‖（φ=arg⟨Z[s],Z[s+N]⟩）。
     増幅期窓（H⊥/H < 1e-6 のフレーム）と飽和後窓（onset+50 以降）で評価。
 (2) 増幅期の回転数 x_early=φ(0)/2π（時計1回転あたり、単位=回転）と窓内ドリフト。
 (3) スペクトル同定: 初期生成子 H=i·K の固有対 (w_k, V_k)。N歩（=時計1回転）での
     モード位相は e^{-i2π w_k} なので予測回転数は (−w_k) mod 1。占有 c_k=|V_k†Z0|² の
     最大モードの予測と実測 x_early の円距離、占有モード間の予測回転数のばらつき。
 (4) 飽和後の回転数 x_late とその最良有理近似 p/q（q≤48、残差 <1e-9 回転でロック判定）、
     判定 q での厳密周期テスト ‖Z[s+qN]−Z[s]‖/‖Z[s]‖。
 (5) 終状態の等振幅偏差（‖Z‖/√M 基準の相対偏差）。
 (6) H⊥/H プレートー（onset+50..500 の平均・標準偏差）と最良有理近似（q≤200）。

出力: rotation_lock_table_v1.csv（N ごと1行）、rotation_lock_analysis_v1.json（詳細）、
      fig_rotation_lock_N3_N40.png（回転数・ロック・プレートーの3面図）。"""
import csv
import hashlib
import json
import math
import os
import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(BASE, '..', 'N3_N40_stage123_sweep_20260905')
RES = os.path.join(PKG, 'results')
CONTROL_STEPS = 10
EARLY_HPERP_MAX = 1e-6
SAT_MARGIN = 50
QMAX_ROT = 48
QMAX_PLATEAU = 200
LOCK_TOL_TURNS = 1e-9

# ---- エンジン（run_N3_N40_stage123_v1.py 14-28行の逐語コピー） ----
def edges(N):
    a,b=np.triu_indices(N,k=1); return a.astype(np.int64),b.astype(np.int64)
def adjacency(N):
    ea,eb=edges(N); M=len(ea); A=np.zeros((M,M),dtype=np.float64)
    for e in range(M):
        share=(ea==ea[e])|(ea==eb[e])|(eb==ea[e])|(eb==eb[e]); share[e]=False; A[e,share]=1.0
    return A
def H_of(z,A):
    H=A*(np.conj(z)[:,None]*z[None,:]); np.fill_diagonal(H,0.0); return H.astype(np.complex128,copy=False)
def one_step(z,A,den):
    # 段3の最小変更（唯一の力学変更点）: 位相のみ生成子 Ĥ の虚部だけを取る H=i·K（K=sin(Δθ) 実反対称）。
    # exp(-iΔτ·iK)=exp(Δτ·K) の実直交回転となり、Z^T Z（零閉塞）と ‖Z‖ を厳密保存する。
    H=H_of(np.exp(1j*np.angle(z)),A); H=(1j*np.imag(H)).astype(np.complex128,copy=False)
    w,V=np.linalg.eigh(H); phase=np.exp(-1j*np.float64(2.0*math.pi/den)*w)
    return (V@(phase*(V.conj().T@z))).astype(np.complex128,copy=False)
# ---- コピーここまで ----

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def circdist(a, b):
    """回転数（単位=回転）の円距離"""
    return abs((a - b + 0.5) % 1.0 - 0.5)

def best_rational(x, qmax):
    x = x % 1.0
    best = (None, None, 2.0)
    for q in range(1, qmax + 1):
        p = round(x * q)
        res = abs(x - p / q)
        if res < best[2] - 1e-18:
            best = (p % q if q > 1 else p, q, res)
    return best  # (p, q, residual_turns)

# ---- 入力ゲート: SHA256 照合 ----
ledger = {}
with open(os.path.join(PKG, 'SHA256SUMS.txt')) as f:
    for line in f:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]
for N in range(3, 41):
    rel = f'results/hm_N{N}_den_{N}_states_500.npz'
    actual = sha256(os.path.join(PKG, rel))
    assert ledger.get(rel) == actual, f'INPUT GATE FAIL: {rel}'
print('INPUT GATE: 38/38 npz SHA256 match ledger', flush=True)

# ---- 正本 CSV 引用: onset と H⊥/H 時系列 ----
onset = {}
with open(os.path.join(RES, 'summary_64bit_with124_N3_N40.csv')) as f:
    for row in csv.DictReader(f):
        if row['series'] == 'N':
            onset[int(row['N'])] = int(row['onset_gt_0.05'])
hperp = {N: np.full(501, np.nan) for N in range(3, 41)}
with open(os.path.join(RES, 'timeseries_64bit_with124_N3_N40.csv')) as f:
    for row in csv.DictReader(f):
        if row['series'] == 'N':
            hperp[int(row['N'])][int(row['step'])] = float(row['Hperp_frac'])

# ---- 本体 ----
rows = []
detail = {}
for N in range(3, 41):
    d = np.load(os.path.join(RES, f'hm_N{N}_den_{N}_states_500.npz'))
    assert int(d['denominator']) == N and int(d['steps']) == 500
    Z = np.asarray(d['Z'], dtype=np.complex128)
    M = N * (N - 1) // 2
    assert Z.shape == (501, M)
    A = adjacency(N)

    # 対照テスト: 10歩 bit 一致
    z = Z[0].copy()
    for t in range(1, CONTROL_STEPS + 1):
        z = one_step(z, A, N)
        assert np.array_equal(z, Z[t]), f'CONTROL FAIL: N={N} step {t}'

    # 剛体回転テスト（全フレーム）
    frames = list(range(0, 501 - N, N))
    phi = np.empty(len(frames)); rres = np.empty(len(frames))
    for i, s in enumerate(frames):
        ip = np.vdot(Z[s], Z[s + N])
        ph = float(np.angle(ip)); phi[i] = ph
        nrm = float(np.linalg.norm(Z[s]))
        rres[i] = float(np.linalg.norm(Z[s + N] - np.exp(1j * ph) * Z[s])) / nrm

    early_idx = [i for i, s in enumerate(frames) if hperp[N][s + N] < EARLY_HPERP_MAX]
    x_early = (phi[0] / (2 * math.pi)) % 1.0
    r_early_max = float(max(rres[i] for i in early_idx)) if early_idx else float('nan')
    drift_early = float(max(circdist((phi[i] / (2 * math.pi)) % 1.0, x_early) for i in early_idx)) if early_idx else float('nan')

    # スペクトル同定（初期生成子）
    H = H_of(np.exp(1j * np.angle(Z[0])), A); H = (1j * np.imag(H)).astype(np.complex128, copy=False)
    w, V = np.linalg.eigh(H)
    occ = np.abs(V.conj().T @ Z[0]) ** 2
    occ = occ / occ.sum()
    occ_idx = np.flatnonzero(occ > 1e-12)
    kdom = int(occ_idx[np.argmax(occ[occ_idx])])
    pred = [(-w[k]) % 1.0 for k in occ_idx]
    d_dom = circdist(x_early, (-w[kdom]) % 1.0)
    spread_occ = float(max(circdist(a, pred[0]) for a in pred)) if len(pred) > 1 else 0.0

    # 飽和後窓
    sc = onset[N]
    late_idx = [i for i, s in enumerate(frames) if sc >= 0 and s >= sc + SAT_MARGIN]
    if len(late_idx) >= 3:
        r_late_max = float(max(rres[i] for i in late_idx))
        x_late = (phi[late_idx[-1]] / (2 * math.pi)) % 1.0
        spread_late = float(max(circdist((phi[i] / (2 * math.pi)) % 1.0, x_late) for i in late_idx))
        lp, lq, lres = best_rational(x_late, QMAX_ROT)
        locked = bool(lres < LOCK_TOL_TURNS)
        # 厳密周期テスト（判定 q）
        s0 = None
        for i in reversed(late_idx):
            if frames[i] + lq * N <= 500:
                s0 = frames[i]; break
        r_period = (float(np.linalg.norm(Z[s0 + lq * N] - Z[s0]) / np.linalg.norm(Z[s0]))
                    if s0 is not None else float('nan'))
        # プレートー
        plat = hperp[N][sc + SAT_MARGIN:501]
        plat_mean = float(np.nanmean(plat)); plat_std = float(np.nanstd(plat))
        pp, pq, pres = best_rational(plat_mean, QMAX_PLATEAU)
    else:
        r_late_max = x_late = spread_late = lres = r_period = plat_mean = plat_std = pres = float('nan')
        lp = lq = pp = pq = None; locked = False

    # 終状態の等振幅偏差
    amps = np.abs(Z[500]); target = float(np.linalg.norm(Z[500])) / math.sqrt(M)
    eq_dev = float(np.max(np.abs(amps - target)) / target)

    rows.append([N, M, sc, len(early_idx), r_early_max, x_early, drift_early,
                 len(occ_idx), float(w[kdom]), float(occ[kdom]), float(d_dom), spread_occ,
                 len(late_idx), r_late_max, x_late, spread_late,
                 lp, lq, lres, locked, r_period, eq_dev, plat_mean, plat_std, pp, pq, pres])
    detail[N] = {
        'occupied_modes': [{'w': float(w[k]), 'occ': float(occ[k]), 'pred_rot_turns': float((-w[k]) % 1.0)}
                           for k in occ_idx[np.argsort(-occ[occ_idx])][:12]],
        'x_early_turns': float(x_early), 'x_late_turns': (None if isinstance(x_late, float) and math.isnan(x_late) else float(x_late)),
    }
    print(f'done N={N}', flush=True)

HEADER = ['N', 'M', 'onset', 'n_early_frames', 'r_early_max', 'x_early_turns', 'drift_early_turns',
          'n_occupied', 'w_dominant', 'occ_dominant', 'd_dom_turns', 'spread_occ_turns',
          'n_late_frames', 'r_late_max', 'x_late_turns', 'spread_late_turns',
          'lock_p', 'lock_q', 'lock_residual_turns', 'locked_1e-9', 'period_residual',
          'eqmod_rel_dev', 'plateau_mean', 'plateau_std', 'plat_p', 'plat_q', 'plat_residual']
with open(os.path.join(BASE, 'rotation_lock_table_v1.csv'), 'w', newline='') as f:
    wr = csv.writer(f); wr.writerow(HEADER); wr.writerows(rows)
with open(os.path.join(BASE, 'rotation_lock_analysis_v1.json'), 'w') as f:
    json.dump({'params': {'control_steps': CONTROL_STEPS, 'early_hperp_max': EARLY_HPERP_MAX,
                          'sat_margin': SAT_MARGIN, 'qmax_rot': QMAX_ROT, 'qmax_plateau': QMAX_PLATEAU,
                          'lock_tol_turns': LOCK_TOL_TURNS},
               'table_header': HEADER,
               'table': rows, 'detail': detail}, f, indent=2, default=str)

# ---- 図: 3面 ----
tab = {r[0]: r for r in rows}
Ns = list(range(3, 41))
odd = [N for N in Ns if N % 2 == 1]; even = [N for N in Ns if N % 2 == 0]
fig, axs = plt.subplots(3, 1, figsize=(12, 14))
ax = axs[0]
for grp, c, lab in ((odd, 'tab:red', 'odd N'), (even, 'tab:blue', 'even N')):
    ax.plot(grp, [tab[N][5] for N in grp], 'o', color=c, label=lab)
ax.set_ylabel('x_early (turns/clock turn)'); ax.set_title('Ramp rotation number (rigid-fit phase at s=0)')
ax.grid(alpha=.3); ax.legend()
ax = axs[1]
for grp, c, lab in ((odd, 'tab:red', 'odd N'), (even, 'tab:blue', 'even N')):
    xs = [N for N in grp if not (isinstance(tab[N][14], float) and math.isnan(tab[N][14]))]
    ax.plot(xs, [tab[N][14] for N in xs], 'o', color=c, label=lab)
    for N in xs:
        if tab[N][19]:
            ax.annotate(f'{tab[N][16]}/{tab[N][17]}', (N, tab[N][14]), textcoords='offset points',
                        xytext=(0, 6), fontsize=7, ha='center')
ax.set_ylabel('x_late (turns/clock turn)'); ax.set_title('Saturated rotation number (locked rationals annotated)')
ax.grid(alpha=.3); ax.legend()
ax = axs[2]
for grp, c, lab in ((odd, 'tab:red', 'odd N'), (even, 'tab:blue', 'even N')):
    xs = [N for N in grp if not (isinstance(tab[N][22], float) and math.isnan(tab[N][22]))]
    ax.plot(xs, [tab[N][22] for N in xs], 'o', color=c, label=lab)
ax.set_ylabel('Hperp/H plateau mean'); ax.set_xlabel('N'); ax.set_title('Saturation plateau')
ax.grid(alpha=.3); ax.legend()
fig.suptitle('Rotation-number locks across N=3..40 (den=N series, readout only)', y=.995)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_rotation_lock_N3_N40.png'), dpi=180)
plt.close(fig)
print('ALL DONE')
