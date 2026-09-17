#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""実行前構造診断（v3 指示書 §4 の8項目を生値で記録する。判定はしない）。

方針（木原・ChatGPT 指示 2026-09-17）:
  監査するのはコード。判定するのは実験後のデータ。
  K 反対称残差・H エルミート残差・H0/Q2 の1step変化・U(1)残差・NaN/Inf 等は
  すべて「停止条件ではなく観測値」。本スクリプトは診断値を
  pre_run_structural_audit.json に生値のまま保存するだけで、
  PASS/FAIL 判定・閾値・実験開始ゲートを持たない（常に exit 0）。

記録項目（ランダム状態と実初期値の双方）:
  1. max(abs(A - A.T))
  2. max(abs(K + K.T))
  3. max(abs(H - H.conj().T))
  4. 大域 U(1) 変換前後の K の最大差
  5. z_e=0 を含む人工状態での該当相互作用成分の値
  6. 一 step 前後の H0 = Z^dagger Z
  7. 一 step 前後の Q2 = Z^T Z
  8. NaN/Inf の有無

K の補正・対称化・正規化は行わない（観測のみ）。

改版記録: v1 は数値閾値による PASS/FAIL と exit 1 ゲートを持っていたが、
2026-09-17 の指示で撤去（閾値緩和ではなくゲート自体を廃止）。
初回実行時の停止原因は K+K^T ~1e-16（numpy SIMD/FMA の丸め非対称）に対する
「厳密0」期待の過剰さであり、コード・物理式の異常ではなかった。

使い方: python3 structural_audit_v1.py <run_id> [<run_id> ...]
  例:   python3 structural_audit_v1.py L8_ma1_mb2_den8 L8_ma1_mb3_den8 L8_ma2_mb4_den8
"""
import os

THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
for _v in THREAD_VARS:
    os.environ[_v] = '1'

import json                                                    # noqa: E402
import platform                                                # noqa: E402
import sys                                                     # noqa: E402
import time                                                    # noqa: E402

import numpy as np                                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from interaction_kernel_theory_v1 import adjacency, K_of, one_step  # noqa: E402
from generator_twowave_v1 import generate                           # noqa: E402
from gen_manifest import build_runs                                 # noqa: E402

AUDIT_SEED = 20260917          # ランダム状態の再現用（事前固定）
U1_PHI = 0.7343218946          # 大域 U(1) 検査の位相（事前固定・非有理角）


def audit_state(tag, z, A, den):
    """一状態に対する §4 の全項目を生値で記録する（判定しない・補正しない）。"""
    K = K_of(z, A)
    H = (1j * K).astype(np.complex128, copy=False)
    scaleK = float(np.max(np.abs(K))) or 1.0

    a_asym = float(np.max(np.abs(A - A.T)))
    k_anti = float(np.max(np.abs(K + K.T)))
    h_herm = float(np.max(np.abs(H - H.conj().T)))

    K2 = K_of(np.exp(1j * U1_PHI) * z, A)
    u1_diff = float(np.max(np.abs(K2 - K)))

    zz = z.copy()
    e0 = 0
    zz[e0] = 0.0
    Kz = K_of(zz, A)
    zero_comp = float(max(np.max(np.abs(Kz[e0, :])), np.max(np.abs(Kz[:, e0]))))

    H0_before = float(np.vdot(z, z).real)
    Q2_before = complex(np.sum(z ** 2))
    z1 = one_step(z, A, den)
    H0_after = float(np.vdot(z1, z1).real)
    Q2_after = complex(np.sum(z1 ** 2))
    h0_drift = abs(H0_after - H0_before) / H0_before
    q2_drift = abs(Q2_after - Q2_before) / max(H0_before, abs(Q2_before))

    finite = bool(np.all(np.isfinite(z.view(np.float64)))
                  and np.all(np.isfinite(K))
                  and np.all(np.isfinite(z1.view(np.float64))))

    return {
        'tag': tag, 'M': int(len(z)), 'den': den,
        'A_asym_max': a_asym,
        'K_antisym_max': k_anti,
        'K_scale_max_abs': scaleK,
        'H_nonherm_max': h_herm,
        'U1_K_maxdiff': u1_diff,
        'U1_K_maxdiff_rel': u1_diff / scaleK,
        'zero_comp_K_max': zero_comp,
        'H0_before': H0_before, 'H0_after': H0_after, 'H0_drift_rel': h0_drift,
        'Q2_before_re': Q2_before.real, 'Q2_before_im': Q2_before.imag,
        'Q2_after_re': Q2_after.real, 'Q2_after_im': Q2_after.imag,
        'Q2_drift_rel': q2_drift,
        'finite': finite,
    }


def main():
    run_ids = sys.argv[1:]
    if not run_ids:
        print('usage: python3 structural_audit_v1.py <run_id> [<run_id> ...]')
        sys.exit(2)
    rows = {r['run_id']: r for r in build_runs()}
    for rid in run_ids:
        if rid not in rows:
            print(f'unknown run_id (manifest 99 に無い): {rid}')
            sys.exit(2)

    results = []
    rng = np.random.default_rng(AUDIT_SEED)
    done_L = set()
    for rid in run_ids:
        row = rows[rid]
        L, ma, mb, den = row['L'], row['ma'], row['mb'], row['den']
        A = adjacency(L)
        z0 = generate(L, ma, mb)
        results.append(audit_state(f'initial_{rid}', z0, A, den))
        if L not in done_L:
            done_L.add(L)
            M = L * (L - 1) // 2
            zr = (rng.standard_normal(M) + 1j * rng.standard_normal(M)).astype(np.complex128)
            results.append(audit_state(f'random_L{L}_seed{AUDIT_SEED}', zr, A, den))

    report = {
        'timestamp_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'mode': 'diagnostic-only（判定・ゲートなし。数値はすべて観測値）',
        'audit_seed': AUDIT_SEED, 'u1_phi': U1_PHI,
        'python': platform.python_version(), 'numpy': np.__version__,
        'run_ids': run_ids,
        'results': results,
    }
    out = os.path.join(HERE, 'pre_run_structural_audit.json')
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=1, ensure_ascii=False)
    for r in results:
        print(f"{r['tag']:>34}: K_anti={r['K_antisym_max']:.3e} "
              f"U1={r['U1_K_maxdiff_rel']:.3e} dH0={r['H0_drift_rel']:.3e} "
              f"dQ2={r['Q2_drift_rel']:.3e} finite={r['finite']}")
    print(f'structural diagnostics recorded ({len(results)} states) -> {out}')


if __name__ == '__main__':
    main()
