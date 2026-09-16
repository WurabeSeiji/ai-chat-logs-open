#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""二波広域サーベイ実行（収集仕様書 v1・APPROVED）。

順序（仕様 §0 事前登録手順）: 本コードと manifest はコミット済みであること。
  1. 互換ゲート（L=12/(2,3)/den=12 と L=6/(2,3)/den=6、t=0..500、scaled gate）
     不合格なら 99 走行を開始せず停止（exit 1）。
  2. 局所積分器監査（代表 3 走行: step1 dense/sparse、step10 shadow、strict one-step 差）
  3. 99 runs 実行（T=4096 全状態保存、毎 step H/Q2/gen_norm のみ）＋ QA
収集中の科学解析は行わない。QA FAIL でもデータは保存し qa_status で記録する。
"""
import csv
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time

import numpy as np
import scipy

HERE = os.path.dirname(os.path.abspath(__file__))
PHASEA = os.path.abspath(os.path.join(HERE, '..', 'PhaseA_実装一式_v1_20260915'))
FROZEN = os.path.abspath(os.path.join(HERE, '..', '動的モード再編成_低ランク凝縮_20260915',
                                      'results', 'states'))
sys.path.insert(0, PHASEA)
sys.path.insert(0, HERE)

from kernel_canonical_dense import adjacency, one_step        # noqa: E402
from kernel_sparse_matrixfree import SparseKernel             # noqa: E402
from generator_phaseA import generate                         # noqa: E402
from gen_manifest import build_runs, T                        # noqa: E402

RUNS_DIR = os.path.join(HERE, 'runs')
os.makedirs(RUNS_DIR, exist_ok=True)

GATE = 1e-12
QA_TOL = 1e-11
LANCZOS_TOL = 1e-14
AUDIT_REPS = [(12, 2, 3, 12), (24, 1, 3, 24), (40, 2, 4, 40)]

COMMIT_SHA = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
                            text=True, cwd=HERE).stdout.strip()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def compat_gate():
    report = []
    ok_all = True
    for L, ma, mb, den, frozen_name in ((12, 2, 3, 12, 'L12_den12_AB_states.npz'),
                                        (6, 2, 3, 6, 'L6_den6_AB_states.npz')):
        ref = np.load(os.path.join(FROZEN, frozen_name))['Z']
        sk = SparseKernel(L)
        z = generate(L, ma, mb)
        worst = 0.0
        for t in range(501):
            d = np.abs(z - ref[t])
            s = np.maximum(1.0, np.abs(ref[t]))
            worst = max(worst, float(np.max(d / s)))
            if t < 500:
                z = sk.one_step(z, den, tol=LANCZOS_TOL)
        ok = worst <= GATE
        ok_all &= ok
        report.append({'gate': f'L{L}_({ma},{mb})_den{den}', 'worst_scaled_diff': worst,
                       'pass': bool(ok)})
        print(f'GATE L{L}: worst={worst:.3e} {"PASS" if ok else "FAIL"}', flush=True)
    return ok_all, report


def strict_one_step(sk, z, den):
    """厳格参照: tol=0（早期終了は厳密0のみ）。breakdown guard は凍結カーネル内蔵
    （Krylov ノルム機械精度消失時に終了 = genuine/happy breakdown、記録して許容）。"""
    return sk.one_step(z, den, tol=0.0)


def audit_local():
    out = []
    for L, ma, mb, den in AUDIT_REPS:
        A = adjacency(L)
        sk = SparseKernel(L)
        z0 = generate(L, ma, mb)
        n0 = np.linalg.norm(z0)
        d1 = one_step(z0, A, den)
        s1 = sk.one_step(z0, den, tol=LANCZOS_TOL)
        strict1 = strict_one_step(sk, z0, den)
        zd, zs = z0.copy(), z0.copy()
        shadow = []
        for t in range(10):
            zd = one_step(zd, A, den)
            zs = sk.one_step(zs, den, tol=LANCZOS_TOL)
            shadow.append(float(np.linalg.norm(zs - zd) / n0))
        out.append({'rep': f'L{L}_({ma},{mb})_den{den}',
                    'step1_dense_vs_sparse': float(np.linalg.norm(s1 - d1) / n0),
                    'step1_sparse_vs_strict': float(np.linalg.norm(s1 - strict1) / n0),
                    'shadow10_dense_vs_sparse': shadow})
        print(f'AUDIT {out[-1]["rep"]}: step1 dvs={out[-1]["step1_dense_vs_sparse"]:.3e} '
              f'strict={out[-1]["step1_sparse_vs_strict"]:.3e} '
              f'shadow10={shadow[-1]:.3e}', flush=True)
    return out


def do_run(row):
    run_id = row['run_id']
    rdir = os.path.join(RUNS_DIR, run_id)
    os.makedirs(rdir, exist_ok=True)
    L, ma, mb, den = row['L'], row['ma'], row['mb'], row['den']
    sk = SparseKernel(L)
    z0 = generate(L, ma, mb)
    M = len(z0)
    t_start = time.time()
    states = np.empty((T + 1, M), np.complex128)
    inv = np.empty((T + 1, 4), np.float64)
    z = z0.copy()
    for t in range(T + 1):
        states[t] = z
        u = np.exp(1j * np.angle(z))
        inv[t] = [float(np.vdot(z, z).real),
                  float(np.sum(z ** 2).real), float(np.sum(z ** 2).imag),
                  float(np.linalg.norm(sk.K_apply(u, z)) / np.linalg.norm(z))]
        if t < T:
            z = sk.one_step(z, den, tol=LANCZOS_TOL)
    wall = time.time() - t_start

    np.savez(os.path.join(rdir, 'states.npz'), Z=states,
             L=np.int64(L), ma=np.int64(ma), mb=np.int64(mb),
             den=np.int64(den), T=np.int64(T))
    with open(os.path.join(rdir, 'invariants.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['t', 'H', 'Q2_re', 'Q2_im', 'generator_norm'])
        w.writerows([[t] + [repr(x) for x in inv[t]] for t in range(T + 1)])

    # QA（科学分析ではない）
    H0 = inv[0, 0]
    Q20 = complex(inv[0, 1], inv[0, 2])
    dH = float(np.max(np.abs(inv[:, 0] - H0) / H0))
    Q2t = inv[:, 1] + 1j * inv[:, 2]
    dQ2 = float(np.max(np.abs(Q2t - Q20)) / max(H0, abs(Q20)))
    finite = bool(np.all(np.isfinite(states.view(np.float64))))
    reasons = []
    if not finite:
        reasons.append('NaN/Inf in states')
    if dH > QA_TOL:
        reasons.append(f'H drift {dH:.3e} > {QA_TOL}')
    if dQ2 > QA_TOL:
        reasons.append(f'Q2 drift {dQ2:.3e} > {QA_TOL}')
    if states.shape[0] != T + 1:
        reasons.append('missing states')
    qa_status = 'PASS' if not reasons else 'FAIL'

    meta = dict(row)
    meta.update(dict(dtau=2.0 * math.pi / den, lanczos_tol=LANCZOS_TOL,
                     initial_condition='z_jk(0)=exp(-i*pi/(2L))*U^(ma*(k-j))'
                                       '+exp(+i*pi/(2L))*U^(mb*(k-j)), U=exp(2*pi*i/L), j<k',
                     code_commit_sha=COMMIT_SHA,
                     python=platform.python_version(), numpy=np.__version__,
                     scipy=scipy.__version__, wall_time_sec=wall,
                     qa_status=qa_status, qa_reason='; '.join(reasons),
                     dH_rel_max=dH, dQ2_rel_max=dQ2))
    with open(os.path.join(rdir, 'metadata.json'), 'w', encoding='utf-8') as f:
        json.dump(meta, f, indent=1, ensure_ascii=False)
    return dict(run_id=run_id, qa_status=qa_status, qa_reason='; '.join(reasons),
                dH_rel_max=dH, dQ2_rel_max=dQ2, wall_time_sec=wall,
                n_states=int(states.shape[0]),
                sha256_states=sha256_file(os.path.join(rdir, 'states.npz')))


def main():
    ok, gate_report = compat_gate()
    audit = audit_local()
    with open(os.path.join(HERE, 'pre_run_audit.json'), 'w', encoding='utf-8') as f:
        json.dump({'compat_gate': gate_report, 'integrator_local_audit': audit,
                   'commit_sha': COMMIT_SHA}, f, indent=1)
    if not ok:
        print('COMPAT GATE FAILED — 99 走行を開始せず停止（仕様 §8）')
        sys.exit(1)

    rows = build_runs()
    qa = []
    for i, row in enumerate(rows):
        r = do_run(row)
        qa.append(r)
        print(f'[{i+1:2d}/99] {r["run_id"]:>22}: {r["qa_status"]} '
              f'dH={r["dH_rel_max"]:.1e} dQ2={r["dQ2_rel_max"]:.1e} '
              f'{r["wall_time_sec"]:.1f}s', flush=True)

    n_pass = sum(1 for r in qa if r['qa_status'] == 'PASS')
    with open(os.path.join(HERE, 'qa_summary.json'), 'w', encoding='utf-8') as f:
        json.dump({'n_runs': len(qa), 'n_pass': n_pass, 'n_fail': len(qa) - n_pass,
                   'qa_tol': QA_TOL, 'commit_sha': COMMIT_SHA, 'runs': qa}, f, indent=1)

    # SHA256SUMS（全 run の全ファイル＋トップレベル）
    entries = []
    for row in rows:
        rdir = os.path.join(RUNS_DIR, row['run_id'])
        for name in ('states.npz', 'metadata.json', 'invariants.csv'):
            p = os.path.join(rdir, name)
            entries.append(f'{sha256_file(p)}  runs/{row["run_id"]}/{name}')
    for name in ('master_manifest.csv', 'qa_summary.json', 'pre_run_audit.json',
                 'gen_manifest.py', 'run_survey.py',
                 '二波完全関係力学_生軌道データベース収集仕様書_v1_20260915.md'):
        p = os.path.join(HERE, name)
        entries.append(f'{sha256_file(p)}  {name}')
    with open(os.path.join(HERE, 'SHA256SUMS.txt'), 'w') as f:
        f.write('\n'.join(entries) + '\n')

    complete = (len(qa) == 99 and all(r['n_states'] == T + 1 for r in qa))
    print(f'\nruns={len(qa)} PASS={n_pass} FAIL={len(qa)-n_pass}')
    print('DATA COLLECTION COMPLETE' if complete else 'INCOMPLETE — 完了宣言しない')


if __name__ == '__main__':
    main()
