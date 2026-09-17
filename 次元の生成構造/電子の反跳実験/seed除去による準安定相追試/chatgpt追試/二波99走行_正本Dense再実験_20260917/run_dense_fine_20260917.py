#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""細分時間刻み走行ランナー den=200L（2026-09-17 木原指示・fine 系列）。

背景: 実測診断で、現行 den=L は生成子スペクトル（Δτ·ω_max ≈ 9〜12 rad/step、
K の成長でさらに悪化）に対して約2桁のサンプリング不足と判明。木原指示により
den = 200·L（初期 Δτ·ω₀ ≈ 0.05 rad、K 成長 6 倍でも ~0.3 rad）の細分系列を
低位 L で実行する。

- 力学・カーネル・初期値生成器・manifest は一切変更しない（同一ファイルを参照）。
  変更は時間刻み den_fine = 200·L と step 数 T_FINE = 200·T = 819200 のみ
  （物理 τ 範囲は旧走行と同一。旧 τ 格子は 200 step ごとに厳密に一致）。
- 全 T_FINE+1 状態を間引かず保存（states.npz、complex128）。
- invariants.csv は毎 step（列は本系列既定の6列）。
- run_id は L{L}_ma{ma}_mb{mb}_den{200L}、保存先は runs_fine/。
- QA は機械的確認のみ（状態数・NaN/Inf・metadata・SHA256）。drift・反対称残差は
  生値記録・閾値分類なし。実行ゲートなし。人為的停止なし。
- 各 run 完了直後に fine 版図化4セット（初期複素 v4 / 終了複素 / インフレーション /
  残差時系列、plot_scripts_run/plot_*_fine_20260917.py、存在 run のみ図化）を自動実行し、
  図を常に最新の走行集合へ更新する（手動図化に依存しない＝再現性確保、2026-09-18 追加）。

使い方: python3 run_dense_fine_20260917.py <base_run_id> [<base_run_id> ...]
  base_run_id は 99-manifest の run_id（例 L8_ma1_mb2_den8）。den・T のみ置換される。
"""
import os

THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
for _v in THREAD_VARS:
    os.environ[_v] = '1'

import csv                                                     # noqa: E402
import hashlib                                                 # noqa: E402
import json                                                    # noqa: E402
import math                                                    # noqa: E402
import platform                                                # noqa: E402
import subprocess                                              # noqa: E402
import sys                                                     # noqa: E402
import time                                                    # noqa: E402

import numpy as np                                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from interaction_kernel_theory_v1 import adjacency, K_of, one_step  # noqa: E402
from generator_twowave_v1 import generate                           # noqa: E402
from gen_manifest import build_runs, T                              # noqa: E402

DEN_MULT = 200
T_FINE = T * DEN_MULT            # 819200（物理 τ 範囲は旧走行と同一）
RUNS_DIR = os.path.join(HERE, 'runs_fine')
QA_PATH = os.path.join(HERE, 'qa_summary_fine.json')
PROGRESS_LOG = os.path.join(HERE, 'run_progress_20260917.log')
CODE_FILES = ('interaction_kernel_theory_v1.py', 'generator_twowave_v1.py',
              'gen_manifest.py', 'run_dense_fine_20260917.py')
# fine 版図化4セット（存在 run のみ図化する版）。粗い系列と同様、各 run 完了直後に
# 自動実行して図を必ず最新の走行集合に更新する（手動図化に依存しない＝再現性確保）。
PLOT_SCRIPTS = ('plot_initial_states_fine_20260917.py',
                'plot_final_states_fine_20260917.py',
                'plot_inflation_hperp_fine_20260917.py',
                'plot_residual_timeseries_fine_20260917.py')


def log_line(msg):
    line = f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())} UTC] {msg}'
    print(line, flush=True)
    with open(PROGRESS_LOG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def run_plots_after(rid):
    """fine 図化4セットを順に実行（read-only 後処理）。失敗しても走行は続行する。"""
    pdir = os.path.join(HERE, 'plot_scripts_run')
    logdir = os.path.join(HERE, 'plot_logs')
    os.makedirs(logdir, exist_ok=True)
    plog = os.path.join(logdir, f'auto_after_fine_{rid}.log')
    ok_all = True
    with open(plog, 'w', encoding='utf-8') as f:
        for name in PLOT_SCRIPTS:
            t0 = time.time()
            r = subprocess.run([sys.executable, name], cwd=pdir,
                               stdout=f, stderr=subprocess.STDOUT)
            dt = time.time() - t0
            ok = (r.returncode == 0)
            ok_all &= ok
            log_line(f'  figure {name}: {"ok" if ok else f"FAIL(exit {r.returncode})"} '
                     f'{dt:.1f}s')
    return ok_all


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def do_run(row):
    L, ma, mb = row['L'], row['ma'], row['mb']
    den = L * DEN_MULT
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    rdir = os.path.join(RUNS_DIR, rid)
    if os.path.exists(os.path.join(rdir, 'states.npz')):
        raise SystemExit(f'{rid}: states.npz が既に存在。上書きしない')
    os.makedirs(rdir, exist_ok=True)
    A = adjacency(L)
    z0 = generate(L, ma, mb)
    z0_sha = hashlib.sha256(z0.tobytes()).hexdigest()
    M = len(z0)
    t_start = time.time()
    states = np.empty((T_FINE + 1, M), np.complex128)
    inv = np.empty((T_FINE + 1, 5), np.float64)
    z = z0.copy()
    for t in range(T_FINE + 1):
        states[t] = z
        K = K_of(z, A)
        inv[t] = [float(np.vdot(z, z).real),
                  float(np.sum(z ** 2).real), float(np.sum(z ** 2).imag),
                  float(np.linalg.norm(K @ z) / np.linalg.norm(z)),
                  float(np.max(np.abs(K + K.T)))]
        if t < T_FINE:
            z = one_step(z, A, den)
    wall = time.time() - t_start

    np.savez(os.path.join(rdir, 'states.npz'), Z=states,
             L=np.int64(L), ma=np.int64(ma), mb=np.int64(mb),
             den=np.int64(den), T=np.int64(T_FINE))
    with open(os.path.join(rdir, 'invariants.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['t', 'H', 'Q2_re', 'Q2_im', 'generator_norm', 'antisym_residual'])
        w.writerows([[t] + [repr(x) for x in inv[t]] for t in range(T_FINE + 1)])

    H0 = inv[0, 0]
    Q20 = complex(inv[0, 1], inv[0, 2])
    dH = float(np.max(np.abs(inv[:, 0] - H0) / H0))
    Q2t = inv[:, 1] + 1j * inv[:, 2]
    dQ2 = float(np.max(np.abs(Q2t - Q20)) / max(H0, abs(Q20)))
    antisym_max = float(np.max(inv[:, 4]))
    finite = bool(np.all(np.isfinite(states.view(np.float64))))
    reasons = []
    if not finite:
        reasons.append('NaN/Inf in states')
    if states.shape[0] != T_FINE + 1:
        reasons.append('missing states')
    qa_status = 'PASS' if not reasons else 'FAIL'

    sha_states = sha256_file(os.path.join(rdir, 'states.npz'))
    meta = dict(row)
    meta.update(dict(
        run_id=rid, den=den, T=T_FINE, den_mult=DEN_MULT,
        base_run_id=row['run_id'], base_den=row['L'],
        dtau=2.0 * math.pi / den,
        series='fine_den200L',
        algorithm='full dense K -> H=iK -> np.linalg.eigh 全固有分解 -> '
                  'スペクトル指数写像（v3 §2 正本、近似なし。時間刻みのみ den=200L）',
        interaction='K_ef = A_ef * Im(conj(z_e) z_f)（振幅込み・係数1・補正なし）',
        code_commit_sha=subprocess.run(['git', 'rev-parse', 'HEAD'],
                                       capture_output=True, text=True,
                                       cwd=HERE).stdout.strip(),
        code_sha256={n: sha256_file(os.path.join(HERE, n)) for n in CODE_FILES},
        z0_sha256=z0_sha, sha256_states=sha_states,
        python=platform.python_version(), numpy=np.__version__,
        thread_env={v: os.environ.get(v) for v in THREAD_VARS},
        wall_time_sec=wall, qa_status=qa_status, qa_reason='; '.join(reasons),
        dH_rel_max=dH, dQ2_rel_max=dQ2, antisym_residual_max=antisym_max))
    with open(os.path.join(rdir, 'metadata.json'), 'w', encoding='utf-8') as f:
        json.dump(meta, f, indent=1, ensure_ascii=False)
    with open(os.path.join(rdir, 'SHA256SUMS.txt'), 'w') as f:
        f.write(f'{sha_states}  states.npz\n'
                f'{sha256_file(os.path.join(rdir, "invariants.csv"))}  invariants.csv\n'
                f'{sha256_file(os.path.join(rdir, "metadata.json"))}  metadata.json\n')
    return dict(run_id=rid, qa_status=qa_status, qa_reason='; '.join(reasons),
                dH_rel_max=dH, dQ2_rel_max=dQ2, antisym_residual_max=antisym_max,
                wall_time_sec=wall, n_states=int(states.shape[0]),
                z0_sha256=z0_sha, sha256_states=sha_states)


def main():
    args = sys.argv[1:]
    if not args:
        print('usage: python3 run_dense_fine_20260917.py <base_run_id> [<base_run_id> ...] '
              '| --upto-L <N>')
        sys.exit(2)
    rows = build_runs()
    by_id = {r['run_id']: r for r in rows}
    if args[0] == '--upto-L':
        if len(args) != 2 or not args[1].isdigit():
            print('usage: python3 run_dense_fine_20260917.py --upto-L <N>')
            sys.exit(2)
        lmax = int(args[1])
        # fine では den=200L に固定されるため、(L,ma,mb) が同じ den40 対照は主系列と
        # 同一 fine run_id に潰れる。fine run_id で重複排除し、既存 states.npz はスキップ。
        selected, seen = [], set()
        for r in rows:
            if r['L'] > lmax:
                continue
            frid = f'L{r["L"]}_ma{r["ma"]}_mb{r["mb"]}_den{r["L"] * DEN_MULT}'
            if frid in seen or os.path.exists(os.path.join(RUNS_DIR, frid, 'states.npz')):
                continue
            seen.add(frid)
            selected.append(r)
        if not selected:
            print(f'L<={lmax} の未走行 fine run はありません')
            sys.exit(0)
        print(f'--upto-L {lmax}: 未走行 {len(selected)} fine run を実行'
              '（走行済み・den40対照の重複はスキップ）')
    else:
        for rid in args:
            if rid not in by_id:
                print(f'unknown base run_id: {rid}')
                sys.exit(2)
        selected = [by_id[rid] for rid in args]
    qa_all = {}
    if os.path.exists(QA_PATH):
        with open(QA_PATH, encoding='utf-8') as f:
            qa_all = {r['run_id']: r for r in json.load(f).get('runs', [])}
    log_line(f'fine batch start (den=200L, T={T_FINE}): '
             + ', '.join(r['run_id'] for r in selected))
    for i, row in enumerate(selected):
        log_line(f'[fine {i + 1}/{len(selected)}] {row["run_id"]} -> den={row["L"] * DEN_MULT} start '
                 f'(M={row["M"]}, T={T_FINE})')
        r = do_run(row)
        qa_all[r['run_id']] = r
        log_line(f'[fine {i + 1}/{len(selected)}] {r["run_id"]}: {r["qa_status"]} '
                 f'dH={r["dH_rel_max"]:.1e} dQ2={r["dQ2_rel_max"]:.1e} '
                 f'antisym={r["antisym_residual_max"]:.1e} {r["wall_time_sec"]:.1f}s')
        runs_list = sorted(qa_all.values(), key=lambda x: x['run_id'])
        with open(QA_PATH, 'w', encoding='utf-8') as f:
            json.dump({'n_runs': len(runs_list),
                       'n_pass': sum(1 for x in runs_list if x['qa_status'] == 'PASS'),
                       'n_fail': sum(1 for x in runs_list if x['qa_status'] == 'FAIL'),
                       'series': 'fine_den200L',
                       'runs': runs_list}, f, indent=1, ensure_ascii=False)
        run_plots_after(r['run_id'])
    log_line(f'fine batch done: {len(selected)} run(s)')


if __name__ == '__main__':
    main()
