#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""二波走行・理論正本 Dense 再実験ランナー（v3 指示書準拠、2026-09-17）。

力学は interaction_kernel_theory_v1.one_step のみを使用する（外部 import なし）。
毎 step: full dense K → H=iK → np.linalg.eigh 全固有分解 → スペクトル指数写像。
Sparse / Lanczos / tol / 打切り / 正規化 / 丸め / 補正は一切ない。

実行対象は run_id 引数で明示指定する（段階実行のため。無指定では走らない）。
99 走行全体は --all99 を明示した場合のみ。

実行ゲートは持たない（木原・ChatGPT 指示 2026-09-17: 監査するのはコード、
判定するのは実験後のデータ。物理量・閉塞量・保存量による人為的停止をしない）。
保存形式は旧サーベイと互換:
  runs/<run_id>/states.npz   … Z (4097×M complex128), L, ma, mb, den, T
  runs/<run_id>/invariants.csv … t,H,Q2_re,Q2_im,generator_norm,antisym_residual
                                 （旧5列＋末尾に反対称残差1列を追加した上位互換）
  runs/<run_id>/metadata.json
QA は機械的確認のみ（4097 状態・NaN/Inf・metadata・SHA256）。H/Q2 drift と
反対称残差は生値で記録し、閾値分類しない（v3 §10）。FAIL でも states を削除しない。

使い方: python3 run_dense_rerun.py <run_id> [<run_id> ...]
        python3 run_dense_rerun.py --upto-L <N>   … manifest の L<=N の未走行 run を全て実行
        python3 run_dense_rerun.py --all99

2026-09-17 木原指示による追加（力学は無変更）:
 - 各 run の完了直後に図化4セット（初期複素 v4 / 終了複素 / インフレーション /
   残差時系列、plot_scripts_run/ の存在runのみ図化版）を subprocess で自動実行
 - 経過ログ run_progress_20260917.log（UTC時刻付き）に run 開始/完了/図化の
   所要時間を追記（stdout にも同内容）。図化失敗は記録して続行（走行は止めない）
"""
import os

THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
for _v in THREAD_VARS:
    os.environ[_v] = '1'

import contextlib                                              # noqa: E402
import csv                                                     # noqa: E402
import hashlib                                                 # noqa: E402
import io                                                      # noqa: E402
import json                                                    # noqa: E402
import math                                                    # noqa: E402
import platform                                                # noqa: E402
import subprocess                                              # noqa: E402
import sys                                                     # noqa: E402
import time                                                    # noqa: E402

import numpy as np                                             # noqa: E402
try:
    import scipy                                               # noqa: E402
    SCIPY_VERSION = scipy.__version__
except ImportError:
    SCIPY_VERSION = None

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from interaction_kernel_theory_v1 import adjacency, K_of, one_step  # noqa: E402
from generator_twowave_v1 import generate                           # noqa: E402
from gen_manifest import build_runs, T                              # noqa: E402

RUNS_DIR = os.path.join(HERE, 'runs')
FP_PATH = os.path.join(HERE, 'environment_fingerprint.json')
QA_PATH = os.path.join(HERE, 'qa_summary.json')
PROGRESS_LOG = os.path.join(HERE, 'run_progress_20260917.log')
CODE_FILES = ('interaction_kernel_theory_v1.py', 'generator_twowave_v1.py',
              'gen_manifest.py', 'run_dense_rerun.py', 'structural_audit_v1.py')
PLOT_SCRIPTS = ('plot_initial_states_all99_v4_20260916.py',
                'plot_final_states_all99_20260916.py',
                'plot_inflation_hperp_all99_20260916.py',
                'plot_residual_timeseries_all99_20260916.py')


def log_line(msg):
    line = f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())} UTC] {msg}'
    print(line, flush=True)
    with open(PROGRESS_LOG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def run_plots_after(rid):
    """図化4セットを順に実行（read-only 後処理）。失敗しても走行は続行する。"""
    pdir = os.path.join(HERE, 'plot_scripts_run')
    logdir = os.path.join(HERE, 'plot_logs')
    os.makedirs(logdir, exist_ok=True)
    plog = os.path.join(logdir, f'auto_after_{rid}.log')
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


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def git_sha():
    try:
        return subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True,
                              text=True, cwd=HERE).stdout.strip()
    except Exception:
        return ''


def cpu_brand():
    if platform.system() == 'Darwin':
        try:
            out = subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'],
                                 capture_output=True, text=True).stdout.strip()
            if out:
                return out
        except Exception:
            pass
    return platform.processor() or platform.machine()


def blas_config():
    try:
        return np.show_config(mode='dicts')
    except TypeError:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            np.show_config()
        return buf.getvalue()


def write_fingerprint(argv):
    rec = {
        'timestamp_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'argv': argv,
        'os': platform.platform(),
        'kernel': platform.uname().release,
        'cpu': cpu_brand(),
        'machine': platform.machine(),
        'python': platform.python_version(),
        'numpy': np.__version__,
        'scipy': SCIPY_VERSION,
        'blas_lapack': blas_config(),
        'thread_env': {v: os.environ.get(v) for v in THREAD_VARS},
        'git_commit_sha': git_sha(),
        'code_sha256': {name: sha256_file(os.path.join(HERE, name))
                        for name in CODE_FILES},
    }
    hist = []
    if os.path.exists(FP_PATH):
        with open(FP_PATH, encoding='utf-8') as f:
            hist = json.load(f)
    hist.append(rec)
    with open(FP_PATH, 'w', encoding='utf-8') as f:
        json.dump(hist, f, indent=1, ensure_ascii=False)
    return rec


def do_run(row, fp):
    run_id = row['run_id']
    rdir = os.path.join(RUNS_DIR, run_id)
    if os.path.exists(os.path.join(rdir, 'states.npz')):
        raise SystemExit(f'{run_id}: states.npz が既に存在。上書きしない（手動で退避してから再実行）')
    os.makedirs(rdir, exist_ok=True)
    L, ma, mb, den = row['L'], row['ma'], row['mb'], row['den']
    A = adjacency(L)
    z0 = generate(L, ma, mb)
    z0_sha = sha256_bytes(z0.tobytes())
    M = len(z0)
    t_start = time.time()
    states = np.empty((T + 1, M), np.complex128)
    inv = np.empty((T + 1, 5), np.float64)
    z = z0.copy()
    for t in range(T + 1):
        states[t] = z
        K = K_of(z, A)              # 診断用の読出し（状態は書き換えない）
        inv[t] = [float(np.vdot(z, z).real),
                  float(np.sum(z ** 2).real), float(np.sum(z ** 2).imag),
                  float(np.linalg.norm(K @ z) / np.linalg.norm(z)),
                  float(np.max(np.abs(K + K.T)))]
        if t < T:
            z = one_step(z, A, den)
    wall = time.time() - t_start

    np.savez(os.path.join(rdir, 'states.npz'), Z=states,
             L=np.int64(L), ma=np.int64(ma), mb=np.int64(mb),
             den=np.int64(den), T=np.int64(T))
    with open(os.path.join(rdir, 'invariants.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['t', 'H', 'Q2_re', 'Q2_im', 'generator_norm', 'antisym_residual'])
        w.writerows([[t] + [repr(x) for x in inv[t]] for t in range(T + 1)])

    # QA（機械的確認のみ。drift/残差は生値記録・閾値分類しない）
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
    if states.shape[0] != T + 1:
        reasons.append('missing states')
    qa_status = 'PASS' if not reasons else 'FAIL'

    sha_states = sha256_file(os.path.join(rdir, 'states.npz'))
    sha_inv = sha256_file(os.path.join(rdir, 'invariants.csv'))
    meta = dict(row)
    meta.update(dict(
        dtau=2.0 * math.pi / den,
        algorithm='full dense K -> H=iK -> np.linalg.eigh 全固有分解 -> '
                  'スペクトル指数写像（v3 指示書 §2 正本、近似なし）',
        interaction='K_ef = A_ef * Im(conj(z_e) z_f)（振幅込み・係数1・補正なし）',
        initial_condition='z_jk(0)=exp(-i*pi/(2L))*U^(ma*(k-j))'
                          '+exp(+i*pi/(2L))*U^(mb*(k-j)), U=exp(2*pi*i/L), j<k',
        code_commit_sha=fp['git_commit_sha'],
        kernel_sha256=fp['code_sha256']['interaction_kernel_theory_v1.py'],
        generator_sha256=fp['code_sha256']['generator_twowave_v1.py'],
        manifest_sha256=fp['code_sha256']['gen_manifest.py'],
        runner_sha256=fp['code_sha256']['run_dense_rerun.py'],
        z0_sha256=z0_sha, sha256_states=sha_states, sha256_invariants=sha_inv,
        python=fp['python'], numpy=fp['numpy'], scipy=fp['scipy'],
        thread_env=fp['thread_env'], wall_time_sec=wall,
        qa_status=qa_status, qa_reason='; '.join(reasons),
        dH_rel_max=dH, dQ2_rel_max=dQ2, antisym_residual_max=antisym_max))
    with open(os.path.join(rdir, 'metadata.json'), 'w', encoding='utf-8') as f:
        json.dump(meta, f, indent=1, ensure_ascii=False)
    sums = [f'{sha_states}  states.npz', f'{sha_inv}  invariants.csv',
            f'{sha256_file(os.path.join(rdir, "metadata.json"))}  metadata.json']
    with open(os.path.join(rdir, 'SHA256SUMS.txt'), 'w') as f:
        f.write('\n'.join(sums) + '\n')
    return dict(run_id=run_id, qa_status=qa_status, qa_reason='; '.join(reasons),
                dH_rel_max=dH, dQ2_rel_max=dQ2, antisym_residual_max=antisym_max,
                wall_time_sec=wall, n_states=int(states.shape[0]),
                z0_sha256=z0_sha, sha256_states=sha_states)


def main():
    args = sys.argv[1:]
    if not args:
        print('usage: python3 run_dense_rerun.py <run_id> [<run_id> ...] '
              '| --upto-L <N> | --all99')
        sys.exit(2)
    rows = build_runs()
    by_id = {r['run_id']: r for r in rows}
    if args == ['--all99']:
        selected = rows
    elif args[0] == '--upto-L':
        if len(args) != 2 or not args[1].isdigit():
            print('usage: python3 run_dense_rerun.py --upto-L <N>')
            sys.exit(2)
        lmax = int(args[1])
        selected = [r for r in rows if r['L'] <= lmax and
                    not os.path.exists(os.path.join(RUNS_DIR, r['run_id'], 'states.npz'))]
        if not selected:
            print(f'L<={lmax} の未走行 run はありません')
            sys.exit(0)
        print(f'--upto-L {lmax}: 未走行 {len(selected)} run を実行（走行済みはスキップ）')
    else:
        for rid in args:
            if rid not in by_id:
                print(f'unknown run_id (manifest 99 に無い): {rid}')
                sys.exit(2)
        selected = [by_id[rid] for rid in args]

    fp = write_fingerprint(sys.argv)
    qa_all = {}
    if os.path.exists(QA_PATH):
        with open(QA_PATH, encoding='utf-8') as f:
            qa_all = {r['run_id']: r for r in json.load(f).get('runs', [])}
    log_line(f'batch start: {len(selected)} run(s) — '
             + ', '.join(r['run_id'] for r in selected))
    for i, row in enumerate(selected):
        log_line(f'[{i + 1:2d}/{len(selected)}] {row["run_id"]} start '
                 f'(M={row["M"]}, T={row["T"]})')
        r = do_run(row, fp)
        qa_all[r['run_id']] = r
        log_line(f'[{i + 1:2d}/{len(selected)}] {r["run_id"]}: {r["qa_status"]} '
                 f'dH={r["dH_rel_max"]:.1e} dQ2={r["dQ2_rel_max"]:.1e} '
                 f'antisym={r["antisym_residual_max"]:.1e} {r["wall_time_sec"]:.1f}s')
        runs_list = sorted(qa_all.values(), key=lambda x: x['run_id'])
        with open(QA_PATH, 'w', encoding='utf-8') as f:
            json.dump({'n_runs': len(runs_list),
                       'n_pass': sum(1 for x in runs_list if x['qa_status'] == 'PASS'),
                       'n_fail': sum(1 for x in runs_list if x['qa_status'] == 'FAIL'),
                       'qa_rule': '機械的確認のみ（4097 states / NaN/Inf / metadata / SHA256）。'
                                  'drift・反対称残差は生値記録、閾値分類なし',
                       'runs': runs_list}, f, indent=1, ensure_ascii=False)
        run_plots_after(r['run_id'])
    log_line(f'batch done: {len(selected)} run(s)')


if __name__ == '__main__':
    main()
