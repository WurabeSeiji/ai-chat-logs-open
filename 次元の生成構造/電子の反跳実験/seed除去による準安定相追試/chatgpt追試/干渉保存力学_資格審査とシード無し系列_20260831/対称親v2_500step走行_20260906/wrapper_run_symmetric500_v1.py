#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対称親 v2 での 500step スイープ——専用ラッパー（物理プログラム無変更・旧データ不可侵）。

木原指示（2026-09-06）: 最新系の 500step 実験を「データだけ変更・物理プログラム一切不変・
過去データ破壊なし・データと図のフォルダを分離・専用ラッパーで実行」。

方式:
 - 物理の正本 `../N3_N40_stage123_sweep_20260905/run_N3_N40_stage123_v1.py` は読み取りのみ
   （冒頭で SHA256 を厳密照合）。ソース中の**データパス定数2行（PARENT_DIR / OUT）だけ**を
   文字列置換し（置換前後の差分行数=2 を検証）、exec で実行する。物理・手順（STEPS=500、
   分母 N−2..N+2,124、全状態保存、図・CSV・metadata 出力）は正本のまま。
 - 親ファイル名はループ内埋め込みのため触らず、ステージング側でファイル名を合わせる:
     parents_control/  = 旧親のコピー（bit 照合）→ 対照走行用
     parents_symmetric/= 対称親 v2 のコピーを旧名 parent_static_N{N:05d}_makeparent_20260905.npz
                         に改名配置（対応表と SHA を manifest に記録、v2 の SHA 台帳と照合）
 - 手順: [1] 対照走行（旧親→ results_control/）→ 既存正本 results/ と 228 npz の全配列 bit 一致・
   CSV 2本の byte 一致を検証（不一致なら中断） [2] 本走行（対称親 v2 → results/）。
 - 旧フォルダへは一切書き込まない。出力は本フォルダ配下のみ。
使い方: python3 wrapper_run_symmetric500_v1.py [--control-only|--skip-control]
"""
import hashlib
import json
import os
import shutil
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
SERIES = os.path.dirname(BASE)
ORIG_PKG = os.path.join(SERIES, 'N3_N40_stage123_sweep_20260905')
ORIG_PROG = os.path.join(ORIG_PKG, 'run_N3_N40_stage123_v1.py')
ORIG_RESULTS = os.path.join(ORIG_PKG, 'results')
ORIG_PARENTS = os.path.join(ORIG_PKG, 'parents')
SYM_PKG = os.path.join(SERIES, '最も対称性の高い初期値_20260906')
SYM_PARENTS = os.path.join(SYM_PKG, 'parents_symmetric')

PROG_SHA = '1abf2353fee2e4f56f05e7a6f149fd086885136beb61ab571b48a56b09691567'  # 正本走行プログラム（無変更の保証）
OLD_PD = "PARENT_DIR=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','parents')"
OLD_OUT = "OUT=os.path.join(ROOT,'N3_N40_stage123_sweep_20260905','results')"
NPZ_NAME = 'parent_static_N{N:05d}_makeparent_20260905.npz'   # 正本が読むファイル名（触らない）
N_MAX = None   # 本走行（指示書 §6: N=3..6 合格につき全 N へ拡張）


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def load_ledger(path):
    led = {}
    for line in open(path):
        p = line.split()
        if len(p) == 2:
            led[p[1]] = p[0]
    return led


def patched_source(parent_dir, out_dir):
    src = open(ORIG_PROG, encoding='utf-8').read()
    assert hashlib.sha256(src.encode()).hexdigest() == PROG_SHA, '正本プログラムのSHAが一致しない（無変更保証の破れ）'
    assert src.count(OLD_PD) == 1 and src.count(OLD_OUT) == 1, 'パス定数行が想定と異なる'
    new = src.replace(OLD_PD, f'PARENT_DIR={parent_dir!r}').replace(OLD_OUT, f'OUT={out_dir!r}')
    if N_MAX is not None:  # 試走（木原指示 2026-09-06: N=6 まで）: ループ上限のみ変更
        old_loop = 'for N in range(3,41):'
        assert new.count(old_loop) == 1
        new = new.replace(old_loop, f'for N in range(3,{N_MAX + 1}):')
    diff = sum(1 for a, b in zip(src.splitlines(), new.splitlines()) if a != b)
    expected = 2 if N_MAX is None else 3
    assert diff == expected and len(src.splitlines()) == len(new.splitlines()), f'置換行数が想定外: {diff}'
    return new


def run_patched(parent_dir, out_dir, tag):
    os.makedirs(out_dir, exist_ok=True)
    print(f'=== 走行 [{tag}] parents={os.path.basename(parent_dir)} → {os.path.basename(out_dir)} ===', flush=True)
    g = {'__name__': f'stage123_{tag}', '__file__': ORIG_PROG}
    exec(compile(patched_source(parent_dir, out_dir), ORIG_PROG, 'exec'), g)


def stage_control():
    d = os.path.join(BASE, 'parents_control')
    os.makedirs(d, exist_ok=True)
    for N in range(3, 41):
        name = NPZ_NAME.format(N=N)
        src = os.path.join(ORIG_PARENTS, name)
        dst = os.path.join(d, name)
        shutil.copyfile(src, dst)
        assert sha(src) == sha(dst)
    print('対照ステージング: 旧親38本 bit 照合済み', flush=True)
    return d


def stage_symmetric():
    led = load_ledger(os.path.join(SYM_PKG, 'SHA256SUMS.txt'))
    d = os.path.join(BASE, 'parents_symmetric_staged')
    os.makedirs(d, exist_ok=True)
    manifest = {}
    for N in range(3, 41):
        src_name = f'parent_symmetric_N{N:05d}_v2.npz'
        src = os.path.join(SYM_PARENTS, src_name)
        h = sha(src)
        assert led[f'parents_symmetric/{src_name}'] == h, f'v2 SHA台帳照合失敗: {src_name}'
        dst_name = NPZ_NAME.format(N=N)
        dst = os.path.join(d, dst_name)
        shutil.copyfile(src, dst)
        assert sha(dst) == h
        manifest[dst_name] = {'source': f'最も対称性の高い初期値_20260906/parents_symmetric/{src_name}', 'sha256': h}
    with open(os.path.join(BASE, 'parents_symmetric_staged_manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    print('対称親ステージング: v2 38本を旧名で配置（SHA台帳照合済み・対応表 manifest 保存）', flush=True)
    return d


def verify_control(out_dir):
    report = {'npz_checked': 0, 'npz_mismatch': [], 'csv': {}}
    for fn in sorted(os.listdir(ORIG_RESULTS)):
        if fn.endswith('.npz'):
            a = np.load(os.path.join(ORIG_RESULTS, fn))
            b = np.load(os.path.join(out_dir, fn))
            same = all(np.array_equal(a[k], b[k]) for k in a.files) and set(a.files) == set(b.files)
            report['npz_checked'] += 1
            if not same:
                report['npz_mismatch'].append(fn)
    for fn in ('timeseries_64bit_with124_N3_N40.csv', 'summary_64bit_with124_N3_N40.csv'):
        same = open(os.path.join(ORIG_RESULTS, fn), 'rb').read() == open(os.path.join(out_dir, fn), 'rb').read()
        report['csv'][fn] = 'byte一致' if same else '不一致'
    ok = (not report['npz_mismatch']) and all(v == 'byte一致' for v in report['csv'].values())
    report['CONTROL'] = 'PASS' if ok else 'FAIL'
    with open(os.path.join(BASE, 'control_report_v1.json'), 'w') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(f"対照検証: npz {report['npz_checked']} 本中 不一致 {len(report['npz_mismatch'])} / CSV {report['csv']} → {report['CONTROL']}", flush=True)
    return ok


if __name__ == '__main__':
    # v3 改修（2026-09-06 木原承認: 旧走行は残さないので修正可）:
    # 対照段・ステージング段を撤去。parents_symmetric_staged/ に配置済みの親
    # （make_parents_v3_centroid_zero_v1.py が書いた重心ゼロ系列 v3）をそのまま読み、
    # 正本プログラム（SHA照合・パス2行置換のみ）で本走行する。
    ps = os.path.join(BASE, 'parents_symmetric_staged')
    run_patched(ps, os.path.join(BASE, 'results'), 'symmetric_v3')
    print('WRAPPER ALL DONE', flush=True)
