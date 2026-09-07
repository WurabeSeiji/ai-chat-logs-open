#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_parent 型データの 3D 図化（指示書 ClaudeCode_3D図化指示_makeparent_selectedN_20260907.md 準拠）。

基準 3D 図: N16_original_note_source_den16_step0_500_3D.html（ChatGPT 提供・仕様の正本）。
本プログラムは基準 HTML を**テンプレートとして読み込み**、レイアウト・軸・アニメーション構造
（PC1×PC2×H⊥/H の3D軌跡＋下段H⊥/H時系列＋Play/Pause＋スライダー、501フレーム）を一切変えず、
以下だけを差し替える（指示書の変更方針 1〜4）:
  1. 入力: full_N3_N40_sweep/states/hm_N{N}_den_{N}_states_500.npz（den=N のみ）
  2. 対象 N = [3, 4, 5, 17, 40]
  3. 出力先: fig3d_makeparent_selectedN/
  4. ファイル名・表題に makeparent と N を含める

量の定義（基準から逆同定し N=16 対照で検証）:
  - p,q = 正本 run_N3_N40_stage123_v1.plane(z0) と同一式（Re z0 正規化、Im z0 を直交化）
  - cp(t)=p·z(t), cq(t)=q·z(t)（複素射影係数）。V=[Re cp, Im cp, Re cq, Im cq] の中心化 PCA 上位2軸が
    PC1/PC2（基準と最大差 2e-14。符号は表示規約: 決定論の svd_flip(U基準) を採用。
    基準サンプルは PC2 の符号が逆＝鏡映で、物理内容は不変）
  - z軸・下段 = H⊥/H（正本 metrics() と同一式、表示用に 1e-31 で床クリップ。基準と最大差 8e-17）
使い方: python3 plot3d_makeparent_selectedN_v1.py            # 5本を生成
        python3 plot3d_makeparent_selectedN_v1.py --control  # N=16 を再生成し基準と配列照合"""
import base64
import copy
import json
import os
import sys
import uuid

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(BASE, 'N16_original_note_source_den16_step0_500_3D.html')
STATES = os.path.join(BASE, 'full_N3_N40_sweep', 'states')
OUTDIR = os.path.join(BASE, 'fig3d_makeparent_selectedN')
NS = [3, 4, 5, 17, 40]
FLOOR = 1e-31


def parse_template():
    t = open(TEMPLATE, encoding='utf-8').read()

    def block(tok, oc, cc):
        i = t.find(tok); j = t.find(oc, i); depth = 0
        for k in range(j, len(t)):
            if t[k] == oc:
                depth += 1
            elif t[k] == cc:
                depth -= 1
                if depth == 0:
                    return j, k + 1
    dj, dk = block('Plotly.newPlot', '[', ']')
    lj = t.find('{', dk); depth = 0
    for k in range(lj, len(t)):
        if t[k] == '{':
            depth += 1
        elif t[k] == '}':
            depth -= 1
            if depth == 0:
                lk = k + 1; break
    fj, fk = block('Plotly.addFrames', '[', ']')
    data = json.loads(t[dj:dk].replace('\\u002f', '/'))
    layout = json.loads(t[lj:lk].replace('\\u002f', '/'))
    frames = json.loads(t[fj:fk].replace('\\u002f', '/'))
    return t, (dj, dk), (lj, lk), (fj, fk), data, layout, frames


def b64(arr, dtype):
    return {'dtype': dtype, 'bdata': base64.b64encode(np.asarray(arr).astype(dtype).tobytes()).decode()}


def series(N):
    d = np.load(os.path.join(STATES, f'hm_N{N}_den_{N}_states_500.npz'))
    assert int(d['denominator']) == N and int(d['steps']) == 500
    S = np.asarray(d['Z'], np.complex128)
    z0 = S[0]
    p = z0.real.astype(np.float64).copy(); p /= np.linalg.norm(p)
    q = z0.imag.astype(np.float64).copy(); q -= np.dot(q, p) * p; q /= np.linalg.norm(q)
    T = S.shape[0]
    f = np.empty(T); cp = np.empty(T, complex); cq = np.empty(T, complex)
    for s in range(T):
        z = S[s]
        a = np.dot(p, z); b = np.dot(q, z)
        zp = z - p * a - q * b
        f[s] = np.vdot(zp, zp).real / np.vdot(z, z).real
        cp[s] = a; cq[s] = b
    fz = np.maximum(f, FLOOR)
    V = np.stack([cp.real, cp.imag, cq.real, cq.imag], axis=1)
    W = V - V.mean(0)
    U, sv, Vt = np.linalg.svd(W, full_matrices=False)
    signs = np.sign(U[np.argmax(np.abs(U), axis=0), np.arange(U.shape[1])])  # svd_flip(U基準)・決定論
    PC = (U * signs) * sv
    return PC[:, 0], PC[:, 1], fz


def build(N, title_h1, title_sub, ann_text, out_path, data0, layout0, frames0, shell, spans):
    X, Y, Z = series(N)
    T = len(X)
    data = copy.deepcopy(data0); layout = copy.deepcopy(layout0); frames = copy.deepcopy(frames0)
    # base traces
    data[0]['customdata'] = b64([0], 'i1')
    data[0]['x'] = b64(X[:1], 'f8'); data[0]['y'] = b64(Y[:1], 'f8'); data[0]['z'] = b64(Z[:1], 'f8')
    data[1]['x'] = [X[0]]; data[1]['y'] = [Y[0]]; data[1]['z'] = [Z[0]]; data[1]['text'] = ['step 0']
    data[2]['x'] = b64(np.arange(T), 'i2'); data[2]['y'] = b64(Z, 'f8')
    data[3]['x'] = [0, 0]; data[3]['y'] = [FLOOR, float(Z[0])]
    data[4]['x'] = [0]; data[4]['y'] = [float(Z[0])]; data[4]['text'] = ['step 0']
    # annotations（下段の題のみ差替え）
    layout['annotations'][1]['text'] = ann_text
    # frames
    for k, fr in enumerate(frames):
        fd = fr['data']
        fd[0]['customdata'] = b64(np.arange(k + 1), 'i2')
        fd[0]['x'] = b64(X[:k + 1], 'f8'); fd[0]['y'] = b64(Y[:k + 1], 'f8'); fd[0]['z'] = b64(Z[:k + 1], 'f8')
        fd[1]['x'] = [X[k]]; fd[1]['y'] = [Y[k]]; fd[1]['z'] = [Z[k]]; fd[1]['text'] = [f'step {k}']
        fd[2]['x'] = [k, k]; fd[2]['y'] = [FLOOR, float(Z[k])]
        fd[3]['x'] = [k]; fd[3]['y'] = [float(Z[k])]; fd[3]['text'] = [f'step {k}']
    # HTML 組み立て（シェルは基準のまま、uuid・タイトルのみ差替え）
    t = shell
    (dj, dk), (lj, lk), (fj, fk) = spans
    enc = lambda o: json.dumps(o, separators=(',', ':'), ensure_ascii=False).replace('/', '\\u002f')
    t = t[:fj] + enc(frames) + t[fk:]
    t = t[:lj] + enc(layout) + t[lk:]
    t = t[:dj] + enc(data) + t[dk:]
    old_id = 'c6d019db-3da4-48ad-bc8b-1d2fde4c20b1'
    t = t.replace(old_id, str(uuid.uuid4()))
    t = t.replace('<title>N16 original note source</title>', f'<title>{title_h1}</title>')
    t = t.replace('<h1>N=16 original note source / den=16 / step 0→500</h1>', f'<h1>{title_h1}</h1>')
    t = t.replace('Source: hm_N16_den_16_states_500.npz materialized directly from Google Drive. No dynamics rerun.',
                  title_sub)
    open(out_path, 'w', encoding='utf-8').write(t)
    return X, Y, Z


def main():
    shell, sp_d, sp_l, sp_f, data0, layout0, frames0 = parse_template()
    spans = (sp_d, sp_l, sp_f)
    os.makedirs(OUTDIR, exist_ok=True)
    if '--control' in sys.argv:
        out = os.path.join(OUTDIR, 'CONTROL_N16_regenerated_3D.html')
        X, Y, Z = build(16, 'N=16 original note source / den=16 / step 0→500',
                        'Source: hm_N16_den_16_states_500.npz materialized directly from Google Drive. No dynamics rerun.',
                        'Original N=16 inflation curve Hperp/H', out, data0, layout0, frames0, shell, spans)
        # 基準の最終フレームと照合
        ref = frames0[-1]['data'][0]
        dec = lambda o: np.frombuffer(base64.b64decode(o['bdata']), dtype=o['dtype'])
        rX, rY, rZ = dec(ref['x']), dec(ref['y']), dec(ref['z'])
        print('CONTROL N16: PC1 max|Δ|=%.2e  PC2(符号込) max|Δ|=%.2e  PC2(鏡映) max|Δ|=%.2e  z max|Δ|=%.2e'
              % (np.max(np.abs(X - rX)), np.max(np.abs(Y - rY)), np.max(np.abs(-Y - rY)), np.max(np.abs(Z - rZ))))
        return
    for N in NS:
        name = f'N{N}_makeparent_den{N}_step0_500_3D.html'
        build(N, f'N={N} make_parent / den={N} / step 0→500',
              f'Source: full_N3_N40_sweep/states/hm_N{N}_den_{N}_states_500.npz (make_parent full-N sweep, bit-identical to canonical). No dynamics rerun.',
              f'make_parent N={N} inflation curve Hperp/H',
              os.path.join(OUTDIR, name), data0, layout0, frames0, shell, spans)
        print('wrote', name)
    print('ALL DONE')


if __name__ == '__main__':
    main()
