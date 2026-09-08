#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=64・2000歩の 3D 図（PC1×PC2×H⊥/H 軌跡＋下段時系列＋Play/スライダー、対話HTML）。

基準: ../N16_original_note_source_den16_step0_500_3D.html（実験14 で仕様確定済みのテンプレート）。
量の定義は plot3d_makeparent_selectedN_v1.py と同一:
  p,q=正本 plane(z0)、cp=p·z, cq=q·z、V=[Re cp,Im cp,Re cq,Im cq] の中心化PCA上位2軸（svd_flip U基準）、
  z軸=正本 metrics() の H⊥/H（1e-31 床クリップ）。
2000歩対応の変更（テンプレート構造 501 フレームは不変）:
  - フレームは stride=4 のサンプル（step 0,4,…,2000 → 501 フレーム。軌跡線は毎歩の全点を保持）
  - 2D x軸 range [0,500]→[0,2000]、スライダーのラベルを実 step 値へ
  - タイトル・注釈・出典を N=64 用に差替え
入力: results/hm_N64_den_64_states_2000.npz
出力: fig3d_N64_makeparent_den64_step0_2000_3D_light.html（軽量版: 3D軌跡線も stride=4 に間引き、ファイル約5MB）"""
import base64
import copy
import json
import os
import uuid

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(BASE, '..', 'N16_original_note_source_den16_step0_500_3D.html')
NPZ = os.path.join(BASE, 'results', 'hm_N64_den_64_states_2000.npz')
OUT = os.path.join(BASE, 'fig3d_N64_makeparent_den64_step0_2000_3D_light.html')
FLOOR = 1e-31
STRIDE = 4  # 2000/4=500 → フレーム 501 本（テンプレートと同数）


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
    return (t, (dj, dk), (lj, lk), (fj, fk),
            json.loads(t[dj:dk].replace('\\u002f', '/')),
            json.loads(t[lj:lk].replace('\\u002f', '/')),
            json.loads(t[fj:fk].replace('\\u002f', '/')))


def b64(arr, dtype):
    return {'dtype': dtype, 'bdata': base64.b64encode(np.asarray(arr).astype(dtype).tobytes()).decode()}


def series():
    d = np.load(NPZ)
    assert int(d['denominator']) == 64 and int(d['steps']) == 2000
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
    signs = np.sign(U[np.argmax(np.abs(U), axis=0), np.arange(U.shape[1])])
    PC = (U * signs) * sv
    return PC[:, 0], PC[:, 1], fz


def main():
    shell, (dj, dk), (lj, lk), (fj, fk), data0, layout0, frames0 = parse_template()
    X, Y, Z = series()
    T = len(X)                      # 2001
    steps = list(range(0, T, STRIDE))
    assert len(steps) == len(frames0) == 501
    data = copy.deepcopy(data0); layout = copy.deepcopy(layout0); frames = copy.deepcopy(frames0)
    data[0]['customdata'] = b64([0], 'i1')
    data[0]['x'] = b64(X[:1], 'f8'); data[0]['y'] = b64(Y[:1], 'f8'); data[0]['z'] = b64(Z[:1], 'f8')
    data[1]['x'] = [X[0]]; data[1]['y'] = [Y[0]]; data[1]['z'] = [Z[0]]; data[1]['text'] = ['step 0']
    data[2]['x'] = b64(np.arange(T), 'i2'); data[2]['y'] = b64(Z, 'f8')
    data[3]['x'] = [0, 0]; data[3]['y'] = [FLOOR, float(Z[0])]
    data[4]['x'] = [0]; data[4]['y'] = [float(Z[0])]; data[4]['text'] = ['step 0']
    layout['annotations'][1]['text'] = 'make_parent N=64 inflation curve Hperp/H (2000 steps)'
    layout['xaxis']['range'] = [0, 2000]
    for j, st in enumerate(steps):
        layout['sliders'][0]['steps'][j]['label'] = str(st)
        fd = frames[j]['data']
        idx = np.arange(0, st + 1, STRIDE)
        fd[0]['customdata'] = b64(idx, 'i2')
        fd[0]['x'] = b64(X[idx], 'f8'); fd[0]['y'] = b64(Y[idx], 'f8'); fd[0]['z'] = b64(Z[idx], 'f8')
        fd[1]['x'] = [X[st]]; fd[1]['y'] = [Y[st]]; fd[1]['z'] = [Z[st]]; fd[1]['text'] = [f'step {st}']
        fd[2]['x'] = [st, st]; fd[2]['y'] = [FLOOR, float(Z[st])]
        fd[3]['x'] = [st]; fd[3]['y'] = [float(Z[st])]; fd[3]['text'] = [f'step {st}']
    enc = lambda o: json.dumps(o, separators=(',', ':'), ensure_ascii=False).replace('/', '\\u002f')
    t = shell
    t = t[:fj] + enc(frames) + t[fk:]
    t = t[:lj] + enc(layout) + t[lk:]
    t = t[:dj] + enc(data) + t[dk:]
    t = t.replace('c6d019db-3da4-48ad-bc8b-1d2fde4c20b1', str(uuid.uuid4()))
    h1 = 'N=64 make_parent / den=64 / step 0→2000'
    t = t.replace('<title>N16 original note source</title>', f'<title>{h1}</title>')
    t = t.replace('<h1>N=16 original note source / den=16 / step 0→500</h1>', f'<h1>{h1}</h1>')
    t = t.replace('Source: hm_N16_den_16_states_500.npz materialized directly from Google Drive. No dynamics rerun.',
                  'Source: results/hm_N64_den_64_states_2000.npz (N=64, den=64, 2000 steps; frames AND trajectory sampled every 4 steps — light version). No dynamics rerun.')
    open(OUT, 'w', encoding='utf-8').write(t)
    print('saved:', OUT)


if __name__ == '__main__':
    main()
