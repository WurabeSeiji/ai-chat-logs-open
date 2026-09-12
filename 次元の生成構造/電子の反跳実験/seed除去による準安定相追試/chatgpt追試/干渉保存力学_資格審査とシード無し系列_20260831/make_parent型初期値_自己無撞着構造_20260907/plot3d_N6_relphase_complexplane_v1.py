#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""N=6 相対位相複素平面図（読出し変更のみ・物理無変更）v1.1。

木原指示（2026-09-12）:
  1) 縦軸の表示サイズを4倍（scene.aspectratio z=4。データ無変更）
  2) 赤い◯（現在ステップカーソル）は廃止
  3) xy平面＝複素平面そのもの。M=15波を各ステップ15点として配置（PC集約の廃止）
  4) 同じ波どうしだけを結ぶ。15本の完全に独立した線分トレース（波をまたぐ連結は禁止）
  5) 表示位相 θ'_m(τ) = θ_m(τ) − θ0(τ)。θ0(τ)は元の赤い◯の回転位相。事前計算して配列格納。
     振幅はそのまま。表示座標 w_m(τ) = z_m(τ)·exp(−i θ0(τ))（読出しのみの回転）
  6) 軌跡は線のみ・透明度50%。現在ステップ位置だけ各波と同色のごく小さな不透明◯（v1.2）

θ0 の向き（v1 の誤りの訂正）:
  series() の PCA は決定論 svd_flip を使い、その PC2 は元図（ChatGPT基準テンプレート）の
  赤い◯に対して鏡映（正本 docstring に記録）。鏡映平面の atan2 は波と逆向きに回るため、
  v1 はθ0を引いたのに回転が2倍になった。元の赤い◯の向きは THETA0_SIGN=-1。
  客観判定を組み込む: 床（onset前）は剛体旋回だから、正しいθ0では床のθ'は静止する。
  本プログラムは床区間のθ'変動 < 1e-2 rad を検証し、破れたら停止する。

物理・データ・走行には一切触れない:
  - 状態は full_N3_N40_sweep/states/hm_N6_den_6_states_500.npz を読むだけ
  - θ0(τ)・H⊥/H は検証済み正本 plot3d_makeparent_selectedN_v1.series(6) の出力から取る
  - HTMLの器（レイアウト・スライダー・フレーム名）は基準テンプレートを流用

出力: fig3d_makeparent_selectedN/N6_relphase_complexplane_den6_step0_500_3D.html
既存 N6_makeparent_den6_step0_500_3D.html は変更しない。
"""
import base64
import copy
import importlib.util
import json
import os
import uuid

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, 'fig3d_makeparent_selectedN',
                   'N6_relphase_complexplane_den6_step0_500_3D.html')
N = 6
THETA0_SIGN = -1          # 元の赤い◯（基準テンプレート）の回転の向き
ONSET_LEVEL = 1e-3        # 床区間の判定: H⊥/H がこの値を超える直前まで
FLOOR_STATIC_TOL = 1e-2   # 床区間のθ'変動許容 [rad]（剛体旋回の静止判定）

_spec = importlib.util.spec_from_file_location(
    'orig', os.path.join(BASE, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(orig)


def b64(arr, dtype):
    return {'dtype': dtype,
            'bdata': base64.b64encode(np.asarray(arr).astype(dtype).tobytes()).decode()}


def main():
    # --- 読出し量の準備（物理式は正本から取得。再実装なし） ---
    X, Y, fz = orig.series(N)                       # PC1, PC2, H⊥/H（検証済み正本）
    theta0 = THETA0_SIGN * np.arctan2(Y, X)         # 元の赤い◯の位相。事前計算・配列格納
    d = np.load(os.path.join(orig.STATES, f'hm_N{N}_den_{N}_states_500.npz'))
    S = np.asarray(d['Z'], np.complex128)           # (501, 15) 状態そのもの
    assert S.shape[0] == len(theta0)
    W = S * np.exp(-1j * theta0)[:, None]           # 表示のみの回転。振幅不変
    T, M = W.shape
    rmax = float(np.abs(W).max()) * 1.05

    # --- 床静止の自動検証（θ0の向きが正しいことの客観判定） ---
    onset = int(np.argmax(fz > ONSET_LEVEL))
    pre = np.unwrap(np.angle(W[:max(onset - 10, 2)]), axis=0)
    floor_var = float(np.max(pre.max(0) - pre.min(0)))
    assert floor_var < FLOOR_STATIC_TOL, (
        f'床区間のθ\'変動 {floor_var:.3e} rad — θ0の向きが誤り（剛体旋回が消えていない）')

    # --- 器はテンプレート流用 ---
    shell, sp_d, sp_l, sp_f, data0, layout0, frames0 = orig.parse_template()
    layout = copy.deepcopy(layout0)

    # トレース構成: [0..14]=15本の独立ファイバー（線のみ・透明度50%）,
    #               [15]=現在位置15点（同色・ごく小・不透明）, [16]=下段曲線, [17]=縦線, [18]=stepマーカー
    palette = ['#636efa', '#EF553B', '#00cc96', '#ab63fa', '#FFA15A',
               '#19d3f3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    color = [palette[e % len(palette)] for e in range(M)]
    fiber_proto = copy.deepcopy(data0[0])
    fiber_proto.pop('customdata', None)
    fiber_proto.pop('hovertemplate', None)
    fiber_proto['mode'] = 'lines+markers'           # 全τの位置に小さな◯を残す
    fiber_proto['opacity'] = 0.5                    # 軌跡（その他）＝透明度50%
    data = []
    for e in range(M):
        tr = copy.deepcopy(fiber_proto)
        tr['line'] = {'color': color[e], 'width': 2}
        tr['marker'] = {'size': 0.5, 'color': color[e]}
        tr['x'] = b64(W[:1, e].real, 'f4')
        tr['y'] = b64(W[:1, e].imag, 'f4')
        tr['z'] = b64(fz[:1], 'f4')
        data.append(tr)
    tips = copy.deepcopy(fiber_proto)
    tips['mode'] = 'markers'
    tips['opacity'] = 1                             # ◯の位置＝透明度0%
    tips['marker'] = {'size': 0.5, 'color': color}  # 従来の1/4のごく小さな◯・各波と同色
    tips.pop('line', None)
    tips['x'] = W[0].real.tolist(); tips['y'] = W[0].imag.tolist(); tips['z'] = [float(fz[0])] * M
    data.append(tips)
    lower_curve = copy.deepcopy(data0[2])
    lower_curve['x'] = b64(np.arange(T), 'i2'); lower_curve['y'] = b64(fz, 'f8')
    vline = copy.deepcopy(data0[3])
    vline['x'] = [0, 0]; vline['y'] = [orig.FLOOR, float(fz[0])]
    smark = copy.deepcopy(data0[4])
    smark['x'] = [0]; smark['y'] = [float(fz[0])]; smark['text'] = ['step 0']
    data += [lower_curve, vline, smark]

    # scene: 複素平面（等方・対称レンジ）＋縦4倍
    layout['scene']['xaxis']['title']['text'] = 'Re w  (phase rel. to θ0)'
    layout['scene']['yaxis']['title']['text'] = 'Im w  (phase rel. to θ0)'
    layout['scene']['xaxis']['range'] = [-rmax, rmax]
    layout['scene']['yaxis']['range'] = [-rmax, rmax]
    layout['scene']['aspectmode'] = 'manual'
    layout['scene']['aspectratio'] = {'x': 1, 'y': 1, 'z': 2.5}
    layout['annotations'][0]['text'] = ('M=15 independent fibers, phase relative to collective '
                                        'rotation θ0, drawn up to current step')
    layout['annotations'][1]['text'] = 'make_parent N=6 inflation curve Hperp/H'
    if 'title' in layout:
        layout['title']['text'] = 'N=6 make_parent / den=6 / relative-phase complex plane / step 0→500'

    # frames: 各ステップkで各ファイバーをk+1点まで＋現在位置15点。フレーム名はテンプレートを踏襲
    frames = []
    for k, fr0 in enumerate(frames0):
        fd = []
        for e in range(M):
            fd.append({'type': 'scatter3d',
                       'x': b64(W[:k + 1, e].real, 'f4'),
                       'y': b64(W[:k + 1, e].imag, 'f4'),
                       'z': b64(fz[:k + 1], 'f4')})
        fd.append({'type': 'scatter3d',
                   'x': W[k].real.tolist(), 'y': W[k].imag.tolist(), 'z': [float(fz[k])] * M})
        fd.append({'type': 'scatter', 'x': [k, k], 'y': [orig.FLOOR, float(fz[k])]})
        fd.append({'type': 'scatter', 'x': [k], 'y': [float(fz[k])], 'text': [f'step {k}']})
        frames.append({'data': fd, 'name': fr0['name'],
                       'traces': list(range(M)) + [M, M + 2, M + 3]})

    # HTML組み立て（正本と同じ手順）
    t = shell
    (dj, dk), (lj, lk), (fj, fk) = sp_d, sp_l, sp_f
    enc = lambda o: json.dumps(o, separators=(',', ':'), ensure_ascii=False).replace('/', '\\u002f')
    t = t[:fj] + enc(frames) + t[fk:]
    t = t[:lj] + enc(layout) + t[lk:]
    t = t[:dj] + enc(data) + t[dk:]
    t = t.replace('c6d019db-3da4-48ad-bc8b-1d2fde4c20b1', str(uuid.uuid4()))
    title = 'N=6 make_parent / den=6 / relative-phase complex plane / step 0→500'
    t = t.replace('<title>N16 original note source</title>', f'<title>{title}</title>')
    t = t.replace('<h1>N=16 original note source / den=16 / step 0→500</h1>', f'<h1>{title}</h1>')
    t = t.replace('Source: hm_N16_den_16_states_500.npz materialized directly from Google Drive. No dynamics rerun.',
                  'Source: full_N3_N40_sweep/states/hm_N6_den_6_states_500.npz. Display-only transform: '
                  'w_m(t)=z_m(t)*exp(-i*theta0(t)), theta0 = collective rotation phase of the original '
                  'red-marker (floor-static verified). No dynamics rerun, no physics change.')
    open(OUT, 'w', encoding='utf-8').write(t)
    print(f'床静止検証: onset=step{onset}, 床区間θ\'変動 {floor_var:.3e} rad < {FLOOR_STATIC_TOL}')
    print('wrote', os.path.relpath(OUT, BASE), f'({os.path.getsize(OUT)/1e6:.1f} MB)')


if __name__ == '__main__':
    main()
