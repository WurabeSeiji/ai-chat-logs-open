#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""相対位相複素平面図（読出し変更のみ・物理無変更）v2 — データ1回埋め込み方式。

v1 の問題: plotly のアニメーションフレームは各フレームが step0..k の全点を再収録するため、
データが Σ(k+1)≈251倍に重複し、N=40 で約1.6GB になった（ブラウザ・GitHub とも不可）。

v2 の方式（木原指示 2026-09-12「最初からデータはすべて渡して表示ON/OFF制御だけ」）:
  - 状態データ W（M波×T step）と H⊥/H を base64(float32) で HTML に**1回だけ**埋め込む。
  - plotly の frames は使わない。自前の range スライダー＋Play が、ブラウザ内で
    埋め込み配列を step0..k に切り出し Plotly.restyle するだけ（表示範囲のON/OFF制御）。
  - データ重複ゼロ。サイズ ≈ M·T·2·4B·(4/3)。N=6≈0.1MB、N=40≈5MB。

図の仕様（v1 で確定・不変）:
  1) θ'_m(τ)=θ_m(τ)−θ0(τ)、θ0=元の赤丸の集団回転位相（THETA0_SIGN=-1、床静止で自動検証）
  2) 表示座標 w_m(τ)=z_m(τ)·exp(−iθ0(τ))。振幅そのまま。物理・データ・走行は無変更
  3) M本の独立ファイバー（波をまたぐ連結なし）。線+全τの小◯（size0.5・各波同色・opacity0.5）
  4) 現在ステップ位置＝各波同色の不透明小◯（size0.5・opacity1）
  5) 縦軸 aspectratio z=2.5。下段に H⊥/H インフレーション曲線＋現在step縦線・マーカー
  θ0・H⊥/H・PC は検証済み正本 plot3d_makeparent_selectedN_v1.series(N) から取得（式の再実装なし）。

使い方: python3 plot3d_relphase_complexplane_v2.py 6      # N=6
        python3 plot3d_relphase_complexplane_v2.py 40     # N=40
出力: fig3d_makeparent_selectedN/N{N}_relphase_complexplane_den{N}_step0_500_3D.html
"""
import base64
import copy
import importlib.util
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
ONSET_LEVEL = 1e-3
# 床静止判定（θ0の向きの客観決定）: 正しい向きは剛体回転を除去して床が静止する。
#  ・絶対基準 FLOOR_STATIC_ABS: 床区間のθ'変動がこの値未満＝正味回転なし（make_parent近似床の
#    残差ゆらぎ数度は許容。誤った向きは同区間で数十rad＝多回転になる）
#  ・比基準 FLOOR_SIGN_RATIO: 逆向きが選んだ向きよりこの倍率以上悪い＝符号が一意に確定
FLOOR_STATIC_ABS = 0.2
FLOOR_SIGN_RATIO = 20.0
PALETTE = ['#636efa', '#EF553B', '#00cc96', '#ab63fa', '#FFA15A',
           '#19d3f3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']

_spec = importlib.util.spec_from_file_location(
    'orig', os.path.join(BASE, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(orig)


def b64f32(arr):
    return base64.b64encode(np.asarray(arr, np.float32).tobytes()).decode()


def build(N):
    # --- 読出し量（物理式は正本から取得・再実装なし） ---
    X, Y, fz = orig.series(N)
    d = np.load(os.path.join(orig.STATES, f'hm_N{N}_den_{N}_states_500.npz'))
    S = np.asarray(d['Z'], np.complex128)
    assert S.shape[0] == len(X)
    onset = int(np.argmax(fz > ONSET_LEVEL))

    # --- θ0の向きは床静止で客観決定（PCA軸向きの任意性のみ。物理ではない） ---
    # 床は相対平衡＝剛体旋回なので、正しいθ0では床区間のθ'が静止する。
    # 二つの向き ±atan2(PC2,PC1) を床静止基準で採点し、静止する方を採る。
    def floor_var_for(sign):
        th0 = sign * np.arctan2(Y, X)
        Wt = S * np.exp(-1j * th0)[:, None]
        pre = np.unwrap(np.angle(Wt[:max(onset - 10, 2)]), axis=0)
        return float(np.max(pre.max(0) - pre.min(0))), th0, Wt
    cand = {s: floor_var_for(s) for s in (+1, -1)}
    sign = min(cand, key=lambda s: cand[s][0])
    floor_var, theta0, W = cand[sign]
    other = max(cand[s][0] for s in cand)
    assert floor_var < FLOOR_STATIC_ABS, (
        f'選んだ向きでも床が静止しない（{floor_var:.3e} rad ≥ {FLOOR_STATIC_ABS}）'
        f'— 床が相対平衡でない可能性')
    assert other > FLOOR_SIGN_RATIO * floor_var, (
        f'両向きの床変動が拮抗（{floor_var:.3e} vs {other:.3e} rad）— θ0の向きが一意に定まらない')
    T, M = W.shape
    rmax = float(np.abs(W).max()) * 1.05

    color = [PALETTE[e % len(PALETTE)] for e in range(M)]
    FLOOR = orig.FLOOR

    # --- レイアウトはテンプレート流用（frames/slider/updatemenus は使わないので除去） ---
    _, _, _, _, _, layout0, _ = orig.parse_template()
    layout = copy.deepcopy(layout0)
    layout.pop('sliders', None)
    layout.pop('updatemenus', None)
    layout['scene']['xaxis']['title']['text'] = 'Re w  (phase rel. to θ0)'
    layout['scene']['yaxis']['title']['text'] = 'Im w  (phase rel. to θ0)'
    layout['scene']['xaxis']['range'] = [-rmax, rmax]
    layout['scene']['yaxis']['range'] = [-rmax, rmax]
    layout['scene']['aspectmode'] = 'manual'
    layout['scene']['aspectratio'] = {'x': 1, 'y': 1, 'z': 2.5}
    layout['annotations'][0]['text'] = (f'M={M} independent fibers, phase relative to collective '
                                        'rotation θ0, drawn up to current step')
    layout['annotations'][1]['text'] = f'make_parent N={N} inflation curve Hperp/H'
    layout.pop('width', None)                # スマホ対応: 固定幅を除きレスポンシブに
    layout.pop('height', None)
    layout['autosize'] = True
    if 'title' in layout:
        layout['title']['text'] = f'N={N} make_parent / den={N} / relative-phase complex plane'

    # --- 初期表示は最終step。トレースは空で作り、読込直後の setStep(k0) が base64 から充填 ---
    #     （軌跡データを JSON と base64 に二重保持しない＝N=40 のサイズ半減） ---
    k0 = T - 1
    data = []
    for e in range(M):                       # [0..M-1] 15/780本の独立ファイバー
        data.append({'type': 'scatter3d', 'mode': 'lines+markers', 'opacity': 0.5,
                     'line': {'color': color[e], 'width': 2},
                     'marker': {'size': 0.5, 'color': color[e]},
                     'x': [], 'y': [], 'z': [], 'showlegend': False,
                     'hovertemplate': f'wave {e}<extra></extra>'})
    data.append({'type': 'scatter3d', 'mode': 'markers', 'opacity': 1,   # [M] 現在位置
                 'marker': {'size': 0.5, 'color': color},
                 'x': [], 'y': [], 'z': [], 'showlegend': False, 'hoverinfo': 'skip'})
    data.append({'type': 'scatter', 'mode': 'lines', 'xaxis': 'x', 'yaxis': 'y',  # [M+1] 下段曲線
                 'line': {'color': '#2ca02c'},
                 'x': list(range(T)), 'y': fz.tolist(), 'showlegend': False, 'hoverinfo': 'skip'})
    data.append({'type': 'scatter', 'mode': 'lines', 'xaxis': 'x', 'yaxis': 'y',  # [M+2] 現在step縦線
                 'line': {'color': '#888', 'dash': 'dot'},
                 'x': [k0, k0], 'y': [FLOOR, float(fz[k0])], 'showlegend': False, 'hoverinfo': 'skip'})
    data.append({'type': 'scatter', 'mode': 'markers+text', 'xaxis': 'x', 'yaxis': 'y',  # [M+3]
                 'marker': {'size': 9, 'color': '#ff7f0e'},
                 'x': [k0], 'y': [float(fz[k0])], 'text': [f'step {k0}'],
                 'textposition': 'top center', 'showlegend': False, 'hoverinfo': 'skip'})

    title = f'N={N} make_parent / den={N} / relative-phase complex plane / step 0→500'
    fig = {'data': data, 'layout': layout, 'config': {'responsive': True}}
    figjson = json.dumps(fig, separators=(',', ':'), ensure_ascii=False)

    # 埋め込みデータ（1回だけ・float32 base64）
    WR = b64f32(np.ascontiguousarray(W.real.T))   # row-major [e*T + t]
    WI = b64f32(np.ascontiguousarray(W.imag.T))
    FZb = b64f32(fz)

    # plotly本体をインライン埋め込み（スマホ・オフラインで開けるよう自己完結化）
    plotly_src = open(os.path.join(BASE, '_lib_plotly-3.3.1.min.js'), encoding='utf-8').read()

    html = f"""<!doctype html><html lang="ja"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>body{{margin:0;font:14px system-ui;-webkit-text-size-adjust:100%}}
h1{{font-size:16px;margin:8px 12px}} p.src{{font-size:11px;margin:0 12px;color:#666}}
#ctl{{padding:8px 12px;display:flex;gap:10px;align-items:center;flex-wrap:wrap}}
#sld{{flex:1;min-width:160px}} button{{font-size:15px;padding:6px 14px}}
#gd{{width:100%;height:78vh;min-height:420px}}</style></head><body>
<h1>{title}</h1>
<p class="src">hm_N{N}_den_{N}_states_500.npz / display-only w=z·exp(-iθ0), floor-static {floor_var:.1e} rad. No physics change.</p>
<div id="ctl"><button id="play">▶ Play</button><input id="sld" type="range" min="0" max="{T-1}" value="{k0}"><span id="lbl">step {k0}</span></div>
<div id="gd"></div>
<script>{plotly_src}</script>
<script>
const M={M}, T={T}, FLOOR={FLOOR};
function d64(s){{const b=atob(s),u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);return new Float32Array(u.buffer);}}
const WR=d64("{WR}"), WI=d64("{WI}"), FZ=d64("{FZb}");
const gd=document.getElementById('gd');
Plotly.newPlot(gd, {figjson}).then(()=>{{ setStep({k0}); }});
function setStep(k){{
  const xs=[],ys=[],zs=[];
  for(let e=0;e<M;e++){{xs.push(WR.subarray(e*T,e*T+k+1));ys.push(WI.subarray(e*T,e*T+k+1));zs.push(FZ.subarray(0,k+1));}}
  const idx=[];for(let e=0;e<M;e++)idx.push(e);
  Plotly.restyle(gd,{{x:xs,y:ys,z:zs}},idx);
  const tx=new Float32Array(M),ty=new Float32Array(M),tz=new Float32Array(M);
  for(let e=0;e<M;e++){{tx[e]=WR[e*T+k];ty[e]=WI[e*T+k];tz[e]=FZ[k];}}
  Plotly.restyle(gd,{{x:[tx],y:[ty],z:[tz]}},[M]);
  Plotly.restyle(gd,{{x:[[k,k]],y:[[FLOOR,FZ[k]]]}},[M+2]);
  Plotly.restyle(gd,{{x:[[k]],y:[[FZ[k]]],text:[['step '+k]]}},[M+3]);
  document.getElementById('lbl').textContent='step '+k;
}}
const sld=document.getElementById('sld'), btn=document.getElementById('play');
sld.addEventListener('input',()=>setStep(+sld.value));
let timer=null;
btn.addEventListener('click',()=>{{
  if(timer){{clearInterval(timer);timer=null;btn.textContent='▶ Play';return;}}
  btn.textContent='⏸ Pause';
  timer=setInterval(()=>{{let k=(+sld.value+1);if(k>=T){{clearInterval(timer);timer=null;btn.textContent='▶ Play';return;}}sld.value=k;setStep(k);}},50);
}});
</script></body></html>"""
    out = os.path.join(BASE, 'fig3d_makeparent_selectedN',
                       f'N{N}_relphase_complexplane_den{N}_step0_500_3D.html')
    open(out, 'w', encoding='utf-8').write(html)
    print(f'N={N}: θ0向き={sign:+d}(床静止で自動決定) 床変動{floor_var:.3e}rad(逆向き{other:.3e}) '
          f'onset=step{onset} / M={M} T={T} / {os.path.getsize(out)/1e6:.2f} MB / {os.path.relpath(out, BASE)}')
    return out, W, fz


if __name__ == '__main__':
    build(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
