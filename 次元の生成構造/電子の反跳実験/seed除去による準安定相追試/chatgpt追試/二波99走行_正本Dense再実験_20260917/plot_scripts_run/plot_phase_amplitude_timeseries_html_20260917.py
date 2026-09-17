#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""位相・振幅時系列の対話HTML・1 run = 1 図（2026-09-17 木原指示 v3）。

左軸=位相 arg z_e(t) ∈ (−π, π]（実線）、右軸=正規化振幅 |z|/max（点マーカー列
＝点線表現）、波毎1色（turbo 28色）。木原指示による機能:

 v2: 右側の凡例メニューからのみ選択するモーダル・ハイライト
   - 凡例クリック = その波の位相・振幅ペアを強調、他を薄く。選択は別の波を
     選ぶまで保持、同じ波を再クリックで解除。凡例の既定動作は無効化。
 v3: 全 4096 step 対応＋横スクロール＋ズーム率
   - データは t=0..4096 の全 4097 状態（WebGL Scattergl、間引きなし）
   - 横軸はドラッグでスクロール（パン。縦軸は固定）
   - ズーム率ボタン x1(全4096) / x10 / x50 / x100 / x200 / x400（窓幅=4096/率）。
     切替時は現在のスクロール位置（左端）を保持。初期値 x400・位置 step 0。
 v4: 振幅の点列を同色・半透明の点線（SVG dash='dot'、hover なし）でも接続。
     表示モードボタン「位相 ON/OFF」「振幅 ON/OFF」を追加（初期値: 両方 ON。
     位相 OFF 中は凡例を振幅側に付け替えて選択 UI を維持）。
 v5: 99走行パターンの選択メニュー（ドロップダウン）を追加。1 run = 1 HTML の
     構成のまま、選択すると同フォルダ内の該当 run の HTML へ遷移する
     （全99を1ファイルに埋め込むと数百MBになるため）。未走行 run は
     「（未走行）」表示で、選択すると通知のみ。

データは states.npz の read-only 参照のみ。1 run = 1 HTML（存在する run だけ）。

使い方: python3 plot_phase_amplitude_timeseries_html_20260917.py
"""
import json
import os
import sys

import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import plotly.graph_objects as go

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, SURVEY)
from gen_manifest import MAIN_L, MAIN_PAIRS, GCD_L, GCD_PAIRS, CTRL_L, CTRL_PAIRS, CTRL_DEN, T  # noqa: E402

ZOOMS = [1, 10, 50, 100, 200, 400]
ZOOM_INIT = 400            # 初期ズーム率（窓幅 4096/400 ≈ 10.24 step、位置 step 0）

SELECT_ZOOM_JS = """
var gd = document.getElementById('{plot_id}');
var TMAX = %(tmax)d;
var RUN_EXISTS = %(run_exists_json)s;   // 99走行 rid → states.npz 生成済みか
var CURRENT_RID = '%(current_rid)s';
var selected = null;               // モーダル選択中の波（meta 値）。null = 非選択
var phaseOn = true, ampOn = true;  // 表示モード（初期値: 両方 ON）
function cls(t){
  if(t.yaxis === 'y2'){ return (t.mode === 'lines') ? 'conn' : 'amp'; }
  return 'phase';
}
function apply(){
  var op = [], w = [], ms = [];
  gd.data.forEach(function(t){
    var c = cls(t);
    var sel = (selected !== null && t.meta === selected);
    var dim = (selected !== null && t.meta !== selected);
    if(c === 'phase'){
      op.push(dim ? 0.10 : (sel ? 1.0 : 0.85));
      w.push(sel ? 3.0 : (dim ? 0.8 : 1.4));
      ms.push(3.2);
    } else if(c === 'amp'){
      op.push(dim ? 0.10 : (sel ? 1.0 : 0.85));
      w.push(1.0);
      ms.push(sel ? 5.5 : (dim ? 2.0 : 3.2));
    } else {  // conn: 振幅を繋ぐ半透明点線
      op.push(dim ? 0.05 : (sel ? 0.8 : 0.35));
      w.push(sel ? 2.2 : 1.0);
      ms.push(3.2);
    }
  });
  Plotly.restyle(gd, {'opacity': op, 'line.width': w, 'marker.size': ms});
}
function applyVis(){
  var vis = gd.data.map(function(t){
    return (cls(t) === 'phase') ? phaseOn : ampOn;
  });
  var sl = gd.data.map(function(t){
    var c = cls(t);
    if(c === 'phase'){ return phaseOn; }
    if(c === 'amp'){ return phaseOn ? false : ampOn; }  // 位相OFF時は凡例を振幅側に
    return false;
  });
  Plotly.restyle(gd, {'visible': vis, 'showlegend': sl});
}
gd.on('plotly_legendclick', function(d){
  var m = gd.data[d.curveNumber].meta;
  selected = (selected === m) ? null : m;   // 同じ波の再クリックで解除
  apply();
  return false;   // 凡例の既定動作（表示切替）を無効化
});
gd.on('plotly_legenddoubleclick', function(d){ return false; });
gd.on('plotly_buttonclicked', function(d){
  var lbl = d.button.label;
  if(lbl.charAt(0) === 'x'){
    var zoom = parseInt(lbl.slice(1));
    var width = TMAX / zoom;
    var left = gd.layout.xaxis.range ? gd.layout.xaxis.range[0] : 0;
    if(left < 0){ left = 0; }
    if(left > TMAX - width){ left = Math.max(0, TMAX - width); }
    Plotly.relayout(gd, {'xaxis.range': [left, left + width]});  // 位置保持で率のみ変更
  } else if(lbl.indexOf('位相') === 0){
    phaseOn = !phaseOn; applyVis();
  } else if(lbl.indexOf('振幅') === 0){
    ampOn = !ampOn; applyVis();
  } else if(lbl.charAt(0) === 'L'){
    var rid = lbl.replace('（未走行）', '').trim();
    if(!RUN_EXISTS[rid]){
      alert(rid + ' はまだ走行していません（states.npz なし）');
      return;
    }
    if(rid !== CURRENT_RID){
      window.location.href = 'phase_amplitude_full_' + rid + '_dense_rerun_20260917.html';
    }
  }
});
"""


def have(L, ma, mb, den):
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    return os.path.exists(os.path.join(SURVEY, 'runs', rid, 'states.npz'))


def all_keys():
    keys = []
    for L in MAIN_L:
        for ma, mb in MAIN_PAIRS:
            keys.append((L, ma, mb, L))
    for L in GCD_L:
        for ma, mb in GCD_PAIRS:
            keys.append((L, ma, mb, L))
    for L in CTRL_L:
        for ma, mb in CTRL_PAIRS:
            keys.append((L, ma, mb, CTRL_DEN))
    return keys


def main():
    keys = all_keys()
    rid_all = [f'L{L}_ma{ma}_mb{mb}_den{den}' for (L, ma, mb, den) in keys]
    exists = {r: have(L, ma, mb, den)
              for r, (L, ma, mb, den) in zip(rid_all, keys)}
    run_buttons = [dict(label=(r if exists[r] else r + '（未走行）'), method='skip')
                   for r in rid_all]
    found = 0
    for (L, ma, mb, den) in all_keys():
        if not have(L, ma, mb, den):
            continue
        found += 1
        rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
        with np.load(os.path.join(SURVEY, 'runs', rid, 'states.npz')) as d:
            Z = np.array(d['Z'])
            assert Z.shape[0] == T + 1
        M = Z.shape[1]
        ja, jb = np.triu_indices(L, k=1)
        amax_run = float(np.max(np.abs(Z)))          # 全体最大振幅（全4097ステップ）
        tt = np.arange(T + 1)
        cols = [mcolors.to_hex(cm.turbo(x)) for x in np.linspace(0.02, 0.98, M)]

        fig = go.Figure()
        for e in range(M):
            name = f'({ja[e]},{jb[e]})'
            ph = np.angle(Z[:, e])
            amp = np.abs(Z[:, e])
            cd = np.stack([ph / np.pi, amp, amp / amax_run], axis=-1)
            fig.add_trace(go.Scattergl(
                x=tt, y=ph, mode='lines',
                line=dict(color=cols[e], width=1.4),
                opacity=0.85, name=name, legendgroup=f'e{e}', meta=e,
                customdata=cd,
                hovertemplate=(f'辺 {name}　位相<br>t=%{{x}}<br>'
                               'arg z=%{y:.4f} rad（%{customdata[0]:.3f}π）'
                               '<extra></extra>')))
            fig.add_trace(go.Scattergl(
                x=tt, y=amp / amax_run, mode='markers', yaxis='y2',
                marker=dict(color=cols[e], size=3.2),
                opacity=0.85, name=name + ' |z|', legendgroup=f'e{e}',
                showlegend=False, meta=e, customdata=cd,
                hovertemplate=(f'辺 {name}　振幅<br>t=%{{x}}<br>'
                               '|z|=%{customdata[1]:.4f}（正規化 %{customdata[2]:.4f}）'
                               '<extra></extra>')))
            fig.add_trace(go.Scatter(          # 振幅を繋ぐ半透明の点線（SVG、hover なし）
                x=tt, y=amp / amax_run, mode='lines', yaxis='y2',
                line=dict(color=cols[e], width=1.0, dash='dot'),
                opacity=0.35, name=name + ' |z| 点線', legendgroup=f'e{e}',
                showlegend=False, meta=e, hoverinfo='skip'))

        pi = float(np.pi)
        width0 = T / ZOOM_INIT
        fig.update_layout(
            title=dict(text=(f'{rid} — 位相・振幅の時系列 t=0..{T}（全ステップ）<br>'
                             '<sup>実線=位相（左軸）・点列＋半透明点線=正規化振幅（右軸）。'
                             '凡例クリック=波を選択保持（再クリックで解除）。'
                             'ドラッグ=横スクロール、ボタン=ズーム率（窓幅 4096/率、位置保持）'
                             f'・位相/振幅の表示ON/OFF（初期 両方ON）。初期 x{ZOOM_INIT}・'
                             f'step 0 起点。振幅正規化: 全体最大 {amax_run:.4f}</sup>'), x=0.5),
            template='plotly_white', width=1180, height=780,
            dragmode='pan',
            xaxis=dict(title='step', range=[0, width0], zeroline=False),
            yaxis=dict(title='位相 arg z_e(t)（実線）',
                       range=[-pi * 1.05, pi * 1.05], fixedrange=True,
                       tickvals=[-pi, -pi / 2, 0, pi / 2, pi],
                       ticktext=['−π', '−π/2', '0', 'π/2', 'π'],
                       zeroline=True, zerolinecolor='lightgray'),
            yaxis2=dict(title=f'振幅 |z_e(t)| / {amax_run:.4f}（点列）',
                        overlaying='y', side='right', range=[0, 1.05],
                        fixedrange=True, showgrid=False),
            legend=dict(font=dict(size=9), itemsizing='constant',
                        title=dict(text='辺 (j,k)')),
            updatemenus=[
                dict(type='buttons', direction='left', x=0.0, y=1.10, xanchor='left',
                     active=ZOOMS.index(ZOOM_INIT),
                     buttons=[dict(label=f'x{z}', method='skip') for z in ZOOMS]),
                dict(type='buttons', direction='left', x=0.72, y=1.10, xanchor='left',
                     buttons=[dict(label='位相 ON/OFF', method='skip'),
                              dict(label='振幅 ON/OFF', method='skip')]),
                dict(type='dropdown', direction='down', x=0.44, y=1.10,
                     xanchor='left', active=rid_all.index(rid),
                     buttons=run_buttons),
            ],
            hovermode='closest')

        outdir = os.path.join(SURVEY, 'plots', 'phase_amplitude')
        os.makedirs(outdir, exist_ok=True)
        out = os.path.join(outdir,
                           f'phase_amplitude_full_{rid}_dense_rerun_20260917.html')
        fig.write_html(out, include_plotlyjs=True, full_html=True,
                       post_script=SELECT_ZOOM_JS % {
                           'tmax': T,
                           'run_exists_json': json.dumps(exists),
                           'current_rid': rid},
                       config={'scrollZoom': False, 'displaylogo': False})
        print(f'{rid}: 全体最大振幅 {amax_run:.6f} -> {os.path.basename(out)} '
              f'({os.path.getsize(out) / 1e6:.1f} MB)')
    if found == 0:
        print('states.npz が存在する run が無い — 図化せず終了')
    else:
        print(f'存在する run: {found}/99')


if __name__ == '__main__':
    main()
