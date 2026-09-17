#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""位相・振幅時系列の対話ビューア HTML（2026-09-17 木原指示 v6・動的ロード方式）。

v5 までの「実行のたびに存在する全 run を検出し run ごとの巨大 HTML を再生成」を
廃止（木原指示）。構成を次の2部に分離する:

 1. ビューア（単一 HTML、run データ非埋め込み・plotly.js 埋め込み）
    plots/phase_amplitude/phase_amplitude_viewer_dense_rerun_20260917.html
    - 99走行ドロップダウンから run を選択すると、data/data_<rid>.js を
      <script> タグ注入で**動的ロード**して描画（file:// でも動作）
    - データが無い run を選ぶと**エラーメッセージ**（未走行または未生成）を表示
    - 機能は v5 と同一: 実線=位相（左軸 −π..π）・点列＋半透明点線=正規化振幅
      （右軸）、凡例クリックのモーダル波選択（保持・再クリック解除）、
      ドラッグ横スクロール、ズーム率 x1..x400（位置保持、初期 x400・step 0）、
      位相/振幅 表示ON/OFF（初期 両ON）
 2. run データ生成（本スクリプトの毎回の仕事）
    plots/phase_amplitude/data/data_<rid>.js … 存在する run について
    位相・正規化振幅を全 4097 step・全辺ぶん JSON で書き出す（丸めなし・
    repr 最短往復精度）。states.npz より新しい既存データはスキップ（増分生成）。

v7（2026-09-17 木原指示）: 位相の振幅カットオフ選択を追加。正規化振幅が
指定値未満の点では位相を非表示（null で線を途切れさせる）。既定=カットオフ無し、
選択肢 0.1 / 0.2 / 0.3 / 0.5 / 0.7 / 1.0。振幅表示・hover・選択には影響しない。
run 切替時も選択中のカットオフを維持適用。

データは states.npz の read-only 参照のみ。

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

VIEWER_JS = """
var gd = document.getElementById('{plot_id}');
var TMAX = __TMAX__;
var COLORS = __COLORS_JSON__;      // M（辺数）→ turbo 色配列
var INITIAL = '__INITIAL__';
var current = null;                // 表示中の run_id
var curData = null;                // 表示中の run データ（カットオフ再適用用）
var selected = null;               // モーダル選択中の波（meta 値）。null = 非選択
var phaseOn = true, ampOn = true;  // 表示モード（初期値: 両方 ON）
var cutoff = null;                 // 位相の振幅カットオフ（正規化値）。null = カットオフ無し
function applyCutoff(){
  // 正規化振幅 < cutoff の点は位相を null にして非表示（線は途切れる）。
  // 対象は位相トレース（各辺の 3e 番目）のみ。振幅表示には影響しない。
  if(!curData){ return; }
  var M = curData.edges.length;
  var idx = [], ys = [];
  for(var e = 0; e < M; e++){
    var ph = curData.phase[e], an = curData.ampn[e];
    var y;
    if(cutoff === null){ y = ph; }
    else {
      y = new Array(ph.length);
      for(var i = 0; i < ph.length; i++){ y[i] = (an[i] >= cutoff) ? ph[i] : null; }
    }
    ys.push(y); idx.push(3 * e);
  }
  Plotly.restyle(gd, {'y': ys}, idx);
}
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
    } else {
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
    if(c === 'amp'){ return phaseOn ? false : ampOn; }
    return false;
  });
  Plotly.restyle(gd, {'visible': vis, 'showlegend': sl});
}
function buildAndRender(d){
  current = d.rid; curData = d; selected = null;
  var M = d.edges.length;
  var cols = COLORS[String(M)];
  var tt = new Array(d.T + 1);
  for(var i = 0; i <= d.T; i++){ tt[i] = i; }
  var traces = [];
  for(var e = 0; e < M; e++){
    var name = '(' + d.edges[e][0] + ',' + d.edges[e][1] + ')';
    var ph = d.phase[e], an = d.ampn[e];
    var cd = new Array(ph.length);
    for(var i2 = 0; i2 < ph.length; i2++){
      cd[i2] = [ph[i2] / Math.PI, an[i2] * d.amax, an[i2]];
    }
    traces.push({type: 'scattergl', x: tt, y: ph, mode: 'lines',
      line: {color: cols[e], width: 1.4}, opacity: 0.85,
      name: name, legendgroup: 'e' + e, meta: e, customdata: cd,
      hovertemplate: '辺 ' + name + '　位相<br>t=%{x}<br>' +
        'arg z=%{y:.4f} rad（%{customdata[0]:.3f}π）<extra></extra>'});
    traces.push({type: 'scattergl', x: tt, y: an, mode: 'markers', yaxis: 'y2',
      marker: {color: cols[e], size: 3.2}, opacity: 0.85,
      name: name + ' |z|', legendgroup: 'e' + e, showlegend: false, meta: e,
      customdata: cd,
      hovertemplate: '辺 ' + name + '　振幅<br>t=%{x}<br>' +
        '|z|=%{customdata[1]:.4f}（正規化 %{customdata[2]:.4f}）<extra></extra>'});
    traces.push({type: 'scatter', x: tt, y: an, mode: 'lines', yaxis: 'y2',
      line: {color: cols[e], width: 1.0, dash: 'dot'}, opacity: 0.35,
      name: name + ' |z| 点線', legendgroup: 'e' + e, showlegend: false,
      meta: e, hoverinfo: 'skip'});
  }
  Plotly.react(gd, traces, gd.layout);
  Plotly.relayout(gd, {'title.text':
    d.rid + ' — 位相・振幅の時系列 t=0..' + d.T + '（全ステップ）<br>' +
    '<sup>実線=位相（左軸）・点列＋半透明点線=正規化振幅（右軸）。' +
    '凡例クリック=波を選択保持（再クリックで解除）。ドラッグ=横スクロール、' +
    'ボタン=ズーム率（位置保持）・位相/振幅 表示ON/OFF。' +
    '振幅正規化: 全体最大 ' + d.amax.toFixed(4) + '</sup>'});
  applyVis();
  apply();
  applyCutoff();   // 選択中のカットオフを新しい run にも適用
}
function loadRun(rid){
  if(rid === current){ return; }
  window.__RUNDATA = null;
  var s = document.createElement('script');
  s.src = 'data/data_' + rid + '.js';
  s.onload = function(){
    if(window.__RUNDATA && window.__RUNDATA.rid === rid){
      buildAndRender(window.__RUNDATA);
    } else {
      alert(rid + ' のデータ形式が不正です');
    }
  };
  s.onerror = function(){
    alert(rid + ' のデータがありません（未走行、またはデータ未生成）。\\n' +
          'plot_phase_amplitude_timeseries_html_20260917.py を実行すると生成されます。');
  };
  document.head.appendChild(s);
}
gd.on('plotly_legendclick', function(d){
  var m = gd.data[d.curveNumber].meta;
  selected = (selected === m) ? null : m;
  apply();
  return false;
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
    Plotly.relayout(gd, {'xaxis.range': [left, left + width]});
  } else if(lbl.indexOf('位相') === 0){
    phaseOn = !phaseOn; applyVis();
  } else if(lbl.indexOf('振幅') === 0){
    ampOn = !ampOn; applyVis();
  } else if(lbl.indexOf('カットオフ') === 0){
    var v = lbl.replace('カットオフ', '').trim();
    cutoff = (v === '無し') ? null : parseFloat(v);
    applyCutoff();
  } else if(lbl.charAt(0) === 'L'){
    loadRun(lbl);
  }
});
loadRun(INITIAL);
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


def write_run_data(key, datadir):
    """存在する run の位相・正規化振幅データ JS を書き出す（増分・丸めなし）。"""
    L, ma, mb, den = key
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    src = os.path.join(SURVEY, 'runs', rid, 'states.npz')
    dst = os.path.join(datadir, f'data_{rid}.js')
    if os.path.exists(dst) and os.path.getmtime(dst) > os.path.getmtime(src):
        return rid, os.path.getsize(dst), 'skip（states.npz より新しい）'
    with np.load(src) as d:
        Z = np.array(d['Z'])
        assert Z.shape[0] == T + 1
    ja, jb = np.triu_indices(L, k=1)
    amax_run = float(np.max(np.abs(Z)))
    payload = {
        'rid': rid, 'L': L, 'T': T, 'amax': amax_run,
        'edges': [[int(ja[e]), int(jb[e])] for e in range(Z.shape[1])],
        'phase': [np.angle(Z[:, e]).tolist() for e in range(Z.shape[1])],
        'ampn': [(np.abs(Z[:, e]) / amax_run).tolist() for e in range(Z.shape[1])],
    }
    with open(dst, 'w', encoding='utf-8') as f:
        f.write('window.__RUNDATA = ')
        json.dump(payload, f, ensure_ascii=False)
        f.write(';\n')
    return rid, os.path.getsize(dst), 'generated'


def write_viewer(outdir, rid_all, initial_rid):
    """run データ非埋め込みのビューア HTML を書き出す。"""
    colors = {}
    for L in sorted({k[0] for k in all_keys()}):
        M = L * (L - 1) // 2
        colors[str(M)] = [mcolors.to_hex(cm.turbo(x))
                          for x in np.linspace(0.02, 0.98, M)]
    pi = float(np.pi)
    width0 = T / ZOOM_INIT
    fig = go.Figure()
    fig.update_layout(
        title=dict(text='run を選択してください（右上のドロップダウン）', x=0.5),
        template='plotly_white', width=1180, height=780,
        dragmode='pan',
        xaxis=dict(title='step', range=[0, width0], zeroline=False),
        yaxis=dict(title='位相 arg z_e(t)（実線）',
                   range=[-pi * 1.05, pi * 1.05], fixedrange=True,
                   tickvals=[-pi, -pi / 2, 0, pi / 2, pi],
                   ticktext=['−π', '−π/2', '0', 'π/2', 'π'],
                   zeroline=True, zerolinecolor='lightgray'),
        yaxis2=dict(title='振幅 |z_e(t)| / 全体最大（点列＋半透明点線）',
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
            dict(type='dropdown', direction='down', x=0.90, y=1.10,
                 xanchor='left', active=0,
                 buttons=[dict(label='カットオフ無し', method='skip')] +
                         [dict(label=f'カットオフ {c}', method='skip')
                          for c in ('0.1', '0.2', '0.3', '0.5', '0.7', '1.0')]),
            dict(type='dropdown', direction='down', x=0.44, y=1.10,
                 xanchor='left', active=rid_all.index(initial_rid),
                 buttons=[dict(label=r, method='skip') for r in rid_all]),
        ],
        hovermode='closest')
    js = (VIEWER_JS
          .replace('__TMAX__', str(T))
          .replace('__COLORS_JSON__', json.dumps(colors))
          .replace('__INITIAL__', initial_rid))
    out = os.path.join(outdir, 'phase_amplitude_viewer_dense_rerun_20260917.html')
    fig.write_html(out, include_plotlyjs=True, full_html=True, post_script=js,
                   config={'scrollZoom': False, 'displaylogo': False})
    return out


def main():
    keys = all_keys()
    rid_all = [f'L{L}_ma{ma}_mb{mb}_den{den}' for (L, ma, mb, den) in keys]
    outdir = os.path.join(SURVEY, 'plots', 'phase_amplitude')
    datadir = os.path.join(outdir, 'data')
    os.makedirs(datadir, exist_ok=True)

    existing = [(k, r) for k, r in zip(keys, rid_all)
                if have(k[0], k[1], k[2], k[3])]
    for k, r in existing:
        rid, size, status = write_run_data(k, datadir)
        print(f'data {rid}: {status} ({size / 1e6:.1f} MB)')
    if not existing:
        print('states.npz が存在する run が無い — データ生成なし（ビューアのみ更新）')
    initial = existing[0][1] if existing else rid_all[0]
    out = write_viewer(outdir, rid_all, initial)
    print(f'viewer -> {os.path.basename(out)} ({os.path.getsize(out) / 1e6:.1f} MB, '
          f'初期表示 run: {initial}, データ生成済み {len(existing)}/99)')


if __name__ == '__main__':
    main()
