#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fine 系列（den=200L・T=819200）用 位相・振幅対話ビューア（2026-09-17 木原指示）。

粗い系列用ビューア v7（plot_phase_amplitude_timeseries_html_20260917.py）の
チャンク動的ロード改造版。fine 系列は 1 run の表示データが約 870MB になり
一括ロード不能のため、時間窓 CHUNK=8192 step ごとに分割して書き出し、
選択中の窓だけを <script> 注入で動的ロードする。

機能（v7 と同一思想）:
 - 実線=位相（左軸 −π..π）・点列＋半透明点線=正規化振幅（右軸、run 全体の
   実測最大で正規化）・波毎1色
 - 凡例クリックのモーダル波選択（保持・再クリック解除）
 - ドラッグ横スクロール（窓内）、ズーム率 x1..x400（窓幅=CHUNK/率、位置保持）
 - 位相/振幅 表示 ON/OFF、位相の振幅カットオフ（無し/0.1/0.2/0.3/0.5/0.7/1.0）
 - fine run 選択ドロップダウン＋時間窓ドロップダウン（100 窓）＋前後窓ボタン
 - データが無い窓/run を選ぶとエラーメッセージ

出力: plots/phase_amplitude_fine/phase_amplitude_fine_viewer_20260917.html
      plots/phase_amplitude_fine/data/data_<rid>_c<ci>.js（増分生成・丸めなし）
データは runs_fine/*/states.npz の read-only 参照のみ。

使い方: python3 plot_phase_amplitude_fine_html_20260917.py
"""
import glob
import json
import os
import sys

import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import plotly.graph_objects as go

HERE = os.path.dirname(os.path.abspath(__file__))
SURVEY = os.path.abspath(os.path.join(HERE, '..'))

CHUNK = 8192
ZOOMS = [1, 10, 50, 100, 200, 400]

VIEWER_JS = """
var gd = document.getElementById('{plot_id}');
var CHUNK = __CHUNK__;
var NCHUNK = __NCHUNK__;
var COLORS = __COLORS_JSON__;
var INITIAL = '__INITIAL__';
var current = null, curChunk = 0, curData = null;
var selected = null;
var phaseOn = true, ampOn = true;
var cutoff = null;
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
      w.push(sel ? 3.0 : (dim ? 0.8 : 1.4)); ms.push(3.2);
    } else if(c === 'amp'){
      op.push(dim ? 0.10 : (sel ? 1.0 : 0.85));
      w.push(1.0); ms.push(sel ? 5.5 : (dim ? 2.0 : 3.2));
    } else {
      op.push(dim ? 0.05 : (sel ? 0.8 : 0.35));
      w.push(sel ? 2.2 : 1.0); ms.push(3.2);
    }
  });
  Plotly.restyle(gd, {'opacity': op, 'line.width': w, 'marker.size': ms});
}
function applyVis(){
  var vis = gd.data.map(function(t){ return (cls(t) === 'phase') ? phaseOn : ampOn; });
  var sl = gd.data.map(function(t){
    var c = cls(t);
    if(c === 'phase'){ return phaseOn; }
    if(c === 'amp'){ return phaseOn ? false : ampOn; }
    return false;
  });
  Plotly.restyle(gd, {'visible': vis, 'showlegend': sl});
}
function applyCutoff(){
  if(!curData){ return; }
  var M = curData.edges.length, idx = [], ys = [];
  for(var e = 0; e < M; e++){
    var ph = curData.phase[e], an = curData.ampn[e], y;
    if(cutoff === null){ y = ph; }
    else {
      y = new Array(ph.length);
      for(var i = 0; i < ph.length; i++){ y[i] = (an[i] >= cutoff) ? ph[i] : null; }
    }
    ys.push(y); idx.push(3 * e);
  }
  Plotly.restyle(gd, {'y': ys}, idx);
}
function buildAndRender(d){
  current = d.rid; curChunk = d.chunk; curData = d; selected = null;
  var M = d.edges.length;
  var cols = COLORS[String(M)];
  var n = d.phase[0].length;
  var tt = new Array(n);
  for(var i = 0; i < n; i++){ tt[i] = d.t0 + i; }
  var traces = [];
  for(var e = 0; e < M; e++){
    var name = '(' + d.edges[e][0] + ',' + d.edges[e][1] + ')';
    var ph = d.phase[e], an = d.ampn[e];
    var cd = new Array(n);
    for(var i2 = 0; i2 < n; i2++){ cd[i2] = [ph[i2] / Math.PI, an[i2] * d.amax, an[i2]]; }
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
  Plotly.relayout(gd, {'xaxis.range': [d.t0, d.t0 + n - 1], 'title.text':
    d.rid + ' — 位相・振幅（fine den=200L）　窓 t=' + d.t0 + '..' + (d.t0 + n - 1) +
    '（' + (d.chunk + 1) + '/' + NCHUNK + '）<br>' +
    '<sup>実線=位相（左軸）・点列＋半透明点線=正規化振幅（右軸）。' +
    '凡例クリック=波選択保持、ドラッグ=窓内スクロール、x?=ズーム率、' +
    '◀/▶=前後の窓。振幅正規化: run 全体最大 ' + d.amax.toFixed(4) + '</sup>'});
  applyVis(); apply(); applyCutoff();
}
function loadChunk(rid, ci){
  if(ci < 0 || ci >= NCHUNK){ return; }
  window.__RUNDATA = null;
  var s = document.createElement('script');
  s.src = 'data/data_' + rid + '_c' + ci + '.js';
  s.onload = function(){
    if(window.__RUNDATA && window.__RUNDATA.rid === rid && window.__RUNDATA.chunk === ci){
      buildAndRender(window.__RUNDATA);
    } else { alert(rid + ' c' + ci + ' のデータ形式が不正です'); }
  };
  s.onerror = function(){
    alert(rid + ' の窓 c' + ci + ' のデータがありません（未走行またはデータ未生成）。');
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
    var width = CHUNK / zoom;
    var lo = current ? curData.t0 : 0;
    var hi = lo + (curData ? curData.phase[0].length - 1 : CHUNK);
    var left = gd.layout.xaxis.range ? gd.layout.xaxis.range[0] : lo;
    if(left < lo){ left = lo; }
    if(left > hi - width){ left = Math.max(lo, hi - width); }
    Plotly.relayout(gd, {'xaxis.range': [left, left + width]});
  } else if(lbl.indexOf('位相') === 0){ phaseOn = !phaseOn; applyVis(); }
  else if(lbl.indexOf('振幅') === 0){ ampOn = !ampOn; applyVis(); }
  else if(lbl.indexOf('カットオフ') === 0){
    var v = lbl.replace('カットオフ', '').trim();
    cutoff = (v === '無し') ? null : parseFloat(v);
    applyCutoff();
  }
  else if(lbl === '◀ 前の窓'){ loadChunk(current, curChunk - 1); }
  else if(lbl === '次の窓 ▶'){ loadChunk(current, curChunk + 1); }
  else if(lbl.indexOf('窓 t=') === 0){
    loadChunk(current, parseInt(lbl.replace('窓 t=', '')) / CHUNK);
  }
  else if(lbl.charAt(0) === 'L'){ loadChunk(lbl, 0); }
});
loadChunk(INITIAL, 0);
"""


def write_run_chunks(src, datadir):
    """1 つの fine run の全チャンク JS を書き出す（増分・丸めなし）。"""
    with np.load(src) as d:
        Z = np.array(d['Z'])
        L = int(d['L'])
        ma, mb, den = int(d['ma']), int(d['mb']), int(d['den'])
        TF = int(d['T'])
    rid = f'L{L}_ma{ma}_mb{mb}_den{den}'
    ja, jb = np.triu_indices(L, k=1)
    edges = [[int(ja[e]), int(jb[e])] for e in range(Z.shape[1])]
    amax = float(np.max(np.abs(Z)))
    nchunk = (TF + CHUNK) // CHUNK          # 819201 状態 → 101 番目は端数
    made = 0
    for ci in range(nchunk):
        t0 = ci * CHUNK
        t1 = min(t0 + CHUNK, TF + 1)
        if t0 >= t1:
            continue
        dst = os.path.join(datadir, f'data_{rid}_c{ci}.js')
        if os.path.exists(dst) and os.path.getmtime(dst) > os.path.getmtime(src):
            continue
        seg = Z[t0:t1]
        payload = {
            'rid': rid, 'chunk': ci, 't0': t0, 'amax': amax, 'edges': edges,
            'phase': [np.angle(seg[:, e]).tolist() for e in range(seg.shape[1])],
            'ampn': [(np.abs(seg[:, e]) / amax).tolist() for e in range(seg.shape[1])],
        }
        with open(dst, 'w', encoding='utf-8') as f:
            f.write('window.__RUNDATA = ')
            json.dump(payload, f, ensure_ascii=False)
            f.write(';\n')
        made += 1
    return rid, L, TF, nchunk, made


def main():
    outdir = os.path.join(SURVEY, 'plots', 'phase_amplitude_fine')
    datadir = os.path.join(outdir, 'data')
    os.makedirs(datadir, exist_ok=True)
    srcs = sorted(glob.glob(os.path.join(SURVEY, 'runs_fine', '*', 'states.npz')))
    if not srcs:
        print('runs_fine に states.npz が無い — 終了')
        sys.exit(1)
    rids, Ls, TF, nchunk = [], set(), None, None
    for src in srcs:
        rid, L, TF, nchunk, made = write_run_chunks(src, datadir)
        rids.append(rid)
        Ls.add(L)
        print(f'{rid}: chunks {nchunk}（新規生成 {made}）')

    colors = {str(L * (L - 1) // 2):
              [mcolors.to_hex(cm.turbo(x))
               for x in np.linspace(0.02, 0.98, L * (L - 1) // 2)]
              for L in sorted(Ls)}
    pi = float(np.pi)
    fig = go.Figure()
    fig.update_layout(
        title=dict(text='fine run と時間窓を選択してください', x=0.5),
        template='plotly_white', width=1180, height=800,
        dragmode='pan',
        xaxis=dict(title='step', range=[0, CHUNK], zeroline=False),
        yaxis=dict(title='位相 arg z_e(t)（実線）',
                   range=[-pi * 1.05, pi * 1.05], fixedrange=True,
                   tickvals=[-pi, -pi / 2, 0, pi / 2, pi],
                   ticktext=['−π', '−π/2', '0', 'π/2', 'π'],
                   zeroline=True, zerolinecolor='lightgray'),
        yaxis2=dict(title='振幅 |z_e(t)| / run 全体最大（点列＋半透明点線）',
                    overlaying='y', side='right', range=[0, 1.05],
                    fixedrange=True, showgrid=False),
        legend=dict(font=dict(size=9), itemsizing='constant',
                    title=dict(text='辺 (j,k)')),
        updatemenus=[
            dict(type='buttons', direction='left', x=0.0, y=1.10, xanchor='left',
                 buttons=[dict(label=f'x{z}', method='skip') for z in ZOOMS]),
            dict(type='buttons', direction='left', x=0.30, y=1.10, xanchor='left',
                 buttons=[dict(label='◀ 前の窓', method='skip'),
                          dict(label='次の窓 ▶', method='skip')]),
            dict(type='dropdown', direction='down', x=0.46, y=1.10, xanchor='left',
                 buttons=[dict(label=f'窓 t={ci * CHUNK}', method='skip')
                          for ci in range(nchunk)]),
            dict(type='buttons', direction='left', x=0.62, y=1.10, xanchor='left',
                 buttons=[dict(label='位相 ON/OFF', method='skip'),
                          dict(label='振幅 ON/OFF', method='skip')]),
            dict(type='dropdown', direction='down', x=0.86, y=1.10, xanchor='left',
                 buttons=[dict(label='カットオフ無し', method='skip')] +
                         [dict(label=f'カットオフ {c}', method='skip')
                          for c in ('0.1', '0.2', '0.3', '0.5', '0.7', '1.0')]),
            dict(type='dropdown', direction='down', x=0.14, y=1.10, xanchor='left',
                 buttons=[dict(label=r, method='skip') for r in rids]),
        ],
        hovermode='closest')
    js = (VIEWER_JS
          .replace('__CHUNK__', str(CHUNK))
          .replace('__NCHUNK__', str(nchunk))
          .replace('__COLORS_JSON__', json.dumps(colors))
          .replace('__INITIAL__', rids[0]))
    out = os.path.join(outdir, 'phase_amplitude_fine_viewer_20260917.html')
    fig.write_html(out, include_plotlyjs=True, full_html=True, post_script=js,
                   config={'scrollZoom': False, 'displaylogo': False})
    print(f'viewer -> {os.path.relpath(out, SURVEY)} '
          f'({os.path.getsize(out) / 1e6:.1f} MB, run {len(rids)} 本, 窓 {nchunk})')


if __name__ == '__main__':
    main()
