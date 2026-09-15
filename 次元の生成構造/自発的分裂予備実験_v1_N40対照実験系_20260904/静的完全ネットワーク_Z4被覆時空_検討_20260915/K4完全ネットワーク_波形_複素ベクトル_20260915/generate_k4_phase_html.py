#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate self-contained HTML overlays using the ORIGINAL saved PNG figures.

Top: exact K4_complete_network_m{2,3}.png figure, with only a red current-state ring added.
Bottom: exact waves_base_harmonic_m{2,3}.png figure, with only a red phase line and red sample dots added.

The base figures are not redrawn or replaced.
"""
from pathlib import Path
import base64, math

OUT = Path('/mnt/data')

CFG = {
    2: {
        'net_size': (1643, 1743),
        'nodes': [(499,1133),(1004,1190),(482,424),(1401,814)], # S0..S3 on original PNG
    },
    3: {
        'net_size': (1643, 1750),
        'nodes': [(1043,716),(490,787),(627,1075),(1210,994)], # S0..S3 on original PNG
    },
}

def data_uri(path: Path) -> str:
    b64 = base64.b64encode(path.read_bytes()).decode('ascii')
    return 'data:image/png;base64,' + b64

# Pixel calibration of the ORIGINAL waveform PNGs (2176 x 1180)
# Matplotlib grid/ticks in saved original image:
# theta=0 at x≈1165.5; 2 rad corresponds to ≈286.17 px.
WAVE_W, WAVE_H = 2176, 1180
X0 = 1165.5
PX_PER_RAD = (2024.0 - 307.0) / 12.0
Y0 = 567.5
PX_PER_AMP = 454.0
WAVE_TOP = 70
WAVE_BOTTOM = 1065

HTML = r'''<!doctype html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root{--red:#dc2626;--border:#cbd5e1;--text:#111827;--muted:#475569}
*{box-sizing:border-box}
body{margin:0;background:#f8fafc;color:var(--text);font-family:system-ui,-apple-system,"Segoe UI","Noto Sans JP",sans-serif}
main{max-width:1180px;margin:0 auto;padding:18px}
h1{font-size:22px;margin:0 0 8px}
.note{color:var(--muted);font-size:14px;margin:0 0 12px}
.controls{display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:0 0 14px}
button{font:inherit;min-height:42px;padding:8px 14px;border:1px solid var(--border);border-radius:9px;background:white;cursor:pointer}
button:hover{background:#f1f5f9}
.badge{border:1px solid var(--border);border-radius:999px;background:white;padding:7px 10px;font-variant-numeric:tabular-nums}
.panel{background:white;border:1px solid var(--border);border-radius:14px;padding:12px;margin-bottom:14px}
.panel h2{font-size:16px;margin:0 0 8px}
.figure{position:relative;width:100%;overflow:hidden}
.figure img{display:block;width:100%;height:auto}
.figure svg{position:absolute;inset:0;width:100%;height:100%;pointer-events:none}
.legend{font-size:13px;color:var(--muted);margin-top:7px;line-height:1.45}
</style>
</head>
<body>
<main>
<h1>__TITLE__</h1>
<div class="note">保存済みの元図をそのまま背景に使い、赤い現在位置だけを重ねています。位相は k=0,1,2,3 の4状態を離散的に移動します。</div>
<div class="controls">
  <button id="play">▶ Play</button>
  <button id="pause">⏸ Pause</button>
  <button id="prev">← Step Back</button>
  <button id="next">Step Forward →</button>
  <span class="badge" id="state">k=0 / θ=0</span>
  <span class="badge" id="amp"></span>
</div>

<section class="panel">
<h2>元の完全ネットワーク図 + 現在状態</h2>
<div class="figure">
  <img src="__NET_IMG__" alt="original K4 network">
  <svg viewBox="0 0 __NET_W__ __NET_H__" preserveAspectRatio="none">
    <circle id="netMarker" cx="0" cy="0" r="27" fill="none" stroke="#dc2626" stroke-width="8"/>
  </svg>
</div>
</section>

<section class="panel">
<h2>元の基底波・倍音図 + 離散位相マーカー</h2>
<div class="figure">
  <img src="__WAVE_IMG__" alt="original base and harmonic waveform">
  <svg viewBox="0 0 2176 1180" preserveAspectRatio="none">
    <line id="phaseLine" x1="0" y1="70" x2="0" y2="1065" stroke="#dc2626" stroke-width="5" stroke-dasharray="14 11"/>
    <circle id="baseDot" cx="0" cy="0" r="13" fill="#dc2626" stroke="white" stroke-width="4"/>
    <circle id="harmDot" cx="0" cy="0" r="13" fill="#dc2626" stroke="white" stroke-width="4"/>
  </svg>
</div>
<div class="legend">赤縦点線 = 現在の離散位相 θ=kπ/2、赤丸 = その位相での元図2波の振幅。背景の波形そのものは変更していません。</div>
</section>
</main>
<script>
(() => {
  const m=__M__;
  const nodes=__NODES__;
  const x0=__X0__, pxPerRad=__PPR__, y0=__Y0__, pxPerAmp=__PPA__;
  let k=0, timer=null, playing=false;
  const interval=900;
  const $=id=>document.getElementById(id);
  function render(){
    const theta=k*Math.PI/2;
    const b=Math.cos(theta);
    const h=Math.cos(m*theta);
    const x=x0+pxPerRad*theta;
    const yb=y0-pxPerAmp*b;
    const yh=y0-pxPerAmp*h;
    $('netMarker').setAttribute('cx',nodes[k][0]);
    $('netMarker').setAttribute('cy',nodes[k][1]);
    $('phaseLine').setAttribute('x1',x); $('phaseLine').setAttribute('x2',x);
    $('baseDot').setAttribute('cx',x); $('baseDot').setAttribute('cy',yb);
    $('harmDot').setAttribute('cx',x); $('harmDot').setAttribute('cy',yh);
    const labs=['0','π/2','π','3π/2'];
    $('state').textContent=`k=${k} / θ=${labs[k]}`;
    $('amp').textContent=`cosθ=${clean(b)} / cos(${m}θ)=${clean(h)}`;
  }
  function clean(v){ return Math.abs(v)<1e-10?'0':String(Math.round(v*1000)/1000); }
  function step(d){k=(k+d+4)%4;render()}
  function play(){if(playing)return;playing=true;timer=setInterval(()=>step(1),interval)}
  function pause(){playing=false;if(timer){clearInterval(timer);timer=null}}
  $('play').onclick=play; $('pause').onclick=pause;
  $('prev').onclick=()=>{pause();step(-1)}; $('next').onclick=()=>{pause();step(1)};
  render(); play();
})();
</script>
</body>
</html>'''

for m,cfg in CFG.items():
    net = data_uri(OUT/f'K4_complete_network_m{m}.png')
    wave = data_uri(OUT/f'waves_base_harmonic_m{m}.png')
    w,h = cfg['net_size']
    s = HTML
    repl = {
        '__TITLE__': f'K4 original network + original waveforms (m={m})',
        '__NET_IMG__': net,
        '__WAVE_IMG__': wave,
        '__NET_W__': str(w), '__NET_H__': str(h),
        '__M__': str(m), '__NODES__': repr(cfg['nodes']),
        '__X0__': repr(X0), '__PPR__': repr(PX_PER_RAD), '__Y0__': repr(Y0), '__PPA__': repr(PX_PER_AMP),
    }
    for a,b in repl.items(): s=s.replace(a,b)
    (OUT/f'K4_network_wave_discrete_m{m}.html').write_text(s,encoding='utf-8')

index = '''<!doctype html><html lang="ja"><head><meta charset="utf-8"><title>K4 corrected HTML</title></head><body style="font-family:sans-serif;padding:24px"><h1>K4 corrected HTML</h1><p><a href="K4_network_wave_discrete_m2.html">m=2</a></p><p><a href="K4_network_wave_discrete_m3.html">m=3</a></p></body></html>'''
(OUT/'K4_network_wave_index.html').write_text(index,encoding='utf-8')
print('written corrected HTMLs using original PNGs')
