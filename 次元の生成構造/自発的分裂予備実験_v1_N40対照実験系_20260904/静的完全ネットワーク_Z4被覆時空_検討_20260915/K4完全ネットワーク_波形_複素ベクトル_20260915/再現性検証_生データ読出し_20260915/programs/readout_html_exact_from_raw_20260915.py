#!/usr/bin/env python3
from pathlib import Path
import json,csv,pandas as pd
HERE=Path(__file__).resolve().parent
RAW=HERE.parent/'K4_reproducibility_raw_20260915'
OUT=HERE/'reproduced_outputs'; OUT.mkdir(parents=True,exist_ok=True)
TEMPLATE=Path('/mnt/data/K4_harmonic_discrete_phase_animation_20260915.html')
html=TEMPLATE.read_text(encoding='utf-8')
# raw state data, preserving signed zero strings then float JSON values (sign is irrelevant in JS after Math.round in original UI)
data={}
for m in (2,3):
    rows=[]
    with open(RAW/'complex_states_animation_exact_i.csv',newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if int(r['m'])==m and r['direction']=='forward':
                rows.append({k:(int(v) if k in ('m','step','k') else (v if k=='direction' else float(v))) for k,v in r.items()})
    rows.sort(key=lambda r:r['step']); data[f'm{m}_states']=rows
    w=pd.read_csv(RAW/f'wave_samples_html_800_m{m}.csv',float_precision='round_trip')
    data[f'm{m}_waves']={c:w[c].tolist() for c in ['base','harmonic','composite']}
D=json.dumps(data,separators=(',',':'))
# inject raw data
html=html.replace("(() => {\n  const root", "(() => {\n  const RAW_DATA="+D+";\n  const root")
# replace path function with array readout preserving same x-index and toFixed behavior
start=html.index('  function pathFor(fn, center, amp){')
end=html.index('\n\n  const basePath=',start)
newfun='''  function pathForValues(values, center, amp){\n    const n=values.length, pts=[];\n    for(let i=0;i<n;i++){\n      const x=W.x0+(W.x1-W.x0)*i/(n-1);\n      const y=center-amp*values[i];\n      pts.push(`${i?'L':'M'}${x.toFixed(2)},${y.toFixed(2)}`);\n    }\n    return pts.join(' ');\n  }'''
html=html[:start]+newfun+html[end:]
# replace cexp with raw-state lookup
start=html.index('  function cexp(q){')
end=html.index('\n  function fmt(z){',start)
html=html[:start]+'''  function stateRow(){ return RAW_DATA[`m${m}_states`][k]; }\n  function cexp(which){\n    const r=stateRow();\n    return which==='base' ? {re:r.base_re,im:r.base_im} : {re:r.harmonic_re,im:r.harmonic_im};\n  }\n'''+html[end:]
# replace static curves
old='''  function renderStaticForM(){\n    basePath.setAttribute('d',pathFor(th=>Math.cos(th),W.yBase,W.amp));\n    harmPath.setAttribute('d',pathFor(th=>Math.cos(m*th),W.yHarm,W.amp));\n    sumPath.setAttribute('d',pathFor(th=>Math.cos(th)+Math.cos(m*th),W.ySum,W.amp*0.72));\n  }'''
new='''  function renderStaticForM(){\n    const w=RAW_DATA[`m${m}_waves`];\n    basePath.setAttribute('d',pathForValues(w.base,W.yBase,W.amp));\n    harmPath.setAttribute('d',pathForValues(w.harmonic,W.yHarm,W.amp));\n    sumPath.setAttribute('d',pathForValues(w.composite,W.ySum,W.amp*0.72));\n  }'''
assert old in html; html=html.replace(old,new)
# replace state extraction and discrete wave marker values
html=html.replace('''    const theta=k*Math.PI/2;\n    const z1=cexp(k);\n    const zm=cexp((m*k)%4);''','''    const r=stateRow();\n    const theta=r.theta_rad;\n    const z1=cexp('base');\n    const zm=cexp('harmonic');''')
html=html.replace('''    const y1=W.yBase-W.amp*Math.cos(theta);\n    const ym=W.yHarm-W.amp*Math.cos(m*theta);\n    const ys=W.ySum-W.amp*0.72*(Math.cos(theta)+Math.cos(m*theta));''','''    const y1=W.yBase-W.amp*r.base_re;\n    const ym=W.yHarm-W.amp*r.harmonic_re;\n    const ys=W.ySum-W.amp*0.72*(r.sum_re);''')
out=OUT/'K4_harmonic_discrete_phase_animation_from_raw_exact_layout_20260915.html'
out.write_text(html,encoding='utf-8')
print(out)
