#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""エイリアス巻数による collapse 再検定（読み出しのみ）。
仮説の精密化（2026-09-05）: 時間の量子化度は q = w_dom/den（1目盛りの巻数）。
流れ領域は q≪1 で素の巻数が時計になるが、ストロボ領域 den≈N は q≈1 なので、
実効時計は最近整数を除いたエイリアス巻数 q_eff = |q − round(q)| のはず。
検定: winding_rates_table_v1.csv の全走行について、
 turns/decade(素) と turns_alias/decade = steps/decade × q_eff の変動係数を
 den 系列別に比較する。ストロボ領域で collapse が復活すれば精密化を支持。
出力: check_aliased_winding_v1.json"""
import csv
import json
import os
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(BASE, 'winding_rates_table_v1.csv'))))

def stats(vals):
    a = np.array(vals)
    return {'n': len(a), 'mean': float(a.mean()), 'std': float(a.std()),
            'cv': float(a.std() / a.mean()) if a.mean() != 0 else None,
            'min': float(a.min()), 'max': float(a.max())}

out = {}
groups = {
    'den=124': lambda r: r['den'] == '124',
    'den=N': lambda r: int(r['den']) == int(r['N']),
    'den=N±1,2': lambda r: r['den'] != '124' and int(r['den']) != int(r['N']),
    'all': lambda r: True,
}
for gname, sel in groups.items():
    sub = [r for r in rows if sel(r)]
    spd = [float(r['steps_per_decade']) for r in sub]
    q = [float(r['w_dom']) / int(r['den']) for r in sub]
    q_eff = [abs(x - round(x)) for x in q]
    t_raw = [s * x for s, x in zip(spd, q)]
    t_alias = [s * x for s, x in zip(spd, q_eff)]
    out[gname] = {
        'q_range': [float(min(q)), float(max(q))],
        'q_eff_range': [float(min(q_eff)), float(max(q_eff))],
        'steps_per_decade': stats(spd),
        'turns_per_decade_raw': stats(t_raw),
        'turns_per_decade_aliased': stats(t_alias),
    }
with open(os.path.join(BASE, 'check_aliased_winding_v1.json'), 'w') as f:
    json.dump(out, f, indent=2)
for g, v in out.items():
    print(f"{g}: q∈[{v['q_range'][0]:.3f},{v['q_range'][1]:.3f}] "
          f"CV steps={v['steps_per_decade']['cv']:.3f} "
          f"raw_turns={v['turns_per_decade_raw']['cv']:.3f} "
          f"alias_turns={v['turns_per_decade_aliased']['cv']:.3f} "
          f"(alias mean {v['turns_per_decade_aliased']['mean']:.3f})")
print('ALL DONE')
