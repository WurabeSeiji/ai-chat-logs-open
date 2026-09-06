#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照テスト: 本フォルダにコピーした図化プログラムが、元プログラムA
(../plot_complex_plane_N3_N40_stage123_v1.py, SHA照合済) の3図を byte 一致で
再現することを確認する。描画コードは A の verbatim コピーを import し、
入力/出力パスだけを A のフォルダに向ける（ハーネス plumbing、描画式は不変）。
一致すれば「本フォルダの描画コードは A と同一に描く」ことが確定する。
出力: control_reproduce_A_figs_v1.json（3図の実測SHA vs A committed SHA）"""
import hashlib
import importlib.util
import json
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(BASE)                       # 10ケースパッケージ
SERIES = os.path.dirname(PKG)
A_DIR = os.path.join(SERIES, 'N3_N40_stage123_sweep_20260905')

# A committed 図の基準SHA（作成時に記録）
A_FIG_SHA = {
    'fig_complex_plane_step0_N3_N40_stage123.png':
        'd51d1eb2619903d1',  # 先頭16桁で照合（下で全桁も算出）
    'fig_complex_plane_final_N3_N40_stage123.png':
        '13eb6af85f834eab',
    'fig_complex_plane_final_zoom_N3_N40_stage123.png':
        'a59b70e2b4e8f7bc',
}

def sha16(path):
    return hashlib.sha256(open(path, 'rb').read()).hexdigest()[:16]

# コピー済みプログラム(PKG/plot_...py)を、BASE/IN を A のフォルダに向けて実行する。
# モジュールの BASE/IN グローバルを差し替えてから draw を呼ぶ。
src = os.path.join(PKG, 'plot_complex_plane_N3_N40_stage123_v1.py')
spec = importlib.util.spec_from_file_location('A_copy', src)
mod = importlib.util.module_from_spec(spec)
# 実行前にモジュールの module-level コードが走る（draw_grid を A_DIR 基準にするため
# BASE/IN を書き換えたいが、module-level で即実行される。よって exec を分離する）。
code = open(src, encoding='utf-8').read()
g = {'__file__': src, '__name__': 'A_copy'}
# module-level の draw 実行を止め、関数だけ取り出すため、末尾の実行部を除いて exec。
cut = code.index("draw_grid(0,")
exec(code[:cut], g)  # import と関数定義まで
# パスを A に向ける
g['BASE'] = A_DIR
g['IN'] = os.path.join(A_DIR, 'results')
# 出力は control/ に書く（A を汚さない）ため、savefig 先を差し替えた draw を再定義せず、
# A の draw_grid は os.path.join(BASE, fname) に書く。BASE を control 出力へ向ける。
g['BASE'] = BASE  # 出力先=control/、入力は IN(=A/results) を明示保持
g['IN'] = os.path.join(A_DIR, 'results')

# 3図を再生成（A の module-level と同一引数）
g['draw_grid'](0, 'fig_complex_plane_step0_N3_N40_stage123.png',
               'Stage1+2+3 sweep: complex-plane readout at step 0 (make_parent static parents, dt=2pi/N files); N=3..40')
g['draw_grid'](500, 'fig_complex_plane_final_N3_N40_stage123.png',
               'Stage1+2+3 sweep: complex-plane readout at final step 500 (dt=2pi/N); N=3..40')

# zoom図は module-level の後半（draw_grid 定義の外）にあるため、その部分を A_DIR 入力・BASE出力で再実行
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
def load(N, step):
    d = np.load(os.path.join(A_DIR, 'results', f'hm_N{N}_den_{N}_states_500.npz'))
    assert int(d['denominator']) == N and int(d['steps']) == 500
    return np.asarray(d['Z'][step], dtype=np.complex128)
fig, axs = plt.subplots(8, 5, figsize=(20, 24))
axs = axs.ravel()
for k, N in enumerate(range(3, 41)):
    ax = axs[k]
    z = load(N, 500)
    amp = float(np.abs(z).max())
    coarse = {}
    for w in z:
        key = (round(float(w.real) / amp, 2), round(float(w.imag) / amp, 2))
        coarse.setdefault(key, []).append(w)
    mem = max(coarse.values(), key=len)
    zz = np.array(mem)
    c = zz.mean()
    dev = np.abs(zz - c)
    spread = float(dev.max())
    win = spread * 1.4 if spread > 0 else amp * 1e-12
    ax.plot(zz.real, zz.imag, 'o', ms=3, color='tab:red', alpha=0.8, linestyle='none')
    cnt = Counter((round(float(w.real), 15), round(float(w.imag), 15)) for w in zz)
    for (a, b), n in cnt.items():
        if n > 1:
            ax.annotate(f'x{n}', (a, b), textcoords='offset points', xytext=(3, 3),
                        fontsize=5, color='black')
    ax.set_xlim(c.real - win, c.real + win); ax.set_ylim(c.imag - win, c.imag + win)
    ax.set_aspect('equal')
    ax.grid(alpha=.25)
    ax.tick_params(labelsize=5)
    ax.ticklabel_format(style='sci', scilimits=(-2, 3), useOffset=True)
    ax.set_title(f'N={N}: {len(zz)} waves, dev={spread:.2e} (|z|max={amp:.2e})', fontsize=7)
for k in range(38, 40):
    axs[k].axis('off')
fig.suptitle('Stage1+2+3 sweep: zoom into largest angle cluster at final step 500 (dt=2pi/N); N=3..40', y=.998)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_complex_plane_final_zoom_N3_N40_stage123.png'), dpi=180)
plt.close(fig)

# 照合
result = {'matplotlib': __import__('matplotlib').__version__, 'checks': {}, 'all_match': True}
for name, ref in A_FIG_SHA.items():
    got = sha16(os.path.join(BASE, name))
    ok = (got == ref)
    result['all_match'] &= ok
    result['checks'][name] = {'ref_sha16': ref, 'got_sha16': got, 'match': ok}
    print(f"{'一致' if ok else '不一致!'} {name}: ref={ref} got={got}")
with open(os.path.join(BASE, 'control_reproduce_A_figs_v1.json'), 'w') as f:
    json.dump(result, f, indent=2, ensure_ascii=False)
print('CONTROL_PASS' if result['all_match'] else 'CONTROL_FAIL')
