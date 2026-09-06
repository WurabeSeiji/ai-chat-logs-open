#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""step0 と完了時（step10000）の複素平面読出し図・10ケース（読み出しのみ・新規走行なし）。
plot_complex_plane_N3_N40_stage123_v1.py（本フォルダに SHA 一致で忠実コピー・対照テスト済）
の per-panel 描画コード（原点からの線分・全体振幅スケール 1.15r・実値目盛・12桁丸め重複の
x本数表記）を verbatim 流用。変更は3点のみ:
 (1) 入力を 500歩スイープ(states_500) から 10000歩スイープ(states_10000) に差し替え（完了時=真の収束状態）。
 (2) ループを N=3..40 の 8×5 グリッドから、代表10ケース N=[3,4,5,6,7,8,10,12,20,40] に置換。
 (3) レイアウトを各N「step0 ｜ 完了時」の1行2列に再構成（対照テスト用: step0と完了時の差を並置）。
入力: ../N3_N40_long10000_20260905/results/hm_N{N}_den_{N}_states_10000.npz（SHA台帳照合）。
出力: fig_step0_final_N{N}.png × 10 と、一覧用 fig_step0_final_contact_10cases.png。"""
import hashlib
import os
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
SERIES = os.path.dirname(BASE)
LONG = os.path.join(SERIES, 'N3_N40_long10000_20260905')
IN = os.path.join(LONG, 'results')

NS = [3, 4, 5, 6, 7, 8, 10, 12, 20, 40]
FINAL = 10000

# 入力SHA台帳（読み出しゲート）
ledger = {}
with open(os.path.join(LONG, 'SHA256SUMS.txt')) as fh:
    for line in fh:
        parts = line.split()
        if len(parts) == 2:
            ledger[parts[1]] = parts[0]


def load(N, step):
    rel = f'results/hm_N{N}_den_{N}_states_10000.npz'
    path = os.path.join(LONG, rel)
    h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    assert ledger[rel] == h, f'INPUT GATE FAIL: {rel}'
    d = np.load(path)
    assert int(d['denominator']) == N and int(d['steps']) == FINAL
    return np.asarray(d['Z'][step], dtype=np.complex128)


def draw_panel(ax, z, title):
    """A の per-panel 描画コードの verbatim 流用（描画式は不変）。"""
    for w in z:
        ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.7, alpha=0.6)
    ax.plot(z.real, z.imag, 'o', ms=2.5, color='tab:red', alpha=0.85, linestyle='none')
    cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in z)
    for (a, b), c in cnt.items():
        if c > 1:
            ax.annotate(f'x{c}', (a, b), textcoords='offset points', xytext=(3, 3),
                        fontsize=6, color='black')
    r = float(np.abs(z).max())
    lim = r * 1.15 if r > 0 else 1.0
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.grid(alpha=.25)
    ax.tick_params(labelsize=7)
    ax.ticklabel_format(style='sci', scilimits=(-2, 3))
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('Re z', fontsize=8); ax.set_ylabel('Im z', fontsize=8)


# (A) 各Nの個別図 10枚: [step0 ｜ 完了時]
for N in NS:
    M = N * (N - 1) // 2
    z0 = load(N, 0)
    zf = load(N, FINAL)
    assert z0.size == M and zf.size == M
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.6))
    draw_panel(axs[0], z0, f'N={N} (M={M}) — step 0')
    draw_panel(axs[1], zf, f'N={N} (M={M}) — final step {FINAL}')
    fig.suptitle(f'Stage1+2+3 (dt=2pi/N): complex-plane readout, N={N}  [step 0 vs final]', y=.99)
    fig.tight_layout()
    fig.savefig(os.path.join(BASE, f'fig_step0_final_N{N}.png'), dpi=170)
    plt.close(fig)
    print(f'N={N}: |z0|max={np.abs(z0).max():.3e} |zf|max={np.abs(zf).max():.3e} 保存 fig_step0_final_N{N}.png', flush=True)

# (B) 一覧用コンタクトシート 10行×2列（step0｜完了時）
fig, axs = plt.subplots(10, 2, figsize=(11, 52))
for i, N in enumerate(NS):
    draw_panel(axs[i, 0], load(N, 0), f'N={N} — step 0')
    draw_panel(axs[i, 1], load(N, FINAL), f'N={N} — final {FINAL}')
fig.suptitle('Stage1+2+3 (dt=2pi/N): step 0 vs final step 10000 — 10 cases (existing data, no new run)', y=.999)
fig.tight_layout()
fig.savefig(os.path.join(BASE, 'fig_step0_final_contact_10cases.png'), dpi=140)
plt.close(fig)
print('contact sheet 保存 fig_step0_final_contact_10cases.png')
print('ALL DONE')
