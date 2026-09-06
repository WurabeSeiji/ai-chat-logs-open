#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v3 系列則（仕様書 §2v3〜§2v3d・§3 A/B案）の複素平面図（2026-09-06 木原指示）。

(1) fig_v3_centroid_zero_series_N3_N40.png —— 重心 Σz=0 が厳密な全38系（8×5 グリッド）
    N=3: B案 120°三つ組 / N=4,5: B案（§2v3b）/ mod4=0: §2v3 / mod4=2: §2v3b / mod4=1: §2v3c / mod4=3: §2v3d
(2) fig_v3_centroid_nonzero_exceptions.png —— 重心が非ゼロの例外（A案 T型 N=3,4,5、1×3）

状態は仕様書の閉形式から生成し、描画前に4条件（‖Z‖²=1・Σz・Σa²−Σb²・Σab）を機械検査して
タイトルに |Σz| を表示する。図様式は plot_complex_plane_N3_N40_stage123_v1.py と同一
（原点からの線分・赤点・12桁丸め重複の x本数表記・等アスペクト）。"""
import math
import os
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = ['Hiragino Sans', 'DejaVu Sans']  # 日本語ラベル対応（macOS）

BASE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(BASE, 'figures'), exist_ok=True)


def build_v3(N):
    """仕様書 §2v3〜§2v3d の重心ゼロ設計（N=3 は B案三つ組、N=4,5 は B案）"""
    M = N * (N - 1) // 2
    if N == 3:  # B案: 120° 三つ組
        r = 1 / math.sqrt(3)
        return np.array([r * np.exp(1j * 2 * math.pi * k / 3) for k in range(3)]), 'B案 120°三つ組'
    m = M % 4
    if m == 0:  # §2v3: 等振幅1値
        q = M // 4; r = 1 / math.sqrt(M)
        Z = [r] * q + [-r] * q + [1j * r] * q + [-1j * r] * q
        return np.array(Z, complex), '1振幅 M/4×4位置'
    if m == 2:  # §2v3b: 両偶数分割・2振幅
        nA = M // 2 + 1; nB = M // 2 - 1
        rA = math.sqrt(1 / (2 * nA)); rB = math.sqrt(1 / (2 * nB))
        Z = [rA] * (nA // 2) + [-rA] * (nA // 2) + [1j * rB] * (nB // 2) + [-1j * rB] * (nB // 2)
        return np.array(Z, complex), f'2振幅 ({nA},{nB})分割'
    # m==1 / m==3: §2v3c / §2v3d 整数トリオ入り
    T = (M - 3) // 2
    q2 = 2 if m == 1 else 1
    q1 = (T - 7) // 2 if m == 1 else T // 2 - 1
    p1 = q1 + 5 if m == 1 else q1 + 1
    tot = (6 + 2 * p1) + (2 * q1 + 8 * q2)
    r = 1 / math.sqrt(tot)
    Z = [2 * r, -r, -r] + [r] * p1 + [-r] * p1
    Z = [x + 0j for x in Z] + [1j * r] * q1 + [-1j * r] * q1 + [2j * r] * q2 + [-2j * r] * q2
    return np.array(Z, complex), f'2振幅 トリオ+対 (q2={q2})'


def build_A(N):
    """§3 A案（T型・旧親と同一構成の正準ゲージ版）"""
    if N == 3:
        return np.array([1 / math.sqrt(2), 0.5j, -0.5j])
    if N == 4:
        r = 1 / (2 * math.sqrt(2))
        return np.array([r, r, r, r, 0.5j, -0.5j])
    if N == 5:
        r42 = math.sqrt(42); r12 = math.sqrt(12)
        return np.array([2 / r42, 2 / r42, 2 / r42, 3 / r42,
                         1j / r12, 1j / r12, 1j / r12, -1j / r12, -1j / r12, -1j / r12])
    raise ValueError(N)


def verify(Z):
    a, b = Z.real, Z.imag
    return dict(norm=float(np.vdot(Z, Z).real), cz=float(abs(Z.sum())),
                pw=float(a @ a - b @ b), ab=float(a @ b))


def panel(ax, Z, title):
    for w in Z:
        ax.plot([0.0, w.real], [0.0, w.imag], color='tab:blue', linewidth=0.7, alpha=0.6)
    ax.plot(Z.real, Z.imag, 'o', ms=2.5, color='tab:red', alpha=0.85, linestyle='none')
    cnt = Counter((round(float(w.real), 12), round(float(w.imag), 12)) for w in Z)
    for (x, y), c in cnt.items():
        if c > 1:
            ax.annotate(f'x{c}', (x, y), textcoords='offset points', xytext=(3, 3), fontsize=5)
    r = float(np.abs(Z).max()); lim = r * 1.2
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
    ax.axhline(0, color='gray', linewidth=0.5, alpha=0.5); ax.axvline(0, color='gray', linewidth=0.5, alpha=0.5)
    ax.grid(alpha=.25); ax.tick_params(labelsize=6); ax.ticklabel_format(style='sci', scilimits=(-2, 3))
    ax.set_title(title, fontsize=8)


if __name__ == '__main__':
    # (1) 重心ゼロ系列 38系
    fig, axs = plt.subplots(8, 5, figsize=(20, 24))
    axs = axs.ravel()
    for k, N in enumerate(range(3, 41)):
        Z, kind = build_v3(N)
        v = verify(Z)
        assert abs(v['norm'] - 1) < 1e-12 and v['cz'] < 1e-12 and abs(v['pw']) < 1e-12 and abs(v['ab']) < 1e-12, (N, v)
        M = N * (N - 1) // 2
        panel(axs[k], Z, f"N={N} (M={M}, mod4={M % 4}) {kind}\n|Σz|={v['cz']:.1e}")
    for k in range(38, 40):
        axs[k].axis('off')
    fig.suptitle('v3 centroid-zero series (all four conditions exact): N=3..40 — N=3: plan B triplet; N=4,5: plan B; rules §2v3–§2v3d', y=.998)
    fig.tight_layout()
    out1 = os.path.join(BASE, 'figures', 'fig_v3_centroid_zero_series_N3_N40.png')
    fig.savefig(out1, dpi=180); plt.close(fig)

    # (2) 重心非ゼロの例外（A案 T型）
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.6))
    for ax, N in zip(axs, (3, 4, 5)):
        Z = build_A(N)
        v = verify(Z)
        assert abs(v['norm'] - 1) < 1e-12 and abs(v['pw']) < 1e-12 and abs(v['ab']) < 1e-12
        M = N * (N - 1) // 2
        panel(ax, Z, f"N={N} (M={M}) A案 T型\nΣz={Z.sum().real:+.6f}（重心非ゼロ・点火実績あり）")
    fig.suptitle('v3 centroid-nonzero exceptions (plan A, T-type = old-parent structure in canonical gauge): N=3,4,5', y=1.0)
    fig.tight_layout()
    out2 = os.path.join(BASE, 'figures', 'fig_v3_centroid_nonzero_exceptions.png')
    fig.savefig(out2, dpi=180); plt.close(fig)
    print('saved:', out1); print('saved:', out2); print('ALL DONE')
