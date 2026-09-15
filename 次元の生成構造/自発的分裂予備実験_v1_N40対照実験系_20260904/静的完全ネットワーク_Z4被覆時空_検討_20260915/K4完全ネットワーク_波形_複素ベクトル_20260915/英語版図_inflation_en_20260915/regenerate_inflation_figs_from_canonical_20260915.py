#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""図1・図2・図4（inflation 系再掲図）の正本データからの再構成 v2（2026-09-15）。

背景:
  K4 論文が再掲する inflation_N40_phase_step0_cross / step500_ring /
  relation_wave_3D の 3 図は、生成スクリプトが保存されていない。
  既存図は共回転フレーム表示 w(τ)=z(τ)·exp(−iθ0(τ)) である（正本
  plot3d_relphase_complexplane_v2.py と同一のゲージ。v1 再構成の生フレーム
  表示では 3D が螺旋円筒になり既存図と不一致 → 本 v2 で共回転に修正）。

正本（読み取りのみ・書き込まない）:
  (1) hm_N40_den_40_states_500.npz : Z 履歴 (501 step × 780 波)
  (2) plot3d_makeparent_selectedN_v1.series(40) : θ0 成分 (X,Y) と H⊥/H 系列 fz
      （式の再実装なし。θ0 の符号は正本 v2 と同一の床静止基準で自動決定）

手順（規約: コピー→対照→最小変更）:
  step A: 共回転 w で日本語対照版 3 図を再構成し、既存図と目視・数値突合（contrast_ja/）。
  step B: ラベル文字列のみ英語の英語版 3 図を生成（en/）。

既知の差異（隠さず記録する）:
  既存図1のタイトルは「H⊥/H = 1.00e-31」だが、正本の N=40 den=40 step0
  実測値は 3.550e-32 である。本再構成は正本実測値を表示する。
"""
import importlib.util
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
GD = ('/Users/kiharahanakira/Library/CloudStorage/'
      'GoogleDrive-kihara.noriaki@gmail.com/マイドライブ/OneDrive/GitHub/'
      'ai-chat-logs-open')
CANON = os.path.join(GD, '次元の生成構造', '電子の反跳実験', 'seed除去による準安定相追試',
                     'chatgpt追試', '干渉保存力学_資格審査とシード無し系列_20260831',
                     'make_parent型初期値_自己無撞着構造_20260907')

_spec = importlib.util.spec_from_file_location(
    'orig', os.path.join(CANON, 'plot3d_makeparent_selectedN_v1.py'))
orig = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(orig)

N = 40
ONSET_LEVEL = 1e-3          # 正本 plot3d_relphase_complexplane_v2.py と同値
FLOOR_STATIC_ABS = 0.2
FLOOR_SIGN_RATIO = 20.0

# --- 読出し量（正本から取得・再実装なし） ---
X, Y, fz = orig.series(N)
d = np.load(os.path.join(orig.STATES, f'hm_N{N}_den_{N}_states_500.npz'))
S = np.asarray(d['Z'], np.complex128)
assert S.shape == (501, 780) and S.shape[0] == len(X)
onset = int(np.argmax(fz > ONSET_LEVEL))

# --- θ0 の向き決定（正本 v2 の床静止基準をそのままコピー） ---
def floor_var_for(sign):
    th0 = sign * np.arctan2(Y, X)
    Wt = S * np.exp(-1j * th0)[:, None]
    pre = np.unwrap(np.angle(Wt[:max(onset - 10, 2)]), axis=0)
    return float(np.max(pre.max(0) - pre.min(0))), th0, Wt

cand = {s: floor_var_for(s) for s in (+1, -1)}
sign = min(cand, key=lambda s: cand[s][0])
floor_var, theta0, W = cand[sign]
other = max(cand[s][0] for s in cand)
assert floor_var < FLOOR_STATIC_ABS
assert other > FLOOR_SIGN_RATIO * floor_var
H0, H500 = float(fz[0]), float(fz[500])
print(f'θ0向き={sign:+d} 床変動={floor_var:.3e}rad(逆向き{other:.3e}) '
      f'H⊥/H: step0={H0:.3e}, step500={H500:.3e}')

JA_FONT = ['Hiragino Sans', 'sans-serif']

LABELS = {
    'ja': dict(
        t1=f'N=40・M=780  床（step 0）＝十字の波配置\nH⊥/H = {H0:.2e}',
        t2=f'N=40・M=780  終端（step 500）＝放射状のリング\nH⊥/H = {H500:.2e}',
        t4=(f'関係波の膨張＝インフレーション（N=40・M=780本の波）\n'
            f'step 500/500　H⊥/H = {H500:.2e}'),
        z4='log₁₀(H⊥/H) ← 膨張（インフレーション）',
    ),
    'en': dict(
        t1=(f'N=40, M=780  floor (step 0) = cross-shaped wave configuration\n'
            f'H⊥/H = {H0:.2e}'),
        t2=(f'N=40, M=780  terminal (step 500) = radial ring\n'
            f'H⊥/H = {H500:.2e}'),
        t4=(f'Expansion of relation waves = inflation (N=40, M=780 waves)\n'
            f'step 500/500   H⊥/H = {H500:.2e}'),
        z4='log₁₀(H⊥/H) ← expansion (inflation)',
    ),
}

CMAP = plt.get_cmap('tab20', 20)
COLS = [CMAP(e % 20) for e in range(780)]


def scatter_fig(w, title, out, lang):
    fig, ax = plt.subplots(figsize=(8.6, 8.6))
    ax.scatter(w.real, w.imag, s=28, c=COLS, alpha=0.9, linewidths=0)
    ax.axhline(0, color='gray', linewidth=0.6, alpha=0.6)
    ax.axvline(0, color='gray', linewidth=0.6, alpha=0.6)
    ax.set_xlim(-0.043, 0.043)
    ax.set_ylim(-0.043, 0.043)
    ax.set_aspect('equal')
    ax.set_xlabel('Re w', fontsize=13)
    ax.set_ylabel('Im w', fontsize=13)
    kw = {'fontfamily': JA_FONT} if lang == 'ja' else {}
    ax.set_title(title, fontsize=13, **kw)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    plt.close(fig)


def traj3d_fig(title, zlabel, out, lang):
    logh = np.log10(fz)
    fig = plt.figure(figsize=(7.2, 10.0))
    ax = fig.add_subplot(111, projection='3d')
    for e in range(780):
        ax.plot(W[:, e].imag, W[:, e].real, logh,
                color=COLS[e], linewidth=0.35, alpha=0.5)
    ax.scatter(W[500, :].imag, W[500, :].real, np.full(780, logh[-1]),
               s=12, c=COLS, alpha=0.9)
    ax.set_xlabel('Im w', fontsize=10)
    ax.set_ylabel('Re w', fontsize=10)
    kw = {'fontfamily': JA_FONT} if lang == 'ja' else {}
    ax.set_zlabel(zlabel, fontsize=10, **kw)
    fig.suptitle(title, fontsize=12, y=0.97, **kw)
    ax.set_box_aspect((1, 1, 2.2))
    ax.view_init(elev=18, azim=-60)
    fig.subplots_adjust(left=0.02, right=0.92, top=0.90, bottom=0.02)
    fig.savefig(out, dpi=140)
    plt.close(fig)


for lang, sub in (('ja', 'contrast_ja'), ('en', 'en')):
    odir = os.path.join(HERE, sub)
    os.makedirs(odir, exist_ok=True)
    L = LABELS[lang]
    suf = '' if lang == 'ja' else '_en'
    scatter_fig(W[0], L['t1'],
                os.path.join(odir, f'inflation_N40_phase_step0_cross{suf}_20260915.png'), lang)
    scatter_fig(W[500], L['t2'],
                os.path.join(odir, f'inflation_N40_phase_step500_ring{suf}_20260915.png'), lang)
    traj3d_fig(L['t4'], L['z4'],
               os.path.join(odir, f'inflation_N40_relation_wave_3D_step500{suf}_20260915.png'), lang)
    print(f'{lang}: 3 figures written to {sub}/')

print('ALL DONE')
