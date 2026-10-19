#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""図「関係波の膨張＝インフレーション（N=40・M=780本の波）」の読み出しの忠実コピー。

出典（読み取りのみ・式の再実装なし。下の二関数は本文を写したもの）：
  series()        ← 次元の生成構造/電子の反跳実験/seed除去による準安定相追試/chatgpt追試/
                     干渉保存力学_資格審査とシード無し系列_20260831/make_parent型初期値_自己無撞着構造_20260907/
                     plot3d_makeparent_selectedN_v1.py の series(N)
  corotating_W()  ← 次元の生成構造/自発的分裂予備実験_v1_N40対照実験系_20260904/静的完全ネットワーク_Z4被覆時空_検討_20260915/
                     K4完全ネットワーク_波形_複素ベクトル_20260915/英語版図_inflation_en_20260915/
                     regenerate_inflation_figs_from_canonical_20260915.py の floor_var_for と θ0 向き決定

変更点：入力を「ファイル名固定の N」から「npz のパス」に変えただけ（den を変えた走行を読むため）。
"""
import numpy as np

FLOOR = 1e-31
ONSET_LEVEL = 1e-3          # 正本 plot3d_relphase_complexplane_v2.py と同値
FLOOR_STATIC_ABS = 0.2
FLOOR_SIGN_RATIO = 20.0


def series_from_states(S):
    # plot3d_makeparent_selectedN_v1.series(N) の本文（np.load 部分を除く）
    z0 = S[0]
    p = z0.real.astype(np.float64).copy(); p /= np.linalg.norm(p)
    q = z0.imag.astype(np.float64).copy(); q -= np.dot(q, p) * p; q /= np.linalg.norm(q)
    T = S.shape[0]
    f = np.empty(T); cp = np.empty(T, complex); cq = np.empty(T, complex)
    for s in range(T):
        z = S[s]
        a = np.dot(p, z); b = np.dot(q, z)
        zp = z - p * a - q * b
        f[s] = np.vdot(zp, zp).real / np.vdot(z, z).real
        cp[s] = a; cq[s] = b
    fz = np.maximum(f, FLOOR)
    V = np.stack([cp.real, cp.imag, cq.real, cq.imag], axis=1)
    W = V - V.mean(0)
    U, sv, Vt = np.linalg.svd(W, full_matrices=False)
    signs = np.sign(U[np.argmax(np.abs(U), axis=0), np.arange(U.shape[1])])  # svd_flip(U基準)・決定論
    PC = (U * signs) * sv
    return PC[:, 0], PC[:, 1], fz


def corotating_W(S, X, Y, fz, strict=True):
    # regenerate_inflation_figs_from_canonical_20260915.py の θ0 向き決定をそのまま
    onset = int(np.argmax(fz > ONSET_LEVEL))

    def floor_var_for(sign):
        th0 = sign * np.arctan2(Y, X)
        Wt = S * np.exp(-1j * th0)[:, None]
        pre = np.unwrap(np.angle(Wt[:max(onset - 10, 2)]), axis=0)
        return float(np.max(pre.max(0) - pre.min(0))), th0, Wt

    cand = {s: floor_var_for(s) for s in (+1, -1)}
    sign = min(cand, key=lambda s: cand[s][0])
    floor_var, theta0, W = cand[sign]
    other = max(cand[s][0] for s in cand)
    ok = (floor_var < FLOOR_STATIC_ABS) and (other > FLOOR_SIGN_RATIO * floor_var)
    if strict:
        assert ok, (floor_var, other)
    return W, theta0, sign, floor_var, other, onset, ok
