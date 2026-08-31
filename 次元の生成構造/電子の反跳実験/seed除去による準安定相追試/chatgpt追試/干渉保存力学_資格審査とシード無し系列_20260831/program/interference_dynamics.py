#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""干渉保存フレームの統一相互作用関数と統一読出し。
フレーム宣言：本パッケージの力学は「干渉保存力学」（新フレーム）。監査 20260831 の帰結 1 に従い明示宣言する。
  read : H_ef = A_ef · conj(z_e) z_f（頂点共有対の複素干渉積の全体。無損失 read。係数自由度なし）
  step : z' = exp(−i Δ H(z)) z（凍結エルミート生成子のユニタリ厳密指数。1 step 1 回・固定 Δ・レジスタなし・IF なし）
現行フレーム（振幅込み実反対称 K・exp(ΔK)）との関係：−iH = K − iP、K = Im H（実反対称＝現行生成子）、
P = Re H（実対称＝これまで捨てていた干渉直交成分）。P を落とすと現行力学に厳密に一致する（資格審査 V0.4）。
保存則：‖z‖² は厳密保存（ユニタリ。連続流 ż=−iH(z)z も z†Hz∈R により保存）。Σz² は保存しない（力学量＝反射チャネル）。
頂点形：(Hz)_e = conj(z_e)·(A z²)_e = conj(z_e)·(S_i+S_j−2z_e²)。エラスティック項は両直交成分間で相殺し、
近傍 z²（第 2 倍音成分）が駆動する位相共役ポンプだけが残る。
相対平衡：H(v)v = μv（μ 実）⇔ z_t = e^{−iΔμt}v が厳密軌道（H は大域位相不変・スケール 2 次同次）。
定理（数値確認済み・資格審査 V0.5）：局所閉塞（S_i≡0）＋等モジュラー ⇒ H(v)v = −2r²v。
N=3 Z3 配置は S_i≠0 だが固有ベクトルで μ = −r²。"""
import math
import numpy as np

L = 124
DELTA = 2 * math.pi / L

def hermitian_H(z, A):
    """無損失 read：H_ef = A_ef conj(z_e) z_f（エルミート）。"""
    return A * (np.conj(z)[:, None] * z[None, :])

def unified_interference_step(z, A, delta=DELTA):
    """統一相互作用（新フレーム）：z' = exp(−i δ H(z)) z。固有分解による厳密ユニタリの 1 回書き込み。"""
    H = hermitian_H(z, A)
    dev = float(np.linalg.norm(H - H.conj().T))
    if dev > 1e-10 * max(1.0, float(np.linalg.norm(H))):
        raise RuntimeError(f"hermiticity failure: {dev}")
    w, V = np.linalg.eigh(H)
    return V @ (np.exp(-1j * delta * w) * (V.conj().T @ z))

def current_frame_step(z, A, delta=DELTA):
    """対照用：現行フレーム（P=Re H を捨てた力学）z' = exp(δK)z、K=Im H。資格審査 V0.4 でのみ使用。"""
    K = A * np.imag(np.conj(z)[:, None] * z[None, :])
    w, V = np.linalg.eigh(1j * K)
    return V @ (np.exp(-1j * delta * w) * (V.conj().T @ z))

def unified_readout(z, A, edges_list):
    """統一読出し（選択なし・全量）：ノルム・閉塞（大域/局所）・新フレーム自己無撞着（μ・残差）・PR・振幅統計。"""
    H2 = float(np.vdot(z, z).real)
    Hz = hermitian_H(z, A) @ z
    mu = float((np.vdot(z, Hz) / np.vdot(z, z)).real)
    res = float(np.linalg.norm(Hz - mu * z) / np.linalg.norm(z))
    d2 = z * z
    glob = float(abs(d2.sum()) / H2)
    N = max(max(e) for e in edges_list) + 1
    S = np.zeros(N, complex)
    for k, (i, j) in enumerate(edges_list):
        S[i] += d2[k]; S[j] += d2[k]
    loc = float(abs(S).max() / H2)
    a2 = np.abs(z) ** 2
    pr = float((a2.sum() ** 2) / (a2 ** 2).sum())
    return dict(H_total=H2, mu_new=mu, residual_new=res, global_closure=glob,
                local_closure=loc, PR=pr,
                amp_min=float(np.abs(z).min()), amp_max=float(np.abs(z).max()),
                amp_std=float(np.abs(z).std()))
