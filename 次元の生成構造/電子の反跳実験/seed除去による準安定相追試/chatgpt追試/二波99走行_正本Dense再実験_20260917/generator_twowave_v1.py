#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase A 初期条件生成器（仕様書 v1.2 §4 準拠）。

生成式（等振幅・半位相刻み、全辺非零は偶数 L で厳密保証）:
  U = exp(2πi/L), c_a = exp(-iπ/(2L)), c_b = exp(+iπ/(2L))
  z_{jk}(0) = c_a U^{m_a Δ} + c_b U^{m_b Δ},  Δ = k - j (j<k, 正本の triu 列挙)

規律: 生成時ラベル・モード番号は provenance であり、力学カーネルへは
{z_e} のみを渡す（dynamics metadata ≠ analysis metadata）。

物理コピー元: PhaseA_実装一式_v1_20260915/generator_phaseA.py
変更は import 元のみ（外部フォルダ import 禁止のため、本実験フォルダ内の
interaction_kernel_theory_v1 から edges を取る）。数式・生成手順は無変更。
"""
import numpy as np

from interaction_kernel_theory_v1 import edges


def generate(L, m_a=2, m_b=3):
    """匿名状態ベクトル Z0 のみを返す（力学へ渡してよい唯一の出力）。"""
    ea, eb = edges(L)
    delta = (eb - ea).astype(np.float64)
    U = np.exp(2j * np.pi / L)
    c_a = np.exp(-1j * np.pi / (2 * L))
    c_b = np.exp(+1j * np.pi / (2 * L))
    z0 = c_a * U ** (m_a * delta) + c_b * U ** (m_b * delta)
    return z0.astype(np.complex128)


def analysis_basis(L, m_a=2, m_b=3):
    """解析側専用（provenance）。Gram 補正射影の B 行列。力学へ渡さない。"""
    ea, eb = edges(L)
    delta = (eb - ea).astype(np.float64)
    U = np.exp(2j * np.pi / L)
    B = np.stack([U ** (m_a * delta), U ** (m_b * delta)], axis=1)
    return B
