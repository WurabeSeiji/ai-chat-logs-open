#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""疎 / matrix-free カーネル実装（仕様書 v1.2 §11 準拠）。

密参照（kernel_canonical_dense.py = 正本コピー）と同一の写像
  F(Z) = exp(Δτ K(Z)) Z,   K_ef = A_ef sin(θ_f − θ_e)
を、密 M×M 行列を形成せずに計算する。

数学的同一性（実装の根拠、物理式の変更なし）:
  (K z)_e = Σ_f A_ef Im(conj(u_e) u_f) z_f          (u = exp(iθ))
          = [conj(u)⊙(A(u⊙z)) − u⊙(A(conj(u)⊙z))] / (2i)
  A y = B(B^T y) − 2y   （B: M×P の辺-頂点接続 0/1 行列。
        (BB^T)_ef = 共有頂点数 = 2(e=f) / 1(共有) / 0（非共有））
  exp(ΔτK) z = exp(−iΔτH) z,  H = iK は Hermitian
        → Hermitian Lanczos による Krylov 近似（許容誤差 tol を明示）。

辺列挙は正本と同一の np.triu_indices（kernel_canonical_dense.edges を import）。
"""
import math

import numpy as np
from scipy.sparse import csr_matrix

from kernel_canonical_dense import edges


class SparseKernel:
    def __init__(self, P):
        self.P = P
        self.ea, self.eb = edges(P)
        self.M = len(self.ea)
        rows = np.concatenate([np.arange(self.M), np.arange(self.M)])
        cols = np.concatenate([self.ea, self.eb])
        data = np.ones(2 * self.M)
        self.B = csr_matrix((data, (rows, cols)), shape=(self.M, self.P))
        self.BT = self.B.T.tocsr()

    def A_apply(self, y):
        return self.B @ (self.BT @ y) - 2.0 * y

    def K_apply(self, u, x):
        return (np.conj(u) * self.A_apply(u * x) - u * self.A_apply(np.conj(u) * x)) / 2j

    def one_step(self, z, den, tol=1e-14, mmax=None):
        """exp(Δτ K(z)) z を Hermitian Lanczos（H=iK）で計算。Δτ=2π/den。"""
        dtau = 2.0 * math.pi / den
        u = np.exp(1j * np.angle(z))

        def H_apply(x):
            return 1j * self.K_apply(u, x)

        M = self.M
        if mmax is None:
            mmax = M
        beta0 = np.linalg.norm(z)
        if beta0 == 0.0:
            return z.copy()
        V = np.empty((mmax + 1, M), dtype=np.complex128)
        alpha = np.zeros(mmax, dtype=np.float64)
        beta = np.zeros(mmax + 1, dtype=np.float64)
        V[0] = z / beta0
        prev = None
        y_final = None
        m_used = mmax
        for m in range(mmax):
            w = H_apply(V[m])
            if m > 0:
                w = w - beta[m] * V[m - 1]
            alpha[m] = np.vdot(V[m], w).real
            w = w - alpha[m] * V[m]
            # 完全再直交化（M が小さいので堅牢性優先）
            for k in range(m + 1):
                w = w - np.vdot(V[k], w) * V[k]
            beta[m + 1] = np.linalg.norm(w)
            T = np.diag(alpha[:m + 1])
            if m > 0:
                T += np.diag(beta[1:m + 1], 1) + np.diag(beta[1:m + 1], -1)
            wT, VT = np.linalg.eigh(T)
            e1 = VT.conj().T[:, 0]
            y = VT @ (np.exp(-1j * dtau * wT) * e1) * beta0
            if prev is not None:
                diff = np.linalg.norm(y[:len(prev)] - np.append(prev, 0.0)[:len(prev)] if len(prev) < len(y) else y - prev[:len(y)])
                # 収束: Krylov 次元追加による結果変化が tol 以下
                full_prev = np.zeros(m + 1, dtype=np.complex128)
                full_prev[:len(prev)] = prev
                diff = np.linalg.norm(y - full_prev)
                if diff <= tol * beta0:
                    y_final = y
                    m_used = m + 1
                    break
            prev = y
            if beta[m + 1] <= 1e-300:
                y_final = y
                m_used = m + 1
                break
            V[m + 1] = w / beta[m + 1]
        if y_final is None:
            y_final = prev
            m_used = mmax
        out = (V[:m_used].T @ y_final).astype(np.complex128, copy=False)
        return out
