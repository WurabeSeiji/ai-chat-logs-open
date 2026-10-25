#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""密行列版の自己無撞着親生成器（2026-10-25）。

旧エンジン run_n_scaling_lowrank_v1.py の make_parent / zero_closure_kernel_seed と同じ手順を、
走行本体 run_subset_stage123_v1.py と同じ H_of から作る密行列 K（式18）で行う。
目的：辺集合が完全グラフの上三角（M=N(N−1)/2）でない場合、特に対角込みの全順序対（M=N²）でも親を作るため。
旧エンジンの解析的グラム行列（set_theta）は「頂点対ごとに共有辺 1 本・自己対なし」の完全グラフ専用で使えない。

手順の対応（旧 → 本版）
  set_theta(θ)                         → K_of(θ, A) = Im(H_of(e^{iθ}, A))        … 式18 の K。A は走行本体と同じ adjacency（0/1）
  eig(J G), idx=argmin(ev.imag) (λ=−iσ_max) → eigh(iK) の最大固有値 σ_max の固有ベクトル  … K v = −iσ_max v ⟺ iK v = +σ_max v
  v = W y, v ← v/‖v‖                   → eigh の固有ベクトル（単位ノルム）
  θ ← arg((1−β)e^{iθ} + β e^{i arg v}), β=0.5 → 同じ
  10 反復ごとに set_theta(arg v) → 残差 μ=Re v†(iKv), r=‖iKv−μv‖ < tol で停止、restarts=3、最良残差を採用 → 同じ
  sigma_spectrum: JG の固有値の正の虚部（>1e-12）降順 → eigh(iK) の正の固有値（>1e-12）降順
  zero_closure_kernel_seed: span(W) の直交補へ射影 → K の値域（SVD、相対しきい 1e-10）の直交補へ射影。u, w の乱数の順・直交化・g=(u+iw)/√2 は同じ
  Z0 = (v + δ g)/‖v + δ g‖, δ=1e-15      → 同じ
  乱数: default_rng(40260722 + 1000N + seed)、消費順「θ⁰ 一様 m 個（リスタート毎）→ ξ₁ 正規 m 個 → ξ₂ 正規 m 個」 → 同じ

避けられない相違：固有ベクトルの大域位相と固有値計算の経路が旧エンジンと違うので反復の経路は一致せず、収束先は同値な別解になり得る。
これは control_dense_parents_vs_engine_v1.py で現構造（上三角）について測る。
"""
import math
import os
import sys
import time

import numpy as np

SEED_BASE = 40260722


def edges_triu(N):
    """原本の辺集合：完全グラフの上三角 i<j。M = N(N−1)/2。"""
    a, b = np.triu_indices(N, k=1)
    return a.astype(np.int64), b.astype(np.int64)


def edges_all(N):
    """走行1 の辺集合：対角込みの全順序対 (i, j)。M = N²。順序は i 外側・j 内側。"""
    a, b = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    return a.ravel().astype(np.int64), b.ravel().astype(np.int64)


def adjacency(ea, eb):
    """走行本体 run_subset_stage123_v1.py の adjacency と同じ規則：頂点を共有すれば 1、それ以外 0、対角 0。"""
    M = len(ea)
    A = np.zeros((M, M), dtype=np.float64)
    for e in range(M):
        share = (ea == ea[e]) | (ea == eb[e]) | (eb == ea[e]) | (eb == eb[e])
        share[e] = False
        A[e, share] = 1.0
    return A


def H_of(z, A):
    """走行本体の H_of と同一（式17）。"""
    H = A * (np.conj(z)[:, None] * z[None, :]); np.fill_diagonal(H, 0.0); return H.astype(np.complex128, copy=False)


def K_of(theta, A):
    """式18 の K = Im Ĥ（実反対称、M×M）。走行本体の one_step が使う 1j*np.imag(H) の K と同一。"""
    return np.imag(H_of(np.exp(1j * theta), A))


def eigenmode(theta, A):
    """iK（エルミート）の最大固有値 σ_max の固有ベクトル（単位ノルム）と全固有値。"""
    w, V = np.linalg.eigh(1j * K_of(theta, A))
    idx = int(np.argmax(w))
    v = V[:, idx].astype(np.complex128)
    return v / np.linalg.norm(v), w


def eigenmode_residual(v, A):
    """旧 _eigenmode_residual と同じ：μ = Re v†(iKv)、r = ‖iKv − μv‖。K は arg v から作る。"""
    kv = K_of(np.angle(v), A) @ v
    mu = float(np.real(np.conj(v) @ (1j * kv)))
    return float(np.linalg.norm(1j * kv - mu * v))


def sigma_spectrum(theta, A):
    w = np.linalg.eigvalsh(1j * K_of(theta, A))
    return np.sort(w[w > 1e-12])[::-1]


def make_parent_dense(A, rng, iters=400, beta=0.5, tol=1e-8, restarts=3, progress=None):
    best = (None, np.inf, None)
    m = A.shape[0]
    for _ in range(restarts):
        theta = rng.uniform(0.0, 2.0 * np.pi, m)
        v = None
        for it in range(iters):
            v, _ = eigenmode(theta, A)
            theta_new = np.angle(v)
            mix = (1.0 - beta) * np.exp(1j * theta) + beta * np.exp(1j * theta_new)
            theta = np.angle(mix)
            if it % 10 == 9:
                res_now = eigenmode_residual(v, A)
                if progress:
                    progress(f"親構成 iter={it+1} 残差={res_now:.2e}")
                if res_now < tol:
                    break
        residual = eigenmode_residual(v, A)
        if residual < best[1]:
            best = (v, residual, sigma_spectrum(np.angle(v), A))
        if residual < tol:
            break
    return best


def zero_closure_seed_dense(theta, A, rng):
    """K(θ) の値域の直交補内の実正規直交対 (u, w) から g = (u + i w)/√2。gᵀg = 0 厳密、g ⊥ range(K)。"""
    K = K_of(theta, A)
    U, s, _ = np.linalg.svd(K)
    Q = U[:, s > 1e-10 * s[0]]
    m = A.shape[0]

    def project_out(g):
        return g - Q @ (Q.T @ g)

    u = project_out(rng.normal(size=m)); u = u / np.linalg.norm(u)
    w = project_out(rng.normal(size=m)); w = w - (w @ u) * u; w = w / np.linalg.norm(w)
    return (u + 1j * w) / math.sqrt(2.0), Q


def build_parent(N, edges_fn, seed=0, delta=1e-15, tol=1e-12, iters=1200, progress=None):
    ea, eb = edges_fn(N)
    A = adjacency(ea, eb)
    rng = np.random.default_rng(SEED_BASE + 1000 * N + seed)
    t0 = time.time()
    v, residual, sig = make_parent_dense(A, rng, iters=iters, tol=tol, progress=progress)
    g, Q = zero_closure_seed_dense(np.angle(v), A, rng)
    Z = v + delta * g
    Z = Z / np.linalg.norm(Z)
    return dict(v=v, g=g, Z0=Z, sigma=sig, residual=np.float64(residual), n=np.int64(N), M=np.int64(len(ea)),
                seed=np.int64(seed), delta=np.float64(delta), tol=np.float64(tol), iters=np.int64(iters),
                rank_K=np.int64(Q.shape[1]), elapsed_s=np.float64(time.time() - t0), ea=ea, eb=eb)


if __name__ == '__main__':
    # 使い方: python parents_dense_v1.py {triu|all} OUTDIR N1 N2 ...
    structure, outdir = sys.argv[1], sys.argv[2]
    Ns = [int(x) for x in sys.argv[3:]]
    fn = {'triu': edges_triu, 'all': edges_all}[structure]
    os.makedirs(outdir, exist_ok=True)
    for N in Ns:
        d = build_parent(N, fn, progress=lambda msg: print(msg, file=sys.stderr, flush=True))
        path = os.path.join(outdir, f'parent_dense_{structure}_N{N:05d}_20261025.npz')
        np.savez_compressed(path, structure=structure, **d)
        status = 'ok' if d['residual'] < 1e-8 else 'NOT_CONVERGED'
        print(f"N={N}: saved M={int(d['M'])} residual={float(d['residual']):.3e} planes={len(d['sigma'])} rank_K={int(d['rank_K'])} "
              f"{float(d['elapsed_s']):.1f}s {status}", flush=True)
    print('DENSE PARENTS DONE')
