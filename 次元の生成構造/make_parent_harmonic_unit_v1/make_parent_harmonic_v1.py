#!/usr/bin/env python3
"""make_parent（倍音対応）単体 v1 — 設計書 ../make_parent(倍音対応).md に従う。

倍音レジスタ拡張状態 Z ∈ C^{M×H} の初期の海を、段ごとの自己無撞着
円偏波閉包（現行 make_parent と同一の不動点反復、枝選択のみパラメータ化）
の等振幅重ね合わせとして生成する。

- 閉塞 Σ_{e,n} Z² = 0 は構造から恒等（射影・補正なし）
- スペクトル分布は等振幅（無名等振幅公理）。n 依存重みは導入しない
- 段の周波数 n·ω₀ はレジスタのメタデータ（構造のみ。力学は範囲外）
- 【禁止】インフレーションプログラム（abl.build_init 系)への組込みは指示があるまで行わない

反復の複製元: run_n_scaling_lowrank_v1.py 158-190行（make_parent）。
枝選択 argmin(Im ev)（'−'）を argmax（'+'）にも開いた点のみが差分。
"""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
ENGINE_FILE = (REPO / "時間軸Q軸とフェルミオンの生成構造" / "検証_対照実験"
               / "第5論文原本_自発的分裂予備実験_v1" / "run_n_scaling_lowrank_v1.py")

_spec = importlib.util.spec_from_file_location("eng_mph", ENGINE_FILE)
eng = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = eng
_spec.loader.exec_module(eng)

ENGINE_SHA256 = hashlib.sha256(ENGINE_FILE.read_bytes()).hexdigest()


def _eigenmode_mu_residual(sys_lr, v):
    """μ = Re v†(iKv) と残差 ‖iKv − μv‖（エンジン 151-155 行と同一の定義）。"""
    kv = sys_lr.kmatvec(v)
    mu = float(np.real(np.conj(v) @ (1j * kv)))
    res = float(np.linalg.norm(1j * kv - mu * v))
    return mu, res


def self_consistent_mode(n_vertices, rng, branch, iters=1200, beta=0.5,
                          tol=1e-12, restarts=3):
    """現行 make_parent（エンジン 158-190 行）の複製。枝選択のみパラメータ化。

    branch: '-' → argmin(Im ev)（現行 make_parent と同一）
            '+' → argmax(Im ev)（鏡像枝）
    戻り値: (v, residual, sigma1, mu, restarts_used)
    """
    sys_lr = eng.LowRankSystem(n_vertices)
    pick = np.argmin if branch == "-" else np.argmax
    best = (None, np.inf)
    used = 0
    for r in range(restarts):
        used = r + 1
        theta = rng.uniform(0.0, 2.0 * np.pi, sys_lr.m)
        v = None
        for it in range(iters):
            sys_lr.set_theta(theta)
            ev, EV = np.linalg.eig(sys_lr.J @ sys_lr.G)
            idx = int(pick(ev.imag))
            v = sys_lr.w(EV[:, idx].astype(complex))
            v = v / np.linalg.norm(v)
            theta_new = np.angle(v)
            mix = (1.0 - beta) * np.exp(1j * theta) + beta * np.exp(1j * theta_new)
            theta = np.angle(mix)
            if it % 10 == 9:
                sys_lr.set_theta(np.angle(v))
                _, res_now = _eigenmode_mu_residual(sys_lr, v)
                if res_now < tol:
                    break
        sys_lr.set_theta(np.angle(v))
        mu, residual = _eigenmode_mu_residual(sys_lr, v)
        if residual < best[1]:
            best = (v, residual)
        if residual < tol:
            break
    v, residual = best
    sys_lr.set_theta(np.angle(v))
    mu, residual = _eigenmode_mu_residual(sys_lr, v)
    sigma1 = float(sys_lr.sigma_spectrum()[0])
    return v, residual, sigma1, mu, used


def top_plane_occupancy(n_vertices, v):
    """v の、自段 K(angle(v)) の選択固有値近傍平面への占有率。"""
    sys_lr = eng.LowRankSystem(n_vertices)
    sys_lr.set_theta(np.angle(v))
    ev, EV = np.linalg.eig(sys_lr.J @ sys_lr.G)
    sig1 = float(np.max(np.abs(ev.imag)))
    cols = [i for i in range(len(ev)) if abs(abs(ev[i].imag) - sig1) < 1e-9 * max(sig1, 1.0)]
    U = np.column_stack([sys_lr.w(EV[:, i].astype(complex)) for i in cols])
    Q, _ = np.linalg.qr(U)
    proj = Q @ (np.conj(Q.T) @ v)
    return float(np.real(np.conj(proj) @ proj) / np.real(np.conj(v) @ v))


def make_parent_harmonic(n_vertices, H, seed, iters=1200, beta=0.5,
                          tol=1e-12, restarts=3, force_branch=None):
    """設計書 §2 の生成器。

    force_branch: None → 段ごとに rng 硬貨投げ（規約、消費順は硬貨→uniform）
                  '+'/'-' → 全段強制（T8 鏡像テスト用）
    戻り値: (Z, info)  Z: (M, H) complex, ‖Z‖_F = 1
    """
    rng = np.random.default_rng(seed)
    m = n_vertices * (n_vertices - 1) // 2
    Z = np.zeros((m, H), dtype=complex)
    levels = []
    for h in range(1, H + 1):
        coin = int(rng.integers(0, 2))          # 段ごとの硬貨は常に消費（消費契約を一定に保つ）
        branch = force_branch if force_branch is not None else ("-" if coin == 0 else "+")
        v, residual, sigma1, mu, used = self_consistent_mode(
            n_vertices, rng, branch, iters=iters, beta=beta, tol=tol, restarts=restarts)
        Z[:, h - 1] = v / np.sqrt(H)
        x, y = v.real, v.imag
        levels.append({
            "n": h, "branch": branch, "mu": mu, "sigma1": sigma1,
            "residual": residual, "restarts_used": used,
            "closure_abs": abs(complex(v @ v)),
            "min_amp": float(np.min(np.abs(v))),
            "rank_re_im": int(np.linalg.matrix_rank(np.column_stack([x, y]), tol=1e-10)),
            "norm_re": float(np.linalg.norm(x)), "norm_im": float(np.linalg.norm(y)),
            "re_dot_im": float(x @ y),
            "occupancy_top_plane": top_plane_occupancy(n_vertices, v),
        })
    info = {"n_vertices": n_vertices, "H": H, "seed": seed,
            "iters": iters, "beta": beta, "tol": tol, "restarts": restarts,
            "force_branch": force_branch,
            "engine_sha256": ENGINE_SHA256,
            "total_closure_abs": abs(complex(np.sum(Z * Z))),
            "frobenius_norm": float(np.linalg.norm(Z)),
            "level_amp": [float(np.linalg.norm(Z[:, k])) for k in range(H)],
            "levels": levels}
    return Z, info
