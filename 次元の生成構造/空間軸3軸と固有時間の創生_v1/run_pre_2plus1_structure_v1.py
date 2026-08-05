#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 P1–P4: 三方向準安定構造の 2＋1 仮説（探索的・判定は記述）

仮説（設計方針_v1.md）: 実体=独立固有ベクトル2本（非直交可）＋従属な回転軸1本。
読出し=固有時計との内積が正規直交3軸を製造する。

P1 占有構造: 準安定窓 [2000,4000] の軌道SVD（親平面除去後）——有意方向の本数。
P2 独立性と非直交性: 一段写像の接線行列の固有分解——占有部分空間に重なる
   固有方向の構成（複素対=平面か実固有=軸か）、物理方向対のGram。
P3 回転軸の従属性: 生成子 K_ef=sin(θ_f−θ_e)（隣接辺・反対称）を占有3D部分空間へ
   射影した反対称3×3の回転軸 â と、上位2方向の外積 n̂ の一致度。回転符号も記録。
P4 読出し直交化: 親平面位相を固有時計とし、時計内積（復調）読出し座標の
   相関行列を、実体側Gramと対比する。

使い方: python3 run_pre_2plus1_structure_v1.py
"""
from __future__ import annotations

import importlib.util
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
ABL = (REPO / "次元の生成構造" / "第8論文_二段階seed除去による準安定相の因果分離"
       / "code" / "run_preliminary_seed_ablation_v1.py")
spec = importlib.util.spec_from_file_location("abl_sp3", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)

N = 5
T_END = 4000
WIN = (2000, 4000)
SAMPLE_EVERY = 5


def edge_adjacency(n):
    edges = list(itertools.combinations(range(n), 2))
    m = len(edges)
    adj = np.zeros((m, m), bool)
    for i, (a, b) in enumerate(edges):
        for j, (c, d) in enumerate(edges):
            if i != j and len({a, b} & {c, d}) == 1:
                adj[i, j] = True
    return adj


def build_K(theta, adj):
    """生成子 K_ef = sin(θ_f − θ_e)（隣接のみ・反対称）。"""
    K = np.sin(theta[None, :] - theta[:, None]) * adj
    return K


def main() -> None:
    t0 = time.time()
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = abl.build_init(N, True)
    M = sys_lr.m
    adj = edge_adjacency(N)

    # ---- 準安定窓まで走行、サンプル収集 ----
    samples = []
    for t in range(T_END):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if WIN[0] <= t < WIN[1] and (t % SAMPLE_EVERY == 0):
            samples.append(Z.copy())
    S = np.array(samples)                      # (ns, M) complex
    ns = S.shape[0]
    print(f"走行完了 T={T_END}, サンプル {ns} 本（窓 {WIN}）")

    # 親平面（複素射影）を除去
    def perp(Zc):
        return Zc - p * (p @ Zc) - q * (q @ Zc)
    Sp = np.array([perp(z) for z in S])
    pow_perp = float(np.mean(np.sum(np.abs(Sp) ** 2, axis=1)))
    pow_tot = float(np.mean(np.sum(np.abs(S) ** 2, axis=1)))
    print(f"親平面外パワー比 = {pow_perp/pow_tot:.4f}")

    # ---- P1: 実20次元でのSVD占有構造 ----
    X = np.hstack([Sp.real, Sp.imag])          # (ns, 2M)
    Xc = X - X.mean(axis=0)
    U_, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    sv_rel = sv / sv[0]
    n_sig = int(np.sum(sv_rel > 0.05))
    print("[P1] 特異値（上位8, /σ1）:", np.round(sv_rel[:8], 4).tolist())
    print(f"     有意方向（>5%）= {n_sig} 本")
    U3 = Vt[:3].T                               # (2M, 3) 占有正規直交基底

    # ---- P2: 接線写像の固有分解（実20×20・有限差分） ----
    Zstar = S[-1].copy()
    wp_fix = wp.copy()
    def step_real(x):
        Zc = x[:M] + 1j * x[M:]
        Zn, _ = abl.evolve(sys_lr, Zc.copy(), wp_fix.copy())
        return np.concatenate([Zn.real, Zn.imag])
    x0 = np.concatenate([Zstar.real, Zstar.imag])
    f0 = step_real(x0)
    eps = 1e-7
    J = np.zeros((2 * M, 2 * M))
    for k in range(2 * M):
        dx = np.zeros(2 * M); dx[k] = eps
        J[:, k] = (step_real(x0 + dx) - f0) / eps
    ev, EV = np.linalg.eig(J)
    # 占有部分空間に重なる固有方向を抽出
    rows = []
    for i in range(len(ev)):
        vec = EV[:, i]
        ov = float(np.linalg.norm(U3.T @ vec.real) ** 2 + np.linalg.norm(U3.T @ vec.imag) ** 2) \
             / float(np.linalg.norm(vec.real) ** 2 + np.linalg.norm(vec.imag) ** 2 + 1e-300)
        rows.append((i, complex(ev[i]), ov))
    rows.sort(key=lambda r: -r[2])
    print("[P2] 占有3D部分空間への重なり上位6固有値:")
    picked = []
    for i, lam, ov in rows[:6]:
        print(f"     λ={lam.real:+.6f}{lam.imag:+.6f}i |λ|={abs(lam):.6f} 重なり={ov:.3f}")
        picked.append((i, lam, ov))
    # 物理方向: 重なり最大の複素固有対 → 平面 (d1,d2)=Re,Im
    top = picked[0]
    vec = EV[:, top[0]]
    d1 = vec.real / np.linalg.norm(vec.real)
    d2 = vec.imag / np.linalg.norm(vec.imag)
    g12 = float(d1 @ d2)
    print(f"     物理方向対（最上位複素固有対のRe/Im）: Gram off-diag d1·d2 = {g12:+.4f}"
          f"（非直交なら≠0）")

    # ---- P3: 生成子射影の回転軸と外積の一致 ----
    theta = np.angle(Zstar)
    K = build_K(theta, adj)                    # M×M 反対称
    K20 = np.zeros((2 * M, 2 * M))
    K20[:M, :M] = K
    K20[M:, M:] = K
    A3 = U3.T @ K20 @ U3
    A3 = 0.5 * (A3 - A3.T)                     # 反対称化（数値）
    w_vec = np.array([A3[2, 1], A3[0, 2], A3[1, 0]])   # A3 x = w × x
    omega = float(np.linalg.norm(w_vec))
    axis = w_vec / max(omega, 1e-300)
    # 物理2方向を占有座標に落として外積
    c1 = U3.T @ d1; c1 /= np.linalg.norm(c1)
    c2 = U3.T @ d2; c2 /= np.linalg.norm(c2)
    n_vec = np.cross(c1, c2)
    n_hat = n_vec / max(np.linalg.norm(n_vec), 1e-300)
    align = float(abs(n_hat @ axis))
    sign = float(np.sign(n_hat @ axis) * np.sign(omega))
    print(f"[P3] 生成子射影の回転レート |ω|={omega:.4e}")
    print(f"     回転軸 â と 外積 n̂=d1×d2 の一致度 |n̂·â| = {align:.4f}")
    print(f"     符号（n̂·â の向き）= {sign:+.0f}")

    # ---- P4: 固有時計との内積読出し ----
    cp = S @ p
    cq = S @ q
    phi = np.unwrap(np.angle(cp + 1j * cq))    # 親平面位相 = 固有時計
    omega_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])
    # 実体方向への複素射影（d1,d2,軸方向）
    d1c = d1[:M] + 1j * d1[M:]; d1c /= np.linalg.norm(d1c)
    d2c = d2[:M] + 1j * d2[M:]; d2c /= np.linalg.norm(d2c)
    ax20 = U3 @ axis
    axc = ax20[:M] + 1j * ax20[M:]; axc /= np.linalg.norm(axc)
    dirs = {"d1": d1c, "d2": d2c, "axis": axc}
    cser = {k: Sp @ np.conj(d) for k, d in dirs.items()}
    # 各座標の支配周波数（時計単位）
    freqs = {}
    for k, c in cser.items():
        F = np.fft.fft(c - c.mean())
        pk = int(np.argmax(np.abs(F)))
        f_per_step = pk / ns if pk <= ns // 2 else (pk - ns) / ns
        freqs[k] = float(2 * np.pi * f_per_step / (omega_clock * SAMPLE_EVERY / SAMPLE_EVERY))
    print(f"[P4] 固有時計レート ω_clock={omega_clock:.4f}/サンプル  "
          f"各座標の支配周波数/時計 = " + ", ".join(f"{k}:{v:+.3f}" for k, v in freqs.items()))
    # 生Gram（実体側・実20次元）
    Draw = np.stack([d1, d2, ax20], axis=1)
    Graw = Draw.T @ Draw
    # 時計復調読出し: r_k(t) = c_k(t) e^{-i n_k φ(t)}（n_k=支配周波数の整数丸め）
    Cread = np.zeros((3, 3), complex)
    keys = list(cser.keys())
    rser = {}
    for k in keys:
        n_k = round(freqs[k])
        rser[k] = cser[k] * np.exp(-1j * n_k * phi)
    for i, ki in enumerate(keys):
        for j, kj in enumerate(keys):
            Cread[i, j] = np.mean(rser[ki] * np.conj(rser[kj]))
    Dn = np.sqrt(np.real(np.diag(Cread)))
    Cn = np.abs(Cread) / np.outer(Dn, Dn)
    print("     実体側Gram（|off-diag|最大）=",
          round(float(np.max(np.abs(Graw - np.diag(np.diag(Graw))))), 4))
    print("     読出し相関行列 |C_ij| =")
    for i in range(3):
        print("       ", np.round(Cn[i], 4).tolist())
    offmax = float(np.max(np.abs(Cn - np.diag(np.diag(Cn)))))
    print(f"     読出し off-diag 最大 = {offmax:.4f}")

    out = {"N": N, "T_END": T_END, "WIN": list(WIN), "n_samples": ns,
           "perp_power_ratio": pow_perp / pow_tot,
           "P1": {"sv_rel_top8": [float(x) for x in sv_rel[:8]], "n_sig_5pct": n_sig},
           "P2": {"top_eigs": [{"lam_re": l.real, "lam_im": l.imag, "abs": abs(l), "overlap": o}
                                for _, l, o in picked],
                   "gram_d1d2": g12},
           "P3": {"omega_gen": omega, "align_axis_cross": align, "sign": sign},
           "P4": {"omega_clock": omega_clock, "dominant_freq_per_clock": freqs,
                   "graw_offdiag_max": float(np.max(np.abs(Graw - np.diag(np.diag(Graw))))),
                   "readout_corr_abs": [[float(Cn[i, j]) for j in range(3)] for i in range(3)],
                   "readout_offdiag_max": offmax},
           "runtime_sec": time.time() - t0}
    (HERE / "pre_2plus1_structure_result_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
