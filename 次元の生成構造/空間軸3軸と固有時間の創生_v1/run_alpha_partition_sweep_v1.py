#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""α分配掃引 v1: 読出しの 0.7/0.3 分配と α 恒等式の cos²θ/sin²θ の照合

動機（2026-08-06 木原指摘）:
  α接点論文の恒等式 1−R = sin²θ = √(4πα) = e_HL = 0.302811（U^n=I ロック角）
  に対し、本論文の未言及の実測 2 件が近い:
    - 回転固有モードの占有3D部分空間への重なり（可読率）≈ 0.70
    - z軸読出しと x,y 読出しの相関 = 0.2986/0.2987
  偏差は逆符号でほぼ同大（±1.4%）。N掃引で 0.6972/0.3028 への収束を検定する。

判定（事前固定・記述的）:
  H_α: N↑ で overlap → cos²θ_α=0.697189, z相関 → sin²θ_α=0.302811（|偏差|単調減少）
  H_finite: どちらかが 0.5 や 0 へ流れる／単調でない → 有限N効果・α非関連
  * どちらでも結果は記録する（仮説→検証→棄却も記録の規約）。

使い方: python3 run_alpha_partition_sweep_v1.py
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
spec1 = importlib.util.spec_from_file_location("pre1_aps", HERE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl
edge_adjacency, build_K = pre1.edge_adjacency, pre1.build_K

T_END = 4000
WIN = (2000, 4000)
SAMPLE_EVERY = 5

ALPHA_INV = 137.035999084
E_HL = float(np.sqrt(4 * np.pi / ALPHA_INV))      # √(4πα) = sin²θ_α
COS2 = 1.0 - E_HL                                  # cos²θ_α


def analyze(n: int) -> dict:
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    M = sys_lr.m
    adj = edge_adjacency(n)
    samples = []
    for t in range(T_END):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if WIN[0] <= t < WIN[1] and (t % SAMPLE_EVERY == 0):
            samples.append(Z.copy())
    S = np.array(samples)
    ns = S.shape[0]

    def perp(Zc):
        return Zc - p * (p @ Zc) - q * (q @ Zc)
    Sp = np.array([perp(z) for z in S])

    X = np.hstack([Sp.real, Sp.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    U3 = Vt[:3].T

    # 接線写像と回転固有対（重なり=可読率）
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
    best_i, best_ov = -1, -1.0
    for i in range(len(ev)):
        if ev[i].imag <= 0:
            continue
        vec = EV[:, i]
        ov = (np.linalg.norm(U3.T @ vec.real) ** 2 + np.linalg.norm(U3.T @ vec.imag) ** 2) \
             / (np.linalg.norm(vec.real) ** 2 + np.linalg.norm(vec.imag) ** 2 + 1e-300)
        if ov > best_ov:
            best_ov, best_i = float(ov), i

    vec = EV[:, best_i]
    d1 = vec.real / np.linalg.norm(vec.real)
    d2 = vec.imag / np.linalg.norm(vec.imag)

    # 回転軸（生成子射影）と外積一致
    theta = np.angle(Zstar)
    K = build_K(theta, adj)
    K20 = np.zeros((2 * M, 2 * M))
    K20[:M, :M] = K
    K20[M:, M:] = K
    A3 = U3.T @ K20 @ U3
    A3 = 0.5 * (A3 - A3.T)
    w_vec = np.array([A3[2, 1], A3[0, 2], A3[1, 0]])
    omega = float(np.linalg.norm(w_vec))
    axis = w_vec / max(omega, 1e-300)
    c1 = U3.T @ d1; c1 /= np.linalg.norm(c1)
    c2 = U3.T @ d2; c2 /= np.linalg.norm(c2)
    n_vec = np.cross(c1, c2)
    n_hat = n_vec / max(np.linalg.norm(n_vec), 1e-300)
    align = float(abs(n_hat @ axis))

    # 固有時計と復調読出し相関（pre1 P4 と同一手順）
    cp = S @ p
    cq = S @ q
    phi = np.unwrap(np.angle(cp + 1j * cq))
    d1c = d1[:M] + 1j * d1[M:]; d1c /= np.linalg.norm(d1c)
    d2c = d2[:M] + 1j * d2[M:]; d2c /= np.linalg.norm(d2c)
    ax20 = U3 @ axis
    axc = ax20[:M] + 1j * ax20[M:]; axc /= np.linalg.norm(axc)
    dirs = {"d1": d1c, "d2": d2c, "axis": axc}
    cser = {k: Sp @ np.conj(d) for k, d in dirs.items()}
    freqs = {}
    omega_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])
    for k, c in cser.items():
        F = np.fft.fft(c - c.mean())
        pk = int(np.argmax(np.abs(F)))
        f_per_step = pk / ns if pk <= ns // 2 else (pk - ns) / ns
        freqs[k] = float(2 * np.pi * f_per_step / omega_clock)
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

    return {"N": n, "M": M, "overlap": best_ov,
            "lam": {"re": float(ev[best_i].real), "im": float(ev[best_i].imag),
                     "abs": float(abs(ev[best_i]))},
            "align_axis_cross": align,
            "corr_xy": float(Cn[0, 1]), "corr_zx": float(Cn[2, 0]),
            "corr_zy": float(Cn[2, 1]),
            "dev_overlap_vs_cos2": float(best_ov - COS2),
            "dev_zcorr_vs_sin2": float(0.5 * (Cn[2, 0] + Cn[2, 1]) - E_HL)}


def main() -> None:
    t0 = time.time()
    print(f"α定数: √(4πα)=sin²θ_α={E_HL:.6f}  cos²θ_α={COS2:.6f}  (α⁻¹={ALPHA_INV})")
    out = {"alpha_inv": ALPHA_INV, "sin2_alpha": E_HL, "cos2_alpha": COS2,
           "T_END": T_END, "WIN": list(WIN), "SAMPLE_EVERY": SAMPLE_EVERY,
           "sweep": []}
    for n in (4, 5, 6, 7, 8, 9, 10, 11, 12):
        t1 = time.time()
        r = analyze(n)
        r["runtime_sec"] = time.time() - t1
        out["sweep"].append(r)
        zc = 0.5 * (r["corr_zx"] + r["corr_zy"])
        print(f"N={n:2d} M={r['M']:2d}: overlap={r['overlap']:.4f} "
              f"(vs cos²θ_α {r['dev_overlap_vs_cos2']:+.4f})  "
              f"z相関={zc:.4f} (vs sin²θ_α {r['dev_zcorr_vs_sin2']:+.4f})  "
              f"xy相関={r['corr_xy']:.4f}  align={r['align_axis_cross']:.4f}  "
              f"[{r['runtime_sec']:.0f}s]")
        out["runtime_sec"] = time.time() - t0
        (HERE / "alpha_partition_sweep_result_v1.json").write_text(
            json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → alpha_partition_sweep_result_v1.json")


if __name__ == "__main__":
    main()
