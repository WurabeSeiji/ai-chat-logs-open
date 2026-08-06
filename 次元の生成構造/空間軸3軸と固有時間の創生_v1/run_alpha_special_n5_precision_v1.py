#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""α特殊解精密検定 v1: N=5 の z相関は有理番地 sin²(23π/124) と厳密一致するか

背景（2026-08-06 木原指摘の再定式化）:
  - α恒等式は無理数ではなく整数比番地: 1−R = sin²(23π/124) = 0.302822（U^248=I の第23根）
  - N=5 は結晶学的制限（許容位数1,2,3,4,6）の外にある最小禁止位数=特殊解。
    第二α根は5倍細分レジスタ（620=5×124）に住む——5の接続。
  - よって掃引のN→∞収束は誤った検定。正しい検定は N=5 単独の精密測定。

方法: N=5 を T=22000 まで走行し、準安定窓 [2000,22000] を 2000 ステップ×10 窓に
分割。各窓で独立に固有対→軸→復調読出し→z相関を測定し、平均±標準誤差を出す。
（固有対と軸は窓ごとに再推定＝系統誤差も窓間分散に反映される。）

判定（事前固定）:
  番地格子: A_m = sin²(mπ/124)。隣接番地間隔は 23 近傍で約 0.012（4%）。
  H_exact23: |mean − A_23| ≤ 2×SE かつ A_23 が最近傍番地 → α番地 23 と同定
  H_other: 別の番地 A_m が最近傍かつ 2SE 内 → その番地を報告（m を記録）
  H_none: どの番地とも 2SE で不整合、または窓間で非定常 → 有理番地仮説を棄却

使い方: python3 run_alpha_special_n5_precision_v1.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec1 = importlib.util.spec_from_file_location("pre1_n5p", HERE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl
edge_adjacency, build_K = pre1.edge_adjacency, pre1.build_K

N = 5
T_END = 22000
WIN_LEN = 2000
WIN_START = 2000
SAMPLE_EVERY = 5

A23 = float(np.sin(23 * np.pi / 124) ** 2)      # 0.302822 = √(4πα) 番地
C23 = 1.0 - A23


def analyze_window(S, p, q, sys_lr, wp_fix, adj):
    """1窓ぶんのサンプル列から z相関・overlap を測る（pre1 P2-P4 と同一手順）。"""
    M = S.shape[1]
    ns = S.shape[0]

    def perp(Zc):
        return Zc - p * (p @ Zc) - q * (q @ Zc)
    Sp = np.array([perp(z) for z in S])
    X = np.hstack([Sp.real, Sp.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    U3 = Vt[:3].T

    Zstar = S[-1].copy()

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

    cp = S @ p
    cq = S @ q
    phi = np.unwrap(np.angle(cp + 1j * cq))
    omega_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])
    d1c = d1[:M] + 1j * d1[M:]; d1c /= np.linalg.norm(d1c)
    d2c = d2[:M] + 1j * d2[M:]; d2c /= np.linalg.norm(d2c)
    ax20 = U3 @ axis
    axc = ax20[:M] + 1j * ax20[M:]; axc /= np.linalg.norm(axc)
    dirs = {"d1": d1c, "d2": d2c, "axis": axc}
    cser = {k: Sp @ np.conj(d) for k, d in dirs.items()}
    freqs = {}
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
    return {"corr_xy": float(Cn[0, 1]), "corr_zx": float(Cn[2, 0]),
            "corr_zy": float(Cn[2, 1]), "overlap": best_ov,
            "omega_clock": omega_clock}


def main() -> None:
    t0 = time.time()
    print(f"番地格子: A_23=sin²(23π/124)={A23:.6f}  cos²={C23:.6f}")
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(N, True)
    adj = edge_adjacency(N)

    windows = []
    samples = []
    win_id = 0
    results = []
    for t in range(T_END):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if t >= WIN_START and (t % SAMPLE_EVERY == 0):
            samples.append(Z.copy())
        if t >= WIN_START and (t - WIN_START + 1) % WIN_LEN == 0 and samples:
            S = np.array(samples)
            r = analyze_window(S, p, q, sys_lr, wp.copy(), adj)
            r["window"] = [WIN_START + win_id * WIN_LEN, WIN_START + (win_id + 1) * WIN_LEN]
            zc = 0.5 * (r["corr_zx"] + r["corr_zy"])
            r["zcorr"] = zc
            results.append(r)
            print(f"窓{win_id:2d} {r['window']}: z相関={zc:.5f}  overlap={r['overlap']:.4f}  "
                  f"xy={r['corr_xy']:.4f}")
            samples = []
            win_id += 1

    zc = np.array([r["zcorr"] for r in results])
    ov = np.array([r["overlap"] for r in results])
    zc_mean, zc_se = float(zc.mean()), float(zc.std(ddof=1) / np.sqrt(len(zc)))
    ov_mean, ov_se = float(ov.mean()), float(ov.std(ddof=1) / np.sqrt(len(ov)))

    # 最近傍番地
    ms = np.arange(1, 62)
    grid = np.sin(ms * np.pi / 124) ** 2
    m_best = int(ms[np.argmin(np.abs(grid - zc_mean))])
    A_best = float(np.sin(m_best * np.pi / 124) ** 2)
    dev23 = zc_mean - A23
    h_exact23 = bool(m_best == 23 and abs(dev23) <= 2 * zc_se)
    h_other = bool(m_best != 23 and abs(zc_mean - A_best) <= 2 * zc_se)

    print("\n==== 集計 ====")
    print(f"z相関 = {zc_mean:.5f} ± {zc_se:.5f} (SE, n={len(zc)})")
    print(f"  対 A_23={A23:.6f}: 偏差 {dev23:+.5f} ({dev23/zc_se:+.1f}σ)")
    print(f"  最近傍番地 m={m_best}: A_{m_best}={A_best:.6f} "
          f"(偏差 {zc_mean-A_best:+.5f}, {(zc_mean-A_best)/zc_se:+.1f}σ)")
    print(f"overlap = {ov_mean:.5f} ± {ov_se:.5f}  対 cos²(23π/124)={C23:.6f}: "
          f"偏差 {ov_mean-C23:+.5f}")
    verdict = "H_exact23" if h_exact23 else ("H_other" if h_other else "H_none")
    print(f"判定: {verdict}")

    out = {"N": N, "T_END": T_END, "WIN_LEN": WIN_LEN, "SAMPLE_EVERY": SAMPLE_EVERY,
           "A23": A23, "C23": C23, "windows": results,
           "zcorr_mean": zc_mean, "zcorr_se": zc_se,
           "overlap_mean": ov_mean, "overlap_se": ov_se,
           "m_best": m_best, "A_best": A_best,
           "dev23": dev23, "dev23_sigma": dev23 / zc_se if zc_se > 0 else None,
           "verdict": verdict, "runtime_sec": time.time() - t0}
    (HERE / "alpha_special_n5_precision_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → alpha_special_n5_precision_result_v1.json")


if __name__ == "__main__":
    main()
