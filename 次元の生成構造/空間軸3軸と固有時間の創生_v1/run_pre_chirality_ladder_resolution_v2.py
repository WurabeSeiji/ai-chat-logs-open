#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 v2: キラリティ・平面梯子・分解能掃引・quadrature読出し初点灯

v1 の発見（回転軸=外積 0.9874・複素固有対1つ=平面・特異値の対構造）を受けた
4系統の検証。ソフト期待（事前記録・探索的、判定は記述）:

C1 キラリティ: 初期状態の複素共役（時間反転＝逆枝）で、P3の回転符号が反転する。
   一致度の大きさは不変。
C2 梯子: 上位の平面対の回転周波数の比（対時計）を実測。調和（1,2,3,…）か
   低調波（1,1/2,…）かを記述。
C3 分解能: 次元の曖昧さ 1−|n̂·â| は N とともに減少する（木原予言:
   低分解能では次元さえ不確定を含む。0.9874 は N=5 の分解能の帰結）。
C4 quadrature読出し: x=Re, y=Im（固有時計復調の直交位相成分）, z=軸射影で
   読み出すと、状態は (x,y) 平面上のゆっくり歳差する円として現れる（記述）。

使い方: python3 run_pre_chirality_ladder_resolution_v2.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
spec1 = importlib.util.spec_from_file_location("pre1", HERE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl
edge_adjacency, build_K = pre1.edge_adjacency, pre1.build_K

T_END = 4000
WIN = (2000, 4000)
SAMPLE_EVERY = 5


def run_and_analyze(n, conjugate=False):
    """走行→占有SVD→接線固有対→回転軸/外積→梯子周波数→quadrature。"""
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    M = sys_lr.m
    if conjugate:
        Z = np.conj(Z)
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
    sv_rel = (sv / sv[0]).tolist()
    U3 = Vt[:3].T

    # 接線写像
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
            best_ov, best_i = ov, i
    lam = complex(ev[best_i])
    vec = EV[:, best_i]
    d1 = vec.real / np.linalg.norm(vec.real)
    d2 = vec.imag / np.linalg.norm(vec.imag)

    # 生成子射影の回転軸と外積
    theta = np.angle(Zstar)
    K = build_K(theta, adj)
    K20 = np.zeros((2 * M, 2 * M))
    K20[:M, :M] = K; K20[M:, M:] = K
    A3 = U3.T @ K20 @ U3
    A3 = 0.5 * (A3 - A3.T)
    w_vec = np.array([A3[2, 1], A3[0, 2], A3[1, 0]])
    omega = float(np.linalg.norm(w_vec))
    axis = w_vec / max(omega, 1e-300)
    c1 = U3.T @ d1; c1 /= np.linalg.norm(c1)
    c2 = U3.T @ d2; c2 /= np.linalg.norm(c2)
    n_hat = np.cross(c1, c2)
    n_hat /= max(np.linalg.norm(n_hat), 1e-300)
    align = float(abs(n_hat @ axis))
    signed = float(n_hat @ axis)

    # 固有時計
    cp = S @ p; cq = S @ q
    phi = np.unwrap(np.angle(cp + 1j * cq))
    om_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])

    # 平面梯子: SVD対 (0,1),(2,3),(4,5),(6,7) の複素射影の支配周波数
    ladder = []
    for k in range(0, 8, 2):
        if k + 1 >= Vt.shape[0]:
            break
        u_a, u_b = Vt[k], Vt[k + 1]
        ca = Xc @ u_a; cb = Xc @ u_b
        c = ca + 1j * cb
        F = np.fft.fft(c - c.mean())
        pk = int(np.argmax(np.abs(F)))
        f = pk / ns if pk <= ns // 2 else (pk - ns) / ns
        w_plane = 2 * np.pi * f
        ladder.append({"pair": k // 2 + 1, "sv_rel": sv_rel[k],
                        "omega_per_sample": float(abs(w_plane)),
                        "ratio_to_clock": float(abs(w_plane) / abs(om_clock))})

    # quadrature読出し（C4）: 平面複素方向へ射影→時計復調→x,y,z
    dc = (d1[:M] + 1j * d1[M:])
    dc = dc / np.linalg.norm(dc)
    cser = Sp @ np.conj(dc)
    r = cser * np.exp(-1j * phi)
    ax20 = U3 @ axis
    axc = (ax20[:M] + 1j * ax20[M:]); axc /= np.linalg.norm(axc)
    zc = Sp @ np.conj(axc)
    zr = zc * np.exp(-1j * phi)
    prec = np.unwrap(np.angle(r))
    prec_rate = float(np.polyfit(np.arange(ns), prec, 1)[0])
    quad = {"r_mean_abs": float(np.mean(np.abs(r))),
            "r_std_abs": float(np.std(np.abs(r))),
            "precession_per_sample": prec_rate,
            "precession_over_clock": prec_rate / om_clock,
            "z_mean_abs": float(np.mean(np.abs(zr)))}

    return {"N": n, "M": M, "conjugate": conjugate,
            "lam": {"re": lam.real, "im": lam.imag, "abs": abs(lam)},
            "overlap": float(best_ov), "sv_rel_top8": sv_rel[:8],
            "omega_gen": omega, "align": align, "signed_align": signed,
            "omega_clock_per_sample": om_clock,
            "ladder": ladder, "quad": quad}


def main() -> None:
    t0 = time.time()
    out = {}

    # ---- C1 キラリティ: N=5 正枝 vs 共役 ----
    print("=== C1 キラリティ検定（N=5, 正 vs 共役初期状態） ===")
    base = run_and_analyze(5, conjugate=False)
    conj = run_and_analyze(5, conjugate=True)
    print(f"  正:   符号つき一致 {base['signed_align']:+.4f}  時計 {base['omega_clock_per_sample']:+.4f}")
    print(f"  共役: 符号つき一致 {conj['signed_align']:+.4f}  時計 {conj['omega_clock_per_sample']:+.4f}")
    flip = np.sign(base["signed_align"]) != np.sign(conj["signed_align"])
    mag_keep = abs(abs(base["align"]) - abs(conj["align"])) < 0.05
    print(f"  → 符号反転={flip}  大きさ保存={mag_keep}")
    out["C1"] = {"base": base, "conj": conj, "sign_flip": bool(flip), "magnitude_kept": bool(mag_keep)}

    # ---- C2 梯子の周波数比（N=5 正枝の結果から） ----
    print("=== C2 平面梯子（N=5）: 各平面の回転数/時計 ===")
    for row in base["ladder"]:
        print(f"  平面{row['pair']}: σ/σ1={row['sv_rel']:.3f} ω/ω_clock={row['ratio_to_clock']:.4f}")
    out["C2"] = base["ladder"]

    # ---- C3 分解能掃引 ----
    print("=== C3 分解能掃引: 次元の曖昧さ 1−|n̂·â| vs N ===")
    res = [{"N": 5, "align": base["align"], "misalign": 1 - base["align"],
            "overlap": base["overlap"]}]
    for n in (4, 6, 8):
        r = run_and_analyze(n, conjugate=False)
        res.append({"N": n, "align": r["align"], "misalign": 1 - r["align"],
                     "overlap": r["overlap"]})
        print(f"  N={n}: |n̂·â|={r['align']:.4f}  曖昧さ={1-r['align']:.4f}  固有対重なり={r['overlap']:.3f}")
        out.setdefault("C3_runs", {})[str(n)] = r
    res.sort(key=lambda x: x["N"])
    print("  掃引まとめ:", [(x["N"], round(x["misalign"], 4)) for x in res])
    mono = all(res[i]["misalign"] >= res[i + 1]["misalign"] - 1e-6 for i in range(len(res) - 1))
    print(f"  単調減少（曖昧さ↓ with N↑）= {mono}")
    out["C3"] = {"sweep": res, "monotone_decrease": bool(mono)}

    # ---- C4 quadrature読出し（N=5） ----
    q = base["quad"]
    print("=== C4 quadrature読出し（N=5） ===")
    print(f"  |r| 平均={q['r_mean_abs']:.4f} ± {q['r_std_abs']:.4f}（円軌道なら std ≪ mean）")
    print(f"  歳差レート/時計 = {q['precession_over_clock']:+.4f}（復調残り＝離調）")
    print(f"  軸方向振幅 |z| 平均 = {q['z_mean_abs']:.4f}")
    out["C4"] = q

    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_chirality_ladder_resolution_result_v2.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
