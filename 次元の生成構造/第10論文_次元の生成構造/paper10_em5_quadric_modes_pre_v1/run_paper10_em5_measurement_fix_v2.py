#!/usr/bin/env python3
"""E-M5 補足v2：測定系の訂正（基準軌道並走差分）——判定基準は v1 と同一

v1 の測定欠陥（診断確定・訂正対象、予言と閾値は不変）:
    (a) Stage 2 のモード座標を「静止参照点 Zs との差」で取ったため、基準点
        自身の永年ドリフト（残差 1.6e-4〜3.8e-4/step オーダー）が混入し、
        unwrap 勾配がほぼゼロ周波数を返した。証拠: Z4 の FFT 基本波
        0.0218 rad/step は線形固有周波数 0.02175 と一致するが、unwrap 測定は
        3.4e-6 を返した（モード座標振幅 1.25e-2 ≫ 励起 1e-3 = ドリフト支配）。
    (b) 「周波数対 7 本（期待 8）」は角度>0 フィルタが実固有値（角度≈0）を
        落としたことによる勘定漏れの可能性——実固有値の本数を明示的に報告する。
    (c) 残差を下げるため緩和を 60000 step に延長（測定点の質の改善であって
        判定基準の変更ではない）。

訂正:
    Stage 2 は基準軌道 Z_b(t)（無摂動）と摂動軌道 Z_p(t) を並走させ、
    差分 δ(t) = frame·(Z_p(t) − Z_b(t)) からモード座標を取る。
    P2/P3/P4 の閾値・判定式は v1 と同一（TOL_P2=1e-8, TOL_P3=1e-6,
    TOL_LOCK=1e-3, 非自明ロック=比≥2）。
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
V1 = HERE / "run_paper10_em5_quadric_modes_pre_v1.py"
spec = importlib.util.spec_from_file_location("em5v1", V1)
v1 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = v1
spec.loader.exec_module(v1)

T_RELAX2 = 60000
T_LONG = v1.T_LONG
EPS_EXC = v1.EPS_EXC


def stage1_v2(A, Z0, label):
    Z = Z0.copy()
    for _ in range(T_RELAX2):
        Z = v1.step(A, Z)
    Zs = Z
    FZ = v1.step(A, Zs)
    phi_star = float(np.angle(np.conj(Zs) @ FZ))
    rot = np.exp(-1j * phi_star)
    residual = float(np.linalg.norm(rot * FZ - Zs))
    m = len(Zs)

    def G_real(xr):
        Zi = xr[:m] + 1j * xr[m:]
        Zo = rot * v1.step(A, Zi)
        return np.concatenate([Zo.real, Zo.imag])

    x0 = np.concatenate([Zs.real, Zs.imag])
    J = np.zeros((2 * m, 2 * m))
    for k in range(2 * m):
        e = np.zeros(2 * m); e[k] = v1.H_FD
        J[:, k] = (G_real(x0 + e) - G_real(x0 - e)) / (2 * v1.H_FD)
    B = v1.tangent_basis(Zs)
    Mred = B.T @ J @ B
    lam, V = np.linalg.eig(Mred)
    order = np.argsort(-np.abs(np.angle(lam)))
    lam, V = lam[order], V[:, order]
    max_dev_circle = float(np.max(np.abs(np.abs(lam) - 1.0)))
    n_real = int(np.sum(np.abs(np.angle(lam)) < 1e-9))
    pos = np.sort(np.abs(np.angle(lam))[np.angle(lam) > 0])
    print(f"  [{label}] Stage1v2: 残差={residual:.3e} max||λ|-1|={max_dev_circle:.3e} "
          f"周波数対={len(pos)}本＋実固有値={n_real}本（期待 対{m-2}）")
    return {"Zs": Zs, "phi_star": phi_star, "residual": residual,
            "lam": lam, "V": V, "B": B, "max_dev_circle": max_dev_circle,
            "freq_pairs": pos, "n_real_eigs": n_real}


def stage2_v2(A, s1, label):
    Zs, phi_star, B, V = s1["Zs"], s1["phi_star"], s1["B"], s1["V"]
    m = len(Zs)
    rot = np.exp(-1j * phi_star)
    exc = B @ np.ones(B.shape[1]); exc = exc / np.linalg.norm(exc)
    Zp = Zs + EPS_EXC * (exc[:m] + 1j * exc[m:])
    Zp = v1.closure_project(Zp)

    n_modes = B.shape[1]
    coords = np.zeros((T_LONG, n_modes), dtype=complex)
    Zb, Zq = Zs.copy(), Zp.copy()
    frame = 1.0 + 0.0j
    for t in range(T_LONG):
        Zb = v1.step(A, Zb)
        Zq = v1.step(A, Zq)
        frame *= rot
        dz = frame * (Zq - Zb)                    # ← 訂正: 基準軌道並走差分
        y = B.T @ np.concatenate([dz.real, dz.imag])
        coords[t] = np.linalg.solve(V, y)

    lam = s1["lam"]
    idx_pos = [i for i in range(len(lam)) if np.angle(lam[i]) > 0]
    freqs, amps = [], []
    for i in idx_pos:
        a = coords[:, i]
        amps.append(float(np.mean(np.abs(a))))
        u = np.unwrap(np.angle(a[max(1, T_LONG - 3000):]))
        freqs.append(abs(float(np.polyfit(np.arange(len(u)), u, 1)[0])))
    fb = np.asarray(freqs)
    fb_valid = fb[fb > 1e-8]
    max_dev, locks, ratio_matrix = v1.em4_lock_stats(fb_valid)
    print(f"  [{label}] Stage2v2: モード周波数={np.array2string(fb, precision=6)}")
    print(f"  [{label}] Stage2v2: 最大比ずれ={max_dev:.3e} 非自明ロック={locks}")
    out = {"mode_freqs_rad_per_step": fb.tolist(), "mode_mean_amps": amps,
           "freq_max_ratio_dev": max_dev, "nontrivial_locks": locks,
           "ratio_matrix": ratio_matrix,
           "linear_freq_pairs": s1["freq_pairs"].tolist()}

    if n_modes == 2:
        sig = coords[:, idx_pos[0]].real
        w = np.hanning(T_LONG)
        F = np.abs(np.fft.rfft(sig * w))
        fgrid = np.fft.rfftfreq(T_LONG)
        peaks = []
        for k in range(2, len(F) - 2):
            if F[k] > F[k - 1] and F[k] > F[k + 1] and F[k] > 1e-9 * F.max():
                d = 0.5 * (F[k - 1] - F[k + 1]) / (F[k - 1] - 2 * F[k] + F[k + 1])
                peaks.append((fgrid[k] + d * (fgrid[1] - fgrid[0]), F[k]))
        peaks.sort(key=lambda x: -x[1])
        top = sorted([p[0] for p in peaks[:3]])
        if len(top) >= 3:
            sp = np.diff(top)
            ladder_dev = float(abs(sp[1] - sp[0]))
        elif len(top) == 2:
            ladder_dev = float(abs(top[1] - 2 * top[0]))
        else:
            ladder_dev = float("nan")
        out["ladder_peaks_cyc_per_step"] = top
        out["ladder_spacing_dev"] = ladder_dev
        print(f"  [{label}] P3v2 ラダー: ピーク={top} 間隔偏差={ladder_dev:.3e}")
    return out


def main() -> None:
    t0 = time.time()
    print("E-M5 補足v2（測定訂正: 並走差分・緩和60000）実行")
    per_g = {}
    for g, label in ((4, "Z4"), (6, "Z6")):
        print(f"\n[g={g}] N={g-1}, M={(g-1)*(g-2)//2}")
        sys_lr, A, Z0, m = v1.lock_config(g)
        s1 = stage1_v2(A, Z0, label)
        rcs = v1.rank_check(s1["Zs"])
        print(f"  [{label}] Stage0(準安定点v2): 物理接次元={rcs['dim_physical_tangent']}"
              f"（期待 {rcs['expected']}）")
        s2 = stage2_v2(A, s1, label)
        per_g[label] = {"m": m, "stage0_meta": rcs,
                        "stage1": {"phi_star": s1["phi_star"], "residual": s1["residual"],
                                    "max_dev_circle": s1["max_dev_circle"],
                                    "freq_pairs": s1["freq_pairs"].tolist(),
                                    "n_freq_pairs": len(s1["freq_pairs"]),
                                    "n_real_eigs": s1["n_real_eigs"],
                                    "expected_pairs": m - 2},
                        "stage2": s2}

    z4, z6 = per_g["Z4"], per_g["Z6"]
    p1 = (z4["stage0_meta"]["dim_physical_tangent"] == 2
          and z6["stage0_meta"]["dim_physical_tangent"] == 16)
    p2 = (z4["stage1"]["max_dev_circle"] < v1.TOL_P2
          and z6["stage1"]["max_dev_circle"] < v1.TOL_P2
          and z6["stage1"]["n_freq_pairs"] == 8
          and z4["stage1"]["n_freq_pairs"] == 1)
    ladder = z4["stage2"].get("ladder_spacing_dev", float("nan"))
    p3 = bool(np.isfinite(ladder) and ladder < v1.TOL_P3)
    p4 = z6["stage2"]["nontrivial_locks"] > 0

    print("\n==== 判定（補足v2・基準は v1 と同一） ====")
    print(f"P1 次元則: {'PASS' if p1 else 'FAIL'}")
    print(f"P2 楕円性(<1e-8, 対1/8): {'PASS' if p2 else 'FAIL'} "
          f"(Z4:{z4['stage1']['max_dev_circle']:.2e}, Z6:{z6['stage1']['max_dev_circle']:.2e}, "
          f"Z6対={z6['stage1']['n_freq_pairs']}+実{z6['stage1']['n_real_eigs']})")
    print(f"P3 等間隔ラダー(<1e-6): {'PASS' if p3 else 'FAIL'} (偏差 {ladder:.3e})")
    print(f"P4 モードレベル整数ロック（v2ゲート）: {'PASS' if p4 else 'FAIL'} "
          f"(非自明ロック {z6['stage2']['nontrivial_locks']})")

    payload = {"experiment": "paper10_em5_measurement_fix_v2",
               "fix": "基準軌道並走差分・緩和60000（判定基準はv1と同一）",
               "params": {"T_RELAX2": T_RELAX2, "T_LONG": T_LONG, "EPS_EXC": EPS_EXC},
               "per_g": per_g,
               "verdicts": {"P1": bool(p1), "P2": bool(p2), "P3": bool(p3), "P4": bool(p4)},
               "runtime_sec": time.time() - t0}
    (HERE / "paper10_em5_result_fix_v2.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({payload['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
