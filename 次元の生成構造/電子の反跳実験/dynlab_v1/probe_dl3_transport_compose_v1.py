#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 追走行 v1 — 輸送の直接検定: 毎步合成輸送 vs 10步直接輸送（判定 M3 の決着）

背景: Davis–Kahan は上界定理であり、ρ>1/2 から「追跡不能」は導けない（査読指摘・正当）。
「毎步必須」を実測事実にするには輸送誤差を直接測る。あわせて M4 の物理検定
（O(3) 最適整列が反射枝を選ぶイベントの有無）と、部分空間ジャンプの計数を行う。

方法（実行前固定）: 物質相後期（15000步以降）で、長さ40窓×10步。
  各窓で、始端枠 V3(n0) から
    合成輸送: 毎步 Procrustes 整列を10回連鎖 → F_comp
    直接輸送: V3(n0+10) を V3(n0) に一回で整列 → F_dir
  同一部分空間上の両者の食い違い回転 R_d = F_dir^T F_comp について
    θ_disc = 回転角（det=+1 のとき）、det(R_d)=−1 は枝不一致イベント
  事前登録判定: θ_disc > 30° または det=−1 の窓が有意（>10%）に存在すれば
  「10步直接輸送は枝を誤る＝毎步輸送が必要」が実測で成立。
  併記: 毎步整列の O(3) 最適解が det=−1 となる步の割合（M4 物理検定）、
        隣接步の上位3部分空間重なり < 0.9 の步（部分空間ジャンプ）の計数。

出力: result_dl3_transport_compose_v1.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNI = HERE.parent.parent / "統一万能関数_v1"
N, DELTA = 16, 0.1
T_SKIP = 15000
N_WIN, W = 40, 10


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    u1 = _load("uni_dl3tc", UNI / "unified_interaction_v1.py")
    eng, _, _ = u1.build_standard_universe(N, DELTA)
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ones = np.ones(N) / np.sqrt(N)

    def V3_now():
        x = eng.C2().sum(axis=(1, 2))
        D2 = np.zeros((N, N))
        D2[ii, jj] = D2[jj, ii] = np.abs(x) ** 2
        B = -0.5 * Jc @ D2 @ Jc
        lam, V = np.linalg.eigh(B)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        lamk, Vk = lam[keep], V[:, keep]
        o = np.argsort(lamk)[::-1]
        return Vk[:, o][:, :3]

    for _ in range(T_SKIP):
        eng.step()

    theta_disc, det_disc = [], []
    o3_reflect_steps = 0
    subspace_jumps = 0
    total_steps = 0
    for _ in range(N_WIN):
        V0 = V3_now()
        F = V0.copy()
        Vp = V0
        for _ in range(W):
            eng.step()
            V = V3_now()
            total_steps += 1
            subspace_jumps += int(np.sum((Vp.T @ V) ** 2) / 3.0 < 0.9)
            U_, S_, Vt_ = np.linalg.svd(F.T @ V)
            R = U_ @ Vt_
            o3_reflect_steps += int(np.linalg.det(R) < 0)
            F = V @ R.T          # 毎步合成輸送
            Vp = V
        V10 = Vp
        U_, S_, Vt_ = np.linalg.svd(V0.T @ V10)
        Fdir = V10 @ (U_ @ Vt_).T   # 10步直接輸送
        Rd = Fdir.T @ F
        d = float(np.linalg.det(Rd))
        det_disc.append(d)
        if d > 0:
            c = (np.trace(Rd) - 1.0) / 2.0
            theta_disc.append(float(np.degrees(np.arccos(np.clip(c, -1, 1)))))
        else:
            theta_disc.append(float("nan"))  # 枝不一致（回転角は定義しない）

    theta = np.array(theta_disc)
    dets = np.array(det_disc)
    n_branch = int(np.sum(dets < 0))
    n_theta30 = int(np.nansum(theta > 30.0))
    n_bad = n_branch + n_theta30
    res = {
        "config": {"N": N, "delta": DELTA, "engine": "unified_interaction_v1",
                   "t_skip": T_SKIP, "n_windows": N_WIN, "window_len": W,
                   "note": "run_dl23_matter_v1 と同一の決定的軌道の再走行"},
        "windows_branch_mismatch": n_branch,
        "windows_theta_gt30": n_theta30,
        "windows_bad_total": n_bad,
        "frac_bad": n_bad / N_WIN,
        "theta_disc_deg": {"median": float(np.nanmedian(theta)),
                           "max": float(np.nanmax(theta))},
        "per_step_o3_reflect_frac": o3_reflect_steps / total_steps,
        "per_step_subspace_jump_frac": subspace_jumps / total_steps,
        "verdict_every_step_required": bool(n_bad / N_WIN > 0.10),
        "elapsed_sec": time.time() - t0,
    }
    np.savez_compressed(HERE / "dl3_transport_compose_v1.npz",
                        theta_disc_deg=theta, det_disc=dets)
    (HERE / "result_dl3_transport_compose_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"窓{N_WIN}本: 枝不一致={n_branch}  θ>30°={n_theta30}  不良率={n_bad/N_WIN:.2f}")
    print(f"θ_disc 中央値={np.nanmedian(theta):.1f}°  最大={np.nanmax(theta):.1f}°")
    print(f"毎步O(3)反射選好率={o3_reflect_steps/total_steps:.3f}  "
          f"部分空間ジャンプ率={subspace_jumps/total_steps:.3f}")
    print(f"毎步輸送必須（実測）: {res['verdict_every_step_required']} "
          f"({res['elapsed_sec']:.0f}s)")


if __name__ == "__main__":
    main()
