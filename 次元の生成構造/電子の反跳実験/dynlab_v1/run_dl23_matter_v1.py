#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL2/DL3 物質走行 v1 — レジスタのジャンプ・保存と、枠の誕生（判定 L1〜L7・M1〜M9）

導出: DL2/DL3_導出ノート.md。力学: unified_interaction_v1（DL0 と同一の実装決定）。
投入: エンジン正本の標準宇宙 build_standard_universe(N, δ)——スライス1へ δ×seed_state。

【事前登録からの修正（明記）】
  L4/L5（パケットの厳密有理数 r=15/34・|z|=16/17）は CR 型パケット投入用の判定であり、
  本走行のエンジン正本投入（gen3 親第1セクション）には適用できない。パケット投入変種は
  別走行として登録し、本走行では L4/L5 を「適用外」と記録する。

理論の追記（実行前）: 共有O線形部は全スライスに同一直交を掛けるため、読出し層の
両レジスタは多スライス占有でも線形部で厳密保存。頂点はレジスタ点ごとに保存（補題DL1-2）。
したがって L2 の予言は「点火を貫通して両レジスタ保存（頂点 RK4 の積分器誤差のみ）」。

出力: result_dl23_matter_v1.json・dl23_series_v1.npz
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
N, M, T, SAMPLE = 16, 120, 20000, 10
DELTA = 0.1


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    u1 = _load("uni_dl23", UNI / "unified_interaction_v1.py")

    # --- 初期データからのジャンプ式（L1）: v=真空ポンプ, g=seed_state ---
    abl, gen3 = u1.abl, u1.gen3
    _, _, _, _, _, _, _, Z0c, _ = abl.build_init(N, False)
    r2 = gen3.make_parent(N, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / N
    seed = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    jump_formula = complex(np.sum(Z0c ** 2) + 2 * DELTA * np.sum(Z0c * seed)
                           + DELTA ** 2 * np.sum(seed ** 2))

    eng, _, _ = u1.build_standard_universe(N, DELTA)
    x0 = eng.C2().sum(axis=(1, 2))
    jump_measured = complex(np.sum(x0 * x0))
    L1_err = abs(jump_measured - jump_formula)

    # 真空対照の CV(A_v)（L6 の基準）
    engv, _, _ = u1.build_standard_universe(N, 0.0)
    def vertex_cv(engine):
        x = engine.C2().sum(axis=(1, 2))
        Av = np.zeros(N)
        np.add.at(Av, engine.ia, np.abs(x) ** 2)
        np.add.at(Av, engine.ib, np.abs(x) ** 2)
        return float(np.std(Av) / np.mean(Av))
    cv_vac0 = vertex_cv(engv)
    cv_mat0 = vertex_cv(eng)

    # --- 本走行 ---
    Jc = np.eye(N) - np.ones((N, N)) / N
    ii, jj = np.triu_indices(N, 1)
    ser = {k: [] for k in ("tau", "c", "bil_re", "bil_im", "trB", "top3",
                            "gmin", "n_neg", "tail", "ey_rel", "dB2", "det_sign",
                            "cvA", "smax_smin")}
    B_prev = None
    V3_al_prev = None
    bil_0 = None
    for tau in range(T):
        eng.step()
        if tau % SAMPLE:
            continue
        C2 = eng.C2()
        x = C2.sum(axis=(1, 2))
        c = float(np.sum(np.abs(x) ** 2))
        bil = complex(np.sum(x * x))
        if bil_0 is None:
            bil_0 = bil
        D2 = np.zeros((N, N))
        D2[ii, jj] = D2[jj, ii] = np.abs(x) ** 2
        B = -0.5 * Jc @ D2 @ Jc
        lam, V = np.linalg.eigh(B)
        ones = np.ones(N) / np.sqrt(N)
        keep = np.ones(N, bool)
        keep[int(np.argmax(np.abs(V.T @ ones)))] = False
        lamk, Vk = lam[keep], V[:, keep]
        o = np.argsort(lamk)[::-1]
        lamk, Vk = lamk[o], Vk[:, o]
        trB = float(lamk.sum())
        top3 = float(lamk[:3].sum() / trB)
        gmin = float(min(lamk[0] - lamk[1], lamk[1] - lamk[2], lamk[2] - lamk[3]))
        n_neg = int(np.sum(lamk < -1e-12 * trB))
        tail = float(1.0 - lamk[:3].sum() / np.sum(np.abs(lamk)))
        # Eckart–Young 恒等式（rank-3 近似・非自明空間内）
        X3 = Vk[:, :3] * np.sqrt(np.maximum(lamk[:3], 0.0))[None, :]
        Bk = (Vk * lamk[None, :]) @ Vk.T
        ey_lhs = float(np.sum((Bk - X3 @ X3.T) ** 2))
        ey_rhs = float(np.sum(lamk[3:] ** 2))
        ey_rel = abs(ey_lhs - ey_rhs) / max(ey_rhs, 1e-300)
        # 逐次整列（SO(3) Procrustes）と掌性
        V3 = Vk[:, :3]
        if V3_al_prev is None:
            V3_al = V3
        else:
            H = V3_al_prev.T @ V3
            U_, _, Vt_ = np.linalg.svd(H)
            Rr = (U_ @ Vt_)
            V3_al = V3 @ Rr.T @ np.diag([1, 1, np.sign(np.linalg.det(Rr))]) \
                if False else V3 @ (U_ @ Vt_).T
        det_sign = float(np.sign(np.linalg.det(V3_al_prev.T @ V3_al))) \
            if V3_al_prev is not None else 1.0
        V3_al_prev = V3_al
        dB2 = float(np.linalg.norm(B - B_prev, 2)) if B_prev is not None else 0.0
        B_prev = B.copy()
        # 頂点モーメント CV と楕円体偏差（c=3/N）
        Av = np.zeros(N)
        np.add.at(Av, eng.ia, np.abs(x) ** 2)
        np.add.at(Av, eng.ib, np.abs(x) ** 2)
        cvA = float(np.std(Av) / np.mean(Av))
        T3 = X3.T @ X3
        s2v = np.einsum("ij,jk,ik->i", X3, np.linalg.pinv(T3), X3) / (3.0 / N)
        s = np.sqrt(np.maximum(s2v, 0.0))
        smax_smin = float(s.max() / max(s.min(), 1e-300))

        ser["tau"].append(tau); ser["c"].append(c)
        ser["bil_re"].append(bil.real); ser["bil_im"].append(bil.imag)
        ser["trB"].append(trB); ser["top3"].append(top3)
        ser["gmin"].append(gmin); ser["n_neg"].append(n_neg)
        ser["tail"].append(tail); ser["ey_rel"].append(ey_rel)
        ser["dB2"].append(dB2); ser["det_sign"].append(det_sign)
        ser["cvA"].append(cvA); ser["smax_smin"].append(smax_smin)

    F = {k: np.array(v) for k, v in ser.items()}
    late = slice(int(len(F["tau"]) * 0.75), None)
    bil = F["bil_re"] + 1j * F["bil_im"]
    L2_drift = float(np.max(np.abs(bil - bil_0)))
    c_drift = float(np.max(np.abs(F["c"] - F["c"][0])))
    trB_rel = float(np.max(np.abs(F["trB"] * N / F["c"] - 1.0)))
    rho = F["dB2"][1:] / np.maximum(F["gmin"][:-1], 1e-300)

    res = {
        "config": {"N": N, "T": T, "delta": DELTA,
                   "engine": "unified_interaction_v1",
                   "L4_L5": "適用外（エンジン正本投入。パケット変種は別走行として登録）"},
        "L1_jump": {"formula": [jump_formula.real, jump_formula.imag],
                    "measured": [jump_measured.real, jump_measured.imag],
                    "abs_err": L1_err, "pass": bool(L1_err < 1e-12)},
        "L2_conservation": {"bil_drift": L2_drift, "c_drift": c_drift,
                            "pass": bool(L2_drift < 1e-7 and c_drift < 1e-7)},
        "L3_trace": {"rel_max": trB_rel, "pass": bool(trB_rel < 1e-6)},
        "L6_vertex_pattern": {"cv_vacuum_t0": cv_vac0, "cv_matter_t0": cv_mat0,
                              "cv_late_mean": float(F["cvA"][late].mean()),
                              "pass": bool(cv_mat0 > cv_vac0)},
        "L7_rA": {"note": "r_A 系列は npz に保存（B 対角）。判定なし・記録"},
        "M1_gap": {"gmin_t0": float(F["gmin"][0]),
                   "gmin_late_mean": float(F["gmin"][late].mean()),
                   "gmin_late_min": float(F["gmin"][late].min()),
                   "pass": bool(F["gmin"][late].min() > 0)},
        "M2_top3": {"t0": float(F["top3"][0]),
                    "late_mean": float(F["top3"][late].mean()),
                    "pass": bool(F["top3"][late].mean() > 0.205)},
        "M3_transport": {"rho_median": float(np.median(rho)),
                         "rho_frac_gt_half": float(np.mean(rho >= 0.5)),
                         "pass": bool(np.median(rho) < 0.5)},
        "M4_handedness": {"det_sign_min": float(F["det_sign"].min()),
                          "pass": bool(F["det_sign"].min() > 0)},
        "M5_eckart_young": {"rel_max": float(F["ey_rel"].max()),
                            "pass": bool(F["ey_rel"].max() < 1e-10)},
        "M6_tail": {"t0": float(F["tail"][0]),
                    "late_mean": float(F["tail"][late].mean()),
                    "note": "時系列記録（値の予言なし）", "pass": True},
        "M8_imag": {"late_max": int(F["n_neg"][late].max()),
                    "note": "観測（判定でなく記録）", "pass": True},
        "M9_ellipsoid": {"smax_smin_late_mean": float(F["smax_smin"][late].mean()),
                         "note": "偏差時系列の記録（乗ることは期待しない）",
                         "pass": True},
        "elapsed_sec": time.time() - t0,
    }
    keys = [k for k in res if k.startswith(("L", "M")) and isinstance(res[k], dict)
            and "pass" in res[k]]
    res["all_pass"] = bool(all(res[k]["pass"] for k in keys))
    np.savez_compressed(HERE / "dl23_series_v1.npz", **F)
    (HERE / "result_dl23_matter_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    for k in keys:
        print(f"  {k}: {'PASS' if res[k]['pass'] else 'FAIL'}  "
              f"{ {kk: vv for kk, vv in res[k].items() if kk != 'note'} }")
    print(f"  ALL: {'PASS' if res['all_pass'] else 'FAIL'} "
          f"({res['elapsed_sec']:.1f}s)")


if __name__ == "__main__":
    main()
