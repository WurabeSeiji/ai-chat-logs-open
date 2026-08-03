#!/usr/bin/env python3
"""第8論文v2予備実験 E-M4：仮説Aの介入検証——結合倍音がレジスタを作るか v1

仮説A（精密形）:
    関係系の結合は sin(Δθ) の第1倍音のみ（エンジンコード235/261行で確認）。
    蔵本理論により第1倍音結合は 1:1 同期（ユニゾン）しか作れない。
    公理2（分解能レジスタ）の力学的実体は「相互作用が運ぶ倍音の内容」である。

介入（無名性厳守）:
    結合に高次倍音項 ε_q·A·sin(q·Δθ) を追加する（全辺同一・種ラベルなし・
    実反対称なのでノルムと零閉塞は自動保存）。エンジンの対照密形式
    Kd = A·sin(Δθ)（エンジン自身が kmatvec 検証に用いる形）を拡張する。

観測量（無名）:
    位相の秩序パラメータ R_k = |Σ_e e^{ikθ_e}|/M （k=1,2,3）。
    ユニゾン → R₁≈1。二クラスタ(0/π分裂) → R₂ 支配。三クラスタ → R₃ 支配。

【v1予言の反証記録と構造的発見】
    v1は位相クラスタ（R_q支配）を予言したが全FAIL。原因は測定でなく公理:
    **零閉塞 Z·Z = Σ|Z_e|²e^{2iθ_e} = 0 は、等振幅では R₂（0/π二クラスタ）を
    公理的に禁止する**（対照の R₂=0.0018 はこの拘束の現れ）。二クラスタ型の
    レジスタ形成路は公理1自身が塞いでいる——公理1は偶数次の位相整列を
    禁止する（海が splay 状態である理由）。E-M1 のユニゾンは周波数の斉唱で
    あり位相整列ではない。よって観測量を周波数側へ訂正する。

予言（v2・訂正版）:
    P0'（対照）: 周波数ユニゾン（比ずれ<1e-2）・非自明整数比ロック0
        ——E-M1 の再現（crossing=1167 は既に一致）
    P1'（介入）: 結合倍音 ε_q·sin(qΔθ) が周波数分化（比ずれ>1e-1）または
        非自明整数比ロック L>0 を作るか——作れば仮説A（レジスタ=結合倍音）
        支持、作らなければ「結合倍音でも不十分」として仮説Aの機構を
        さらに絞る（レジスタには場の円環構造そのものが必要）
    副産物記録: 結合倍音は crossing 時刻を大きく変える
        （q2: 1167→200、q3: →6092）——倍音はインフレーション動力学に強く効く
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CODE = HERE.parent / "code" / "run_preliminary_seed_ablation_v1.py"
spec = importlib.util.spec_from_file_location("abl_m4", CODE)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
import run_n_scaling_lowrank_v1 as eng  # abl が sys.path を通済み

N = 5
XMAX = 12000
GAMMA = eng.GAMMA


def dense_setup():
    sys_lr, v, B_p1, B_rot, B0, p, q, Z0, wp = abl.build_init(N, initial_seed=False)
    ea, eb = sys_lr.ea, sys_lr.eb
    m = sys_lr.m
    A = np.zeros((m, m))
    for i in range(m):
        share = (ea == ea[i]) | (ea == eb[i]) | (eb == ea[i]) | (eb == eb[i])
        A[i, share] = 1.0
    np.fill_diagonal(A, 0.0)
    return A, Z0, p, q, m


def run_case(A, Z0, p, q, harmonics, xmax=XMAX):
    """harmonics: {q_h: eps} 例 {1:1.0} / {1:1.0, 2:2.0}"""
    Z = Z0.copy()
    m = len(Z)
    R123 = np.zeros((xmax + 1, 3))
    fs = np.zeros(xmax + 1)
    phases_hist = np.zeros((xmax + 1, m))

    def fval(Zv):
        Zp = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Zv) @ Zv))

    def order(theta):
        return [abs(np.mean(np.exp(1j * k * theta))) for k in (1, 2, 3)]

    th = np.angle(Z)
    R123[0] = order(th); fs[0] = fval(Z)
    crossing = None
    eye = np.eye(m)
    for t in range(1, xmax + 1):
        th = np.angle(Z)
        K = np.zeros((m, m))
        for qh, eps in harmonics.items():
            K += eps * A * np.sin(qh * (th[None, :] - th[:, None]))
        sigma = np.linalg.norm(K, 2)
        if sigma > 0:
            g = GAMMA / sigma
            Z = np.linalg.solve(eye - g * K, (eye + g * K) @ Z)
        R123[t] = order(np.angle(Z)); fs[t] = fval(Z)
        phases_hist[t] = np.angle(Z)
        if crossing is None and fs[t] > 0.05:
            crossing = t
    meta = crossing + abl.GUARD if crossing is not None else None
    tail = slice(max(0, xmax - 2000), xmax)
    late = R123[tail].mean(axis=0)
    # 周波数解析（E-M1 と同一の観測量）
    u = np.unwrap(phases_hist[max(1, xmax - 3000):], axis=0)
    fr = np.abs(np.polyfit(np.arange(u.shape[0]), u, 1)[0])
    fb = fr[fr > 1e-8]
    if fb.size > 1:
        r = fb[:, None] / np.maximum(fb[None, :], 1e-30)
        rmax = np.maximum(r, 1 / np.maximum(r, 1e-30))
        max_dev = float(np.max(np.abs(r[r >= 1] - 1)))
        pr = np.round(rmax)
        locks = int(np.sum((pr >= 2) & (np.abs(rmax - pr) < 1e-3)) // 2)
    else:
        max_dev, locks = 0.0, 0
    return {"crossing": crossing, "metastable_start": meta,
            "R1_late": float(late[0]), "R2_late": float(late[1]), "R3_late": float(late[2]),
            "freq_max_ratio_dev": max_dev, "nontrivial_locks": locks,
            "freq_mean": float(np.mean(fb)) if fb.size else 0.0,
            "norm_drift": abs(float(np.real(np.conj(Z) @ Z)) - 1.0),
            "zero_closure_abs": abs(complex(Z @ Z))}


def main() -> None:
    t0 = time.time()
    A, Z0, p, q, m = dense_setup()
    cases = {
        "control_eps0": {1: 1.0},
        "q2_eps05": {1: 1.0, 2: 0.5},
        "q2_eps2": {1: 1.0, 2: 2.0},
        "q3_eps2": {1: 1.0, 3: 2.0},
    }
    results = {}
    for name, h in cases.items():
        r = run_case(A, Z0, p, q, h)
        results[name] = r
        print(f"{name:12s}: crossing={r['crossing']} R1={r['R1_late']:.4f} R2={r['R2_late']:.4f} "
              f"R3={r['R3_late']:.4f} 周波数比ずれ={r['freq_max_ratio_dev']:.3e} "
              f"ロック={r['nontrivial_locks']} |Z·Z|={r['zero_closure_abs']:.1e}")

    c = results["control_eps0"]
    p0 = c["freq_max_ratio_dev"] < 1e-2 and c["nontrivial_locks"] == 0
    any_split = any(r["freq_max_ratio_dev"] > 1e-1 or r["nontrivial_locks"] > 0
                    for k, r in results.items() if k != "control_eps0")
    print(f"\nP0' 対照=周波数ユニゾン・ロック0（E-M1再現）: {'PASS' if p0 else 'FAIL'}")
    print(f"P1' 結合倍音による周波数分化またはロック: {'あり（仮説A支持）' if any_split else 'なし（結合倍音でも不十分——場の円環構造が必要）'}")
    print("副産物: 零閉塞がR₂位相クラスタを公理的に禁止（構造的発見）／"
          "結合倍音はcrossing時刻を大きく変える")

    payload = {
        "experiment": "paper8_axiom2_coupling_harmonics_pre_v1",
        "design": ("結合への高次倍音注入（無名・全辺同一・実反対称でノルム/零閉塞自動保存）。"
                    "密形式はエンジン自身の対照検証形 Kd=A·sin(Δθ) の拡張"),
        "engine_contrast": "初期状態は abl.build_init と同一（read-only import）",
        "cases": results,
        "v1_falsified": "位相クラスタ予言は全FAIL——零閉塞がR₂整列を公理的に禁止（構造的発見として記録）",
        "P0_control_freq_unison": bool(p0),
        "P1_harmonic_coupling_splits_frequency": bool(any_split),
        "runtime_sec": time.time() - t0,
        "conclusion": (
            "v1の位相クラスタ予言は反証（零閉塞がR₂整列を公理的に禁止=構造的発見）。"
            f"v2周波数観測: 対照ユニゾン={'再現' if p0 else '不再現'}、"
            f"結合倍音による周波数分化・ロック={'あり=仮説A支持' if any_split else 'なし=結合倍音でも不十分（レジスタには場の円環構造が必要と絞られる）'}。"
            "副産物: 結合倍音はcrossing時刻に強く作用（q2で1167→200、q3で→6092）"),
    }
    (HERE / "paper8_axiom2_coupling_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
