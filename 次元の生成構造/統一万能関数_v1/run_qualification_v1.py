#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""統一万能関数の資格審査 v1 — 正本との対照テスト（全通過が使用条件）

Q1 相互作用の忠実性: UnifiedEngine が TB v2 の検証済みエンジン
   （正本 stage3 HairEngine のコピー・genesis 物理を再現済み）と
   ビット単位で同一軌道（N=5・T=200・標準宇宙構成）。
Q2 正本退化: UnifiedEngine(Nη=1) が正本 VertexEngineV2（stage3）と
   ビット単位で同一軌道（N=5・T=200・genesis v1 構成のレジスタ5帯）。
Q3 読出しの正本一致: unified_readout の f₂・r が run_genesis_v1 の
   インライン式と数値一致（<1e-14 相対）——同一軌道上で全步比較。
Q4 受動性: パネル併走の有無で終状態ビット単位同一。
Q5 二体正本の再輸出: collision_step_exact が読み込め、50衝突で
   閉塞ドリフト <1e-12（正本性能）。
使い方: python3 run_qualification_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def main():
    t0 = time.time()
    ui = load("ui_q", HERE / "unified_interaction_v1.py")
    ur = load("ur_q", HERE / "unified_readout_v1.py")
    tb = load("tb_q", HERE.parent / "電子の反跳実験" / "run_tb_genesis_universe_v2.py")

    results = {}

    # Q1: UnifiedEngine ≡ TB v2 エンジン（ビット同一）
    e1, p2, q2 = ui.build_standard_universe(5, 1e-2)
    e2, p2b, q2b = tb.build_universe(5, 1e-2)
    ok = np.array_equal(e1.C, e2.C)
    for t in range(200):
        e1.step(); e2.step()
        if not np.array_equal(e1.C, e2.C):
            ok = False
            break
    results["Q1"] = bool(ok)
    print(f"Q1 相互作用の忠実性（TB v2検証済エンジンとビット同一・T=200）: "
          f"{'通過' if ok else '不成立'}")

    # Q2: Nη=1 退化 ≡ 正本 VertexEngineV2
    n, m = 5, 10
    _, v, _, _, _, _, _, Z0c, wp0 = ui.abl.build_init(n, False)
    r2 = ui.gen3.make_parent(n, seed=2)
    Csec = np.fft.fft(r2.relation_waves, axis=1) / n
    seed_state = Csec[:, 1] / np.linalg.norm(Csec[:, 1])
    C0 = np.zeros((m, 5), complex)
    C0[:, 2] = Z0c
    C0[:, 1] = 1e-2 * seed_state
    eV2 = ui.VertexEngineV2(n, C0.copy(), wp0, vertex_on=True)
    eU = ui.UnifiedEngine(n, C0.copy().reshape(m, 5, 1), wp0, vertex_on=True)
    ok = True
    for t in range(200):
        eV2.step(); eU.step()
        if not np.array_equal(eV2.C, eU.C):
            ok = False
            break
    results["Q2"] = bool(ok)
    print(f"Q2 正本退化（Nη=1 ≡ VertexEngineV2・ビット同一・T=200）: "
          f"{'通過' if ok else '不成立'}")

    # Q3: 読出しの正本一致（f₂・r を genesis v1 のインライン式と全步比較）
    eng, p2, q2 = ui.build_standard_universe(5, 1e-2)
    prev_flat = eng.C.reshape(-1).copy()
    carry_flat, carry_gen = None, None
    worst_f2 = 0.0; worst_r = 0.0
    for t in range(300):
        eng.step()
        C2 = eng.C2()
        pan = ur.g_panel(C2, p2, q2, carry_flat, carry_gen)
        carry_flat = pan["_carry"]["C_flat"]
        carry_gen = pan["_carry"]["c_gen"]
        # 正本インライン式（run_genesis_v1 の定義そのまま）
        Z2 = C2[:, 2, 0]
        Zp = Z2 - p2 * (p2 @ Z2) - q2 * (q2 @ Z2)
        f2_ref = float(np.real(np.conj(Zp) @ Zp)) / float(np.real(np.conj(Z2) @ Z2))
        cf = eng.C.reshape(-1)
        ip = np.vdot(prev_flat, cf)
        ph = ip / abs(ip) if abs(ip) > 0 else 1.0
        r_ref = float(np.linalg.norm(cf - ph * prev_flat))
        prev_flat = cf.copy()
        if f2_ref > 0:
            worst_f2 = max(worst_f2, abs(pan["f2"] - f2_ref) / f2_ref)
        if "r" in pan and r_ref > 0:
            worst_r = max(worst_r, abs(pan["r"] - r_ref) / r_ref)
    ok = worst_f2 < 1e-14 and worst_r < 1e-14
    results["Q3"] = bool(ok)
    print(f"Q3 読出し正本一致（f₂差{worst_f2:.1e}・r差{worst_r:.1e}<1e-14）: "
          f"{'通過' if ok else '不成立'}")

    # Q4: 受動性
    ea, p2, q2 = ui.build_standard_universe(5, 1e-2)
    eb, _, _ = ui.build_standard_universe(5, 1e-2)
    carry_flat, carry_gen = None, None
    for t in range(150):
        ea.step()
        pan = ur.g_panel(ea.C2(), p2, q2, carry_flat, carry_gen)
        carry_flat = pan["_carry"]["C_flat"]; carry_gen = pan["_carry"]["c_gen"]
        eb.step()
    ok = np.array_equal(ea.C, eb.C)
    results["Q4"] = bool(ok)
    print(f"Q4 受動性（パネル有無でビット同一・T=150）: {'通過' if ok else '不成立'}")

    # Q5: 二体正本の再輸出
    params = ui.two_body_base.Params(high_n=63, recursive_collision_count=200)
    sp = ui.two_body_base.build_source_params(params)
    n2, ne2 = sp.chi_grid_n, sp.eta_grid_n
    rng = np.random.default_rng(7)
    a = (rng.normal(size=n2 * ne2) + 1j * rng.normal(size=n2 * ne2)) * 0.1
    b = -1j * a
    pw = float(np.sum(np.abs(a) ** 2 + np.abs(b) ** 2))
    for j in range(50):
        a, b, _ = ui.collision_step_exact(a, b, sp)
    cl1 = complex(np.sum(a * a) + np.sum(b * b))
    drift = abs(cl1) / pw   # b=−ia は閉塞恒等0 → ノルム規格化の絶対評価
    ok = drift < 1e-12
    results["Q5"] = bool(ok)
    print(f"Q5 二体正本再輸出（50衝突の閉塞ドリフト {drift:.1e}<1e-12）: "
          f"{'通過' if ok else '不成立'}")

    allpass = all(results.values())
    print("資格審査: " + ("ALL PASS — 統一関数は使用可" if allpass else "不成立あり"))
    results.update({"all_pass": bool(allpass), "worst_f2": worst_f2,
                    "worst_r": worst_r, "two_body_closure_drift": drift,
                    "runtime_sec": time.time() - t0})
    (HERE / "qualification_result_v1.json").write_text(
        json.dumps(results, indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
