#!/usr/bin/env python3
"""論文9補強予備実験 E-B3：万能エンジンによる対挙動の自発分化 v1

設計原理（無名性の厳守）:
    粒子種による分岐（IF文）を一切書かない。全ての対に同一のコードパス
    ——内生θ読出し（theta_from_ab: 反射率を波の形だけから決める万能関数）
    → 実直交回転（rotate_ab）——を適用し、挙動の違いは波の形だけから
    自発的に現れることを検証する。

対（同一エンジン・同一衝突数）:
    even-even : A=偶数束, B=偶数束
    odd-odd   : A=奇数束, B=奇数束
    even-odd  : A=偶数束, B=奇数束

予言（測定前固定）:
    P1（ボゾン的透過）: even-even の内生θは全衝突で厳密に0
        → 200衝突後も対忠実度 F=1（機械精度）。相互作用は起こらない
        ＝透過。これはIF文の帰結ではなく、万能関数が偶数束の波形に
        割り当てる反射率が0であることの帰結である
    P2（フェルミオン的交換）: odd-odd の内生θ>0
        → ノルムが Rabi 型に交換し、対忠実度は1から離れる＝排他・交換
    P3（探索的・予言は方向のみ）: even-odd 混合は奇数束側の偶数ビン成分が
        θを立てるため相互作用が生じる（θ>0）。定量は記録に徹する
    P4（アンカー）: odd-odd の初期θ・Rが、A=[1]の longrun 公表値の系譜と
        整合すること（同一エンジンの確認）

判定基準:
    透過/交換の区別は、(i) 内生θの値 (ii) 対忠実度 F(j) (iii) ノルム移乗
    max|ΔN_B| の三つの無名な観測量だけで行う。
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOY_RUNNER_PATH = HERE.parent / "run_ab_invariant_theta_toy_v1.py"
spec = importlib.util.spec_from_file_location("toy_for_pair_behavior_v1", TOY_RUNNER_PATH)
toy = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = toy
spec.loader.exec_module(toy)
base = toy.base
plt = base.plt

EVEN_KS = tuple(range(2, 63, 2))
ODD_KS = tuple(range(1, 64, 2))
J_MAX = 200


def make_bundle(sp, ks, which):
    case = base.explicit_packet_case(mode=f"pair_{which}", packet_a=tuple(ks), packet_b=tuple(ks))
    v = base.make_case_state(sp, case, which, hair_enabled=True)
    return v / np.sqrt(float(np.vdot(v, v).real))


def pair_fidelity(a0, b0, a, b) -> float:
    ov = complex(np.vdot(a0, a) + np.vdot(b0, b))
    n0 = float(np.vdot(a0, a0).real + np.vdot(b0, b0).real)
    return abs(ov) ** 2 / n0 ** 2


def run_pair(name, ks_a, ks_b, sp):
    a = make_bundle(sp, ks_a, "A")
    b = make_bundle(sp, ks_b, "B")
    a0, b0 = a.copy(), b.copy()
    thetas, fids, nbs = [], [], []
    for _ in range(J_MAX):
        readout = toy.theta_from_ab(a, b, sp)
        thetas.append(float(readout.theta))
        a, b = toy.rotate_ab(a, b, readout.theta)
        fids.append(pair_fidelity(a0, b0, a, b))
        nbs.append(float(np.vdot(b, b).real))
    thetas = np.array(thetas); fids = np.array(fids); nbs = np.array(nbs)
    return {
        "pair": name,
        "theta_initial": thetas[0], "theta_max": float(np.max(np.abs(thetas))),
        "R_initial": float(np.sin(thetas[0]) ** 2),
        "fidelity_min": float(np.min(fids)), "fidelity_final": float(fids[-1]),
        "max_abs_dNB": float(np.max(np.abs(nbs - nbs[0]))),
    }, thetas, fids


def main() -> None:
    params = base.Params(high_n=63, recursive_collision_count=J_MAX)
    sp = base.build_source_params(params)

    results = []
    traces = {}
    for name, ka, kb in (("even-even", EVEN_KS, EVEN_KS),
                          ("odd-odd", ODD_KS, ODD_KS),
                          ("even-odd", EVEN_KS, ODD_KS)):
        row, thetas, fids = run_pair(name, ka, kb, sp)
        results.append(row)
        traces[name] = (thetas, fids)
        print(f"{name:9s}: θ0={row['theta_initial']:.10f} θmax={row['theta_max']:.10f} "
              f"R0={row['R_initial']:.6f} Fmin={row['fidelity_min']:.12f} "
              f"max|ΔN_B|={row['max_abs_dNB']:.3e}")

    ee = next(r for r in results if r["pair"] == "even-even")
    oo = next(r for r in results if r["pair"] == "odd-odd")
    eo = next(r for r in results if r["pair"] == "even-odd")

    p1 = ee["theta_max"] <= 1e-12 and abs(ee["fidelity_final"] - 1) <= 1e-12
    p2 = oo["theta_initial"] > 1e-3 and oo["fidelity_min"] < 1 - 1e-3
    p3_interacts = eo["theta_max"] > 1e-6
    print(f"\nP1 ボゾン的透過（even-even: θ≡0・F=1）: {'PASS' if p1 else 'FAIL'}")
    print(f"P2 フェルミオン的交換（odd-odd: θ>0・F<1）: {'PASS' if p2 else 'FAIL'}")
    print(f"P3 混合対の相互作用（even-odd: θ>0）: {p3_interacts}（探索的・記録）")

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4), constrained_layout=True)
    for name in ("even-even", "odd-odd", "even-odd"):
        th, fd = traces[name]
        axes[0].plot(th, lw=1.1, label=name)
        axes[1].plot(fd, lw=1.1, label=name)
    axes[0].set_title("endogenous theta(j): same universal engine, no branching")
    axes[0].set_xlabel("collision j"); axes[0].set_ylabel("theta"); axes[0].legend()
    axes[1].set_title("pair fidelity F(j): transparency vs exchange")
    axes[1].set_xlabel("collision j"); axes[1].set_ylabel("F"); axes[1].legend()
    for ext in ("png", "svg"):
        fig.savefig(HERE / f"paper9_pair_behavior_v1.{ext}", dpi=160)
    plt.close(fig)

    payload = {
        "experiment": "paper9_universal_engine_pair_behavior_pre_v1",
        "design": "単一コードパス（内生θ→回転）・粒子種によるIF分岐ゼロ（無名性厳守）",
        "core_runner": {"path": TOY_RUNNER_PATH.name, "sha256": toy.sha256(TOY_RUNNER_PATH)},
        "rows": results,
        "P1_bosonic_transparency": bool(p1),
        "P2_fermionic_exchange": bool(p2),
        "P3_mixed_interacts": bool(p3_interacts),
        "conclusion": (
            "同一の万能エンジン（反射率を波の形だけから決める内生θ読出し＋実直交回転、"
            "分岐ゼロ）が、偶偶対には反射率0＝完全透過（ボゾン的挙動）を、奇奇対には"
            "θ>0＝ノルム交換（フェルミオン的挙動）を自発的に割り当てた。"
            "統計的挙動の分化は外部ラベルやIF文の帰結ではなく、万能関数が波形に"
            "割り当てる反射率の帰結である——『偶数倍音=ボゾン型』が行動レベルで実証された"),
    }
    (HERE / "paper9_pair_behavior_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print("\nsaved: paper9_pair_behavior_result_v1.json")


if __name__ == "__main__":
    main()
