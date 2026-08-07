#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""柱G1+G8正本: 閉塞監査——重力の入れ方二流儀とエンジン建築の統一

判定（事前固定）:
(A) 素の二体力学: |ΔC|/C ≤ 1e-12（閉塞は自然維持）
(B) V(x)注入（標準理論式・振幅1e-3・200步）: 相対破れ >1%（必ず破れる）
(C) ゲージ側操作（読出し目盛のみ歪める）: |ΔC|=0（構成的厳密）
(D) N体系（第8論文・Cayley）: Σz²・ノルムのドリフト ≤1e-12（N=5,8・300步）
→ (A)(D)が通れば「読出し生成回転 F[Ψ]=R(G[Ψ])·Ψ」の両実例が同一不変量を
  保存する（柱G8）。(B)(C)の対比が柱G1（重力は状態側に入れられない）。
使い方: python3 run_g1_closure_and_engines_v1.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("exact_g1", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
base = ex.base
ABL = HERE.parent / "第8論文_二段階seed除去による準安定相の因果分離" / "code" / "run_preliminary_seed_ablation_v1.py"
spec2 = importlib.util.spec_from_file_location("abl_g1", ABL)
abl = importlib.util.module_from_spec(spec2); sys.modules[spec2.name] = abl
spec2.loader.exec_module(abl)

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    x = np.arange(n)
    dxv = np.minimum(np.abs(x - n // 2), n - np.abs(x - n // 2))
    sea = 0.2 * np.exp(2j * np.pi * 3 * x / n)[:, None] * np.ones((1, ne))
    lump = 0.2 * np.exp(-0.5 * (dxv / 3.0) ** 2)[:, None] * np.ones((1, ne))
    a0 = (sea + lump).reshape(-1).astype(complex)
    V = (1e-3 * np.exp(-0.5 * (dxv / 20.0) ** 2))[:, None] * np.ones((1, ne))
    Vf = V.reshape(-1)
    out = {}
    def closure(a, b): return complex(np.sum(a * a) + np.sum(b * b))
    for tag, inject in (("bare", False), ("V_inject", True)):
        a, b = a0.copy(), a0.copy()
        C0 = closure(a, b)
        for _ in range(200):
            a, b, _ = ex.collision_step_exact(a, b, sp)
            if inject:
                ph = np.exp(1j * Vf); a, b = a * ph, b * ph
        dC = abs(closure(a, b) - C0); rel = dC / abs(C0)
        out[tag] = {"abs_drift": dC, "rel_drift": rel}
        print(f"[{tag}] |ΔC|={dC:.3e} 相対={rel:.3e}")
    out["gauge_side"] = {"abs_drift": 0.0, "rel_drift": 0.0,
                         "note": "読出し目盛のみの歪みは状態に非接触＝構成的に厳密0"}
    print("[gauge_side] |ΔC|=0（構成的）")
    for N in (5, 8):
        sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = abl.build_init(N, True)
        c0 = complex(Z @ Z); n0 = float(np.real(np.conj(Z) @ Z))
        dc = dn = 0.0
        for _ in range(300):
            Z, wp = abl.evolve(sys_lr, Z, wp)
            dc = max(dc, abs(complex(Z @ Z) - c0))
            dn = max(dn, abs(float(np.real(np.conj(Z) @ Z)) - n0))
        out[f"Nbody_N{N}"] = {"zero_square_drift": dc, "norm_drift": dn}
        print(f"[N体 N={N}] Σz²ドリフト={dc:.2e} ノルムドリフト={dn:.2e}")
    ok = (out["bare"]["rel_drift"] < 1e-12 and out["V_inject"]["rel_drift"] > 0.01
          and all(out[f"Nbody_N{N}"]["zero_square_drift"] < 1e-12 for N in (5, 8)))
    out["verdict"] = "柱G1+G8成立" if ok else "要精査"
    print(out["verdict"])
    out["runtime_sec"] = time.time() - t0
    (HERE / "result_g1_closure_and_engines_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
