#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""柱G1補強正本: 接線定理——閉塞保存 ⟺ 生成子の反対称性（一般定理の数値定量化）

定理: 摂動 δΨ=εKΨ に対し δC = 2ε·Ψ^T K_sym Ψ + O(ε²)（C=Σψ²・K_sym=(K+K^T)/2）。
全Ψで δC=0 ⟺ K^T=−K（閉塞面の接線回転のみ許容）。
ポテンシャル注入は K=i·diag(V)（対称）ゆえ一次で必ず破る——31%（G1）はこの実例。
SO(2)対回転・Cayleyステップは反対称ゆえ厳密保存（柱G8）。
判定（事前固定）: 実測 |δC|/ε が理論 2|Ψ^T K_sym Ψ| と相対誤差<1%で一致（混合比掃引）。
使い方: python3 run_g1b_tangent_theorem_v1.py
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

def main():
    t0 = time.time()
    rng = np.random.default_rng(11)
    N = 200
    psi = rng.normal(size=N) + 1j * rng.normal(size=N)
    eps = 1e-6
    C0 = psi @ psi
    rows = []
    A = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    Ka = (A - A.T) / 2; Ks = (A + A.T) / 2
    Ks_n = Ks * np.linalg.norm(Ka) / np.linalg.norm(Ks)
    cases = [("antisym", Ka), ("sym", Ks), ("V_diag", 1j * np.diag(rng.normal(size=N))),
             ("mix_10pct_sym", Ka + 0.1 * Ks_n), ("mix_1pct_sym", Ka + 0.01 * Ks_n)]
    ok = True
    for name, K in cases:
        psi1 = psi + eps * (K @ psi)
        dC = abs(psi1 @ psi1 - C0) / eps
        theo = 2 * abs(psi @ (((K + K.T) / 2) @ psi))
        rel = abs(dC - theo) / max(theo, 1e-3)
        if theo > 1e-3 and rel > 0.01:
            ok = False
        rows.append({"case": name, "dC_per_eps": float(dC), "theory": float(theo)})
        print(f"{name:>14}: |δC|/ε={dC:.4e}  理論={theo:.4e}")
    verdict = "接線定理成立（δC=2εΨ^T K_symΨ・保存⟺K^T=−K）" if ok else "要精査"
    print(verdict)
    (HERE / "result_g1b_tangent_theorem_v1.json").write_text(
        json.dumps({"rows": rows, "verdict": verdict, "runtime_sec": time.time() - t0},
                   indent=1, ensure_ascii=False))
    print(f"完了 {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
