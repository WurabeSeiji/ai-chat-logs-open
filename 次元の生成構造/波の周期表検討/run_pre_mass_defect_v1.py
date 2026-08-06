#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""予備実験 P1: 質量欠損——凝縮体は複合粒子か

仮説（事前記録）: 凝縮体=3構成子の束縛状態なら、構成子質量の和 Σm_i と
複合質量 M の間に系統的な差（欠損または超過）が出る。

規約（事前固定）:
  - 構成子 = 占有SVDの上位3平面（(v1,v2),(v3,v4),(v5,v6)）。
  - 各平面の読出し対を集団時計で復調し 2×2 Gram Γ_i を構成。
    質量_i = √max(detΓ_i,0)（絶対）、無次元度 = 4detΓ_i/(trΓ_i)²。
  - 複合の定義2種を併記（どちらも事前固定・後から選ばない）:
      定義A（混合）: Γ_tot = ΣΓ_i → M_A = √detΓ_tot
      定義B（最上位対のみ・空間時間論文の凝縮体質量と同一）: M_B = √detΓ_1
  - 欠損指標: Δ_A = Σm_i − M_A。線形性の基準: Γが共通固有基底なら
    √det は劣加法（√det(Γ1+Γ2) ≥ √detΓ1+√detΓ2 が2×2半正定値で成立）
    ——つまり混合は必ず M_A ≥ Σm_i（数学的超過）。物理的欠損は
    これを上回る/下回る系統偏差として読む。
  - 窓2分割（[2000,3000],[3000,4000]）で安定性を確認。N=5,6,8。

読み（記述的）:
  R1 各構成子の質量が同程度か階層的か（質量階層の第一痕跡）
  R2 M_A/Σm_i の比が窓・Nで安定か（束縛の秩序変数候補）
  R3 平面間Gramの非対角（構成子間コヒーレンス）——束縛の直接証拠候補

使い方: python3 run_pre_mass_defect_v1.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SPACE = HERE.parent / "空間軸3軸と固有時間の創生_v1"
spec1 = importlib.util.spec_from_file_location(
    "pre1_md", SPACE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1)
sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl

T_END = 4000
WINS = ((2000, 3000), (3000, 4000))
SAMPLE_EVERY = 5
N_LIST = (5, 6, 8)


def measure(n: int, win) -> dict:
    sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
    M = sys_lr.m
    samples = []
    for t in range(T_END):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        if win[0] <= t < win[1] and (t % SAMPLE_EVERY == 0):
            samples.append(Z.copy())
    S = np.array(samples)
    ns = S.shape[0]
    Sp = S - np.outer(S @ p, p) - np.outer(S @ q, q)

    X = np.hstack([Sp.real, Sp.imag])
    Xc = X - X.mean(axis=0)
    _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)

    cp = S @ p
    cq = S @ q
    phi = np.unwrap(np.angle(cp + 1j * cq))
    omega_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])

    def demod_dir(d):
        dc = d[:M] + 1j * d[M:]
        dc = dc / np.linalg.norm(dc)
        c = Sp @ np.conj(dc)
        F = np.fft.fft(c - c.mean())
        pk = int(np.argmax(np.abs(F)))
        f = pk / ns if pk <= ns // 2 else (pk - ns) / ns
        n_k = round(2 * np.pi * f / omega_clock) if omega_clock != 0 else 0
        return c * np.exp(-1j * n_k * phi)

    planes = []
    rs = []
    for j in range(3):
        d1, d2 = Vt[2 * j], Vt[2 * j + 1]
        r1, r2 = demod_dir(d1), demod_dir(d2)
        rs.append((r1, r2))
        G11 = float(np.mean(np.abs(r1) ** 2))
        G22 = float(np.mean(np.abs(r2) ** 2))
        G12 = complex(np.mean(r1 * np.conj(r2)))
        detG = G11 * G22 - abs(G12) ** 2
        T = 0.5 * (G11 + G22)
        planes.append({"plane": j + 1,
                        "power": G11 + G22,
                        "mass_abs": float(np.sqrt(max(detG, 0.0))),
                        "mass_deg": float(detG / T ** 2) if T > 0 else 0.0,
                        "sv_pair": [float(sv[2 * j]), float(sv[2 * j + 1])]})

    # 定義A: 混合Gram
    Gt = np.zeros((2, 2))
    for (r1, r2) in rs:
        Gt[0, 0] += float(np.mean(np.abs(r1) ** 2))
        Gt[1, 1] += float(np.mean(np.abs(r2) ** 2))
        g12 = complex(np.mean(r1 * np.conj(r2)))
        Gt[0, 1] += abs(g12)
        Gt[1, 0] += abs(g12)
    detA = Gt[0, 0] * Gt[1, 1] - Gt[0, 1] * Gt[1, 0]
    M_A = float(np.sqrt(max(detA, 0.0)))

    # 構成子間コヒーレンス（平面1の第1チャネルと平面2/3の第1チャネル）
    def xcoh(a, b):
        return float(abs(np.mean(a * np.conj(b)))
                     / max(np.sqrt(np.mean(np.abs(a) ** 2) * np.mean(np.abs(b) ** 2)), 1e-300))
    cross = {"p1p2": xcoh(rs[0][0], rs[1][0]),
             "p1p3": xcoh(rs[0][0], rs[2][0]),
             "p2p3": xcoh(rs[1][0], rs[2][0])}

    sum_mi = float(sum(pl["mass_abs"] for pl in planes))
    return {"N": n, "win": list(win), "omega_clock": omega_clock,
            "planes": planes, "sum_mass": sum_mi,
            "M_A_mix": M_A, "M_B_top": planes[0]["mass_abs"],
            "ratio_A": M_A / sum_mi if sum_mi > 0 else None,
            "cross_coherence": cross}


def main() -> None:
    t0 = time.time()
    out = {"T_END": T_END, "WINS": [list(w) for w in WINS],
           "SAMPLE_EVERY": SAMPLE_EVERY, "N_LIST": list(N_LIST), "runs": []}
    for n in N_LIST:
        for win in WINS:
            r = measure(n, win)
            out["runs"].append(r)
            pls = r["planes"]
            print(f"N={n} 窓{win}: 構成子質量=[" +
                  ", ".join(f"{pl['mass_abs']:.3e}" for pl in pls) + "]")
            print(f"       Σm_i={r['sum_mass']:.3e}  M_A(混合)={r['M_A_mix']:.3e} "
                  f"比={r['ratio_A']:.3f}  M_B(最上位)={r['M_B_top']:.3e}")
            print(f"       構成子間コヒーレンス={ {k: round(v,3) for k,v in r['cross_coherence'].items()} }")
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_mass_defect_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s → pre_mass_defect_result_v1.json")


if __name__ == "__main__":
    main()
