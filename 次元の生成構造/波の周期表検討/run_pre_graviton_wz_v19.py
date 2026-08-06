#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v19: A=グラビトン複合励起（四重極2ω線）／B=離調-質量関係（W/Z検定）

Part A（事前記録）: グラビトン=ℓ=1量子2個の複合（v16でℓ=2線形枠なし→double copy
再解釈）。計量的観測量=読出しの二乗 c²(t) は枠回転で e^{2iβ}（ℓ=2厳密・代数）。
力学的実在の判定 H_gcomp: c² スペクトルに 2ω_clock の線が立ち（|f₂/f₁−2|<0.02）、
線幅比 w₂/w₁ ≤ 3（コヒーレントな四重極波＝複合グラビトンが持続的に存在）。

Part B（事前記録）: 生成子の副固有平面（λ/λ₁<1）は天然の離調種。各固有平面の
実測回転レート比 ρ_i と分散補償Gram質量²_i を測る。
判定 H_WZ: 質量² が離調 δ=|1−ρ| とともに単調増加（相関r>0.7）
——「t巻き（自前時計の離調）⇒質量」＝W/Z有質量の機構検定。

使い方: python3 run_pre_graviton_wz_v19.py
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SPACE = HERE.parent / "空間軸3軸と固有時間の創生_v1"
spec1 = importlib.util.spec_from_file_location("pre1_v19", SPACE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1); sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl
edge_adjacency, build_K = pre1.edge_adjacency, pre1.build_K

T_END = 4000; WIN = (2000, 4000); SE = 5

def linewidth(F, pk):
    h = F[pk] / 2
    lo = pk
    while lo > 0 and F[lo] > h: lo -= 1
    hi = pk
    while hi < len(F) - 1 and F[hi] > h: hi += 1
    return hi - lo

def main():
    t0 = time.time()
    out = {"scan": []}
    for n in (5, 6, 8):
        sys_lr, v, _, _, _, p, q, Z, wp = abl.build_init(n, True)
        M = sys_lr.m
        adj = edge_adjacency(n)
        samples = []
        for t in range(T_END):
            Z, wp = abl.evolve(sys_lr, Z, wp)
            if WIN[0] <= t < WIN[1] and (t % SE == 0):
                samples.append(Z.copy())
        S = np.array(samples); ns = S.shape[0]
        Sp = S - np.outer(S @ p, p) - np.outer(S @ q, q)
        X = np.hstack([Sp.real, Sp.imag]); Xc = X - X.mean(axis=0)
        _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
        d1c = (Vt[0][:M] + 1j * Vt[0][M:]); d1c /= np.linalg.norm(d1c)
        d2c = (Vt[1][:M] + 1j * Vt[1][M:]); d2c /= np.linalg.norm(d2c)
        c1 = Sp @ np.conj(d1c) + 1j * (Sp @ np.conj(d2c))
        # Part A: 基本線 vs 四重極線
        pad = 4
        F1 = np.abs(np.fft.fft(c1 - c1.mean(), n=pad * ns))
        F2 = np.abs(np.fft.fft(c1 ** 2 - np.mean(c1 ** 2), n=pad * ns))
        p1 = int(np.argmax(F1)); p2 = int(np.argmax(F2))
        L = pad * ns
        f1 = p1 / L if p1 <= L // 2 else (p1 - L) / L
        f2 = p2 / L if p2 <= L // 2 else (p2 - L) / L
        ratio = f2 / f1 if f1 != 0 else float("nan")
        w1 = linewidth(F1, p1); w2 = linewidth(F2, p2)
        ok_A = bool(abs(ratio - 2) < 0.02 and w2 <= 3 * w1)
        print(f"N={n} PartA: f₂/f₁={ratio:.4f}  線幅比 w₂/w₁={w2}/{w1}={w2/max(w1,1):.2f} → "
              f"{'四重極線あり' if ok_A else '不成立'}")
        # Part B: 生成子固有平面ごとの離調と質量
        theta = np.angle(S[-1])
        K = build_K(theta, adj)
        A = np.zeros((2 * M, 2 * M)); A[:M, :M] = K; A[M:, M:] = K
        A = 0.5 * (A - A.T)
        ev, EV = np.linalg.eig(A)
        # 時計レート（親平面）
        cp = S @ p; cq = S @ q
        phi = np.unwrap(np.angle(cp + 1j * cq))
        om_clock = float(np.polyfit(np.arange(ns), phi, 1)[0])
        rows = []
        seen = set()
        idx = np.argsort(-np.abs(ev.imag))
        for i in idx[:16]:
            if ev[i].imag <= 1e-9: continue
            lam = float(ev[i].imag)
            key = round(lam, 6)
            if key in seen: continue
            seen.add(key)
            vr = EV[:, i].real; vi = EV[:, i].imag
            if np.linalg.norm(vr) < 1e-12 or np.linalg.norm(vi) < 1e-12: continue
            vr /= np.linalg.norm(vr); vi /= np.linalg.norm(vi)
            r1 = Xc @ vr; r2 = Xc @ vi
            occ = float(np.mean(r1 ** 2 + r2 ** 2))
            c = r1 + 1j * r2
            # 実測回転レート
            rate = float(np.angle(np.sum(c[1:] * np.conj(c[:-1]))))
            rho = rate / om_clock
            # 分散補償Gram質量²（2チャネル=Re/Im投影、各自前時計で復調）
            t_ = np.arange(ns)
            z1 = r1.astype(complex); z2 = r2.astype(complex)
            # 実系列はヒルベルト的に複素化: c と conj(c) の対で代用
            ca = c; cb = np.conj(c)
            wA_ = np.angle(np.sum(ca[1:] * np.conj(ca[:-1])))
            wB_ = np.angle(np.sum(cb[1:] * np.conj(cb[:-1])))
            Ad = ca * np.exp(-1j * wA_ * t_); Bd = cb * np.exp(-1j * wB_ * t_)
            Gaa = float(np.mean(np.abs(Ad) ** 2)); Gbb = float(np.mean(np.abs(Bd) ** 2))
            Gab = complex(np.mean(Ad * np.conj(Bd)))
            T_ = 0.5 * (Gaa + Gbb)
            m2 = float((Gaa * Gbb - abs(Gab) ** 2) / T_ ** 2) if T_ > 0 else 0.0
            if occ < 1e-8: continue
            rows.append({"rho": rho, "delta": abs(1 - abs(rho)), "mass2": m2, "occ": occ})
        rows.sort(key=lambda r: r["delta"])
        ds = [r["delta"] for r in rows]; m2s = [r["mass2"] for r in rows]
        corr = float(np.corrcoef(ds, m2s)[0, 1]) if len(rows) >= 4 else None
        txt = ", ".join(f"δ={r['delta']:.2f}:m²={r['mass2']:.3f}" for r in rows[:6])
        print(f"N={n} PartB: {txt}")
        print(f"    離調-質量相関 r={corr}")
        out["scan"].append({"N": n, "A": {"ratio": ratio, "w1": w1, "w2": w2, "pass": ok_A},
                             "B": {"rows": rows, "corr": corr}})
    okA = all(r["A"]["pass"] for r in out["scan"])
    corrs = [r["B"]["corr"] for r in out["scan"] if r["B"]["corr"] is not None]
    okB = bool(corrs and all(c > 0.7 for c in corrs))
    print(f"\nH_gcomp（複合四重極波）= {okA}   H_WZ（離調⇒質量, r>0.7全N）= {okB}")
    out["H_gcomp"] = okA; out["H_WZ"] = okB
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_graviton_wz_result_v19.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
