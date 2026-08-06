#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v16: グラビトンℓ=2予言の検定——生成子スペクトルにℓ=2枠は存在するか

方法（事前記録）: 凝縮体の回転生成子 K̂ の固有値スペクトル ±iλ_j を全分解し、
占有平面のλ₁で校正した比 λ_j/λ₁ を列挙。ℓ=2枠（比≈2）の存否と、
窓サンプルの各固有平面への占有率を測る。
判定: H_slot: 比∈[1.8,2.2] の固有平面が存在 → グラビトン枠あり（占有率も報告、
未占有なら「枠あり・真空」）。存在しなければ予言は棄却。N=5,6,8。
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SPACE = HERE.parent / "空間軸3軸と固有時間の創生_v1"
spec1 = importlib.util.spec_from_file_location("pre1_v16", SPACE / "run_pre_2plus1_structure_v1.py")
pre1 = importlib.util.module_from_spec(spec1); sys.modules[spec1.name] = pre1
spec1.loader.exec_module(pre1)
abl = pre1.abl
edge_adjacency, build_K = pre1.edge_adjacency, pre1.build_K

T_END = 4000; WIN = (2000, 4000); SE = 5

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
        S = np.array(samples)
        Sp = S - np.outer(S @ p, p) - np.outer(S @ q, q)
        X = np.hstack([Sp.real, Sp.imag]); Xc = X - X.mean(axis=0)
        _, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
        theta = np.angle(S[-1])
        K = build_K(theta, adj)
        A = np.zeros((2 * M, 2 * M)); A[:M, :M] = K; A[M:, M:] = K
        A = 0.5 * (A - A.T)
        ev, EV = np.linalg.eig(A)
        lam = np.abs(ev.imag)
        idx = np.argsort(-lam)
        # 占有平面1の λ 校正: 平面1方向への射影が最大の固有対
        d1 = Vt[0]
        ovl = [abs(np.dot(d1, EV[:, i].real)) + abs(np.dot(d1, EV[:, i].imag)) for i in range(len(ev))]
        lam1 = lam[int(np.argmax(ovl))]
        # 各固有平面の占有率（窓パワー射影）
        Ptot = float(np.sum(Xc ** 2))
        rows = []
        seen = set()
        for i in idx:
            l_ = lam[i] / max(lam1, 1e-300)
            key = round(l_, 3)
            if key in seen or lam[i] < 1e-9: continue
            seen.add(key)
            vec_r, vec_i = EV[:, i].real, EV[:, i].imag
            nr = np.linalg.norm(vec_r); ni = np.linalg.norm(vec_i)
            occ = 0.0
            if nr > 1e-12: occ += float(np.sum((Xc @ (vec_r / nr)) ** 2))
            if ni > 1e-12: occ += float(np.sum((Xc @ (vec_i / ni)) ** 2))
            rows.append({"ratio": float(l_), "occupancy": occ / Ptot})
        rows.sort(key=lambda r: -r["ratio"])
        slot2 = [r for r in rows if 1.8 <= r["ratio"] <= 2.2]
        top = ", ".join(f"{r['ratio']:.3f}({r['occupancy']:.1%})" for r in rows[:8])
        print(f"N={n}: λ/λ₁スペクトル上位 = {top}")
        print(f"    ℓ=2枠: {'あり ' + str([(round(r['ratio'],3), round(r['occupancy'],4)) for r in slot2]) if slot2 else '**なし**'}")
        out["scan"].append({"N": n, "spectrum": rows, "slot2": slot2})
    out["H_slot"] = bool(any(r["slot2"] for r in out["scan"]))
    print(f"H_slot（ℓ=2枠の存在）= {out['H_slot']}")
    out["runtime_sec"] = time.time() - t0
    (HERE / "pre_graviton_slot_result_v16.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.0f}s")

if __name__ == "__main__":
    main()
