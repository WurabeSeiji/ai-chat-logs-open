#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7-F3c — 初期配置アンサンブルによる±力の決着（実行前固定）

F3/F3b の教訓: 決定論系で長時間平均のσ判定は不適（分離が対蹠へ飽和し差が洗われる）。
決着は初期配置アンサンブル: 初期分離5通り × φ∈{0,π}（フェルミオン型対・Z₂シート）で
初期窓（T=2000）の平均分離を測り、**全配置で φ=0 > φ=π（同シート=分離大）が揃うか**
を符号検定で判定する（5/5 なら二項 p=1/32、事前登録）。

出力: result_dl7_f3c_v1.json・dl7_f3c_series_v1.npz
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"
ODD = tuple(range(1, 18, 2))
T = 2000
GEOMS = [(-30.0, 20.0), (-30.0, 30.0), (-30.0, 45.0), (-20.0, 40.0), (-45.0, 30.0)]

def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s)
    sys.modules[n] = m; s.loader.exec_module(m); return m

def main():
    t0 = time.time()
    uni = _load("uni_f3c", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_f3c", EXP / "run_cr0_control_no_theta_v2.py")
    base = uni.two_body_base; step = uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)

    def pair(dega, degb, phi0):
        c = base.explicit_packet_case(
            mode=f"f3c_{dega}_{degb}", packet_a=ODD, packet_b=ODD,
            packet_a_shift=cr0.shift_for_deg(dega, slope, icept),
            packet_b_shift=cr0.shift_for_deg(degb, slope, icept))
        a = base.make_case_state(sp, c, "A", hair_enabled=True)
        b = base.make_case_state(sp, c, "B", hair_enabled=True)
        a /= np.sqrt(np.vdot(a, a).real); b /= np.sqrt(np.vdot(b, b).real)
        return a, np.exp(1j * phi0) * b

    def run(dega, degb, phi0):
        a, b = pair(dega, degb, phi0)
        sep = np.empty(T)
        for t in range(T):
            a, b, _ = step(a, b, sp)
            ta, _ = cr0.circle_position(a, nc, ne)
            tb, _ = cr0.circle_position(b, nc, ne)
            sep[t] = abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb))))))
        return sep

    rows = []
    seps = {}
    for g in GEOMS:
        s0 = run(g[0], g[1], 0.0)
        sP = run(g[0], g[1], np.pi)
        d = float(s0.mean() - sP.mean())
        rows.append({"geom": g, "mean_phi0": float(s0.mean()),
                     "mean_phiPi": float(sP.mean()), "diff": d,
                     "sign": int(np.sign(d))})
        seps[f"{g}"] = (s0, sP)
        print(f"  geom={g}: φ=0 {s0.mean():7.2f}°  π {sP.mean():7.2f}°  差 {d:+8.3f}°")
    signs = [r["sign"] for r in rows]
    n_pos = sum(1 for s in signs if s > 0)
    unanimous = bool(abs(sum(signs)) == len(signs))
    res = {"config": {"T": T, "geoms": GEOMS, "packet_odd": list(ODD)},
           "rows": rows, "n_pos_of_5": n_pos, "unanimous": unanimous,
           "binom_p_two_sided": float(2 * 0.5 ** len(signs)) if unanimous else None,
           "verdict": ("同シート(φ=0)=分離大が全配置で成立——±力の方向が確定"
                       if unanimous and n_pos == 5 else
                       "逆シート=分離大が全配置で成立" if unanimous else
                       "配置により符号が割れる——未決着"),
           "elapsed_sec": time.time() - t0}
    np.savez_compressed(HERE / "dl7_f3c_series_v1.npz",
                        **{f"s0_{i}": seps[f"{g}"][0] for i, g in enumerate(GEOMS)},
                        **{f"sP_{i}": seps[f"{g}"][1] for i, g in enumerate(GEOMS)})
    (HERE / "result_dl7_f3c_v1.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"符号検定: +{n_pos}/5  全会一致={unanimous}  → {res['verdict']}")
    print(f"({res['elapsed_sec']:.0f}s)")

if __name__ == "__main__":
    main()
