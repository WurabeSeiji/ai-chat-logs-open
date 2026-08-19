#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL7-F3b — F3-2 の統計延長走行（T=48000・判定基準は F3 と同一の 3σ・実行前固定）

F3 の分離分岐（差 49.8°・2.86σ）はしきい未達だが差が巨大——判定は動かさず
走行を4倍に延長して決着させる。あわせて Z₂ ロック度と位相系列も全長で記録。
出力: result_dl7_f3b_v1.json・dl7_f3b_series_v1.npz
"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
UNI = EXP.parent / "統一万能関数_v1"
ODD = tuple(range(1, 18, 2))
T = 48000
NBLK = 24

def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s)
    sys.modules[n] = m; s.loader.exec_module(m); return m

def main():
    t0 = time.time()
    G = _load("g5_f3b", UNI / "unified_readout_v5.py")
    uni = _load("uni_f3b", UNI / "unified_interaction_v1.py")
    cr0 = _load("cr0_f3b", EXP / "run_cr0_control_no_theta_v2.py")
    base = uni.two_body_base; step = uni.collision_step_exact
    sp = base.build_source_params(base.Params(high_n=63, recursive_collision_count=200))
    nc, ne = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = cr0.calibrate_shift(sp, nc, ne)
    def fermion(which):
        c = base.explicit_packet_case(
            mode=f"f3b_{which}", packet_a=ODD, packet_b=ODD,
            packet_a_shift=cr0.shift_for_deg(-30.0, slope, icept),
            packet_b_shift=cr0.shift_for_deg(+30.0, slope, icept))
        s = base.make_case_state(sp, c, which, hair_enabled=True)
        return s / np.sqrt(np.vdot(s, s).real)
    def run(phi0):
        a = fermion("A"); b = np.exp(1j * phi0) * fermion("B")
        sep = np.empty(T); ph = np.empty(T)
        for t in range(T):
            a, b, _ = step(a, b, sp)
            ta, _ = cr0.circle_position(a, nc, ne)
            tb, _ = cr0.circle_position(b, nc, ne)
            sep[t] = abs(float(np.degrees(np.angle(np.exp(1j * (ta - tb))))))
            ph[t] = G.g_pair_overlap(a, b)["overlap_phase"]
        return sep, ph
    s0, p0 = run(0.0)
    sP, pP = run(np.pi)
    def blk(x):
        b_ = x.reshape(NBLK, -1).mean(axis=1)
        return float(b_.mean()), float(b_.std(ddof=1) / np.sqrt(NBLK))
    m0, e0 = blk(s0); mP, eP = blk(sP)
    d = m0 - mP; se = float(np.hypot(e0, eP))
    z2 = {"phi0": float(np.mean(np.cos(p0[T//2:]) ** 2)),
          "phiPi": float(np.mean(np.cos(pP[T//2:]) ** 2))}
    res = {"config": {"T": T, "n_blocks": NBLK, "packet_odd": list(ODD)},
           "z2_lock_late": z2,
           "sep_phi0": [m0, e0], "sep_phiPi": [mP, eP],
           "diff": d, "combined_se": se, "sigma": abs(d) / max(se, 1e-300),
           "split_3sigma": bool(abs(d) > 3 * se),
           "direction": ("同シート(φ=0)が大分離=斥力側・逆シートが引力側"
                         if d > 0 else "逆向き"),
           "elapsed_sec": time.time() - t0}
    np.savez_compressed(HERE / "dl7_f3b_series_v1.npz",
                        sep0=s0, sepP=sP, ph0=p0, phP=pP)
    (HERE / "result_dl7_f3b_v1.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"Z₂ロック(後半): φ=0: {z2['phi0']:.3f}  π: {z2['phiPi']:.3f}")
    print(f"分離: φ=0: {m0:.3f}±{e0:.3f}°  π: {mP:.3f}±{eP:.3f}°")
    print(f"差={d:+.3f}° ({res['sigma']:.2f}σ)  3σ分岐={res['split_3sigma']}  {res['direction']}")
    print(f"({res['elapsed_sec']:.0f}s)")

if __name__ == "__main__":
    main()
