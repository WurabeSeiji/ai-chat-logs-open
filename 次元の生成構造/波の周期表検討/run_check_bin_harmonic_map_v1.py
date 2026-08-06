#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""検算スクリプト: bin→χ周波数写像（±1シフト）と束の毛（η巻きm=+1）の実測
（スピン統計定理のラベル整合の根拠。当初ヒアドキュメント→正式化）"""
from __future__ import annotations
import importlib.util, json, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
UIM = HERE.parent / "万能非弾性写像_managed_v1"
spec = importlib.util.spec_from_file_location("chkmap", UIM / "run_ignition_fate_exact_v3.py")
ex = importlib.util.module_from_spec(spec); sys.modules[spec.name] = ex
spec.loader.exec_module(ex)
v1, toy, base = ex.v1, ex.toy, ex.base

def main():
    t0 = time.time()
    params = base.Params(high_n=63, recursive_collision_count=200)
    sp = base.build_source_params(params)
    n, ne = sp.chi_grid_n, sp.eta_grid_n
    ms = np.arange(ne); mm = np.where(ms <= ne // 2, ms, ms - ne)
    out = {"chi_n": int(n), "eta_n": int(ne), "bundles": []}
    for label, bins in (("偶bin束(30,32,34)", (30, 32, 34)),
                          ("奇bin束(29,31,33)", (29, 31, 33)), ("種bin(21)", (21,))):
        a = v1.make_bundle(sp, bins, "A", scale=1.0)
        bz = np.zeros_like(a)
        freqs, power = toy.combined_chi_power(a, bz, sp)
        tot = power.sum()
        even_pow = float(power[(np.abs(freqs) % 2 == 0)].sum() / tot)
        idx = np.argsort(-power)[:6]
        tops = [[int(freqs[i]), float(power[i])] for i in idx if power[i] > 1e-6]
        f = np.fft.fft(a.reshape(n, ne), axis=1, norm="ortho")
        Pm = np.sum(np.abs(f) ** 2, axis=0)
        eta_occ = [[int(mm[i]), float(Pm[i] / Pm.sum())] for i in range(ne)
                    if Pm[i] / Pm.sum() > 1e-6]
        print(f"{label}: χ周波数トップ={tops[:4]} 偶パワー比={even_pow:.4f} η占有={eta_occ}")
        out["bundles"].append({"label": label, "bins": list(bins),
                                 "chi_freq_top": tops, "even_freq_frac": even_pow,
                                 "eta_occupancy": eta_occ})
    # 判定: bin→周波数±1シフト（偶bin→奇周波数・逆も）と毛m=+1
    b0, b1 = out["bundles"][0], out["bundles"][1]
    h_shift = bool(b0["even_freq_frac"] < 1e-6 and b1["even_freq_frac"] > 1 - 1e-6)
    h_hair = all(len(b["eta_occupancy"]) == 1 and b["eta_occupancy"][0][0] == 1
                 for b in out["bundles"])
    print(f"H_shift（bin↔周波数パリティ反転）= {h_shift}  H_hair（毛m=+1厳密）= {h_hair}")
    out["H_shift"] = h_shift; out["H_hair"] = h_hair
    out["runtime_sec"] = time.time() - t0
    (HERE / "check_bin_harmonic_map_result_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"完了 {out['runtime_sec']:.1f}s")

if __name__ == "__main__":
    main()
