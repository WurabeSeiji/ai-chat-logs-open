#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR4: 倍音構成を変えた相対位相の系列（1步刻み・CR3 と同一形式）

CR3 は元の実験条件（A:1〜17 / B:1〜3）を 1步刻みで見た。本スクリプトは
同じ計器立てのまま倍音構成だけを変えて系列にする。

  case17_3 : A=1〜17 / B=1〜3   元の実験条件（非対称）
  case3_3  : A=1〜3  / B=1〜3   対称・少倍音
  case17_17: A=1〜17 / B=1〜17  対称・多倍音

いずれも 1步刻み（間引きなし）。Δθ の周期は 5 步台なので、20 步おきの
標本化ではナイキスト（40 步）を割ってエイリアスする（CR3 で確認済み）。

出力: cr4_relative_<tag>_data_v1.json（各ケース）
図の生成は build_cr4_pages_v1.py が cr_relative_shell_v1.html と組み合わせる。

使い方:
  python3 make_cr4_relative_series_v1.py            # 全ケース
  python3 make_cr4_relative_series_v1.py case3_3    # 個別
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
UNI = HERE.parent / "統一万能関数_v1"


def _load(n, p):
    s = importlib.util.spec_from_file_location(n, p)
    m = importlib.util.module_from_spec(s)
    sys.modules[n] = m
    s.loader.exec_module(m)
    return m


_uni = _load("uni_cr4", UNI / "unified_interaction_v1.py")
K = _load("kin_cr4", UNI / "unified_kinetic_v1.py")
_cr0 = _load("cr0_cr4", HERE / "run_cr0_control_no_theta_v2.py")
_cr1 = _load("cr1_cr4", HERE / "run_cr1_kinetic_feedback_v1.py")

CASES = {
    "case17_3":  {"a": tuple(range(1, 18)), "b": tuple(range(1, 4)),
                  "label": "A:倍音1〜17 / B:倍音1〜3", "note": "元の実験条件（非対称）"},
    "case3_3":   {"a": tuple(range(1, 4)),  "b": tuple(range(1, 4)),
                  "label": "A・B ともに倍音1〜3", "note": "対称・少倍音"},
    "case17_17": {"a": tuple(range(1, 18)), "b": tuple(range(1, 18)),
                  "label": "A・B ともに倍音1〜17", "note": "対称・多倍音"},
}
DEG_A, DEG_B = -30.0, +30.0
T_STEPS = 4000
OMEGA0 = np.pi / 72.0


def run_case(tag):
    cfg = CASES[tag]
    sp = _uni.two_body_base.build_source_params(
        _uni.two_body_base.Params(high_n=63, recursive_collision_count=200))
    n_chi, n_eta = int(sp.chi_grid_n), int(sp.eta_grid_n)
    slope, icept, _ = _cr0.calibrate_shift(sp, n_chi, n_eta)
    case = _uni.two_body_base.explicit_packet_case(
        mode="cr4_" + tag, packet_a=cfg["a"], packet_b=cfg["b"],
        packet_a_shift=_cr0.shift_for_deg(DEG_A, slope, icept),
        packet_b_shift=_cr0.shift_for_deg(DEG_B, slope, icept))
    a = _uni.two_body_base.make_case_state(sp, case, "A", hair_enabled=True)
    b = _uni.two_body_base.make_case_state(sp, case, "B", hair_enabled=True)

    _, z0a = _cr0.circle_position(a, n_chi, n_eta)
    _, z0b = _cr0.circle_position(b, n_chi, n_eta)

    omega, v = OMEGA0, 0.0
    H = {k: [] for k in ("chi", "mid", "pr", "Rp", "r", "v", "omega",
                         "cloRe", "cloIm")}

    for _ in range(T_STEPS):
        pa, _ = _cr0.circle_position(a, n_chi, n_eta)
        pb, _ = _cr0.circle_position(b, n_chi, n_eta)
        chi = float(np.angle(np.exp(1j * (pa - pb))))
        mid = float(np.angle(np.exp(1j * pa) + np.exp(1j * pb)))
        r_now = float(_cr1.toy.theta_from_ab(a, b, sp).reflection_rate)
        acc = -4.0 * np.sin(omega / 2.0) ** 2 * chi
        v += acc
        omega += (1.0 - r_now) * acc          # κ = 1−r

        RpA, _, _, _ = _cr1.cone_components(a, n_chi, n_eta)
        Z = complex(np.sum(a * a) + np.sum(b * b))
        H["chi"].append(round(chi, 8)); H["mid"].append(round(mid, 6))
        H["pr"].append(round(_cr0.participation_ratio(a, n_chi, n_eta), 2))
        H["Rp"].append(float(f"{RpA:.6g}")); H["r"].append(round(r_now, 6))
        H["v"].append(float(f"{v:.6g}")); H["omega"].append(round(omega, 8))
        H["cloRe"].append(float(f"{Z.real:.4g}"))
        H["cloIm"].append(float(f"{Z.imag:.4g}"))

        a = K.k_translate_flat(a, -v, n_chi, n_eta)
        a, b, _ = _uni.collision_step_exact(a, b, sp)

    chi = np.array(H["chi"])
    dchi = np.gradient(chi)

    # スペクトル診断（1步刻みなのでナイキストは 2 步）
    y = chi - chi.mean()
    Y = np.abs(np.fft.rfft(y * np.hanning(len(y)))); f = np.fft.rfftfreq(len(y)); Y[0] = 0
    tot = float(np.sum(Y ** 2))
    idx = np.argsort(Y)[::-1]
    pk = [i for i in idx if 0 < i < len(Y) - 1 and Y[i] > Y[i - 1] and Y[i] > Y[i + 1]][:6]
    spec = [{"period": float(1 / f[i]), "share": float(100 * Y[i] ** 2 / tot)} for i in pk]
    lo = (f > 0) & (1 / np.maximum(f, 1e-12) >= 100)
    share100 = float(100 * np.sum(Y[lo] ** 2) / tot)
    band = (f > 0) & (1 / np.maximum(f, 1e-12) > 129.6) & (1 / np.maximum(f, 1e-12) < 158.4)
    share144 = float(100 * np.max(np.where(band, Y, 0)) ** 2 / tot) if band.sum() else 0.0

    out = {
        "experiment": "cr4_relative_series_v1", "case": tag,
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "config": {"packet_a": list(cfg["a"]), "packet_b": list(cfg["b"]),
                   "label": cfg["label"], "note": cfg["note"],
                   "deg_a": DEG_A, "deg_b": DEG_B, "T": T_STEPS,
                   "every": 1, "omega0": OMEGA0, "kappa": "1-r"},
        "z0": [z0a, z0b],
        "spec": spec, "share_period_ge100": share100, "share_144band": share144,
        "chi": H["chi"], "dchi": [float(f"{x:.6g}") for x in dchi],
        "mid": H["mid"], "pr": H["pr"], "Rp": H["Rp"], "r": H["r"],
        "v": H["v"], "omega": H["omega"],
        "cloRe": H["cloRe"], "cloIm": H["cloIm"],
    }
    p = HERE / f"cr4_relative_{tag}_data_v1.json"
    p.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")),
                 encoding="utf-8")
    return out, p


def main() -> None:
    tags = sys.argv[1:] or list(CASES)
    rows = []
    for tag in tags:
        t0 = time.time()
        print(f"走行中 {tag} ({CASES[tag]['label']}) …", end="", flush=True)
        out, p = run_case(tag)
        chi = np.degrees(np.array(out["chi"]))
        Rp = np.array(out["Rp"], float)
        rows.append((tag, out, chi, Rp))
        print(f" 完了 {time.time()-t0:.1f}s  {p.stat().st_size//1024}KB")

    print()
    print(f"{'ケース':11} {'|z|':>8} {'PR初期':>8} {'Δθ範囲[°]':>18} "
          f"{'R′²変動[%]':>10} {'第1周期':>9} {'占有%':>7} {'144帯%':>8} {'r平均':>8}")
    print("-" * 100)
    for tag, out, chi, Rp in rows:
        pr = out["pr"]
        print(f"{tag:11} {out['z0'][0]:8.6f} {pr[0]:8.2f} "
              f"[{chi.min():+7.2f},{chi.max():+7.2f}] "
              f"{100*Rp.std()/Rp.mean():10.4f} {out['spec'][0]['period']:9.4f} "
              f"{out['spec'][0]['share']:7.2f} {out['share_144band']:8.4f} "
              f"{float(np.mean(out['r'])):8.6f}")
    print("\n※144帯% = 周期 129.6〜158.4 步（144±10%）の帯内最大ピークの占有。"
          "\n  復元力 a≈−ω²Δθ の固有周期は 2π/ω = 144 步。")


if __name__ == "__main__":
    main()
