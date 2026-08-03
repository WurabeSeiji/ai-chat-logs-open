#!/usr/bin/env python3
"""E-M10：倍音海軌道への物質誕生予兆分析 v1

目的（本来の目的への復帰）:
    倍音対応 make_parent の閉包から出発した E-M9 軌道（N=5、対照＋倍音海8段）に、
    物質誕生の予兆観測（E-M8 の確定済み観測系）を適用する。
    (1) 窓別の整数比ロック（E-M4 111-116行と同一基準・比≥2・偏差<1e-3）
    (2) 透明度（隣接対 |sinΔθ|<0.05 の割合）と participation ratio
    (3) 新観測: 各段の到達ユニゾン周波数 ω_n の段間比較——段間比が整数比を
        成すか（倍音レジスタで初めて問える予兆。E-M4 基準で判定）

固定予言:
    P1（E-M8 の再現・対照）: 対照軌道は窓別ロック 0、遷移帯で周波数分化と
        透明度閉止、PR は等分配へ（E-M8 実測の再現）
    P2（開いた問い）: 破れ族の倍音閉包生まれの軌道が、窓別ロック > 0 または
        E-M8 と質的に異なる分光を見せるか
    P3（開いた問い）: 段間ユニゾン周波数比 ω_m/ω_n に非自明整数比（E-M4 基準）
        が現れるか

規約: 種ラベル IF 分岐なし・abl read-only import・SHA-256 記録・反証も記録。
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAPER8 = HERE.parent
REPO = PAPER8.parent.parent
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
MPH = REPO / "次元の生成構造" / "make_parent_harmonic_unit_v1" / "make_parent_harmonic_v1.py"

spec = importlib.util.spec_from_file_location("abl_m10", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("mph_m10", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)
import run_n_scaling_lowrank_v1 as eng

N = 5
H = 8
SEED = 40260801
XMAX = 12000
WIN = 500
STRIDE = 250
TOL_LOCK = 1e-3
FREQ_MIN = 1e-8
TRANSPARENT_EPS = 0.05


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def em4_lock_stats(fb: np.ndarray):
    fb = fb[fb > FREQ_MIN]
    if fb.size > 1:
        r = fb[:, None] / np.maximum(fb[None, :], 1e-30)
        rmax = np.maximum(r, 1 / np.maximum(r, 1e-30))
        max_dev = float(np.max(np.abs(r[r >= 1] - 1)))
        pr = np.round(rmax)
        locks = int(np.sum((pr >= 2) & (np.abs(rmax - pr) < TOL_LOCK)) // 2)
    else:
        max_dev, locks = 0.0, 0
    return max_dev, locks


def adjacency_pairs(sys_lr):
    ea, eb = sys_lr.ea, sys_lr.eb
    m = sys_lr.m
    pairs = []
    for i in range(m):
        share = (ea == ea[i]) | (ea == eb[i]) | (eb == ea[i]) | (eb == eb[i])
        for j in np.nonzero(share)[0]:
            if j > i:
                pairs.append((i, int(j)))
    return np.array(pairs, dtype=np.int64)


def analyze(v0, wp, label):
    sys_lr = abl.LowRankSystem(N)
    sys_lr.set_theta(np.angle(v0))
    m = sys_lr.m
    pairs = adjacency_pairs(sys_lr)
    Z = v0.copy()
    phases = np.zeros((XMAX + 1, m))
    prs = np.zeros(XMAX + 1)
    phases[0] = np.angle(Z); prs[0] = eng.participation_ratio(Z)
    for t in range(1, XMAX + 1):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        phases[t] = np.angle(Z); prs[t] = eng.participation_ratio(Z)

    windows = []
    for s in range(0, XMAX - WIN + 1, STRIDE):
        w = slice(s, s + WIN)
        u = np.unwrap(phases[w], axis=0)
        fr = np.abs(np.polyfit(np.arange(WIN), u, 1)[0])
        max_dev, locks = em4_lock_stats(fr)
        dth = phases[w][:, pairs[:, 0]] - phases[w][:, pairs[:, 1]]
        sabs = np.abs(np.sin(dth))
        windows.append({"center": s + WIN // 2, "locks": locks,
                        "max_ratio_dev": max_dev,
                        "transparent_frac": float((sabs < TRANSPARENT_EPS).mean()),
                        "pr_mean": float(prs[w].mean())})
    # 到達ユニゾン周波数（最終3000step、E-M1/M4 と同一の抽出）
    u = np.unwrap(phases[XMAX - 3000:], axis=0)
    fr = np.abs(np.polyfit(np.arange(u.shape[0]), u, 1)[0])
    fb = fr[fr > FREQ_MIN]
    omega = float(np.mean(fb)) if fb.size else 0.0
    max_locks = max(w["locks"] for w in windows)
    peak_dev = max(w["max_ratio_dev"] for w in windows)
    tmin = min(w["transparent_frac"] for w in windows)
    print(f"  [{label}] 窓ロック最大={max_locks} 分化ピーク={peak_dev:.3f} "
          f"透明度最小={tmin:.3f} 到達ω={omega:.6f}")
    return {"windows": windows, "omega_late": omega,
            "max_window_locks": max_locks, "peak_ratio_dev": peak_dev,
            "min_transparent_frac": tmin}


def main() -> None:
    t0 = time.time()
    print("E-M10 物質誕生予兆分析（倍音海軌道）実行")
    print(f"  import: ABL {sha256(ABL)[:16]}…  MPH {sha256(MPH)[:16]}…")
    results = {"imports": {"abl": sha256(ABL), "mph": sha256(MPH),
                            "engine": mph.ENGINE_SHA256},
               "params": {"N": N, "H": H, "SEED": SEED, "XMAX": XMAX, "WIN": WIN,
                           "STRIDE": STRIDE, "TOL_LOCK": TOL_LOCK}}

    series = {}
    _, v, _, _, _, _, _, Z0, wp0 = abl.build_init(N, False)
    series["control"] = {"family": "control",
                          **analyze(Z0, wp0.copy(), "control")}

    Zh, info = mph.make_parent_harmonic(N, H, SEED, iters=2000, restarts=10, tol=1e-12)
    for h in range(H):
        lv = info["levels"][h]
        fam = "N-1" if abs(lv["sigma1"] - (N - 1)) < 1e-9 else "broken"
        v0 = Zh[:, h] * np.sqrt(H)
        wp = np.random.default_rng(90000 + h).normal(size=len(v0))
        series[f"n{h+1}"] = {"family": fam, "sigma1": lv["sigma1"],
                              **analyze(v0, wp, f"段n={h+1}（{fam}）")}

    # P3: 段間ユニゾン周波数比（E-M4 基準）
    omegas = {k: s["omega_late"] for k, s in series.items() if s["omega_late"] > FREQ_MIN}
    keys = sorted(omegas)
    fb = np.array([omegas[k] for k in keys])
    inter_dev, inter_locks = em4_lock_stats(fb)
    print("\n段間到達ω一覧:", {k: f"{omegas[k]:.6f}" for k in keys})
    print(f"段間比 最大ずれ={inter_dev:.3e} 非自明整数比ロック={inter_locks}")

    p1 = (series["control"]["max_window_locks"] == 0
          and series["control"]["peak_ratio_dev"] > 0.1
          and series["control"]["min_transparent_frac"] < 0.05)
    p2 = any(s["max_window_locks"] > 0 for k, s in series.items()
             if s["family"] == "broken")
    p3 = inter_locks > 0
    print(f"\nP1 対照=E-M8再現（ロック0・分化・透明度閉止）: {'PASS' if p1 else 'FAIL'}")
    print(f"P2 倍音閉包生まれの軌道に窓ロック: {'あり' if p2 else 'なし（予兆なし——反証記録）'}")
    print(f"P3 段間ω比に非自明整数比: {'あり' if p3 else 'なし（全段ユニゾンの再現）'}")

    results["series"] = series
    results["inter_level"] = {"omegas": omegas, "max_dev": inter_dev,
                               "nontrivial_locks": inter_locks}
    results["verdicts"] = {"P1": bool(p1), "P2": bool(p2), "P3": bool(p3)}
    results["runtime_sec"] = time.time() - t0
    (HERE / "paper8_em10_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
