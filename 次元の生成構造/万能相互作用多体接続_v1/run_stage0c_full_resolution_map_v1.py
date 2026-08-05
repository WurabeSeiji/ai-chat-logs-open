#!/usr/bin/env python3
"""GATE-0c: 分解能の全数地図 N=3..144——滑り・周波数比・親安定性のN依存構造

目的（2026-08-05 木原氏指示: N=1から144まで全て調べる）:
    GATE-0b（12点）では gcd(144)共約説・厳密単調説がともに棄却され、
    滑りは r~N^-1.18 の滑らかな減衰に見えた。全数地図でこの読みを検証し、
    算術構造（特定Nの外れ値・周波数の有理数吸着）を探す。
    N=1,2 は力学が退化（M=0 / K=0）するため除外し、その旨を記録する。

観測量（系列ごと・全てO(M)/step）:
    r_tail  = 末尾窓の一段残差中央値（滑り）
    O_tail  = 末尾窓の平面安定度中央値（Δ=10）
    ρ_tail  = 末尾窓の一歩位相進み ω=arg⟨Z(t),Z(t+1)⟩ の中央値 ÷ (π/72)
              ——エンジン時計との周波数比。有理数吸着があれば算術構造の直接証拠
    親系列のみ: crossing・t_launch・r(0) → 平衡分類
              （equilibrium-stable / equilibrium-burst / non-equilibrium）

系列（各N）:
    白色起源親 1本（v3生成器、シード2..9の最初の収束値を記録）
    セクター: N≤60 は全許容 k（2k≢0 mod N）、N>60 は決定的標本12本
        k ∈ {1,2,3,N-1,N-2,⌊N/6⌋,⌊N/4⌋,⌊N/3⌋,⌊2N/5⌋,⌊N/2⌋-1,⌊3N/5⌋,⌊2N/3⌋}
        （許容かつ非重複に整理）

事前固定の解析規則（事後変更禁止）:
    A1 冪則検定: log r_tail 対 log N の回帰（全セクター中央値）。
        GATE-0bの r~N^-1.18 が全域で持続するか。
    A2 外れ値規則: 各Nの log(セクター中央値 r) が移動窓（±3点）中央値から
        ±log(2) 超ずれる N を算術候補として全数列挙（方向問わず）。
    A3 周波数吸着: ρ の小分母有理数 p/q（q≤12）への最近接距離 d(ρ) を全系列で記録。
        d < 1e-4 の系列を「吸着」とし、吸着率のN依存を報告（探索的）。
    A4 親分類地図: N ごとの安定/不安定/非平衡の分布を全数報告（探索的・
        開始様式論文の未解明(i)への一次資料）。

方法: 力学無改変・read-only import（SHA記録）。T=4000、末尾窓1000。
    N単位でプロセス並列（8ワーカー）。乱数シードは GATE-0b と同一規約
    （親wp=91000、セクターwp=92000+k）。

使い方: python3 run_stage0c_full_resolution_map_v1.py
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SERIES_DIR = HERE.parent
ABL = SERIES_DIR / "第8論文_二段階seed除去による準安定相の因果分離" / "code" / "run_preliminary_seed_ablation_v1.py"
GEN3 = SERIES_DIR / "make_parent_white_managed_v1" / "make_parent_white_harmonics_n_only_v3.py"
NS = list(range(3, 145))
FULL_SECTOR_MAX = 60
T_LONG = 4000
DELTA_O = 10
TAIL = 1000
OMEGA_ENGINE = math.pi / 72.0
WORKERS = 8

_mods = {}


def _load():
    if _mods:
        return _mods["abl"], _mods["gen3"]
    import importlib.util
    spec = importlib.util.spec_from_file_location("abl_g0c", ABL)
    abl = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = abl
    spec.loader.exec_module(abl)
    spec2 = importlib.util.spec_from_file_location("gen3_g0c", GEN3)
    gen3 = importlib.util.module_from_spec(spec2)
    sys.modules[spec2.name] = gen3
    spec2.loader.exec_module(gen3)
    _mods["abl"] = abl
    _mods["gen3"] = gen3
    return abl, gen3


def orth_plane(Z):
    p = Z.real / np.linalg.norm(Z.real)
    q = Z.imag - (Z.imag @ p) * p
    nq = np.linalg.norm(q)
    if nq < 1e-300:
        return np.column_stack([p, np.zeros_like(p)])
    return np.column_stack([p, q / nq])


def run_one(abl, n, Z0, wp, track_f):
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(Z0))
    Z = Z0 / np.linalg.norm(Z0)
    if track_f:
        p0 = Z.real / np.linalg.norm(Z.real)
        q0 = Z.imag - (Z.imag @ p0) * p0
        q0 = q0 / np.linalg.norm(q0)
    rs = np.zeros(T_LONG)
    ws = np.zeros(T_LONG)
    fs = np.zeros(T_LONG + 1) if track_f else None
    Os = []
    B_prev = orth_plane(Z)
    r0 = None
    for t in range(1, T_LONG + 1):
        Znew, wp = abl.evolve(sys_lr, Z, wp)
        ip = np.conj(Z) @ Znew
        ph = ip / abs(ip) if abs(ip) > 0 else 1.0
        rr = float(np.linalg.norm(Znew - ph * Z))
        rs[t - 1] = rr
        ws[t - 1] = float(np.angle(ip))
        if r0 is None:
            r0 = rr
        Z = Znew
        if track_f:
            Zp = Z - p0 * (p0 @ Z) - q0 * (q0 @ Z)
            fs[t] = float(np.real(np.conj(Zp) @ Zp))
        if t % DELTA_O == 0:
            B = orth_plane(Z)
            Os.append(float(np.sum((B_prev.T @ B) ** 2) / 2.0))
            B_prev = B
    Os = np.array(Os)
    out = {"r_tail": float(np.median(rs[-TAIL:])),
           "O_tail": float(np.median(Os[-(TAIL // DELTA_O):])),
           "rho_tail": float(np.median(ws[-TAIL:]) / OMEGA_ENGINE),
           "r0": float(r0)}
    if track_f:
        crossing = next((t for t, f in enumerate(fs) if f > 0.05), None)
        below = np.nonzero(fs < 1e-20)[0]
        out["crossing"] = crossing
        out["t_launch"] = int(below.max()) if below.size else 0
        out["f_final"] = float(fs[-1])
    return out


def sector_list(n):
    allowed = [k for k in range(1, n) if (2 * k) % n != 0]
    if n <= FULL_SECTOR_MAX:
        return allowed
    cand = [1, 2, 3, n - 1, n - 2, n // 6, n // 4, n // 3,
            (2 * n) // 5, n // 2 - 1, (3 * n) // 5, (2 * n) // 3]
    out = []
    for k in cand:
        if k in allowed and k not in out:
            out.append(k)
    return out


def process_n(n):
    abl, gen3 = _load()
    t0 = time.time()
    m = n * (n - 1) // 2
    parent = None
    gen_seed = None
    for seed in range(2, 10):
        try:
            r = gen3.make_parent(n, seed=seed)
            parent = r
            gen_seed = seed
            break
        except Exception:
            continue
    if parent is None:
        return n, {"error": "make_parent failed seeds 2..9"}
    C = np.fft.fft(parent.relation_waves, axis=1) / n
    vp = parent.parent_vector / np.linalg.norm(parent.parent_vector)

    pres = run_one(abl, n, vp, np.random.default_rng(91000).normal(size=m), True)
    if pres["r0"] > 1e-8:
        pclass = "non-equilibrium"
    elif pres["crossing"] is not None:
        pclass = "equilibrium-burst"
    else:
        pclass = "equilibrium-stable"

    secs = {}
    for k in sector_list(n):
        Zk = C[:, k] / np.linalg.norm(C[:, k])
        wp = np.random.default_rng(92000 + k).normal(size=m)
        secs[str(k)] = run_one(abl, n, Zk, wp, False)
    med_r = float(np.median([s["r_tail"] for s in secs.values()]))
    med_O = float(np.median([s["O_tail"] for s in secs.values()]))
    return n, {"gen_seed": gen_seed, "n_sectors": len(secs),
                "parent": {**pres, "class": pclass},
                "sector_median_r": med_r, "sector_median_O": med_O,
                "sectors": secs, "runtime_sec": time.time() - t0}


def nearest_rational_dist(x, qmax=12):
    best = 1e9
    bestpq = None
    for q in range(1, qmax + 1):
        p = round(x * q)
        d = abs(x - p / q)
        if d < best:
            best = d
            bestpq = (int(p), q)
    return best, bestpq


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    t0 = time.time()
    print(f"GATE-0c 全数地図 N=3..144  ABL {sha256(ABL)[:16]}…  GEN3 {sha256(GEN3)[:16]}…")
    results = {}
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        for n, res in ex.map(process_n, sorted(NS, reverse=True)):
            results[str(n)] = res
            if "error" in res:
                print(f"  N={n}: {res['error']}")
            else:
                print(f"  N={n:3d}: 親={res['parent']['class']:20s} "
                      f"sec中央値 r={res['sector_median_r']:.2e} O={res['sector_median_O']:.6f} "
                      f"({res['runtime_sec']:.0f}s)", flush=True)

    ok = {int(n): r for n, r in results.items() if "error" not in r}
    ns = np.array(sorted(ok))
    rs = np.array([ok[n]["sector_median_r"] for n in ns])
    A = np.vstack([np.log(ns), np.ones_like(ns, float)]).T
    coef, _, _, _ = np.linalg.lstsq(A, np.log(rs), rcond=None)
    pred = A @ coef
    ss = 1 - np.sum((np.log(rs) - pred) ** 2) / np.sum((np.log(rs) - np.log(rs).mean()) ** 2)

    logr = np.log(rs)
    outliers = []
    for i, n in enumerate(ns):
        lo, hi = max(0, i - 3), min(len(ns), i + 4)
        local = np.median(np.concatenate([logr[lo:i], logr[i + 1:hi]]))
        if abs(logr[i] - local) > math.log(2):
            outliers.append({"N": int(n), "dev_factor": float(math.exp(logr[i] - local))})

    absorbed = []
    for n in ns:
        for k, s in ok[n]["sectors"].items():
            d, pq = nearest_rational_dist(s["rho_tail"])
            if d < 1e-4:
                absorbed.append({"N": int(n), "k": int(k), "rho": s["rho_tail"],
                                  "pq": list(pq)})
    pclasses = {}
    for n in ns:
        pclasses.setdefault(ok[n]["parent"]["class"], []).append(int(n))

    print(f"\nA1 冪則: r ~ N^{coef[0]:.3f} (R²={ss:.4f})")
    print(f"A2 外れ値（±3窓・2倍規則）: {outliers if outliers else 'なし'}")
    print(f"A3 有理数吸着（q≤12, d<1e-4）: {len(absorbed)}系列")
    print(f"A4 親分類: " + ", ".join(f"{k}:{len(v)}本" for k, v in pclasses.items()))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
    plt.rcParams["axes.unicode_minus"] = False
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    axes[0].loglog(ns, rs, "o-", ms=3, lw=0.5, color="#4C78A8")
    axes[0].loglog(ns, np.exp(pred), "k--", lw=0.8,
                    label=f"fit N^{coef[0]:.2f}")
    for o in outliers:
        axes[0].axvline(o["N"], color="red", lw=0.5, alpha=0.5)
    axes[0].set_ylabel("sector median tail r")
    axes[0].legend(fontsize=8)
    Ovals = [ok[n]["sector_median_O"] for n in ns]
    axes[1].plot(ns, Ovals, "o-", ms=3, lw=0.5, color="#4C78A8")
    axes[1].axhline(0.999, color="k", ls=":", lw=0.8)
    axes[1].set_ylabel("sector median tail O")
    cls_color = {"equilibrium-stable": "#2E7D32", "equilibrium-burst": "#7F7F7F",
                 "non-equilibrium": "#E45756"}
    for n in ns:
        c = ok[n]["parent"]["class"]
        axes[2].plot(n, {"equilibrium-stable": 0, "equilibrium-burst": 1,
                          "non-equilibrium": 2}[c], "s", ms=4, color=cls_color[c])
    axes[2].set_yticks([0, 1, 2])
    axes[2].set_yticklabels(["eq-stable", "eq-burst", "non-eq"])
    axes[2].set_ylabel("parent class")
    axes[2].set_xlabel("N")
    axes[0].set_title("GATE-0c: full resolution map N=3..144")
    fig.tight_layout()
    fig.savefig(HERE / "fig_gate0c_full_map.png", dpi=130)
    plt.close(fig)

    out = {"NS": NS, "excluded": {"N=1": "M=0", "N=2": "K=0 (single edge)"},
           "imports": {"abl": sha256(ABL), "generator_v3": sha256(GEN3)},
           "criteria": {"T_LONG": T_LONG, "TAIL": TAIL, "FULL_SECTOR_MAX": FULL_SECTOR_MAX,
                         "rules": ["A1 power-law fit", "A2 outlier ±3-window factor-2",
                                    "A3 rational absorption q<=12 d<1e-4",
                                    "A4 parent class map"]},
           "analysis": {"A1_exponent": float(coef[0]), "A1_r2": float(ss),
                         "A2_outliers": outliers, "A3_absorbed": absorbed,
                         "A4_parent_classes": {k: v for k, v in pclasses.items()}},
           "results": results, "runtime_sec": time.time() - t0}
    (HERE / "gate0c_full_resolution_map_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
