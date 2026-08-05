#!/usr/bin/env python3
"""GATE-0b: 分解能Nの掃引——セクター収束のN選択性の検定

目的（2026-08-05 木原氏仮説）:
    GATE-0でN=5セクターのみが開ゲート署名（末尾 O≥0.999）に収束しなかった。
    木原氏仮説: 系はNに対して対称ではない（U^n=I）。最低解像度は任意ではなく
    力学的に選択される（E8・α問題・結晶学的制限と整合する読み）。

機構候補（実行前に固定する判別予言の根拠）:
    エンジンの一段はCayley変換で厳密に π/72 回転（γ=tan(π/144)）
    ＝144step周期の内部時計。144=2⁴·3²。ロックは共約性を要求するなら、
    Nの素因数が{2,3}に収まるN（144と共約）はロック可能、
    素因数5以上を含むNはロック不能で滑りが持続するはず。
    実測済み2点: N=40(gcd=8)収束 / N=5(gcd=1)滑り——読みと整合。

判別予言（実行前固定・事後変更禁止）:
    H_selective（木原氏仮説）: セクター末尾の滑りはNに非単調で、
        素因数⊂{2,3}群 {6,8,9,12,16,18} は末尾中央値O ≥ 0.999、
        素因数≥5を含む群 {5,7,10,11,13,14} は末尾中央値O < 0.999
        （群内のNの大小に依らず群で分かれる）
    H_artifact（帰無仮説）: 滑りはNとともに単調減少し、素因数構成に依らない
    判定: 各Nの全許容セクターの末尾中央値O（系列毎に中央値→N毎に中央値）を
        両群で比較。群分離かつ非単調ならH_selective、単調ならH_artifact。
    探索的記録（合否外）: 末尾中央値r、セクターkごとの構造（kのmod N位数）、
        gcd(N,144)との順序相関。

方法: GATE-0と同一の計装（力学無改変・O(M)/step）。各Nで v3生成器
    （seed=2、非収束時は3,4,…の最初の収束シードを記録）から許容セクター
    （2k≢0 mod N）を全数注入、T=6000、末尾1000stepのr・O中央値を記録。

使い方: python3 run_stage0b_resolution_sweep_v1.py
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).resolve().parent
SERIES_DIR = HERE.parent
PAPER8 = SERIES_DIR / "第8論文_二段階seed除去による準安定相の因果分離"
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"
GEN3 = SERIES_DIR / "make_parent_white_managed_v1" / "make_parent_white_harmonics_n_only_v3.py"
NS = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18]
SMOOTH = {6, 8, 9, 12, 16, 18}      # 素因数 ⊂ {2,3}（144と共約な群）
ROUGH = {5, 7, 10, 11, 13, 14}      # 素因数 ≥5 を含む群
T_LONG = 6000
DELTA_O = 10
TAIL = 1000
O_LOCK = 0.999

spec = importlib.util.spec_from_file_location("abl_g0b", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("gen3_g0b", GEN3)
gen3 = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = gen3
spec2.loader.exec_module(gen3)


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def orth_plane(Z):
    p = Z.real / np.linalg.norm(Z.real)
    q = Z.imag - (Z.imag @ p) * p
    nq = np.linalg.norm(q)
    if nq < 1e-300:
        return np.column_stack([p, np.zeros_like(p)])
    return np.column_stack([p, q / nq])


def run_series(n, Z0, wp):
    sys_lr = abl.LowRankSystem(n)
    sys_lr.set_theta(np.angle(Z0))
    Z = Z0 / np.linalg.norm(Z0)
    rs = np.zeros(T_LONG)
    Os = []
    B_prev = orth_plane(Z)
    for t in range(1, T_LONG + 1):
        Znew, wp = abl.evolve(sys_lr, Z, wp)
        ip = np.conj(Z) @ Znew
        ph = ip / abs(ip) if abs(ip) > 0 else 1.0
        rs[t - 1] = float(np.linalg.norm(Znew - ph * Z))
        Z = Znew
        if t % DELTA_O == 0:
            B = orth_plane(Z)
            Os.append(float(np.sum((B_prev.T @ B) ** 2) / 2.0))
            B_prev = B
    Os = np.array(Os)
    n_tail_O = TAIL // DELTA_O
    return {"tail_median_r": float(np.median(rs[-TAIL:])),
            "tail_median_O": float(np.median(Os[-n_tail_O:]))}


def make_parent_with_fallback(n):
    for seed in range(2, 10):
        try:
            r = gen3.make_parent(n, seed=seed)
            return r, seed
        except Exception:
            continue
    raise RuntimeError(f"make_parent failed for N={n} seeds 2..9")


def main() -> None:
    t0 = time.time()
    print(f"GATE-0b 分解能掃引 N={NS}  ABL {sha256(ABL)[:16]}…  GEN3 {sha256(GEN3)[:16]}…")
    results = {}
    for n in NS:
        r, used_seed = make_parent_with_fallback(n)
        C = np.fft.fft(r.relation_waves, axis=1) / n
        m = n * (n - 1) // 2
        ks = [k for k in range(n) if (2 * k) % n != 0]
        secs = {}
        for k in ks:
            Zk = C[:, k] / np.linalg.norm(C[:, k])
            wp = np.random.default_rng(92000 + k).normal(size=m)
            secs[f"k{k}"] = run_series(n, Zk, wp)
        med_O = float(np.median([s["tail_median_O"] for s in secs.values()]))
        med_r = float(np.median([s["tail_median_r"] for s in secs.values()]))
        locked = med_O >= O_LOCK
        group = "smooth" if n in SMOOTH else "rough"
        results[str(n)] = {"gen_seed": used_seed, "n_sectors": len(ks),
                            "median_tail_O": med_O, "median_tail_r": med_r,
                            "locked": locked, "group": group,
                            "gcd144": int(np.gcd(n, 144)), "sectors": secs}
        print(f"  N={n:2d} [{group}] gcd(144)={np.gcd(n,144):2d} "
              f"セクター{len(ks):2d}本: 末尾中央値 O={med_O:.6f} r={med_r:.2e} "
              f"{'LOCK' if locked else 'SLIDE'}")

    smooth_lock = [results[str(n)]["locked"] for n in NS if n in SMOOTH]
    rough_lock = [results[str(n)]["locked"] for n in NS if n in ROUGH]
    h_selective = all(smooth_lock) and not any(rough_lock)
    os_by_n = [(n, results[str(n)]["median_tail_O"]) for n in NS]
    monotone = all(os_by_n[i + 1][1] >= os_by_n[i][1] - 1e-6
                   for i in range(len(os_by_n) - 1))
    print(f"\n判定: H_selective（群完全分離）={h_selective}  "
          f"H_artifact（N単調）={monotone}")

    fig, ax = plt.subplots(figsize=(9, 5))
    for n in NS:
        d = results[str(n)]
        c = "#4C78A8" if d["group"] == "smooth" else "#E45756"
        ax.plot(n, d["median_tail_O"], "o", ms=9, color=c)
        for s in d["sectors"].values():
            ax.plot(n, s["tail_median_O"], ".", ms=3, color=c, alpha=0.4)
    ax.axhline(O_LOCK, color="k", ls=":", lw=0.8, label=f"lock threshold O={O_LOCK}")
    ax.plot([], [], "o", color="#4C78A8", label="prime factors ⊂ {2,3} (÷144)")
    ax.plot([], [], "o", color="#E45756", label="contains prime ≥5")
    ax.set_xlabel("N (resolution)")
    ax.set_ylabel("tail median plane stability O")
    ax.set_title("GATE-0b: sector convergence vs resolution N")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "fig_gate0b_resolution_sweep.png", dpi=130)
    plt.close(fig)

    out = {"NS": NS, "SMOOTH": sorted(SMOOTH), "ROUGH": sorted(ROUGH),
           "imports": {"abl": sha256(ABL), "generator_v3": sha256(GEN3)},
           "criteria": {"T_LONG": T_LONG, "TAIL": TAIL, "O_LOCK": O_LOCK,
                         "H_selective": "smooth all locked AND rough none locked",
                         "H_artifact": "median_tail_O monotone in N"},
           "verdict": {"H_selective": bool(h_selective),
                        "H_artifact_monotone": bool(monotone)},
           "results": results, "runtime_sec": time.time() - t0}
    (HERE / "gate0b_resolution_sweep_v1.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({out['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
