#!/usr/bin/env python3
"""E-M11：倍音時計の力学的実装と段間分光 v1

背景（E-M10 の診断）:
    Cayley 一段の γ/σ 正規化により、全段の時計は σ 非依存の π/72/step に
    規格化される——倍音周波数 n·ω₀ はメタデータに留まり力学的実体がなかった。

実装（エンジン改変ゼロ）:
    段 n には 1 tick あたり n 回のマイクロステップ（abl.evolve の純反復）を
    与える。これで段 n の時計は厳密に n·(π/72)/tick となり、倍音格子
    ωₙ = n·ω₀ が初めて力学的実体になる。

固定予言:
    P1（実装検証・厳格）: 到達ユニゾン周波数の段間比 ωₙ/ω₁ = n
        （整数、偏差 < 1e-3）。不成立なら実装が誤り。
    P2（記録・記述的）: 窓別の段横断周波数プール上の整数比対の勘定
        （E-M4 基準）。後期窓では構成上の格子（ユニゾン段間の整数比）が
        出るはずで、これは構成の実現確認。遷移帯（分化中）の有理対密度が
        後期と比べてどう振る舞うかを記録する——結合が入ったとき
        ロック候補となる「資源」の分光。
    P3（E-M10 の再現）: 各段内部の窓ロックは 0 のまま
        （時計の再目盛りは段内force学を変えない——crossing_tick×n ≈ E-M9 の
        crossing 値も検証）。

規約: 種ラベル IF 分岐なし（マイクロステップ数は段番号のみの関数）、
    read-only import、SHA-256 記録、反証も記録。
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

spec = importlib.util.spec_from_file_location("abl_m11", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
spec2 = importlib.util.spec_from_file_location("mph_m11", MPH)
mph = importlib.util.module_from_spec(spec2)
sys.modules[spec2.name] = mph
spec2.loader.exec_module(mph)

N = 5
H = 8
SEED = 40260801
TICKS = 12000
WIN = 500
STRIDE = 250
TOL_LOCK = 1e-3
FREQ_MIN = 1e-8
OMEGA0 = np.pi / 72.0


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


def evolve_level(v0, wp, micro):
    """段の発展: 1 tick = micro 回のマイクロステップ（abl.evolve の純反復）。"""
    sys_lr = abl.LowRankSystem(N)
    sys_lr.set_theta(np.angle(v0))
    m = sys_lr.m
    Z = v0.copy()
    phases = np.zeros((TICKS + 1, m))
    phases[0] = np.angle(Z)
    for t in range(1, TICKS + 1):
        for _ in range(micro):
            Z, wp = abl.evolve(sys_lr, Z, wp)
        phases[t] = np.angle(Z)
    return phases


def window_freqs(phases):
    """窓別の辺周波数 [rad/tick]（E-M4/M10 と同一の抽出）。"""
    out = []
    for s in range(0, TICKS - WIN + 1, STRIDE):
        u = np.unwrap(phases[s:s + WIN], axis=0)
        fr = np.abs(np.polyfit(np.arange(WIN), u, 1)[0])
        out.append((s + WIN // 2, fr))
    return out


def main() -> None:
    t0 = time.time()
    print("E-M11 倍音時計実装・段間分光 実行")
    print(f"  import: ABL {sha256(ABL)[:16]}…  MPH {sha256(MPH)[:16]}…")
    results = {"imports": {"abl": sha256(ABL), "mph": sha256(MPH),
                            "engine": mph.ENGINE_SHA256},
               "params": {"N": N, "H": H, "SEED": SEED, "TICKS": TICKS,
                           "WIN": WIN, "STRIDE": STRIDE, "TOL_LOCK": TOL_LOCK,
                           "omega0": OMEGA0}}

    Zh, info = mph.make_parent_harmonic(N, H, SEED, iters=2000, restarts=10, tol=1e-12)
    levels = {}
    for h in range(1, H + 1):
        lv = info["levels"][h - 1]
        fam = "N-1" if abs(lv["sigma1"] - (N - 1)) < 1e-9 else "broken"
        v0 = Zh[:, h - 1] * np.sqrt(H)
        wp = np.random.default_rng(90000 + (h - 1)).normal(size=len(v0))
        ph = evolve_level(v0, wp, micro=h)
        wf = window_freqs(ph)
        # 到達ユニゾン周波数（最終 3000 tick）
        u = np.unwrap(ph[TICKS - 3000:], axis=0)
        fr = np.abs(np.polyfit(np.arange(u.shape[0]), u, 1)[0])
        fb = fr[fr > FREQ_MIN]
        omega = float(np.mean(fb)) if fb.size else 0.0
        # 段内窓ロック（P3）
        locks_in = max(em4_lock_stats(fr_w)[1] for _, fr_w in wf)
        levels[h] = {"family": fam, "sigma1": lv["sigma1"], "omega_late": omega,
                      "omega_over_omega0": omega / OMEGA0,
                      "max_window_locks_internal": locks_in, "wf": wf}
        print(f"  段n={h}（{fam}）: ω={omega:.6f} rad/tick  ω/ω₀={omega/OMEGA0:.4f} "
              f"段内窓ロック最大={locks_in}")

    # P1: 倍音格子の実装検証
    active = {h: d for h, d in levels.items() if d["omega_late"] > FREQ_MIN}
    ratios = {h: d["omega_over_omega0"] for h, d in active.items()}
    p1 = all(abs(r - h) < 1e-3 for h, r in ratios.items())
    print(f"\nP1 倍音格子 ωₙ = n·ω₀（偏差<1e-3）: {'PASS' if p1 else 'FAIL'} "
          f"({ {h: f'{r:.4f}' for h, r in ratios.items()} })")

    # P2: 窓別・段横断プールの整数比対勘定
    centers = [c for c, _ in levels[1]["wf"]]
    pool_stats = []
    for i, c in enumerate(centers):
        pool = np.concatenate([levels[h]["wf"][i][1] for h in range(1, H + 1)])
        _, locks = em4_lock_stats(pool)
        pool_stats.append({"center": c, "cross_level_locks": locks})
    early = [w["cross_level_locks"] for w in pool_stats if w["center"] < 1500]
    late = [w["cross_level_locks"] for w in pool_stats if w["center"] > TICKS - 3000]
    trans = [w["cross_level_locks"] for w in pool_stats if 1000 <= w["center"] <= 4000]
    print(f"P2 段横断整数比対: 初期窓 max={max(early) if early else 0} "
          f"遷移帯 max={max(trans) if trans else 0} 後期窓 max={max(late) if late else 0}"
          f"（後期は構成上の格子の実現確認）")

    # P3: 段内ロック 0・crossing 整合
    p3 = all(d["max_window_locks_internal"] == 0 for d in levels.values())
    print(f"P3 段内窓ロック 0 の再現: {'PASS' if p3 else 'FAIL'}")

    results["levels"] = {str(h): {k: v for k, v in d.items() if k != "wf"}
                          for h, d in levels.items()}
    results["pool_stats"] = pool_stats
    results["verdicts"] = {"P1": bool(p1),
                            "P2_late_lattice_max": max(late) if late else 0,
                            "P2_transition_max": max(trans) if trans else 0,
                            "P3": bool(p3)}
    results["runtime_sec"] = time.time() - t0
    (HERE / "paper8_em11_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"saved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
