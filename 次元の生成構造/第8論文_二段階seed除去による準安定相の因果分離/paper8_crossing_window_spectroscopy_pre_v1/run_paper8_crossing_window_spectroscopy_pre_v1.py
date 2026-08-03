#!/usr/bin/env python3
"""E-M8：交差窓スペクトル診断——インフレーション停止機構と「見落とされた過渡生成」の検査 v1

背景（2026-08-04 木原氏指摘）:
    これまでの周波数・ロック測定（E-M1/M4/M5）は全て最終 2000-3000 step
    （沈静化後）の窓でのみ行われ、crossing の最中——位相が乱れ θ が開き、
    衝突が起こりうる唯一の窓——を一度も見ていない。標準宇宙論では
    インフレーションは粒子生成（再加熱）によって終わる。この系の停止機構も
    「倍音の重合＝過渡的凝縮」である可能性があり、生成は過渡的に起こって
    いるが生き残れないだけかもしれない。傍証: E-M4 副産物（結合倍音が
    crossing 時刻を桁で動かす: q2 1167→200, q3 →6092）。

対立仮説（実行前固定）:
    H-A（停止＝倍音重合・再加熱型）: crossing 近傍の窓で過渡的な非自明整数比
        ロック（E-M4 基準）が出現する、または振幅が少数辺に凝縮する
        （PR の急落）。かつ最終窓ではロック 0（溶解）。
    H-B（停止＝運動学的飽和、現行記述）: 全窓でロック 0、PR に遷移帯の
        急落なし——停止前後でスペクトル構造に質的変化はない。

判定基準（実行前固定）:
    遷移帯 = 窓中心 ∈ [crossing−500, crossing+GUARD]。後期帯 = 最後の2窓。
    ロック判定 = E-M4 111-116行と同一（比≥2・丸め偏差<1e-3、周波数>1e-8）。
    凝縮判定 = PR_min(遷移帯) / median(PR(crossing前帯)) < 0.5。
    H-A 支持 ⟺ （遷移帯ロック最大値>0 かつ 後期帯ロック=0）または 凝縮判定成立。
    H-B 支持 ⟺ 遷移帯ロック=0 かつ 凝縮判定不成立。
    それ以外は「混在」として生データごと記録する。

観測量（無名・機械的、種ラベルなし）:
    (1) 滑走窓（幅500・刻み250）ごとの辺別瞬間周波数（位相 unwrap の勾配、
        E-M1/M4 と同一の抽出法）→ 窓ごとのロック数・最大比ずれ・周波数広がり
    (2) participation ratio PR(t)（エンジン内蔵 participation_ratio、毎step）
    (3) 透明度: 隣接対（頂点共有）の |sin Δθ| の窓平均と、
        |sin Δθ|<0.05 の対の割合（海の自己透明性の開閉）
    (4) f(t)（既定の拡大座標、crossing 検出は f>0.05 でエンジンと同一）

再現性:
    軌道生成は abl.build_init(n, initial_seed=False)＋abl.evolve の read-only
    import（E-M1 と同一、乱数種は build_init 内で固定 40260722+1000n）。
    import 元3ファイルの SHA-256 を結果 JSON に記録。新しい機構・介入は一切
    注入しない——既存の決定論的走行を、正しい窓で初めて観測するだけである。

対象: N=5（crossing≈1166, E-M1 と同一）主測定、N=40（crossing≈2011, E-M3 と
    同一初期条件）確認。XMAX=12000 共通。
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
ABL = PAPER8 / "code" / "run_preliminary_seed_ablation_v1.py"

spec = importlib.util.spec_from_file_location("abl_m8", ABL)
abl = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = abl
spec.loader.exec_module(abl)
import run_n_scaling_lowrank_v1 as eng  # abl が sys.path を通済み

XMAX = 12000
WIN = 500
STRIDE = 250
TOL_LOCK = 1e-3            # E-M4 と同一
FREQ_MIN = 1e-8            # E-M4 と同一
TRANSPARENT_EPS = 0.05
PR_CONDENSE_RATIO = 0.5
PRE_ZONE_MARGIN = 500      # 遷移帯の開始 = crossing − 500


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def em4_lock_stats(fb: np.ndarray):
    """E-M4 111-116行と同一の判定。"""
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
    """頂点共有の辺対リスト（E-M4 dense_setup と同一の共有判定）。"""
    ea, eb = sys_lr.ea, sys_lr.eb
    m = sys_lr.m
    pairs = []
    for i in range(m):
        share = (ea == ea[i]) | (ea == eb[i]) | (eb == ea[i]) | (eb == eb[i])
        for j in np.nonzero(share)[0]:
            if j > i:
                pairs.append((i, int(j)))
    return np.array(pairs, dtype=np.int64)


def run_case(n: int):
    sys_lr, v, B_p1, B_rot, B0, p, q, Z, wp = abl.build_init(n, initial_seed=False)
    m = sys_lr.m
    pairs = adjacency_pairs(sys_lr)

    def fval(Zv):
        Zp_ = Zv - p * (p @ Zv) - q * (q @ Zv)
        return float(np.real(np.conj(Zp_) @ Zp_)) / float(np.real(np.conj(Zv) @ Zv))

    phases = np.zeros((XMAX + 1, m), dtype=np.float64)
    amps = np.zeros((XMAX + 1, m), dtype=np.float64)
    prs = np.zeros(XMAX + 1)
    fs = np.zeros(XMAX + 1)
    phases[0] = np.angle(Z); amps[0] = np.abs(Z)
    prs[0] = eng.participation_ratio(Z); fs[0] = fval(Z)
    crossing = None
    for t in range(1, XMAX + 1):
        Z, wp = abl.evolve(sys_lr, Z, wp)
        phases[t] = np.angle(Z); amps[t] = np.abs(Z)
        prs[t] = eng.participation_ratio(Z); fs[t] = fval(Z)
        if crossing is None and fs[t] > 0.05:
            crossing = t
    guard_end = (crossing + abl.GUARD) if crossing is not None else None

    # 滑走窓解析
    windows = []
    starts = list(range(0, XMAX - WIN + 1, STRIDE))
    for s in starts:
        w = slice(s, s + WIN)
        center = s + WIN // 2
        u = np.unwrap(phases[w], axis=0)
        fr = np.abs(np.polyfit(np.arange(WIN), u, 1)[0])
        max_dev, locks = em4_lock_stats(fr)
        fb = fr[fr > FREQ_MIN]
        spread = float(np.std(fb) / np.mean(fb)) if fb.size > 1 else 0.0
        dth = phases[w][:, pairs[:, 0]] - phases[w][:, pairs[:, 1]]
        sabs = np.abs(np.sin(dth))
        windows.append({
            "center": center,
            "locks": locks,
            "max_ratio_dev": max_dev,
            "freq_spread_rel": spread,
            "freq_mean": float(np.mean(fb)) if fb.size else 0.0,
            "pr_min": float(prs[w].min()),
            "pr_mean": float(prs[w].mean()),
            "transparency_mean_sin": float(sabs.mean()),
            "transparent_pair_frac": float((sabs < TRANSPARENT_EPS).mean()),
        })

    # 帯域の切り分けと判定（固定基準）
    verdict = {"crossing": crossing, "guard_end": guard_end}
    if crossing is not None:
        t_lo, t_hi = crossing - PRE_ZONE_MARGIN, guard_end
        trans = [w for w in windows if t_lo <= w["center"] <= t_hi]
        pre = [w for w in windows if w["center"] < t_lo]
        late = windows[-2:]
        locks_trans = max((w["locks"] for w in trans), default=0)
        locks_late = max((w["locks"] for w in late), default=0)
        pr_min_trans = min((w["pr_min"] for w in trans), default=float("nan"))
        pr_pre_med = float(np.median([w["pr_mean"] for w in pre])) if pre else float("nan")
        condense = bool(pre and pr_min_trans / pr_pre_med < PR_CONDENSE_RATIO)
        ha = (locks_trans > 0 and locks_late == 0) or condense
        hb = (locks_trans == 0) and (not condense)
        verdict.update({
            "locks_transition_max": locks_trans,
            "locks_late_max": locks_late,
            "pr_min_transition": pr_min_trans,
            "pr_pre_median": pr_pre_med,
            "condensation": condense,
            "H_A_supported": bool(ha),
            "H_B_supported": bool(hb),
            "mixed": bool(not ha and not hb),
        })
    return windows, verdict, m


def print_timeline(label, windows, verdict):
    cr = verdict["crossing"]
    print(f"\n[{label}] crossing={cr} guard_end={verdict['guard_end']}")
    print("  窓中心  ロック 比ずれ最大   周波数広がり  PR_min   透明対割合  |sinΔθ|平均")
    for w in windows:
        zone = ""
        if cr is not None:
            if cr - PRE_ZONE_MARGIN <= w["center"] <= verdict["guard_end"]:
                zone = " ◀遷移帯"
        print(f"  {w['center']:6d}  {w['locks']:4d}  {w['max_ratio_dev']:.3e}  "
              f"{w['freq_spread_rel']:.3e}  {w['pr_min']:8.3f}  "
              f"{w['transparent_pair_frac']:.3f}      {w['transparency_mean_sin']:.4f}{zone}")


def main() -> None:
    t0 = time.time()
    print("E-M8 交差窓スペクトル診断 実行")
    print(f"  import: ABL sha256={sha256(ABL)[:16]}…")
    print(f"  import: ENG sha256={sha256(Path(eng.__file__))[:16]}…")
    results = {"imports": {"abl": sha256(ABL), "engine": sha256(Path(eng.__file__))},
               "params": {"XMAX": XMAX, "WIN": WIN, "STRIDE": STRIDE,
                           "TOL_LOCK": TOL_LOCK, "FREQ_MIN": FREQ_MIN,
                           "TRANSPARENT_EPS": TRANSPARENT_EPS,
                           "PR_CONDENSE_RATIO": PR_CONDENSE_RATIO,
                           "PRE_ZONE_MARGIN": PRE_ZONE_MARGIN}}

    for n in (5, 40):
        windows, verdict, m = run_case(n)
        print_timeline(f"N={n} (M={m})", windows, verdict)
        if "H_A_supported" in verdict:
            tag = ("H-A支持（過渡生成/凝縮あり）" if verdict["H_A_supported"]
                   else "H-B支持（質的変化なし）" if verdict["H_B_supported"] else "混在")
            print(f"  判定: {tag}  遷移帯ロック最大={verdict['locks_transition_max']} "
                  f"後期ロック={verdict['locks_late_max']} 凝縮={verdict['condensation']} "
                  f"(PR遷移min={verdict['pr_min_transition']:.3f} / 前帯中央値={verdict['pr_pre_median']:.3f})")
        results[f"N{n}"] = {"windows": windows, "verdict": verdict, "m": m}

    results["runtime_sec"] = time.time() - t0
    (HERE / "paper8_crossing_window_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print(f"\nsaved ({results['runtime_sec']:.0f}s)")


if __name__ == "__main__":
    main()
