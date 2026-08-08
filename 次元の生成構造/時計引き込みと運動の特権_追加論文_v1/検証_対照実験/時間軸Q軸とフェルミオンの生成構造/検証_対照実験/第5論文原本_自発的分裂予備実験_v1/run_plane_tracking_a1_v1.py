#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 段階A1（摂動なし）：混合準安定状態の回転平面分解・時間追跡・占有量分類。

第二準安定窓で、瞬時生成子 K(arg Z_t) を密行列 eig で2次元回転平面群（σ値でグループ化、
縮退は部分空間）へ分解し、占有量 E_j(t)=|Π_{P_j(t)}Z_t|² を測る。平面は射影重なり
½Tr(Π_j(t)Π_k(t+1)) の最大対応で時間追跡する。摂動は加えない。

報告：平面数／平面別平均占有量スペクトル／支配・中間・低占有分類／追跡重なり／縮退／
平面自体の回転／低占有候補平面数／有効ランク。N=10,20,40。原本エンジンは不変更で import。

使い方: python3 run_plane_tracking_a1_v1.py 10 20 40
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_transverse_stability_v1 import reconstruct_metastable, GUARD, LEARN, VALID

RESULT_DIR = HERE / "plane_tracking_a1_result_v1"
SAMPLE_EVERY = 50          # 準安定窓の抽出間隔


def generator_plane_groups(sys_lr, Z, sig_rel=1e-6):
    """K(arg Z) を密 eig でσ群（縮退部分空間）へ分解。各群 (σ, dim, B) を占有降順で返す。"""
    sys_lr.set_theta(np.angle(Z))
    M = sys_lr.m
    K = np.column_stack([sys_lr.kmatvec(np.eye(M)[:, j]) for j in range(M)])
    w, V = np.linalg.eig(K)
    sig = w.imag
    smax = float(sig.max())
    thr = sig_rel * smax
    raw = {}
    for i in range(M):
        if sig[i] > thr:
            key = round(float(sig[i]), 6)
            raw.setdefault(key, []).extend([V[:, i].real, V[:, i].imag])
    groups = []
    for s, cols in raw.items():
        Q, R = np.linalg.qr(np.column_stack(cols))
        keep = np.abs(np.diag(R)) > 1e-8
        B = Q[:, keep]
        groups.append({"sigma": s, "dim": int(B.shape[1]), "B": B})
    return groups, smax, M


def occupation(B, Z):
    a, b = Z.real, Z.imag
    return float(np.sum((B.T @ a) ** 2) + np.sum((B.T @ b) ** 2))


def overlap(Ba, Bb):
    """½Tr(Π_a Π_b) を min(dim) で正規化。=1 で同一平面群。"""
    F = Ba.T @ Bb
    return float(np.sum(F ** 2) / min(Ba.shape[1], Bb.shape[1]))


def run_a1(n):
    rec = reconstruct_metastable(n)
    sys_lr = rec["sys"]
    p, q = rec["p"], rec["q"]
    t0 = rec["crossing"] + GUARD
    t1 = rec["crossing"] + GUARD + LEARN + VALID
    times = [t for t in range(t0, t1 + 1, SAMPLE_EVERY) if t in rec["Zs"]]
    # 各サンプルで平面群分解・占有量
    snaps = []
    for t in times:
        Z = rec["Zs"][t]
        groups, smax, M = generator_plane_groups(sys_lr, Z)
        for g in groups:
            g["occ"] = occupation(g["B"], Z)
        ker_occ = 1.0 - sum(g["occ"] for g in groups)
        snaps.append({"t": t, "groups": groups, "smax": smax, "ker": max(0.0, ker_occ),
                      "f": p is not None and float(1 - (abs(p @ Z) ** 2 + abs(q @ Z) ** 2) / np.real(np.conj(Z) @ Z))})
    M = sys_lr.m

    # 時間追跡：連続サンプル間で σ 近傍＋射影重なり最大の対応
    first = snaps[0]["groups"]
    tracks = [{"sigma0": g["sigma"], "dim": g["dim"], "occ": [g["occ"]],
               "B_first": g["B"], "B_last": g["B"], "min_overlap": 1.0} for g in first]
    prev = first
    for sn in snaps[1:]:
        cur = sn["groups"]
        used = set()
        for ti, tr in enumerate(tracks):
            # 直近 prev の対応群 = tr に紐付く群（σ近傍で prev から探す）
            # 単純化：tr.B_last と cur 各群の重なり最大を割当
            best_j, best_o = -1, -1.0
            for j, cg in enumerate(cur):
                if j in used:
                    continue
                if abs(cg["dim"] - tr["dim"]) > 0:
                    continue
                o = overlap(tr["B_last"], cg["B"])
                if o > best_o:
                    best_o, best_j = o, j
            if best_j >= 0:
                used.add(best_j)
                tr["occ"].append(cur[best_j]["occ"])
                tr["B_last"] = cur[best_j]["B"]
                tr["min_overlap"] = min(tr["min_overlap"], best_o)
            else:
                tr["occ"].append(np.nan)
        prev = cur

    # 平面別統計
    plane_stats = []
    for tr in tracks:
        occ = np.array(tr["occ"], float)
        occ = occ[~np.isnan(occ)]
        plane_stats.append({
            "sigma0": tr["sigma0"], "dim": tr["dim"],
            "mean_occ": float(np.mean(occ)), "min_occ": float(np.min(occ)),
            "max_occ": float(np.max(occ)),
            "track_min_overlap": float(tr["min_overlap"]),
            "self_rotation_start_end_overlap": overlap(tr["B_first"], tr["B_last"]),
            "n_tracked": int(len(occ)),
        })
    plane_stats.sort(key=lambda x: -x["mean_occ"])

    # 分類（平均占有量のギャップ）：支配/中間/低占有
    occs = np.array([ps["mean_occ"] for ps in plane_stats])
    ker_mean = float(np.mean([sn["ker"] for sn in snaps]))
    # 数値床：最小占有と ker の比較。低占有 = mean_occ < 1e-3、中間 = 1e-3〜0.05、支配 = >0.05
    for ps in plane_stats:
        m = ps["mean_occ"]
        ps["class"] = "dominant" if m > 0.05 else ("intermediate" if m > 1e-3 else "low")

    # 有効ランク：窓のフレーム共分散
    cols = []
    for t in times:
        Z = rec["Zs"][t]
        cols.append(Z.real); cols.append(Z.imag)
    X = np.column_stack(cols)
    s = np.linalg.svd(X, compute_uv=False)
    lam = s ** 2
    eff_rank = float((lam.sum() ** 2) / np.sum(lam ** 2))

    # 帯（band）分析：個別低占有平面は追跡不能なので、σで支配帯と低占有帯へ束ねる。
    # 各サンプルで σ>0.75σmax を支配帯、それ以外(σ>0)を低占有帯としてまとめた射影で占有・追跡。
    def band_bases(groups, smax):
        dom, low = [], []
        for g in groups:
            (dom if g["sigma"] > 0.75 * smax else low).append(g["B"])
        Bd = np.linalg.qr(np.column_stack(dom))[0] if dom else None
        Bl = np.linalg.qr(np.column_stack(low))[0] if low else None
        return Bd, Bl
    band = {"dom_occ": [], "low_occ": [], "low_dim": [], "low_sigma_lo": [], "low_sigma_hi": []}
    Bl_first = Bl_last = Bd_first = Bd_last = None
    for k, sn in enumerate(snaps):
        Z = rec["Zs"][sn["t"]]
        Bd, Bl = band_bases(sn["groups"], sn["smax"])
        band["dom_occ"].append(occupation(Bd, Z))
        band["low_occ"].append(occupation(Bl, Z))
        band["low_dim"].append(int(Bl.shape[1]))
        lows = [g["sigma"] for g in sn["groups"] if g["sigma"] <= 0.75 * sn["smax"]]
        band["low_sigma_lo"].append(min(lows)); band["low_sigma_hi"].append(max(lows))
        if k == 0:
            Bl_first, Bd_first = Bl, Bd
        Bl_last, Bd_last = Bl, Bd
    band_report = {
        "dom_band_mean_occ": float(np.mean(band["dom_occ"])),
        "low_band_mean_occ": float(np.mean(band["low_occ"])),
        "low_band_dim": int(np.round(np.mean(band["low_dim"]))),
        "low_band_sigma_range": [float(np.mean(band["low_sigma_lo"])), float(np.mean(band["low_sigma_hi"]))],
        "low_band_sigma_over_dom": float(np.mean(band["low_sigma_hi"])) / float(np.mean([sn["smax"] for sn in snaps])),
        "dom_band_self_rotation_start_end": overlap(Bd_first, Bd_last),
        "low_band_self_rotation_start_end": overlap(Bl_first, Bl_last),
    }

    n_plane_groups = len(plane_stats)
    n_low = sum(1 for ps in plane_stats if ps["class"] == "low")
    n_inter = sum(1 for ps in plane_stats if ps["class"] == "intermediate")
    n_dom = sum(1 for ps in plane_stats if ps["class"] == "dominant")
    degenerate = [ps for ps in plane_stats if ps["dim"] > 2]

    report = {
        "n": n, "m": M, "crossing": rec["crossing"],
        "window": [t0, t1], "n_samples": len(snaps),
        "n_plane_groups_nonzero_sigma": n_plane_groups,
        "kernel_mean_occupation": ker_mean,
        "effective_rank_frame": eff_rank,
        "classification_counts": {"dominant": n_dom, "intermediate": n_inter, "low": n_low},
        "n_degenerate_groups(dim>2)": len(degenerate),
        "plane_stats": plane_stats,
        "sigma_max_mean": float(np.mean([sn["smax"] for sn in snaps])),
        "f_window_mean": float(np.mean([sn["f"] for sn in snaps])),
        "band": band_report,
    }
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    with open(RESULT_DIR / f"a1_N{n:05d}.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ns", type=int, nargs="+")
    args = ap.parse_args()
    for n in args.ns:
        r = run_a1(n)
        print(f"\n=== N={n} (M={r['m']}) 準安定窓 f≈{r['f_window_mean']:.3f} σmax≈{r['sigma_max_mean']:.2f} ===")
        print(f"回転平面群数(σ>0)={r['n_plane_groups_nonzero_sigma']}  核平均占有={r['kernel_mean_occupation']:.4f}  有効ランク={r['effective_rank_frame']:.2f}")
        print(f"分類: 支配={r['classification_counts']['dominant']} 中間={r['classification_counts']['intermediate']} 低占有={r['classification_counts']['low']}  縮退群(dim>2)={r['n_degenerate_groups(dim>2)']}")
        print(f"{'順位':>3} {'σ':>8} {'dim':>3} {'平均占有':>10} {'追跡重なり':>9} {'自己回転':>8} {'分類':>10}")
        for i, ps in enumerate(r["plane_stats"][:12]):
            print(f"{i+1:>3} {ps['sigma0']:>8.3f} {ps['dim']:>3} {ps['mean_occ']:>10.5f} {ps['track_min_overlap']:>9.4f} {ps['self_rotation_start_end_overlap']:>8.4f} {ps['class']:>10}")
        b = r["band"]
        print(f"帯: 支配帯 σ≈{r['sigma_max_mean']:.2f} 占有={b['dom_band_mean_occ']:.5f} 自己回転={b['dom_band_self_rotation_start_end']:.3f} | "
              f"低占有帯 σ∈[{b['low_band_sigma_range'][0]:.2f},{b['low_band_sigma_range'][1]:.2f}] (σ/σdom={b['low_band_sigma_over_dom']:.3f}) "
              f"dim={b['low_band_dim']} 占有={b['low_band_mean_occ']:.2e} 帯自己回転={b['low_band_self_rotation_start_end']:.3f}")


if __name__ == "__main__":
    main()
