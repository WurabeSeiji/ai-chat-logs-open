#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""α根への最接近値を全走行から取る v1

契機: 全480枚監査（`分析記録_全図面監査_v2.md`）で、N=4 の走行の
「α根からの距離」パネルが 10⁻³ 付近まで落ち込むのを目視した。
既存報告 `結果報告_シード型別δ一括掃引_T42000_v1.md` は
「保存済み node 別 R の全40走行最大は 0.678602632934 で α根に未到達」と述べているが、
それは **N=12・T=42000 の 40 条件**についての値であり、
**N=4 の走行と長時間走行（T=300000）は含まれていない**。

本器は正本 NPZ を全数走査し、保存されている読出しすべてについて
α根 R_α = cos²(23π/124) への最接近を取る。新規走行はしない（read-only）。

走査する読出し（母体が保存しているもの）:
  rec_m_r_mean   node 別 R = scale·sin²θ の node 平均
  rec_m_r_med    同 中央値
  rec_m_r_min    同 最小（node 別の下端）
  rec_m_r_max    同 最大（node 別の上端）  ← 「どれか1つの node がどこまで行ったか」
  rec_m_r_raw    奇数8帯 /(奇数8帯+偶数7帯・ポンプ込み)
  rec_m_r_nopump 奇数8帯 /(奇数8帯+偶数6帯・ポンプ除く)

注意（本器の限界）: これらはいずれも**受動的な診断量**であり、
非線形頂点が実際に使う状態依存の反射率 R_e は保存されていない。
本器が測るのは「保存された読出しがα根にどこまで近づいたか」までである。

出力: result_alpha_root_closest_v1.json（＋標準出力に要約）

使い方: python3 probe_alpha_root_closest_approach_v1.py
"""
from __future__ import annotations
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
R_ALPHA = math.cos(23 * math.pi / 124) ** 2      # 0.6971779275566593
KEYS = ["rec_m_r_mean", "rec_m_r_med", "rec_m_r_min", "rec_m_r_max",
        "rec_m_r_raw", "rec_m_r_nopump"]


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def parse(name: str) -> dict:
    """ファイル名から走行の素性を読む（nsweep_<mode>_[T..]_[d..]_[rep-..]_N<k>_v2.npz）。"""
    s = name[len("nsweep_"):-len("_v2.npz")]
    parts = s.split("_")
    N = int(parts[-1][1:]) if parts[-1].startswith("N") else None
    T, delta, rep = 4000, None, ""
    mode_parts = []
    for p in parts[:-1]:
        if p.startswith("T") and p[1:].isdigit():
            T = int(p[1:])
        elif p.startswith("d") and any(c.isdigit() for c in p[1:]):
            try:
                delta = float(p[1:])
            except ValueError:
                mode_parts.append(p)
        elif p.startswith("rep-"):
            rep = p[4:]
        else:
            mode_parts.append(p)
    return {"mode": "_".join(mode_parts), "T": T, "delta": delta,
            "N": N, "rep": rep}


def main() -> None:
    t0 = time.time()
    files = sorted(HERE.glob("nsweep_*_v2.npz"))
    rows = []
    print(f"走査対象 NPZ: {len(files)} 件")
    for f in files:
        z = np.load(f)
        meta = parse(f.name)
        rec = {"file": f.name, **meta, "readouts": {}}
        for k in KEYS:
            if k not in z.files:
                continue
            v = np.asarray(z[k], float)
            v = v[np.isfinite(v)]
            if v.size == 0:
                continue
            d = np.abs(v - R_ALPHA)
            i = int(np.argmin(d))
            rec["readouts"][k] = {
                "min_abs_dist": float(d[i]),
                "value_at_min": float(v[i]),
                "step_at_min": i + 1,
                "max_value": float(v.max()),
                "crossed": bool(v.max() >= R_ALPHA),
            }
        if rec["readouts"]:
            rows.append(rec)
    out = {"generator": {"script": Path(__file__).name,
                         "sha256": sha256(Path(__file__).resolve())},
           "R_alpha": R_ALPHA,
           "definition": "R_alpha = cos^2(23*pi/124)",
           "n_files": len(files), "rows": rows}

    print(f"\nα根 R_α = {R_ALPHA:.15f}\n")
    for k in KEYS:
        cand = [(r["readouts"][k]["min_abs_dist"], r) for r in rows if k in r["readouts"]]
        if not cand:
            continue
        cand.sort(key=lambda x: x[0])
        best_d, best = cand[0]
        b = best["readouts"][k]
        n_cross = sum(1 for r in rows if k in r["readouts"] and r["readouts"][k]["crossed"])
        gmax = max(r["readouts"][k]["max_value"] for r in rows if k in r["readouts"])
        print(f"=== {k} ===")
        print(f"  最接近 |r−R_α| = {best_d:.6e}   （相対 {best_d/R_ALPHA*100:.4f}%）")
        print(f"    そのときの値 r = {b['value_at_min']:.12f}  第 {b['step_at_min']} 回")
        print(f"    走行: mode={best['mode']} N={best['N']} T={best['T']} "
              f"δ={best['delta']} rep={best['rep'] or '—'}")
        print(f"  全走行での最大値 = {gmax:.12f}")
        print(f"  α根を越えた走行 = {n_cross} / {len(cand)} 件")
        out.setdefault("summary", {})[k] = {
            "closest": {"abs_dist": best_d, "value": b["value_at_min"],
                        "step": b["step_at_min"], "file": best["file"],
                        "mode": best["mode"], "N": best["N"], "T": best["T"],
                        "delta": best["delta"], "rep": best["rep"]},
            "global_max": gmax, "n_crossed": n_cross, "n_runs": len(cand)}
        print()

    # N=12・T=42000 の 40 条件だけに絞った既存報告との照合
    sub = [r for r in rows if r["N"] == 12 and r["T"] == 42000 and not r["rep"]]
    if sub:
        gmax40 = max(r["readouts"]["rec_m_r_max"]["max_value"]
                     for r in sub if "rec_m_r_max" in r["readouts"])
        out["legacy_check_N12_T42000_no_rep"] = {
            "n_runs": len(sub), "max_r_max": gmax40,
            "reported_in_existing_md": 0.678602632934}
        print(f"[既存報告との照合] N=12・T=42000・rep なし {len(sub)} 件の "
              f"node 別 R 最大 = {gmax40:.12f}")
        print(f"                   既存 md の記載値 = 0.678602632934")

    out["runtime_sec"] = time.time() - t0
    (HERE / "result_alpha_root_closest_v1.json").write_text(
        json.dumps(out, indent=1, ensure_ascii=False))
    print(f"\n完了 {out['runtime_sec']:.0f}s → result_alpha_root_closest_v1.json")


if __name__ == "__main__":
    main()
