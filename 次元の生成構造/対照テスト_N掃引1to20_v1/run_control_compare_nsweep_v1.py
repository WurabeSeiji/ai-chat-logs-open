#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照テスト: N掃引 1→20 の再現一致検定（数値＋図）

判定（事前固定）:
  C1 runtime 系を除く全数値フィールドが完全一致（NaN は NaN と一致とみなす）
  C2 不一致があれば最大絶対差・最大相対差を報告（≤1e-12 なら機械精度一致）
  C3 整数フィールド（N など）完全一致
  F1 図 PNG のバイト一致（md5）
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE / "ベースライン_result_tb_nsweep_1to20_v1.json"
CTRL = HERE / "result_tb_nsweep_1to20_v1.json"
OUT = HERE / "control_compare_nsweep_result_v1.json"

SKIP = {"runtime_sec", "runtime", "elapsed_sec", "timestamp", "date"}
INT_KEYS = {"N", "Nn", "Neta", "T", "seed"}


def walk(a, b, path="", acc=None):
    if acc is None:
        acc = {"n_num": 0, "n_exact": 0, "n_nan": 0, "diffs": [], "struct": []}
    if isinstance(a, dict):
        if set(a) != set(b):
            acc["struct"].append((path, "キー集合が異なる",
                                  sorted(set(a) ^ set(b))[:8]))
            return acc
        for k in a:
            if k in SKIP:
                continue
            walk(a[k], b[k], f"{path}/{k}", acc)
    elif isinstance(a, list):
        if len(a) != len(b):
            acc["struct"].append((path, f"長さ {len(a)} vs {len(b)}"))
            return acc
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f"{path}[{i}]", acc)
    elif isinstance(a, bool) or a is None or isinstance(a, str):
        if a != b:
            acc["struct"].append((path, f"{a!r} vs {b!r}"))
    elif isinstance(a, (int, float)):
        acc["n_num"] += 1
        fa, fb = float(a), float(b)
        if math.isnan(fa) and math.isnan(fb):
            acc["n_exact"] += 1
            acc["n_nan"] += 1
        elif fa == fb:
            acc["n_exact"] += 1
        else:
            d = abs(fa - fb)
            acc["diffs"].append({
                "path": path, "base": a, "ctrl": b, "abs": d,
                "rel": d / max(abs(fa), abs(fb), 1e-300),
                "int_key": path.rsplit("/", 1)[-1].split("[")[0] in INT_KEYS})
    else:
        acc["struct"].append((path, f"未対応型 {type(a).__name__}"))
    return acc


def md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest()


def main() -> None:
    acc = walk(json.loads(BASE.read_text()), json.loads(CTRL.read_text()))
    max_abs = max((d["abs"] for d in acc["diffs"]), default=0.0)
    max_rel = max((d["rel"] for d in acc["diffs"]), default=0.0)
    c1 = not acc["diffs"] and not acc["struct"]
    c3 = not any(d["int_key"] for d in acc["diffs"])

    figs = []
    for c in sorted(HERE.glob("fig_nsweep_*.png")):
        b = HERE / f"ベースライン_{c.name}"
        if not b.exists():
            figs.append({"図": c.name, "ベースライン": "なし"})
            continue
        mc, mb = md5(c), md5(b)
        figs.append({"図": c.name, "F1_バイト一致": mc == mb,
                     "md5_対照": mc, "md5_基準": mb,
                     "バイト数": [c.stat().st_size, b.stat().st_size]})
    f1 = all(f.get("F1_バイト一致") for f in figs if "F1_バイト一致" in f)

    out = {"対照テスト": "N掃引 1→20（run_tb_nsweep_1to20_v1.py 無改変）",
           "比較数値フィールド数": acc["n_num"], "完全一致数": acc["n_exact"],
           "うちNaN一致": acc["n_nan"], "差分件数": len(acc["diffs"]),
           "最大絶対差": max_abs, "最大相対差": max_rel,
           "構造差": acc["struct"], "判定C1_完全一致": c1,
           "判定C2_機械精度一致": max_abs <= 1e-12, "判定C3_整数一致": c3,
           "判定F1_全図バイト一致": f1, "差分詳細": acc["diffs"][:40],
           "図の比較": figs}
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False))

    print(f"数値フィールド {acc['n_num']}  完全一致 {acc['n_exact']}"
          f"（NaN一致 {acc['n_nan']}）  差分 {len(acc['diffs'])}")
    print(f"最大絶対差 {max_abs:.3e}  最大相対差 {max_rel:.3e}")
    print(f"C1 完全一致 {c1}  C2 機械精度一致 {max_abs <= 1e-12}  C3 整数一致 {c3}")
    for s in acc["struct"][:10]:
        print("  構造差:", s)
    for d in acc["diffs"][:10]:
        print(f"  差分 {d['path']}: {d['base']} vs {d['ctrl']} abs={d['abs']:.3e}")
    ng = [f["図"] for f in figs if f.get("F1_バイト一致") is False]
    print(f"図 {len(figs)} 枚  全バイト一致 {f1}" + (f"  不一致: {ng}" if ng else ""))
    print(f"→ {OUT.name}")


if __name__ == "__main__":
    main()
