#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照テスト: 番地走査（未使用実験）の再現一致検定

目的: 波の周期表 v1/v2 の図P1 に使われた走査（実験一覧 #1・
run_periodic_address_scan_v1.py）を、本論文フォルダ内で無改変コピーから
再実行し、公開ベースライン JSON とビット一致するかを判定する。

規約（事前固定）:
  - 原本スクリプトは無改変（md5 照合済み）。依存チェーンも無改変コピー。
  - 比較対象は runtime_sec 以外の全数値フィールド（再帰的に全走査）。
  - 判定 C1: 全数値フィールドが完全一致（浮動小数の表現まで同一）。
  - 判定 C2: 一致しない場合、最大絶対差・最大相対差を報告し、
             1e-12 以下なら「機械精度一致」、それを超えれば不一致とする。
  - 判定 C3: 整数フィールド（N, M, n_sig）は完全一致必須。

使い方: python3 run_control_compare_address_scan_v1.py
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCAN = HERE / "次元の生成構造" / "番地走査_v1"
BASE = SCAN / "ベースライン_periodic_address_scan_result_v1.json"
CTRL = SCAN / "periodic_address_scan_result_v1.json"
OUT = HERE / "control_compare_address_scan_result_v1.json"

INT_KEYS = {"N", "M", "n_sig", "T_END", "SAMPLE_EVERY"}
SKIP_KEYS = {"runtime_sec"}


def walk(a, b, path="", acc=None):
    """2つのJSONを再帰比較し、(path, a, b) の差分と一致件数を集める。"""
    if acc is None:
        acc = {"n_num": 0, "n_exact": 0, "diffs": [], "type_mismatch": [],
               "int_mismatch": []}
    if isinstance(a, dict):
        if set(a) != set(b):
            acc["type_mismatch"].append((path, "キー集合が異なる"))
            return acc
        for k in a:
            if k in SKIP_KEYS:
                continue
            walk(a[k], b[k], f"{path}/{k}", acc)
    elif isinstance(a, list):
        if len(a) != len(b):
            acc["type_mismatch"].append((path, f"長さ {len(a)} vs {len(b)}"))
            return acc
        for i, (x, y) in enumerate(zip(a, b)):
            walk(x, y, f"{path}[{i}]", acc)
    elif isinstance(a, bool) or a is None or isinstance(a, str):
        if a != b:
            acc["type_mismatch"].append((path, f"{a!r} vs {b!r}"))
    elif isinstance(a, (int, float)):
        acc["n_num"] += 1
        key = path.rsplit("/", 1)[-1].split("[")[0]
        if a == b:
            acc["n_exact"] += 1
        else:
            d = abs(float(a) - float(b))
            rel = d / max(abs(float(a)), abs(float(b)), 1e-300)
            acc["diffs"].append({"path": path, "base": a, "ctrl": b,
                                 "abs": d, "rel": rel})
            if key in INT_KEYS:
                acc["int_mismatch"].append((path, a, b))
    else:
        acc["type_mismatch"].append((path, f"未対応型 {type(a).__name__}"))
    return acc


def main() -> None:
    base = json.loads(BASE.read_text())
    ctrl = json.loads(CTRL.read_text())
    acc = walk(base, ctrl)

    max_abs = max((d["abs"] for d in acc["diffs"]), default=0.0)
    max_rel = max((d["rel"] for d in acc["diffs"]), default=0.0)
    c1 = (not acc["diffs"]) and (not acc["type_mismatch"])
    c2 = max_abs <= 1e-12
    c3 = not acc["int_mismatch"]

    # 主要観測量の並べ比較（表示用）
    table = []
    for rb, rc in zip(base["scan"], ctrl["scan"]):
        table.append({
            "N": rb["N"],
            "clock_over_pi72": [rb["clock_over_pi72"], rc["clock_over_pi72"]],
            "clock_per_step": [rb["clock_over_pi72"] / base["SAMPLE_EVERY"],
                               rc["clock_over_pi72"] / ctrl["SAMPLE_EVERY"]],
            "n_sig": [rb["n_sig"], rc["n_sig"]],
            "perp_ratio": [rb["partitions"]["perp_ratio"],
                           rc["partitions"]["perp_ratio"]],
            "mass_deg": [rb["mass_deg"], rc["mass_deg"]],
        })

    out = {
        "対照テスト": "番地走査 v1（波の周期表 実験一覧 #1・図P1）",
        "ベースライン": str(BASE.relative_to(HERE)),
        "対照実行": str(CTRL.relative_to(HERE)),
        "比較した数値フィールド数": acc["n_num"],
        "完全一致数": acc["n_exact"],
        "差分件数": len(acc["diffs"]),
        "最大絶対差": max_abs,
        "最大相対差": max_rel,
        "型・キー不一致": acc["type_mismatch"],
        "整数フィールド不一致": acc["int_mismatch"],
        "判定C1_完全一致": c1,
        "判定C2_機械精度一致(1e-12)": c2,
        "判定C3_整数完全一致": c3,
        "差分詳細": acc["diffs"][:50],
        "主要観測量_並べ比較": table,
    }
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False))

    print(f"比較数値フィールド: {acc['n_num']}  完全一致: {acc['n_exact']}  "
          f"差分: {len(acc['diffs'])}")
    print(f"最大絶対差: {max_abs:.3e}  最大相対差: {max_rel:.3e}")
    print(f"C1 完全一致: {c1}   C2 機械精度一致: {c2}   C3 整数一致: {c3}")
    if acc["type_mismatch"]:
        print("型・キー不一致:", acc["type_mismatch"][:10])
    for d in acc["diffs"][:10]:
        print(f"  差分 {d['path']}: base={d['base']} ctrl={d['ctrl']} "
              f"abs={d['abs']:.3e}")
    print(f"→ {OUT.name}")


if __name__ == "__main__":
    main()
