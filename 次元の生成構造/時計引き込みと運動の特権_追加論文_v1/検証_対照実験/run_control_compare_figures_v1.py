#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""対照テスト（図）: 図P1〜P6 の再現一致検定

目的: 数値 JSON の一致（run_control_compare_address_scan_v1.py）に加え、
公開図そのものが対照環境で再現するかを判定する。図P1 の入力は
**対照実行で生成した** periodic_address_scan_result_v1.json であり、
残り5本の JSON はベースラインの read-only import。

規約（事前固定）:
  - 図生成スクリプト make_paper_figures_periodic_v1.py は無改変（md5 照合済み）。
  - 判定 F1: PNG バイト列が公開版と完全一致（md5 一致）。
  - 判定 F2: F1 が不成立の場合、画素配列の最大絶対差と不一致画素率を報告し、
             最大絶対差 0 なら「画素完全一致（メタデータのみ差）」とする。
  - 判定 F3: 画像サイズ（幅・高さ・チャネル数）一致。

使い方: python3 run_control_compare_figures_v1.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
CTRL_DIR = HERE / "次元の生成構造" / "波の周期表検討"
PUB_DIR = HERE.parent.parent / "波の周期表検討"
OUT = HERE / "control_compare_figures_result_v1.json"

FIGS = ["fig_p1_clock_universality_v1.png",
        "fig_p2_stable_vs_resonance_v1.png",
        "fig_p3_charged_lifetime_walk_v1.png",
        "fig_p4_rectification_v1.png",
        "fig_p5_ledger_v1.png",
        "fig_p6_cyclic_conservation_v1.png"]

# 図P1 のみ対照実行の走査結果を入力に使う（他はベースライン import）
CTRL_INPUT = {"fig_p1_clock_universality_v1.png": "対照実行の走査JSON",
              "fig_p2_stable_vs_resonance_v1.png": "ベースラインJSON",
              "fig_p3_charged_lifetime_walk_v1.png": "ベースラインJSON",
              "fig_p4_rectification_v1.png": "ベースラインJSON",
              "fig_p5_ledger_v1.png": "ベースラインJSON",
              "fig_p6_cyclic_conservation_v1.png": "ベースラインJSON"}


def md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest()


def main() -> None:
    rows = []
    for name in FIGS:
        c, p = CTRL_DIR / name, PUB_DIR / name
        if not c.exists() or not p.exists():
            rows.append({"図": name, "状態": "欠落",
                         "対照": c.exists(), "公開": p.exists()})
            continue
        mc, mp = md5(c), md5(p)
        ac = np.asarray(Image.open(c).convert("RGBA")).astype(np.int16)
        ap = np.asarray(Image.open(p).convert("RGBA")).astype(np.int16)
        same_shape = ac.shape == ap.shape
        if same_shape:
            d = np.abs(ac - ap)
            max_abs = int(d.max())
            bad = int((d.max(axis=2) > 0).sum())
            total = int(d.shape[0] * d.shape[1])
        else:
            max_abs, bad, total = -1, -1, -1
        rows.append({
            "図": name,
            "入力": CTRL_INPUT[name],
            "md5_対照": mc, "md5_公開": mp,
            "F1_バイト一致": mc == mp,
            "F3_サイズ一致": same_shape,
            "形状": list(ac.shape),
            "画素_最大絶対差": max_abs,
            "不一致画素数": bad, "総画素数": total,
            "F2_画素完全一致": max_abs == 0,
            "バイト数": [c.stat().st_size, p.stat().st_size],
        })

    f1 = all(r.get("F1_バイト一致") for r in rows)
    f2 = all(r.get("F2_画素完全一致") for r in rows)
    f3 = all(r.get("F3_サイズ一致") for r in rows)
    out = {
        "対照テスト": "図P1〜P6 の再現一致検定",
        "図生成スクリプト": "make_paper_figures_periodic_v1.py（無改変・md5 0e62ba5c0baae7f7de6ce327d85199bf）",
        "判定F1_全図バイト一致": f1,
        "判定F2_全図画素完全一致": f2,
        "判定F3_全図サイズ一致": f3,
        "図ごとの結果": rows,
    }
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False))

    print(f"{'図':40} {'バイト一致':>8} {'画素最大差':>8} {'不一致画素':>10}")
    for r in rows:
        print(f"{r['図']:40} {str(r.get('F1_バイト一致')):>8} "
              f"{r.get('画素_最大絶対差'):>8} "
              f"{r.get('不一致画素数')}/{r.get('総画素数')}")
    print(f"F1 全図バイト一致: {f1}   F2 全図画素完全一致: {f2}   "
          f"F3 サイズ一致: {f3}")
    print(f"→ {OUT.name}")


if __name__ == "__main__":
    main()
