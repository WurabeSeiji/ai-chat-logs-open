#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全走行図のコンタクトシート生成 v1（監査用・論文には載せない）

目的: 本フォルダに 480 枚ある走行ごとの診断図を、人が全部目視できる形に焼く。
families:
  4panel       … f₂ / 次元 / 時間 / 閉塞残差
  mix          … 帯パワーと混合率
  ledger       … 128 セル帳簿と target セル
  summary      … N 掃引の要約
  birth_matrix … 誕生マトリクス
出力: sheet_<family>_<i>.png（監査用。git には入れるが論文には載せない）

使い方: python3 make_contact_sheets_v1.py [families...]
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

plt.rcParams["font.family"] = ["Hiragino Sans", "sans-serif"]

HERE = Path(__file__).resolve().parent
FAMILIES = ["4panel", "mix", "ledger", "summary", "birth_matrix"]
PER_SHEET = 12          # 3 列 × 4 行
NCOL = 3


def short(name: str) -> str:
    s = name.replace("fig_", "").replace("_v2.png", "")
    s = s.replace("fermion_family", "F5").replace("boson_family", "B3")
    s = s.replace("_rep-", " / ")
    return s


def build(family: str) -> int:
    files = sorted(p for p in HERE.glob(f"fig_*_{family}*_v2.png"))
    if not files:
        print(f"  {family}: 0 枚")
        return 0
    n_sheet = 0
    for s in range(0, len(files), PER_SHEET):
        chunk = files[s:s + PER_SHEET]
        nrow = (len(chunk) + NCOL - 1) // NCOL
        fig, axes = plt.subplots(nrow, NCOL, figsize=(NCOL * 6.4, nrow * 4.2))
        axes = axes.ravel() if hasattr(axes, "ravel") else [axes]
        for a in axes:
            a.axis("off")
        for a, f in zip(axes, chunk):
            a.imshow(mpimg.imread(f))
            a.set_title(short(f.name), fontsize=7)
        n_sheet += 1
        out = HERE / f"sheet_{family}_{n_sheet:02d}.png"
        fig.suptitle(f"監査用コンタクトシート — {family} "
                     f"（{s+1}〜{s+len(chunk)} / 全 {len(files)} 枚）", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(out, dpi=95)
        plt.close(fig)
        print(f"  → {out.name}")
    print(f"  {family}: {len(files)} 枚 → {n_sheet} シート")
    return n_sheet


if __name__ == "__main__":
    fams = sys.argv[1:] or FAMILIES
    total = 0
    for fam in fams:
        total += build(fam)
    print(f"合計 {total} シート")
