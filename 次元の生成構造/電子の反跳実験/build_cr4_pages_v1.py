#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR4: 相対位相の図を生成する（テンプレート＋データ → 自己完結HTML）

テンプレート: cr_relative_shell_v1.html（__TITLE__/__EYEBROW__/__SUB__/__NOTE__/__DATA__）
データ:       cr4_relative_<tag>_data_v1.json（make_cr4_relative_series_v1.py が生成）
出力:         cr4_relative_<tag>_v1.html

使い方:
  python3 build_cr4_pages_v1.py            # 全ケース
  python3 build_cr4_pages_v1.py case3_3    # 個別
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHELL = HERE / "cr_relative_shell_v1.html"

TITLES = {
    "case17_3":  "非対称二体の相対位相",
    "case3_3":   "三倍音どうしの相対位相",
    "case17_17": "十七倍音どうしの相対位相",
}


def build(tag: str) -> Path:
    d = json.loads((HERE / f"cr4_relative_{tag}_data_v1.json").read_text(encoding="utf-8"))
    cfg = d["config"]
    sp0 = d["spec"][0]
    M_a, M_b = len(cfg["packet_a"]), len(cfg["packet_b"])
    Rp = d["Rp"]
    rp_const = (min(Rp) == max(Rp))
    r_mean = sum(d["r"]) / len(d["r"])

    sub = (f"  τ 方向の共通回転を取り除き、相対位相 Δθ だけを 1步刻みで見る。"
           f"{cfg['label']}（{cfg['note']}）。<br>\n"
           f"  Δθ は ±60° に閉じたまま、支配周期 {sp0['period']:.4f} 步で振動する。"
           f"復元力が作るはずの 144 步は現れない（帯内占有 {d['share_144band']:.4f}%）。")

    note = f"""  <b>何を見せているか。</b> 位相円を巡る運動から共通回転を差し引くと、
  残るのは相対位相 Δθ だけになる。振り子の振れ角が Δθ で、左右の玉が A と B の
  相対位置、円の大きさが波束の広がりである。構成は {cfg['label']}、
  初期位置 A={cfg['deg_a']}° / B={cfg['deg_b']}°、κ は宣言せず毎步の透過率 1−r を読む
  （実測平均 r = {r_mean:.6f}）。
  <br><br>
  <b>位置の確度。</b> 単一倍音では |ψ|² が周波数 2k しか持たず第1円周モーメントが
  恒等的にゼロになり、位置が定義できない。M 本重ねると確度は
  <code>|z| = 1 − 1/M</code> で立ち上がる。本ケースは A が M={M_a} で
  |z| = {d['z0'][0]:.6f}、B が M={M_b} で |z| = {d['z0'][1]:.6f}。
  波束の広がりは初期 {d['pr'][0]:.2f} セル（円周 512 のうち）。
  <br><br>
  <b>曲率半径 R′²。</b> {"本ケースでは <b>定数</b>で、動径は固定表示になる。R′ が構造を持つには倍音の非対称が要る。" if rp_const else f"本ケースでは <b>{100*(max(Rp)-min(Rp))/ (sum(Rp)/len(Rp)):.1f}% の幅で変動</b>する。"}
  <br><br>
  <b>単振動ではない。</b> 復元力 <code>a = −4sin²(ω/2)·Δθ ≈ −ω²Δθ</code> が作る固有振動の
  周期は <code>2π/ω = 144 步</code>（普遍時計の一周）だが、<b>144±10% の帯の占有は
  {d['share_144band']:.4f}%</b>、周期 100 步以上の総和も {d['share_period_ge100']:.2f}% にとどまる。
  支配は {sp0['period']:.4f} 步（占有 {sp0['share']:.2f}%）。1步の加速度が速度そのものと
  同程度で、144 步かけて振り戻す前に、並進と衝突の写像が持つ速い過程が支配している。
  位相空間の軌道も復元力の単振動ではなく、この速い振動のものである。
  <br><br>
  <b>+1 ボタンと矢印キーで1步ずつ送れる。</b> τ の1步は動力学の最小単位で、
  これより細かい刻みは存在しない（離散が基礎で、連続はその極限）。
  20 步おきに間引くとナイキスト（40 步）を割ってエイリアスするため、
  本データは間引いていない。
"""

    html = (SHELL.read_text(encoding="utf-8")
            .replace("__TITLE__", TITLES[tag])
            .replace("__EYEBROW__", f"CR4 · 1步刻み · {cfg['label']}")
            .replace("__SUB__", sub)
            .replace("__NOTE__", note)
            .replace("__DATA__", json.dumps(d, ensure_ascii=False, separators=(",", ":"))))
    out = HERE / f"cr4_relative_{tag}_v1.html"
    out.write_text(html, encoding="utf-8")
    assert "__" not in html.split("<script>")[0].replace("__DATA__", ""), "未置換あり"
    return out


def main() -> None:
    tags = sys.argv[1:] or list(TITLES)
    for tag in tags:
        p = build(tag)
        print(f"生成: {p.name}  {p.stat().st_size//1024}KB  ({TITLES[tag]})")


if __name__ == "__main__":
    main()
