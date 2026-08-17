#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CR5 の図を作る（テンプレート cr_waveform_shell_v1.html にデータを差し込む）

前提: make_cr5_waveform_series_v1.py が cr5_waveform_<tag>_data_v1.json を出力済み。
出力: cr5_waveform_<tag>_v1.html

使い方: python3 build_cr5_pages_v1.py [tag ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHELL = HERE / "cr_waveform_shell_v1.html"

TAGS = ["case17_3", "case3_3", "case17_17", "even", "even96"]

EYEBROW = "CR5 · 合成波形"

SUB = ("A と B の合成波形を、<b>それぞれの中心位相を 0° に揃えて</b>重ねたもの。"
       "中心は円周第1モーメントの偏角（Δθ を測るのと同じ計器）で決め、"
       "中心での搬送波位相は割り戻してある（割り戻さないと毎步回って形が読めない）。"
       "太線が合成波形、細線が包絡。<b>正規化は A・B・全時刻で共通の一つの係数</b>なので、"
       "高さの差はそのままパワーの差を表す。")

NOTE = ("<b>円の大きさは廃止した。</b>CR4 までは円の半径を A の参加率 PR から作っていたが、"
        "(1) その一つの値を A の円と B の円の両方に使っており B 側は何も表していなかった、"
        "(2) 正規化がその走行の実測範囲なのでケース間で比較できなかった、"
        "(3) PR はパワー由来で、像が 180° 間隔に割れる構成では二つの像を合わせて数えてしまう——"
        "の三点で量として成立していなかった。本版の ○ は固定サイズで、大きさは何も表さない。"
        "<br><br>"
        "広がりをスカラーに潰すのをやめ、波形そのものを出している。"
        "潰さなければ「幅をどう定義するか」の選択が要らず、像が割れていればそれが見える。"
        "<br><br>"
        "1 步刻み（間引きなし）。Δθ の支配周期は π/θ 步で、"
        "θ は AB 合成 χ スペクトルの |k|≥4 かつ偶数のパワー占有率 r = sin²θ から決まる。"
        "r は倍音を数えた有理数しか取らないので、<b>τ の刻みは離散</b>であり細分できない"
        "（probe_harmonic_composition_v1.py）。")

TITLE = {
    "case17_3":  "荷電二体波の合成波形 — 非対称（A:1〜17 / B:1〜3）",
    "case3_3":   "荷電二体波の合成波形 — 対称・少倍音（A・B とも 1〜3）",
    "case17_17": "荷電二体波の合成波形 — 対称・多倍音（A・B とも 1〜17）",
    "even":      "荷電二体波の合成波形 — 偶数優位（1,2,3,4,6,8,10,12,14,16）",
    "even96":    "荷電二体波の合成波形 — 偶数優位を帯域96まで拡大（50本）",
}

EXTRA = {
    "even": ("<br><br><b>この構成では像が二つに割れる。</b>"
             "偶数倍音が奇数倍音を圧倒するため χ の周期が 180° になり、"
             "149.77° と 329.77° に像が立つ（高さ比 0.600、|z| = 0.300）。"
             "偶数倍音のみなら高さ比 1.000 で |z| = 0.000000 になる。"
             "その代わり r = 1/20 まで下がり、1 周期を 13.89 步で読める"
             "（元の実験条件は 5.38 步）。"),
    "even96": ("<br><br><b>偶数優位を帯域 96 まで伸ばした構成（50本）。</b>"
               "r = 1/100 まで下がり、1 周期を 31.36 步（実測 30.77 步）で読める。"
               "元の実験条件 5.38 步の 5.7 倍細かく、1 步あたりの Δθ は 2.40°。"
               "<br><br>"
               "像は 149.77° と 329.77° の 2 つのままだが、高さ比が 0.600 → "
               "<b>0.918</b> に上がり、二像がほぼ等価になった。その結果 |z| は "
               "0.300 → <b>0.060</b> まで落ちている。位置が曖昧になったのではなく、"
               "<b>どちらの像か区別がつかなくなった</b>という意味である。"
               "<br><br>"
               "そのため<b>中心位相の決定が悪条件</b>になっている点に注意が要る。"
               "中心は円周第1モーメントで決めているが、走行中 |z| は 0.0524〜0.0600 "
               "しかない。二像がほぼ等しいとき第1モーメントは差の部分しか拾わないので、"
               "波形の中心合わせはこのケースでは弱い根拠しか持たない。"),
}


def build(tag: str) -> Path:
    shell = SHELL.read_text(encoding="utf-8")
    dp = HERE / f"cr5_waveform_{tag}_data_v1.json"
    data = dp.read_text(encoding="utf-8")
    cfg = json.loads(data)["config"]
    sub = (SUB + f"<br><br>構成: <b>{cfg['label']}</b>（{cfg['note']}）／"
           f"r = {cfg['r0']:.6f}、θ = {cfg['theta0']*180/3.141592653589793:.3f}°、"
           f"予測周期 π/θ = {cfg['period_pred']:.4f} 步／T = {cfg['T']} 步・1步刻み")
    html = (shell.replace("__TITLE__", TITLE[tag])
                 .replace("__EYEBROW__", EYEBROW)
                 .replace("__SUB__", sub)
                 .replace("__NOTE__", NOTE + EXTRA.get(tag, ""))
                 .replace("__DATA__", data))
    p = HERE / f"cr5_waveform_{tag}_v1.html"
    p.write_text(html, encoding="utf-8")
    return p


def main() -> None:
    for tag in (sys.argv[1:] or TAGS):
        p = build(tag)
        print(f"{p.name}  {p.stat().st_size//1024}KB")


if __name__ == "__main__":
    main()
