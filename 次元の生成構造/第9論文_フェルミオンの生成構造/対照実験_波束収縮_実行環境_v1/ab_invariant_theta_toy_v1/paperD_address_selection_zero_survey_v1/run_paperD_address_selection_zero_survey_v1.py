#!/usr/bin/env python3
"""論文D予備実験 D9-零次：住所選択の算術的候補汎関数の全数調査 v1

位置づけ:
    【探索的調査——予言を置かない】選択問題の第五候補（質量平衡選択）の
    設計に先立ち、レジスタ124の全既約住所について、単純な算術的選別汎関数が
    物理住所（流れ規約 39/124、電荷 cos²(39π/124)=0.302822）を選ぶかを
    全数調査する。単純な選別が全て 39 を外すなら、選択は算術のみでは
    決まらず、動力学（ノルムセクター）または親からの時計継承を要する——
    という消去の一歩として記録する。

調査する汎関数（事前登録・全数計算）:
    F1: 数え上げメニューのタング {p/q : q∈{1,2,3,4,6}} への最遠（maximin距離）
    F2: Farey重み付き距離 min_q |m/124 - p/q|·q²（大タングからの実効距離）
    F3: 連分数の最大部分商（小さいほど noble = 捕獲されにくい）
    F4: π/4（数え上げ均衡点）への近さ
"""

from __future__ import annotations

import csv
import json
from fractions import Fraction
from math import gcd, pi
from pathlib import Path

HERE = Path(__file__).resolve().parent
N = 124
PHYSICAL_M = 39  # 流れ規約の物理住所（電荷 cos²(39π/124) = 0.302822）


def continued_fraction(fr: Fraction):
    out = []
    p, q = fr.numerator, fr.denominator
    while q:
        a, r = divmod(p, q)
        out.append(a)
        p, q = q, r
    return out


def main() -> None:
    menu = sorted({Fraction(p, q) for q in (1, 2, 3, 4, 6) for p in range(0, q + 1)})
    addresses = [m for m in range(1, N) if m % 2 == 1 and gcd(m, N) == 1 and m < N // 2 + 10]
    # 流れ規約の物理領域（0<m<62 が独立、対称の重複を避けつつ39近傍を含む）

    rows = []
    for m in addresses:
        x = Fraction(m, N)
        d1 = min(abs(float(x - t)) for t in menu)
        d2 = min(abs(float(x - t)) * t.denominator ** 2 for t in menu if t.denominator > 0)
        cf = continued_fraction(x)
        f3 = max(cf[1:]) if len(cf) > 1 else 0
        f4 = abs(float(x) - 0.25)  # π/4 ⇔ 比 1/4（角/π 単位: m/124 vs 31/124=0.25）
        rows.append({"m": m, "ratio": float(x), "F1_maximin_menu": d1,
                     "F2_farey_weighted": d2, "F3_max_partial_quotient": f3,
                     "F4_dist_balance_point": f4, "cf": str(cf)})

    winners = {
        "F1（メニュー最遠）": max(rows, key=lambda r: r["F1_maximin_menu"])["m"],
        "F2（Farey重み最遠）": max(rows, key=lambda r: r["F2_farey_weighted"])["m"],
        "F3（最小の最大部分商=最noble）": min(rows, key=lambda r: r["F3_max_partial_quotient"])["m"],
        "F4（均衡点に最近接）": min(rows, key=lambda r: r["F4_dist_balance_point"])["m"],
    }
    phys = next(r for r in rows if r["m"] == PHYSICAL_M)

    print("=== D9-零次: 汎関数の勝者 vs 物理住所 39 ===")
    hit = []
    for name, w in winners.items():
        mark = "★物理住所と一致" if w == PHYSICAL_M else ""
        print(f"{name}: m = {w} {mark}")
        hit.append(w == PHYSICAL_M)
    print(f"\n物理住所 39/124 のプロファイル: {phys}")
    ranked1 = sorted(rows, key=lambda r: -r["F1_maximin_menu"])
    rank39_f1 = next(i for i, r in enumerate(ranked1, 1) if r["m"] == PHYSICAL_M)
    print(f"F1 での 39 の順位: {rank39_f1}/{len(rows)}")

    with (HERE / "paperD_address_survey_rows_v1.csv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    payload = {
        "experiment": "paperD_address_selection_zero_survey_v1",
        "kind": "exploratory_survey_no_predictions",
        "register": N, "physical_address_flow": PHYSICAL_M,
        "winners": winners,
        "any_functional_selects_physical": bool(any(hit)),
        "physical_profile": phys,
        "physical_rank_in_F1": rank39_f1,
        "conclusion": (
            "調査した単純算術汎関数（メニュー最遠・Farey重み・noble度・均衡点近接）の"
            f"勝者はそれぞれ {list(winners.values())} であり、物理住所39を選ぶものは"
            f"{'あった' if any(hit) else 'なかった'}。単純な算術選別では素電荷の住所は"
            "決まらない——選択はノルムセクターの動力学（質量平衡・D9本実験）または"
            "親閉鎖からの時計継承を要する、という消去の一歩として記録する"),
    }
    (HERE / "paperD_address_survey_result_v1.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    print("\nsaved: paperD_address_survey_result_v1.json")


if __name__ == "__main__":
    main()
