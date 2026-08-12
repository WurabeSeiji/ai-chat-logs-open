#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""選択層 S v1 — 曖昧さの束から確定値を作る操作（宣言されるべき人為）

  確定値 = S ∘ G

統一場の万能読出し関数 v2（unified_readout_v2.py）は、選択をせずに
重みつきの読み値の束を返す。実験が「確定した一つの値」を必要とするとき、
どの選択を適用したかを **仕様に宣言した上で** ここの選択子を使う。

規約:
  S1 選択子は必ず引数として閾値・基準を受け取り、既定値を仕様書に記録する。
     （G はパラメータフリー・S はパラメータを持つ——この分離が設計の要点）
  S2 選択子は G の束のみを引数に取り、状態 C2 を直接読まない。
  S3 実験の記録には「S の宣言 ＋ G の束」を対で残す。従来の確定値は
     どの S を適用したかを添えて初めて再現可能になる。

正当化: 木原の観測理論——粒子・位置・時間が確定して見えるのは、系が
確定しているからではなく、読出し側が選択しているから。選択を関数として
書き出すことで、「どの選択の下で何が確定して見えるか」が実験の対象になる。
"""
from __future__ import annotations
import numpy as np

# 既定の選択パラメータ（仕様書に記録する値・v1 の暗黙値を明示化したもの）
FLOOR_CONTENT = 1e-280      # 内容パワーの存在床（v1 g_position_1d の暗黙値）
FLOOR_OVERLAP = 1e-30       # 時計の重なり床（v1 g_matter_clock の暗黙値）
FLOOR_CARRIER = 1e-12       # 時計の担体ノルム床（v1 g_matter_clock の暗黙値）


def s_present(bundle, floor=FLOOR_CONTENT):
    """内容の存在判定（連続重み → 真偽）。v1 g_position_1d の present に一致。"""
    return bool(bundle["content_power"] > floor)


def s_position_maxmoment(bundle, floor=FLOOR_CONTENT):
    """位置の被覆選択: 巻き m=1,2 のうち重みの大きい方を採る。
    v1 g_position_1d の被覆判定（if |z2|>|z1|）と厳密に同一の選択。
    戻り値 {"x", "cover", "present"}。"""
    if not s_present(bundle, floor):
        return {"x": None, "cover": None, "present": False}
    w = bundle["pos_weight"]
    if len(w) >= 2 and w[1] > w[0]:
        return {"x": float(bundle["pos_x"][1]), "cover": 2, "present": True}
    return {"x": float(bundle["pos_x"][0]), "cover": 1, "present": True}


def s_position_argmax(bundle, floor=FLOOR_CONTENT):
    """位置の被覆選択（一般形）: 全ての巻き m のうち重み最大を採る。
    Nn が大きい環境で m≥3 の被覆が現れうるため、maxmoment の一般化。"""
    if not s_present(bundle, floor):
        return {"x": None, "cover": None, "present": False}
    i = int(np.argmax(bundle["pos_weight"]))
    return {"x": float(bundle["pos_x"][i]), "cover": int(bundle["pos_m"][i]),
            "present": True}


def s_clock_acquirable(bundle, floor_overlap=FLOOR_OVERLAP,
                       floor_carrier=FLOOR_CARRIER):
    """固有時間の取得可否判定（連続重み → 真偽）。
    v1 g_matter_clock の acquirable と一致。戻り値 {"omega", "acquirable"}。"""
    ok = (bundle["overlap"] > floor_overlap
          and np.sqrt(bundle["carrier_power"]) > floor_carrier)
    return {"omega": float(bundle["phase"]) if ok else None,
            "acquirable": bool(ok)}


def s_support_count(bundle):
    """台集合の指示: セルごとの「厳密非零の参加本数」。
    【警告】これは個数の量ではない。監査で count ≡ M×(power>0) と実測され、
    帳簿の台集合以外の情報を持たない（数値塵も凝縮体も同じ M 本になる）。
    台集合そのものを見たい場合に限り、選択として宣言して使うこと。
    （G の cell_support は台集合の濃度という純粋な構造量であり、それ自体は
    正しい。誤りは v1 でこれを「粒子数＝波の本数」と呼んだことにある。）"""
    return bundle["cell_support"]


def s_occupancy(bundle, eps):
    """占有数 n = セルパワー / ε（ε = 1 局在ドメインあたりのパワー）。
    ε は較正実験（run_tb_epsilon_calibration_v1.py）で決める測定値であり、
    G には置けない（パラメータフリー規約 R2）。ε を宣言して初めて
    「何個」と言える——これが本系における個数の定義である。"""
    return bundle["cell_power"] / eps


def s_unwrap_phase(phases, prev_unwrapped=None):
    """周期量の直線展開（巻数の累積）: 円上の位相列 → 連続量。
    G は円上の値しか返さない（R7）。直線時間は選択層で作られる帳簿である。"""
    ph = np.asarray(phases, float)
    out = np.unwrap(ph)
    if prev_unwrapped is not None:
        out = out - out[0] + prev_unwrapped
    return out
