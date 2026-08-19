#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""DL3 解析 v1 — ホロノミーイベント集合 H と系譜交換集合 L の完全対応表（第6回査読対応）

|H|=18・|L|=18・|H∩L|=17 なら「交換を伴わないホロノミーイベント」が必ず存在する
（集合論的必然）。保存済み系列 dl3_branch_lineage_v1.npz から：
 (1) χ=-1 窓を併合してホロノミーイベント（クラスタ）H を列挙
 (2) 各クラスタの被覆区間（窓 k は步 k..k+W を覆う）に 3↔4 系譜交換が入るか照合
 (3) 三分類 A（交換あり・χ=-1）/ B（交換あり・χ=+1＝静かな系譜交換）/
     C（交換なし・χ=-1＝ホロノミーのみ）の完全な対応表を出力
 (4) C イベントの時刻・窓数・窓内最小 Δ34 を特定

出力: result_dl3_holonomy_exchange_map_v1.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
W = 10


def main():
    D = np.load(HERE / "dl3_branch_lineage_v1.npz")
    tau, lam, ex34, chi = D["tau"], D["lam"], D["ex34"], D["chi"]
    ev = chi < 0
    n_win = len(chi)
    d34 = lam[:, 2] - lam[:, 3]

    # (1) クラスタ併合
    ev_int = ev.astype(int)
    starts = np.where(np.diff(np.concatenate([[0], ev_int])) == 1)[0]
    ends = np.where(np.diff(np.concatenate([ev_int, [0]])) == -1)[0]
    ex_idx = np.where(ex34)[0]  # 交換は步 k→k+1（窓指数系）

    clusters = []
    matched_ex = set()
    for s, e in zip(starts, ends):
        cover = (s, e + W)     # クラスタ窓群が覆う步範囲
        ex_in = [int(k) for k in ex_idx if cover[0] <= k <= cover[1]]
        matched_ex.update(ex_in)
        seg = d34[s:min(e + W + 1, len(d34))]
        clusters.append({
            "tau_range": [int(tau[s]), int(tau[min(e + W, len(tau) - 1)])],
            "n_windows": int(e - s + 1),
            "n_exchanges_in_cover": len(ex_in),
            "exchange_taus": [int(tau[k]) for k in ex_in],
            "min_d34_in_cover": float(seg.min()),
        })

    only_holonomy = [c for c in clusters if c["n_exchanges_in_cover"] == 0]
    multi_ex = [c for c in clusters if c["n_exchanges_in_cover"] >= 2]
    only_exchange = [int(tau[k]) for k in ex_idx if int(k) not in matched_ex]

    res = {
        "config": {"W": W, "source": "dl3_branch_lineage_v1.npz（決定的系列の事後解析）"},
        "n_holonomy_events": len(clusters),
        "n_exchanges": int(ex34.sum()),
        "n_intersection": len(clusters) - len(only_holonomy),
        "A_both": len(clusters) - len(only_holonomy),
        "B_exchange_only_taus": only_exchange,
        "C_holonomy_only": only_holonomy,
        "clusters_with_multiple_exchanges": multi_ex,
        "clusters": clusters,
    }
    (HERE / "result_dl3_holonomy_exchange_map_v1.json").write_text(
        json.dumps(res, indent=1, ensure_ascii=False))
    print(f"H={len(clusters)}  L={int(ex34.sum())}  "
          f"H∩L={res['n_intersection']}")
    print(f"B（交換のみ・静かな系譜交換）: {only_exchange}")
    print(f"C（ホロノミーのみ）: "
          f"{[(c['tau_range'], c['n_windows'], '%.2e' % c['min_d34_in_cover']) for c in only_holonomy]}")
    if multi_ex:
        print(f"複数交換を含むクラスタ: {[(c['tau_range'], c['exchange_taus']) for c in multi_ex]}")


if __name__ == "__main__":
    main()
