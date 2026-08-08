#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""論文7 N=5,40,300 有効方向 比較（§15 表・図10 q3q4・図11 s5/s1）。解釈なし。"""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = Path(__file__).resolve().parent
BASE = CODE.parent
NS = [5, 40, 300]


def load(n):
    return json.load(open(BASE / "diagnostics" / f"N{n:05d}_saturation.json"))


def main():
    cd = BASE / "figures" / "comparison"; cd.mkdir(parents=True, exist_ok=True)
    data = {}
    for n in NS:
        try:
            data[n] = load(n)
        except FileNotFoundError:
            print(f"  N={n} saturation JSON なし")
    ns = [n for n in NS if n in data]

    # 図10: q3,q4 (crossing相対) 各N
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    for n in ns:
        raw = BASE / "raw" / f"N{n:05d}_dimension_saturation_v2" / f"q_svd_N{n:05d}.csv"
        r = list(csv.DictReader(open(raw)))
        rt = np.array([float(x["relative_time"]) for x in r])
        a1.plot(rt, [float(x["q3"]) for x in r], lw=0.7, label=f"N={n}")
        a2.plot(rt, [float(x["q4"]) for x in r], lw=0.7, label=f"N={n}")
    a1.set_title("Fig10a: q3 (all N)"); a1.set_xlabel("time - crossing"); a1.legend()
    a2.set_title("Fig10b: q4 (all N)"); a2.set_xlabel("time - crossing"); a2.legend()
    fig.tight_layout(); fig.savefig(cd / "sat_fig10_q3q4_N5_40_300.png", dpi=130); plt.close(fig)

    # 図11: s5/s1 と 全時間結合 s_j/s_1
    fig, ax = plt.subplots(figsize=(9, 6))
    for n in ns:
        s = np.array(data[n]["global_s"])
        ax.semilogy(np.arange(1, len(s) + 1), np.clip(s / s[0], 1e-20, None), "o-", ms=4, label=f"N={n}")
    ax.axhline(1e-8, color="r", ls="--", lw=0.8, label="1e-8")
    ax.set_xlabel("index j"); ax.set_ylabel("s_j / s_1 (log)")
    ax.set_title("Fig11: global combined s_j/s_1 (all N)"); ax.legend()
    fig.tight_layout(); fig.savefig(cd / "sat_fig11_sj_s1_N5_40_300.png", dpi=130); plt.close(fig)

    # §15 比較表
    occ = {5: (1.000000, 4.55e-14), 40: (1.000000, 8.27e-18), 300: (1.000000, 3.0e-17)}  # v2/paper6由来
    table = {}
    for n in ns:
        d = data[n]
        cl = d["closure"]
        eps_max = max(v["eps_fro_max"] for v in cl.values())
        eta_max = max(v["eta4_max"] for v in cl.values())
        table[n] = {
            "M": d["M"], "gram_rank": d["gram_rank_final"],
            "q3_final": d["q3_final"], "q4_final": d["q4_final"], "rank_q": d["rank_q_final"],
            "s5_over_s1": d["s5_over_s1"], "max_closure_eps_fro": eps_max, "max_eta4": eta_max,
            "stream_dim_primary": d["stream_final_dim_by_tau"].get("1e-10"),
            "dom_occupation": occ.get(n, (None, None))[0],
            "max_nondom_occupation": occ.get(n, (None, None))[1],
        }
    with open(cd / "saturation_comparison_table.json", "w", encoding="utf-8") as fh:
        json.dump(table, fh, indent=2, ensure_ascii=False)
    print("[比較] fig10/11 + saturation_comparison_table.json")
    print(f"{'N':>4} {'M':>7} {'gram_rank':>9} {'q3':>8} {'q4':>8} {'rank_q':>6} {'s5/s1':>10} {'閉鎖max':>9} {'stream次数':>10}")
    for n in ns:
        t = table[n]
        print(f"{n:>4} {t['M']:>7} {t['gram_rank']:>9} {t['q3_final']:>8.4f} {t['q4_final']:>8.4f} "
              f"{t['rank_q']:>6} {t['s5_over_s1']:>10.3e} {t['max_closure_eps_fro']:>9.2e} {t['stream_dim_primary']:>10}")


if __name__ == "__main__":
    main()
