#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第7論文 v2：N=5 vs N=40 比較図（§5 の 18,19,20）＋比較表。解釈なし。"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CODE = Path(__file__).resolve().parent
BASE = CODE.parent
NS = [5, 40]


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def crossing(n):
    return json.load(open(BASE / "diagnostics" / f"N{n:05d}.json"))["crossing"]


def q_series(n):
    r = rows(BASE / "raw" / f"N{n:05d}" / "q_svd.csv")
    t = np.array([int(float(x["time"])) for x in r])
    return t, np.array([float(x["q3"]) for x in r]), np.array([float(x["q4"]) for x in r])


def dom_origin(n):
    r = rows(BASE / "raw" / f"N{n:05d}" / "initial_space_origin.csv")
    byt = defaultdict(list)
    for x in r:
        byt[int(float(x["step"]))].append(x)
    ts = sorted(byt); of = []; onf = []
    for t in ts:
        dom = max(byt[t], key=lambda x: float(x["sigma_representative"]))
        of.append(float(dom["overlap_with_initial_floor_space"]))
        onf.append(float(dom["overlap_with_initial_nonfloor_space"]))
    return np.array(ts), np.array(of), np.array(onf)


def dom_occ(n):
    r = rows(BASE / "raw" / f"N{n:05d}" / "clusters.csv")
    byt = defaultdict(list)
    for x in r:
        byt[int(float(x["time"]))].append(x)
    ts = sorted(byt); occ = []; nondom = []
    for t in ts:
        cs = byt[t]
        dom = max(cs, key=lambda x: float(x["sigma_over_sigma1"]))
        occ.append(float(dom["occupation_fraction"]))
        others = [float(x["occupation_fraction"]) for x in cs if x is not dom]
        nondom.append(max(others) if others else 0.0)
    return np.array(ts), np.array(occ), np.array(nondom)


def main():
    cd = BASE / "figures" / "comparison"; cd.mkdir(parents=True, exist_ok=True)
    # 図18: q3,q4 比較
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    for n in NS:
        t, q3, q4 = q_series(n); cr = crossing(n)
        a1.plot(t - cr, q3, lw=0.7, label=f"N={n}")
        a2.plot(t - cr, q4, lw=0.7, label=f"N={n}")
    a1.set_title("Fig18a: q3 = 3rd singular value of [B0|Bdom]"); a1.set_xlabel("time - crossing"); a1.set_ylabel("q3"); a1.legend()
    a2.set_title("Fig18b: q4 = 4th singular value"); a2.set_xlabel("time - crossing"); a2.set_ylabel("q4"); a2.legend()
    fig.tight_layout(); fig.savefig(cd / "fig18_q3q4_compare.png", dpi=130); plt.close(fig)

    # 図19: 支配平面 origin 比較
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    for n in NS:
        t, of, onf = dom_origin(n); cr = crossing(n)
        a1.plot(t - cr, of, lw=0.7, label=f"N={n}")
        a2.plot(t - cr, onf, lw=0.7, label=f"N={n}")
    a1.set_title("Fig19a: dominant plane O_initial_floor"); a1.set_xlabel("time - crossing"); a1.set_ylabel("O_floor"); a1.legend()
    a2.set_title("Fig19b: dominant plane O_initial_nonfloor"); a2.set_xlabel("time - crossing"); a2.set_ylabel("O_nonfloor"); a2.legend()
    fig.tight_layout(); fig.savefig(cd / "fig19_dominant_origin_compare.png", dpi=130); plt.close(fig)

    # 図20: 瞬時支配占有 比較
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    for n in NS:
        t, occ, nd = dom_occ(n); cr = crossing(n)
        a1.plot(t - cr, occ, lw=0.7, label=f"N={n}")
        a2.semilogy(t - cr, np.clip(nd, 1e-20, None), lw=0.7, label=f"N={n}")
    a1.set_title("Fig20a: instantaneous dominant occupation E_dom/|Z|^2"); a1.set_xlabel("time - crossing"); a1.set_ylabel("E_dom"); a1.legend()
    a2.set_title("Fig20b: max non-dominant occupation (log)"); a2.set_xlabel("time - crossing"); a2.set_ylabel("max non-dom occ"); a2.legend()
    fig.tight_layout(); fig.savefig(cd / "fig20_dominant_occupation_compare.png", dpi=130); plt.close(fig)

    # 比較表（最終時刻）
    comp = {}
    for n in NS:
        end = max(int(float(x["time"])) for x in rows(BASE / "raw" / f"N{n:05d}" / "q_svd.csv"))
        qr = [x for x in rows(BASE / "raw" / f"N{n:05d}" / "q_svd.csv") if int(float(x["time"])) == end][0]
        cl = [x for x in rows(BASE / "raw" / f"N{n:05d}" / "clusters.csv") if int(float(x["time"])) == end]
        dom = max(cl, key=lambda x: float(x["sigma_over_sigma1"]))
        nondom = max((float(x["occupation_fraction"]) for x in cl if x is not dom), default=0.0)
        org = [x for x in rows(BASE / "raw" / f"N{n:05d}" / "initial_space_origin.csv") if int(float(x["step"])) == end]
        domo = max(org, key=lambda x: float(x["sigma_representative"]))
        band = sorted(float(x["sigma_over_sigma1"]) for x in cl)
        comp[n] = {
            "q3": float(qr["q3"]), "q4": float(qr["q4"]),
            "rank_B0_Bdom": int(sum(float(qr[f"q{i}"]) > 1e-8 for i in (1, 2, 3, 4))),
            "dom_O_floor": float(domo["overlap_with_initial_floor_space"]),
            "dom_O_nonfloor": float(domo["overlap_with_initial_nonfloor_space"]),
            "dom_occupation": float(dom["occupation_fraction"]),
            "max_nondom_occupation": nondom,
            "nondom_sigma_band": [band[0], band[-2] if len(band) > 1 else band[0]],
        }
    with open(cd / "comparison_table_N5_N40.json", "w", encoding="utf-8") as fh:
        json.dump(comp, fh, indent=2, ensure_ascii=False)
    print("[comparison] fig18/19/20 + comparison_table_N5_N40.json")
    for n in NS:
        c = comp[n]
        print(f"  N={n}: q3={c['q3']:.4f} q4={c['q4']:.4f} rank={c['rank_B0_Bdom']} "
              f"dom_Ofloor={c['dom_O_floor']:.4f} dom_occ={c['dom_occupation']:.6f} maxnondom={c['max_nondom_occupation']:.2e}")


if __name__ == "__main__":
    main()
