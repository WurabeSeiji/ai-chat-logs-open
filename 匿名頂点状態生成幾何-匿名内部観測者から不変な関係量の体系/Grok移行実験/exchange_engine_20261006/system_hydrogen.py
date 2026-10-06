#!/usr/bin/env python3
"""系の部品：水素。固定一式 ../exchange_rel_20261005/exchange_cascade.py の準位表と A 係数を、エンジンの形式に写すだけ。"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
FIXED = HERE.parent / "exchange_rel_20261005"
sys.path.insert(0, str(FIXED))
from matplotlib import font_manager  # noqa: E402
font_manager.fontManager.addfont = lambda path: None
import exchange_cascade as X  # noqa: E402
import engine as EN  # noqa: E402

RANK = {"e1": 1, "m1": 1, "e2": 2, "gw": 2, "2ph": 0}


def build(gw_weight=None):
    gw_weight = X.W_GW if gw_weight is None else gw_weight
    states = []
    for i, (n, l, se, sp) in enumerate(X.ST):
        j, F = X.jf_of(l, se, sp)
        states.append(dict(label=(n, l, se, sp), name="%d%s_{%d/2} F=%d" % (n, "spdfg"[l], int(2 * j), int(F)),
                           E_eV=X.ETOT[i] * X.EV, J=F, M_kg=(1.0 + X.MU_P + X.ETOT[i]) * EN.ME_SI, n=n, l=l, j=j))
    channels = [(src, dst, kind, RANK[kind], A) for (src, dst, kind, A) in X.edges(gw_weight) if A > 0.0]
    init = X.IDX[(5, 1, 1, -1)]
    return states, channels, init


if __name__ == "__main__":
    states, channels, init = build()
    print("states %d, channels %d, init %s" % (len(states), len(channels), states[init]["name"]))
