#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""初期状態と S の定義を作って、内訳を表示する。状態の更新は実行しない。

  python3 check_definition.py patterns/<パターン>.json
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
import unified_engine as eng
import unified_init as init


def main() -> None:
    pattern = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    print("実行パターン", pattern)
    psi, names, metas, normalization = init.initial_state(pattern)
    kind, arg = init.action_definition(psi.shape[0])
    n = psi.shape[0]
    print("規格化", normalization)
    print("状態の成分の数", n, "/ S の成分の数", n * n)
    for m in metas:
        print("  ", m["case_id"], "lambda_A", m["lambda_A"], "lambda_B", m["lambda_B"], "N", m["N"], "C", m["C"], "D", m["D"])
    print("初期状態で 0 でない成分:")
    for name, v in zip(names, psi):
        if v != 0.0:
            print(f"   {name:<12} {v!r}")
    print("S の成分の種類ごとの数:")
    for k, c in init.describe(kind).items():
        print(f"   {k:<11} {c:>6}")
    nonzero = int((kind != eng.K_ZERO).sum())
    print("値が 0 でない成分になり得る数", nonzero, "/ 常に 0 の成分の数", n * n - nonzero)
    # 交差成分（行と列の関係が違う成分）がすべて 0 の成分であること
    owner = np.array([nm.split(".")[0] for nm in names])
    cross = owner[:, None] != owner[None, :]
    print("交差成分の数", int(cross.sum()), "/ そのうち 0 以外の種類の成分", int((kind[cross] != eng.K_ZERO).sum()))
    # 各成分の関数が読む成分は、その行と同じ関係の成分だけであること
    bad = 0
    for i in range(n):
        for j in range(n):
            if kind[i, j] not in (eng.K_ZERO, eng.K_ONE, eng.K_RATE_DE, eng.K_CONST_H):
                used = {eng.K_IDENT: 1, eng.K_RATE_DP: 4, eng.K_RATE_DT: 2, eng.K_STAGE_HALF: 2, eng.K_STAGE_FULL: 2,
                        eng.K_COMBINE: 4, eng.K_MID: 2, eng.K_ADD: 2, eng.K_ROT_U_RE: 2, eng.K_ROT_U_IM: 2,
                        eng.K_ROT_H_RE: 3, eng.K_ROT_H_IM: 3, eng.K_CLOCK: 2}[int(kind[i, j])]
                for t in range(used):
                    if owner[arg[i, j, t]] != owner[i]:
                        bad += 1
    print("ほかの関係の成分を読む成分の数", bad)
    save_index, saved_names = init.saved_components(names)
    rows = int(pattern["macrosteps"]) + 1
    print("macrostep ごとに保存する項目", len(saved_names), "個:", ", ".join(saved_names))
    print("保存する行数", f"{rows:,}", "/ 列数", 1 + len(saved_names), "/ 圧縮前の大きさ", f"{rows * (1 + len(saved_names)) * 8 / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
