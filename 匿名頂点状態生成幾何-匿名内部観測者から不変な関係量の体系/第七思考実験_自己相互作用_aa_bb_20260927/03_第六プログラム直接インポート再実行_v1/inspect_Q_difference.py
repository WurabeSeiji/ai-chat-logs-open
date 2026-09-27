#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q と t_readout の違いが、どの step で、どれだけの大きさで出ているかを調べる。

今回の self_aa_n1・self_aa_n3 と、第六の保存データ n1_repulsive・n3_repulsive を、全行読んで比べる。読むだけ。
差は「その値の最小刻み（ulp）の何個ぶんか」で数える。
結果は comparison_results/Q_difference_detail.json に書く。
"""
from __future__ import annotations
import json
import numpy as np
import h5py
from wrapper_common import HERE
from compare_with_previous import P6, part_files, paper6_raw_files, Stream, total_rows, BLOCK

COL_Q, COL_T = 7, 11


def inspect(new_case: str, paper6_case: str) -> dict:
    new = part_files(HERE / "cases" / new_case)
    old, kind = paper6_raw_files(P6 / paper6_case)
    n = min(total_rows(new), total_rows(old))
    sa, sb = Stream(new), Stream(old)
    bins = {}
    ulp_hist_q, ulp_hist_t = {}, {}
    first_diff = None
    last_diff_step = None
    examples = []
    done = 0
    while done < n:
        k = min(BLOCK, n - done)
        a, b = sa.take(k), sb.take(k)
        step = a[:, 0].astype(np.int64)
        for col, hist in ((COL_Q, ulp_hist_q), (COL_T, ulp_hist_t)):
            x, y = a[:, col], b[:, col]
            fin = np.isfinite(x) & np.isfinite(y)
            x, y = x[fin], y[fin]
            if x.size == 0:
                continue
            ulp = np.spacing(np.maximum(np.abs(x), np.abs(y)))
            d = np.rint((x - y) / ulp).astype(np.int64)
            vals, counts = np.unique(d, return_counts=True)
            for v, c in zip(vals, counts):
                hist[int(v)] = hist.get(int(v), 0) + int(c)
        qa, qb = a[:, COL_Q], b[:, COL_Q]
        dq = ~((qa == qb) | (np.isnan(qa) & np.isnan(qb)))
        for s0 in range(int(step[0]) // 100_000 * 100_000, int(step[-1]) + 1, 100_000):
            m = (step >= s0) & (step < s0 + 100_000)
            e = bins.setdefault(s0, [0, 0])
            e[0] += int(m.sum())
            e[1] += int((dq & m).sum())
        if dq.any():
            idx = np.nonzero(dq)[0]
            if first_diff is None:
                i = int(idx[0])
                first_diff = {"step": int(step[i]), "Q_new": float(a[i, COL_Q]).hex(), "Q_paper6": float(b[i, COL_Q]).hex(),
                              "Q_new_decimal": repr(float(a[i, COL_Q])), "Q_paper6_decimal": repr(float(b[i, COL_Q]))}
                for j in idx[:5]:
                    examples.append({"step": int(step[j]), "Q_new": float(a[j, COL_Q]).hex(), "Q_paper6": float(b[j, COL_Q]).hex()})
            last_diff_step = int(step[int(idx[-1])])
        done += k
    # 最終行
    with h5py.File(new[-1], "r") as hf:
        last_new = hf["raw_macro_trajectory"][-1]
    with h5py.File(old[-1], "r") as hf:
        last_old = hf["raw_macro_trajectory"][-1]
    out = {
        "compared": new_case + " (this run) vs paper 6 " + paper6_case + " (" + kind + ")",
        "rows_compared": int(n),
        "Q_difference_in_ulp_histogram": {str(k): v for k, v in sorted(ulp_hist_q.items())},
        "t_readout_difference_in_ulp_histogram": {str(k): v for k, v in sorted(ulp_hist_t.items())},
        "first_row_with_Q_difference": first_diff,
        "first_examples": examples,
        "last_step_with_Q_difference": last_diff_step,
        "final_row": {"step": int(last_new[0]), "Q_new": float(last_new[COL_Q]).hex(), "Q_paper6": float(last_old[COL_Q]).hex(),
                      "Q_equal": bool(last_new[COL_Q] == last_old[COL_Q] or (np.isnan(last_new[COL_Q]) and np.isnan(last_old[COL_Q])))},
        "rows_with_Q_difference_per_100000_steps": {str(k): {"rows": v[0], "rows_with_Q_difference": v[1]} for k, v in sorted(bins.items())},
    }
    print(json.dumps({k: v for k, v in out.items() if k != "rows_with_Q_difference_per_100000_steps"}, ensure_ascii=False, indent=1))
    return out


def main() -> None:
    result = {"n1": inspect("self_aa_n1", "n1_repulsive"), "n3": inspect("self_aa_n3", "n3_repulsive")}
    (HERE / "comparison_results" / "Q_difference_detail.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
