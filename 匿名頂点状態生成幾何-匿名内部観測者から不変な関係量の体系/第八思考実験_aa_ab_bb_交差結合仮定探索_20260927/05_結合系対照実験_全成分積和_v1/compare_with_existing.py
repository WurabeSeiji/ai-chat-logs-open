#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""結合系の実行結果を、同じ実験パターンの既存データと全行比べる。

  python3 compare_with_existing.py comparisons/<比較の定義>.json

読むだけで、どのデータも書き換えない。
比べる行は、両方にある macrostep のすべて（先頭から、行数の少ない方の最後まで）。
比べる項目
  状態 10 項目   P, E, U_re, U_im, H_re, H_im, Q, N, C, D
  読み出し 3 項目 t = TAU0*log(Q),  x = P*H_re,  y = P*H_im
    結合系の読み出しは、保存した状態から第六思考実験の readout と同じ式で計算する。
結果は comparison_results/<comparison_id>/ に JSON と Markdown で保存する。
"""
from __future__ import annotations
import json
import math
import sys
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
from numba import njit

import unified_engine as eng

HERE = Path(__file__).resolve().parent
SERIES = HERE.parents[1]
STATE_ITEMS = ("P", "E", "U_re", "U_im", "H_re", "H_im", "Q", "N", "C", "D")
READOUT_ITEMS = ("t_readout", "x_readout", "y_readout")


@njit(cache=True)
def readout(P, H_re, H_im, Q, t, x, y):
    """第六思考実験の readout と同じ式。"""
    for k in range(P.shape[0]):
        t[k] = eng.TAU0 * math.log(Q[k])
        x[k] = P[k] * H_re[k]
        y[k] = P[k] * H_im[k]


def read_parts(folder: Path):
    """raw_macro_partNNN.h5 を番号順に読み、行をつないで返す。"""
    files = sorted(folder.glob("raw_macro_part[0-9][0-9][0-9].h5"))
    if not files:
        raise RuntimeError(f"no raw_macro_partNNN.h5 in {folder}")
    blocks, columns, expected = [], None, 0
    for f in files:
        with h5py.File(f, "r") as hf:
            ds = hf["raw_macro_trajectory"]
            cols = json.loads(hf.attrs["columns_json"])
            if columns is None:
                columns = cols
            elif cols != columns:
                raise RuntimeError(f"column mismatch in {f}")
            if int(hf.attrs["start_step"]) != expected:
                raise RuntimeError(f"{f.name}: start_step {int(hf.attrs['start_step'])} != {expected}")
            blocks.append(ds[:])
            expected += ds.shape[0]
    data = np.concatenate(blocks)
    steps = data[:, columns.index("step")]
    if not np.array_equal(steps, np.arange(data.shape[0], dtype=np.float64)):
        raise RuntimeError(f"steps are not 0,1,2,... in {folder}")
    return data, columns, [f.name for f in files]


def compare_column(a, b):
    """同じ長さの 2 列を比べる。"""
    same = (a == b) | (np.isnan(a) & np.isnan(b))
    differ = ~same
    count = int(differ.sum())
    finite = np.isfinite(a) & np.isfinite(b)
    out = {
        "rows": int(a.shape[0]),
        "differing_values": count,
        "first_differing_step": int(np.argmax(differ)) if count else None,
        "last_differing_step": int(a.shape[0] - 1 - np.argmax(differ[::-1])) if count else None,
        "nonfinite_in_unified": int((~np.isfinite(a)).sum()),
        "nonfinite_in_existing": int((~np.isfinite(b)).sum()),
        "max_abs_difference": 0.0,
        "max_relative_difference": 0.0,
        "difference_in_smallest_steps_min": 0,
        "difference_in_smallest_steps_max": 0,
    }
    m = differ & finite
    if m.any():
        d = a[m] - b[m]
        scale = np.maximum(np.abs(a[m]), np.abs(b[m]))
        out["max_abs_difference"] = float(np.abs(d).max())
        out["max_relative_difference"] = float((np.abs(d) / scale).max())
        ulp = np.rint(d / np.spacing(scale)).astype(np.int64)
        out["difference_in_smallest_steps_min"] = int(ulp.min())
        out["difference_in_smallest_steps_max"] = int(ulp.max())
    return out


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: python3 compare_with_existing.py comparisons/<definition>.json")
    definition = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    cid = definition["comparison_id"]
    out = HERE / "comparison_results" / cid
    out.mkdir(parents=True, exist_ok=True)

    uni_dir = HERE / definition["unified_result"]
    uni, uni_cols, uni_files = read_parts(uni_dir)
    summary = json.loads((uni_dir / "generator_summary.json").read_text(encoding="utf-8"))
    run_def = json.loads((uni_dir / "run_definition.json").read_text(encoding="utf-8"))
    print(f"unified: {uni.shape[0]:,} rows, {len(uni_cols)} columns", flush=True)

    result = {
        "comparison_id": cid,
        "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "unified": {"folder": definition["unified_result"], "rows": int(uni.shape[0]), "files": uni_files,
                    "pattern": run_def["pattern"], "nan_terms_not_added": summary["nan_terms_not_added"]},
        "relations": {},
    }

    for rel, base in definition["baselines"].items():
        ref, ref_cols, ref_files = read_parts(SERIES / base["folder"])
        n = min(uni.shape[0], ref.shape[0])
        print(f"{rel}: existing {ref.shape[0]:,} rows ({base['label']}); compare {n:,} rows", flush=True)
        r = {"existing": dict(base, rows=int(ref.shape[0]), files=ref_files),
             "rows_unified": int(uni.shape[0]), "rows_existing": int(ref.shape[0]), "rows_compared": int(n),
             "unified_rows_without_existing_data": int(max(0, uni.shape[0] - ref.shape[0])),
             "state": {}, "readout": {}}
        col = {name: uni[:n, uni_cols.index(f"{rel}.{name}")] for name in STATE_ITEMS}
        for name in STATE_ITEMS:
            r["state"][name] = compare_column(col[name], ref[:n, ref_cols.index(name)])
        t, x, y = np.empty(n), np.empty(n), np.empty(n)
        readout(np.ascontiguousarray(col["P"]), np.ascontiguousarray(col["H_re"]),
                np.ascontiguousarray(col["H_im"]), np.ascontiguousarray(col["Q"]), t, x, y)
        for name, values in zip(READOUT_ITEMS, (t, x, y)):
            r["readout"][name] = compare_column(values, ref[:n, ref_cols.index(name)])
        last = n - 1
        r["last_compared_row"] = {"step": last,
                                  "unified": {k: float(col[k][last]) for k in STATE_ITEMS},
                                  "existing": {k: float(ref[last, ref_cols.index(k)]) for k in STATE_ITEMS}}
        r["state_items_all_same"] = [k for k in STATE_ITEMS if r["state"][k]["differing_values"] == 0]
        r["state_items_with_difference"] = [k for k in STATE_ITEMS if r["state"][k]["differing_values"] > 0]
        result["relations"][rel] = r

    # 結合系の中での aa と bb
    aa = np.column_stack([uni[:, uni_cols.index(f"aa.{k}")] for k in STATE_ITEMS])
    bb = np.column_stack([uni[:, uni_cols.index(f"bb.{k}")] for k in STATE_ITEMS])
    result["aa_vs_bb_in_unified"] = {k: compare_column(aa[:, i], bb[:, i]) for i, k in enumerate(STATE_ITEMS)}

    (out / "comparison.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    # 分析結果の Markdown
    L = []
    pat = run_def["pattern"]
    L += [f"# 結合系の対照実験と既存データの比較（{cid}）", "",
          f"作成: {result['created']}　プログラム: `compare_with_existing.py`　数値の元: `comparison.json`", "",
          "## 比べたもの", "",
          f"結合系の実行パターン `{pat['pattern_id']}`（n={pat['charge_multiple_n']}、a の符号 {pat['sign_a']:+d}、b の符号 {pat['sign_b']:+d}、"
          f"{pat['macrosteps']:,} macrostep）。保存した行は {uni.shape[0]:,} 行。足さなかった NaN の項は {summary['nan_terms_not_added']:,} 個。", "",
          "| 関係 | 既存データ | 既存データを作ったプログラム | 実行環境 | 既存の行数 | 比べた行数 |", "|---|---|---|---|---:|---:|"]
    for rel, r in result["relations"].items():
        e = r["existing"]
        L.append(f"| {rel} | {e['label']} | {e['program']} | {e['environment']} | {e['rows']:,} | {r['rows_compared']:,} |")
    L += ["", "比べたのは、両方にある macrostep のすべてである。間引いていない。", ""]
    for rel, r in result["relations"].items():
        if r["unified_rows_without_existing_data"]:
            L.append(f"{rel} は、結合系の最後の {r['unified_rows_without_existing_data']:,} 行に対応する既存データがない"
                     f"（既存データは macrostep {r['rows_existing'] - 1:,} で終わる）。")
    L += ["", "## 状態 10 項目の比較", "",
          "| 関係 | 項目 | 比べた行数 | 違う値の個数 | 最初に違う macrostep | 最大の差 | 最大の相対差 | 最小刻みで数えた差 |",
          "|---|---|---:|---:|---:|---:|---:|---|"]
    for rel, r in result["relations"].items():
        for k in STATE_ITEMS:
            c = r["state"][k]
            first = "—" if c["first_differing_step"] is None else f"{c['first_differing_step']:,}"
            ulp = "—" if c["differing_values"] == 0 else f"{c['difference_in_smallest_steps_min']:+d} 〜 {c['difference_in_smallest_steps_max']:+d}"
            L.append(f"| {rel} | {k} | {c['rows']:,} | {c['differing_values']:,} | {first} | "
                     f"{c['max_abs_difference']:.3e} | {c['max_relative_difference']:.3e} | {ulp} |")
    L += ["", "## 読み出し 3 項目の比較", "",
          "結合系の読み出しは、保存した状態から `t = TAU0*log(Q)`、`x = P*H_re`、`y = P*H_im` で計算した。", "",
          "| 関係 | 項目 | 比べた行数 | 違う値の個数 | 最初に違う macrostep | 最大の差 | 最大の相対差 | 最小刻みで数えた差 |",
          "|---|---|---:|---:|---:|---:|---:|---|"]
    for rel, r in result["relations"].items():
        for k in READOUT_ITEMS:
            c = r["readout"][k]
            first = "—" if c["first_differing_step"] is None else f"{c['first_differing_step']:,}"
            ulp = "—" if c["differing_values"] == 0 else f"{c['difference_in_smallest_steps_min']:+d} 〜 {c['difference_in_smallest_steps_max']:+d}"
            L.append(f"| {rel} | {k} | {c['rows']:,} | {c['differing_values']:,} | {first} | "
                     f"{c['max_abs_difference']:.3e} | {c['max_relative_difference']:.3e} | {ulp} |")
    L += ["", "## 関係ごとのまとめ", "",
          "| 関係 | 全行で同じだった状態項目 | 違いがあった状態項目 |", "|---|---|---|"]
    for rel, r in result["relations"].items():
        L.append(f"| {rel} | {'、'.join(r['state_items_all_same']) or 'なし'} | {'、'.join(r['state_items_with_difference']) or 'なし'} |")
    L += ["", "## 比べた最後の行の値", "",
          "| 関係 | macrostep | 項目 | 結合系 | 既存データ |", "|---|---:|---|---:|---:|"]
    for rel, r in result["relations"].items():
        for k in ("P", "Q", "H_re", "H_im"):
            L.append(f"| {rel} | {r['last_compared_row']['step']:,} | {k} | {r['last_compared_row']['unified'][k]!r} | {r['last_compared_row']['existing'][k]!r} |")
    L += ["", "## 結合系の中での aa と bb", "",
          "| 項目 | 行数 | 違う値の個数 |", "|---|---:|---:|"]
    for k in STATE_ITEMS:
        c = result["aa_vs_bb_in_unified"][k]
        L.append(f"| {k} | {c['rows']:,} | {c['differing_values']:,} |")
    L += ["", "## 分かったこと", ""]
    for rel, r in result["relations"].items():
        e = r["existing"]
        if not r["state_items_with_difference"]:
            L.append(f"- {rel}: 状態 10 項目のすべてが、比べた {r['rows_compared']:,} 行の全行で既存データと同じであった。")
        else:
            L.append(f"- {rel}: 状態 10 項目のうち {len(r['state_items_all_same'])} 項目（{'、'.join(r['state_items_all_same'])}）が、"
                     f"比べた {r['rows_compared']:,} 行の全行で既存データと同じであった。"
                     f"違いがあったのは {'、'.join(r['state_items_with_difference'])} である。")
            for k in r["state_items_with_difference"]:
                c = r["state"][k]
                L.append(f"  - {k}: 違う値は {c['differing_values']:,} 個。差は最小刻みで数えて "
                         f"{c['difference_in_smallest_steps_min']:+d} から {c['difference_in_smallest_steps_max']:+d}、"
                         f"相対差の最大は {c['max_relative_difference']:.3e}。既存データの実行環境: {e['environment']}。")
        ro_same = [k for k in READOUT_ITEMS if r["readout"][k]["differing_values"] == 0]
        ro_diff = [k for k in READOUT_ITEMS if r["readout"][k]["differing_values"] > 0]
        L.append(f"  - 読み出し: 同じだった項目 {'、'.join(ro_same) or 'なし'}。違いがあった項目 {'、'.join(ro_diff) or 'なし'}。")
    same_ab = all(result["aa_vs_bb_in_unified"][k]["differing_values"] == 0 for k in STATE_ITEMS)
    L.append(f"- 結合系の中で、aa と bb の状態 10 項目は全行で{'同じであった' if same_ab else '違いがあった'}。")
    L += ["", "## 再現", "", "```sh", f"./run_compare.sh comparisons/{cid}.json", "```", ""]
    (out / "comparison_ja.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L), flush=True)


if __name__ == "__main__":
    main()
