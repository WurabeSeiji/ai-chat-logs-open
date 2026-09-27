#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import csv, json, math
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE / "original_initial_rows_paper6.csv"
NEW = HERE / "generated_initial_rows_G_q0_c1.csv"
OUT_JSON = HERE / "initial_row_verification_result.json"
OUT_MD = HERE / "initial_row_verification_result_ja.md"

FIELDS = ["micro","q_index","P","N","C","D","k1p","k2p","k3p","k4p","dP","Q"]


def load(path):
    with path.open(newline="", encoding="utf-8") as f:
        return {r["case_id"]: r for r in csv.DictReader(f)}

old = load(OLD); new = load(NEW)
results = []
all_pass = True
for case in old:
    o = old[case]; n = new[case]
    diffs = {}
    case_pass = True
    for fld in FIELDS:
        ov = float(o[fld]); nv = float(n[fld])
        diff = nv - ov
        exact = (nv == ov)
        diffs[fld] = {"old": ov, "new": nv, "diff": diff, "exact_equal": exact}
        case_pass &= exact
    all_pass &= case_pass
    results.append({"case_id": case, "pass": case_pass, "fields": diffs})

payload = {
    "purpose": "Compare Paper 6 saved micro=0 rows against rows generated under G=q0=c=1.",
    "all_pass": all_pass,
    "comparison": "exact IEEE-754 float equality after CSV parse",
    "cases": results,
}
OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")

lines = [
    "# 第6論文補講: G=q0=c=1 初期化同値検証", "",
    f"- 判定: **{'PASS' if all_pass else 'FAIL'}**", 
    "- 比較対象: 保存済み5ケースの `raw_first_macro_microsteps.csv` の micro=0 行",
    "- 新規側: `G=q0=c=1`, `q0=1`, `M/q0=20/3`, 等質量 `m_A=m_B=M/2` から初期化",
    "- 比較方法: CSV読込後の各診断列の完全一致", "",
    "| case | result | P | N | C | D | Q |",
    "|---|---|---:|---:|---:|---:|---:|",
]
for rec in results:
    f = rec["fields"]
    lines.append(f"| {rec['case_id']} | {'PASS' if rec['pass'] else 'FAIL'} | {f['P']['new']} | {f['N']['new']} | {f['C']['new']} | {f['D']['new']} | {f['Q']['new']} |")
lines += ["", "全ケースで micro, q_index, P, N, C, D, k1p, k2p, k3p, k4p, dP, Q が保存済み初期rowと完全一致した。" if all_pass else "不一致あり。JSONを参照。"]
OUT_MD.write_text("\n".join(lines)+"\n", encoding="utf-8")
print(json.dumps({"all_pass": all_pass, "json": str(OUT_JSON), "md": str(OUT_MD)}, ensure_ascii=False))
if not all_pass:
    raise SystemExit(1)
