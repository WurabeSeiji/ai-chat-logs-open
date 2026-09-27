#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""今回の出力を、保存済みの出力と比べる。読むだけで、どの出力も書き換えない。

比べる相手
  A. 第六思考実験の保存データ（同じ初期状態のケース）
       self_*_n1  <->  n1_repulsive   (C=0.91, D=0)
       self_*_n3  <->  n3_repulsive   (C=0.19, D=0)
  B. 第七思考実験の旧プログラム run_paper7_self_interaction.py の出力
  C. 今回の aa と bb どうし

結果は comparison_results/ に JSON と Markdown で書く。
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

import h5py
import numpy as np

from wrapper_common import HERE, SELF_CASES

SERIES = HERE.parents[1]
P6 = SERIES / "第五思考実験_乗法状態内部展開と近似次数検証_20260924" / "13_charged8_integer_multiple_5case_strict_numeric_20260925" / "02_CASES"
P7OLD = HERE.parent / "01_独立自己相互作用数値実験_v1"
OUT = HERE / "comparison_results"

PAPER6_CASE = {1: "n1_repulsive", 3: "n3_repulsive"}
OLD7_CASE = {1: ("self_n1", "sampled_trajectory.csv"), 3: ("self_n3_fresh", "sampled_trajectory_prefix100k.csv")}
BLOCK = 250_000


def part_files(case_dir: Path):
    m = json.loads((case_dir / "part_manifest.json").read_text(encoding="utf-8"))
    raw = case_dir / "raw" if (case_dir / "raw").is_dir() else case_dir
    return [raw / p["file"] for p in m["parts"]]


def paper6_raw_files(case_dir: Path):
    """第六の保存データのファイルを選ぶ。

    フォルダにある raw_macro_partNNN.h5 を番号順に並べ、各ファイルが持つ start_step と行数が
    切れ目なくつながっていればそれを使う。なければ、行を保ったまま分割し直した
    raw_macro_drive_partNNN.h5 を使う。part_manifest.json はフォルダの中身と合わない場合が
    あるので、ファイル自身が持つ情報で確かめる。
    """
    raw = case_dir / "raw"
    for pattern, key in (("raw_macro_part[0-9][0-9][0-9].h5", "start_step"),
                         ("raw_macro_drive_part[0-9][0-9][0-9].h5", "drive_global_start_row")):
        files = sorted(raw.glob(pattern))
        if not files:
            continue
        expected = 0
        ok = True
        for f in files:
            with h5py.File(f, "r") as hf:
                if int(hf.attrs[key]) != expected:
                    ok = False
                    break
                expected += hf["raw_macro_trajectory"].shape[0]
        if ok:
            return files, pattern
    raise RuntimeError(f"no contiguous raw files in {raw}")


def rows(files):
    for f in files:
        with h5py.File(f, "r") as hf:
            ds = hf["raw_macro_trajectory"]
            for i in range(0, ds.shape[0], BLOCK):
                yield ds[i:min(i + BLOCK, ds.shape[0])]


def columns_of(files):
    with h5py.File(files[0], "r") as hf:
        return json.loads(hf.attrs["columns_json"])


class Stream:
    """複数ファイルにまたがる行を、先頭から指定した行数ずつ取り出す。"""

    def __init__(self, files):
        self.it = rows(files)
        self.buf = np.empty((0, 15))

    def take(self, n):
        while self.buf.shape[0] < n:
            try:
                self.buf = np.concatenate([self.buf, next(self.it)])
            except StopIteration:
                break
        out, self.buf = self.buf[:n], self.buf[n:]
        return out


def total_rows(files):
    n = 0
    for f in files:
        with h5py.File(f, "r") as hf:
            n += hf["raw_macro_trajectory"].shape[0]
    return n


def compare_rows(files_a, files_b):
    """2 つの全行データを、先頭から同じ行どうしで比べる。"""
    cols = columns_of(files_a)
    rows_a, rows_b = total_rows(files_a), total_rows(files_b)
    n = min(rows_a, rows_b)
    sa, sb = Stream(files_a), Stream(files_b)
    differing = np.zeros(len(cols), dtype=np.int64)
    max_abs = np.zeros(len(cols))
    nonfinite_mismatch = np.zeros(len(cols), dtype=np.int64)
    done = 0
    while done < n:
        k = min(BLOCK, n - done)
        a, b = sa.take(k), sb.take(k)
        fa, fb = np.isfinite(a), np.isfinite(b)
        same = (a == b) | (np.isnan(a) & np.isnan(b))
        differing += (~same).sum(axis=0)
        nonfinite_mismatch += (fa != fb).sum(axis=0)
        with np.errstate(invalid="ignore"):
            d = np.where(fa & fb, np.abs(a - b), 0.0)
        max_abs = np.maximum(max_abs, d.max(axis=0))
        done += k
    return {
        "rows_new": int(rows_a), "rows_other": int(rows_b), "rows_compared": int(done),
        "columns": cols,
        "differing_values_per_column": {c: int(v) for c, v in zip(cols, differing)},
        "max_abs_difference_per_column": {c: float(v) for c, v in zip(cols, max_abs)},
        "finite_vs_nonfinite_mismatch_per_column": {c: int(v) for c, v in zip(cols, nonfinite_mismatch)},
        "differing_values_total": int(differing.sum()),
    }


def compare_vectors(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return {"shape_new": list(a.shape), "shape_other": list(b.shape), "comparable": False}
    same = (a == b) | (np.isnan(a) & np.isnan(b))
    both = np.isfinite(a) & np.isfinite(b)
    return {"comparable": True, "values": int(a.size), "differing_values": int((~same).sum()),
            "max_abs_difference": float(np.where(both, np.abs(a - b), 0.0).max()) if a.size else 0.0}


def compare_text(path_a: Path, path_b: Path):
    if not path_b.exists():
        return {"comparable": False, "reason": "other file not found"}
    a, b = path_a.read_bytes(), path_b.read_bytes()
    return {"comparable": True, "bytes_new": len(a), "bytes_other": len(b), "identical": a == b}


def compare_summary(new: dict, other: dict, keys):
    out = {}
    for k_new, k_other in keys:
        vn, vo = new.get(k_new), other.get(k_other)
        out[k_new] = {"new": vn, "other": vo, "equal": vn == vo}
    return out


def sampled_rows_from_new(files, steps):
    """旧プログラムが保存した step と同じ行を、今回の全行データから取り出す。"""
    wanted = np.asarray(sorted(set(int(s) for s in steps)), dtype=np.int64)
    found = {}
    for block in rows(files):
        st = block[:, 0].astype(np.int64)
        lo, hi = int(st[0]), int(st[-1])
        sel = wanted[(wanted >= lo) & (wanted <= hi)]
        for s in sel:
            found[int(s)] = block[int(s) - lo].copy()
    return found


def compare_old_sampled(files_new, csv_path: Path):
    if not csv_path.exists():
        return {"comparable": False, "reason": "old sampled csv not found"}
    old = np.genfromtxt(csv_path, delimiter=",", names=True)
    cols_new = columns_of(files_new)
    common = [c for c in cols_new if c in old.dtype.names]
    found = sampled_rows_from_new(files_new, old["step"])
    differing = {c: 0 for c in common}
    max_abs = {c: 0.0 for c in common}
    missing = 0
    for rec in old:
        s = int(rec["step"])
        if s not in found:
            missing += 1
            continue
        row = found[s]
        for c in common:
            a, b = float(row[cols_new.index(c)]), float(rec[c])
            same = (a == b) or (math.isnan(a) and math.isnan(b))
            if not same:
                differing[c] += 1
                if math.isfinite(a) and math.isfinite(b):
                    max_abs[c] = max(max_abs[c], abs(a - b))
    return {"comparable": True, "old_rows": int(old.shape[0]), "old_rows_not_in_new": missing,
            "columns_compared": common,
            "columns_only_in_old": [c for c in old.dtype.names if c not in cols_new],
            "differing_values_per_column": differing, "max_abs_difference_per_column": max_abs,
            "differing_values_total": int(sum(differing.values()))}


def main() -> None:
    OUT.mkdir(exist_ok=True)
    result = {"A_vs_paper6": {}, "B_vs_old_paper7": {}, "C_aa_vs_bb": {}}
    new_files = {}
    for cid, n, sa, sb in SELF_CASES:
        cdir = HERE / "cases" / cid
        new_files[cid] = part_files(cdir)
        gen_new = json.loads((cdir / "generator_summary.json").read_text(encoding="utf-8"))

        # A. 第六の保存データ
        p6dir = P6 / PAPER6_CASE[n]
        p6files, p6kind = paper6_raw_files(p6dir)
        gen_p6 = json.loads((p6dir / "generator_summary.json").read_text(encoding="utf-8"))
        a = {"paper6_case": PAPER6_CASE[n], "paper6_raw_files_used": p6kind, "paper6_raw_file_count": len(p6files)}
        a["raw_all_rows"] = compare_rows(new_files[cid], p6files)
        a["final_full_state"] = compare_vectors(np.load(cdir / "final_full_state.npy"), np.load(p6dir / "final_full_state.npy"))
        a["raw_first_macro_microsteps_csv"] = compare_text(cdir / "raw_first_macro_microsteps.csv", p6dir / "raw" / "raw_first_macro_microsteps.csv")
        a["generator_summary"] = compare_summary(gen_new, gen_p6, (
            ("macro_steps", "macro_steps"), ("microsteps", "microsteps"), ("final_P", "final_P"),
            ("t_cross_r20", "t_cross_r20"), ("phi_cross_r20", "phi_cross_r20"),
            ("first_nonfinite_time_step", "first_nonfinite_time_step"),
            ("first_nonfinite_time_P", "first_nonfinite_time_P"),
            ("initial_state", "initial_state"), ("max_identity_drift", "max_identity_drift")))
        rc_new, rc_p6 = cdir / "reference_comparison.json", p6dir / "reference_comparison.json"
        if rc_new.exists() and rc_p6.exists():
            rn, ro = json.loads(rc_new.read_text()), json.loads(rc_p6.read_text())
            a["reference_comparison"] = {k: {"new": rn.get(k), "other": ro.get(k), "equal": rn.get(k) == ro.get(k)}
                                         for k in ("generator_crossing", "reference_crossing", "abs_crossing_errors",
                                                   "same_phase_residuals_full_run", "same_time_residuals_finite_readout_only",
                                                   "identity_drifts", "time_readout")}
        result["A_vs_paper6"][cid] = a

        # B. 第七の旧プログラムの出力
        old_name, old_csv = OLD7_CASE[n]
        odir = P7OLD / "results" / old_name
        b = {"old_case": old_name}
        b["sampled_rows"] = compare_old_sampled(new_files[cid], odir / old_csv)
        b["raw_first_macro_microsteps_csv"] = compare_text(cdir / "raw_first_macro_microsteps.csv", odir / "raw_first_macro_microsteps.csv")
        old_sum = json.loads((odir / "summary.json").read_text(encoding="utf-8"))
        b["old_run_reached_target"] = old_sum.get("reached_target")
        if old_sum.get("reached_target"):
            b["final_full_state"] = compare_vectors(np.load(cdir / "final_full_state.npy"), np.load(odir / "final_full_state.npy"))
            b["summary"] = compare_summary(gen_new, old_sum, (("macro_steps", "macro_steps_executed"), ("final_P", "final_P")))
        else:
            b["note"] = "旧プログラムは途中までしか走らせていないので、最終状態と総 step 数は比べられない"
            b["old_macro_steps_executed"] = old_sum.get("macro_steps_executed")
        result["B_vs_old_paper7"][cid] = b

    # C. aa と bb
    for n in (1, 3):
        aa, bb = f"self_aa_n{n}", f"self_bb_n{n}"
        c = {"raw_all_rows": compare_rows(new_files[aa], new_files[bb]),
             "final_full_state": compare_vectors(np.load(HERE / "cases" / aa / "final_full_state.npy"),
                                                 np.load(HERE / "cases" / bb / "final_full_state.npy")),
             "raw_first_macro_microsteps_csv": compare_text(HERE / "cases" / aa / "raw_first_macro_microsteps.csv",
                                                            HERE / "cases" / bb / "raw_first_macro_microsteps.csv")}
        result["C_aa_vs_bb"][f"n{n}"] = c

    (OUT / "comparison_with_previous.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = ["# 今回の出力と保存済み出力の比較", "",
             "`compare_with_previous.py` の出力。数値はすべて `comparison_with_previous.json` と同じ。", ""]
    lines += ["## A. 第六思考実験の保存データとの比較（全行）", "",
              "| 今回のケース | 第六のケース | 今回の行数 | 第六の行数 | 違う値の個数 | 最大の差 | 最終状態 41 成分の違い | 最初の macrostep の CSV |",
              "|---|---|---:|---:|---:|---:|---:|---|"]
    for cid, a in result["A_vs_paper6"].items():
        r = a["raw_all_rows"]
        lines.append(f"| {cid} | {a['paper6_case']} | {r['rows_new']:,} | {r['rows_other']:,} | {r['differing_values_total']:,} | "
                     f"{max(r['max_abs_difference_per_column'].values()):.3e} | {a['final_full_state'].get('differing_values')} | "
                     f"{'同一' if a['raw_first_macro_microsteps_csv'].get('identical') else '相違'} |")
    lines += ["", "### A の内訳：違いがあった列", "",
              "| 今回のケース | 列 | 違う値の個数 | 最大の差 |", "|---|---|---:|---:|"]
    for cid, a in result["A_vs_paper6"].items():
        r = a["raw_all_rows"]
        hit = [c for c in r["columns"] if r["differing_values_per_column"][c]]
        for c in hit:
            lines.append(f"| {cid} | {c} | {r['differing_values_per_column'][c]:,} | {r['max_abs_difference_per_column'][c]:.3e} |")
        if not hit:
            lines.append(f"| {cid} | なし | 0 | 0 |")
    lines += ["", "違いがなかった列: " + "、".join(
        c for c in next(iter(result["A_vs_paper6"].values()))["raw_all_rows"]["columns"]
        if all(a["raw_all_rows"]["differing_values_per_column"][c] == 0 for a in result["A_vs_paper6"].values())), ""]
    lines += ["### A の内訳：要約ファイルの値", "",
              "| 今回のケース | macro_steps | final_P | t_cross_r20 | phi_cross_r20 | first_nonfinite_time_step |", "|---|---|---|---|---|---|"]
    for cid, a in result["A_vs_paper6"].items():
        g = a["generator_summary"]
        lines.append(f"| {cid} | " + " | ".join(("同じ" if g[k]["equal"] else "違う") + f"（{g[k]['new']}）"
                                                 for k in ("macro_steps", "final_P", "t_cross_r20", "phi_cross_r20", "first_nonfinite_time_step")) + " |")
    lines += ["", "## B. 第七思考実験の旧プログラムの出力との比較", "",
              "| 今回のケース | 旧ケース | 旧の保存行数 | 比べた列の数 | 違う値の個数 | 最大の差 | 最終状態 41 成分の違い | 最初の macrostep の CSV |",
              "|---|---|---:|---:|---:|---:|---|---|"]
    for cid, b in result["B_vs_old_paper7"].items():
        s = b["sampled_rows"]
        fs = b.get("final_full_state", {}).get("differing_values", "比較不能（旧は途中まで）")
        lines.append(f"| {cid} | {b['old_case']} | {s.get('old_rows')} | {len(s.get('columns_compared', []))} | {s.get('differing_values_total')} | "
                     f"{max(s.get('max_abs_difference_per_column', {'-': 0.0}).values()):.3e} | {fs} | "
                     f"{'同一' if b['raw_first_macro_microsteps_csv'].get('identical') else '相違'} |")
    lines += ["", "## C. 今回の aa と bb の比較（全行）", "",
              "| n | aa の行数 | bb の行数 | 違う値の個数 | 最大の差 | 最終状態 41 成分の違い |",
              "|---:|---:|---:|---:|---:|---:|"]
    for k, c in result["C_aa_vs_bb"].items():
        r = c["raw_all_rows"]
        lines.append(f"| {k[1:]} | {r['rows_new']:,} | {r['rows_other']:,} | {r['differing_values_total']:,} | "
                     f"{max(r['max_abs_difference_per_column'].values()):.3e} | {c['final_full_state'].get('differing_values')} |")
    (OUT / "comparison_with_previous_ja.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
