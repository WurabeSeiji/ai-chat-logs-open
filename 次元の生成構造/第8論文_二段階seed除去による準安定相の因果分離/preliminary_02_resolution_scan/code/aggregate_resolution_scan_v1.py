#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第2予備実験 §12 集計CSV と §14 再現性差分。解釈しない（集計・突合のみ）。

summary/all_runs_manifest.csv, all_runs_final_values.csv, all_stop_reasons.csv, all_diagnostics.csv。
（all_fixed_band_regressions.csv は compute_prefixed_regressions_v1 が出力）
再現性: exec1 vs exec2 を突合し、完全一致しない場合のみ diagnostics/reproducibility_diff_<run_id>.csv。
"""
import csv
import json
from pathlib import Path

CODE = Path(__file__).resolve().parent
P2 = CODE.parent
RAW = P2 / "raw"
SUM = P2 / "summary"; SUM.mkdir(exist_ok=True)
DIAG = P2 / "diagnostics"; DIAG.mkdir(exist_ok=True)

CMP_COLS = ["step", "stop_reason", "f_outside", "a_outside", "quantization_l2"]


def cfg_of(d):
    return json.load(open(d / "run_config.json"))


def diag_of(d):
    return json.load(open(d / "run_diagnostics.json"))


def rows_of(d):
    return list(csv.DictReader(open(d / "timeseries.csv")))


def main():
    runs = sorted([d for d in RAW.iterdir() if d.is_dir() and (d / "run_config.json").exists()])

    # manifest
    with open(SUM / "all_runs_manifest.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run_id", "execution_index", "N", "M", "p", "Delta_ref", "Delta_actual",
                    "resolution_operator", "max_step"])
        for d in runs:
            c = cfg_of(d)
            w.writerow([c["run_id"], c["execution_index"], c["N"], c["M"], c["p"], c["Delta_ref"],
                        c["Delta_actual"], c["resolution_operator"], c["max_step"]])

    # final values
    with open(SUM / "all_runs_final_values.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run_id", "N", "p", "Delta_ref", "Delta_actual", "resolution_operator",
                    "stop_step", "stop_reason", "final_f_outside", "final_a_outside", "final_q3", "final_q4",
                    "final_rank_Q", "n_saved_rows"])
        for d in runs:
            c = cfg_of(d); g = diag_of(d); r = rows_of(d)
            last = r[-1] if r else {}
            w.writerow([c["run_id"], c["N"], c["p"], c["Delta_ref"], c["Delta_actual"], c["resolution_operator"],
                        g["stop_step"], g["stop_reason"], g["final_f_outside"], g["final_a_outside"],
                        last.get("q3", ""), last.get("q4", ""), last.get("rank_Q", ""), g["n_saved_rows"]])

    # stop reasons
    with open(SUM / "all_stop_reasons.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["run_id", "N", "p", "Delta_ref", "resolution_operator", "stop_reason", "stop_step"])
        for d in runs:
            c = cfg_of(d); g = diag_of(d)
            w.writerow([c["run_id"], c["N"], c["p"], c["Delta_ref"], c["resolution_operator"],
                        g["stop_reason"], g["stop_step"]])

    # diagnostics
    with open(SUM / "all_diagnostics.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        keys = ["run_id", "stop_reason", "stop_step", "n_saved_rows", "final_Z_sha256",
                "max_quantization_l2", "max_retraction_correction_l2", "max_closure_after_reproj",
                "max_norm_dev_after_reproj", "final_f_outside", "final_a_outside"]
        w.writerow(keys)
        for d in runs:
            g = diag_of(d); w.writerow([g.get(k, "") for k in keys])

    # reproducibility: pair exec1/exec2 by base run_id
    bases = {}
    for d in runs:
        c = cfg_of(d)
        base = c["run_id"].rsplit("_exec", 1)[0]
        bases.setdefault(base, {})[c["execution_index"]] = d
    n_mismatch = 0
    for base, ex in sorted(bases.items()):
        if 1 not in ex or 2 not in ex:
            continue
        d1, d2 = ex[1], ex[2]
        g1, g2 = diag_of(d1), diag_of(d2)
        r1, r2 = rows_of(d1), rows_of(d2)
        diffs = []
        if g1["stop_step"] != g2["stop_step"]:
            diffs.append(["stop_step", "", g1["stop_step"], g2["stop_step"], ""])
        if g1["stop_reason"] != g2["stop_reason"]:
            diffs.append(["stop_reason", "", g1["stop_reason"], g2["stop_reason"], ""])
        if g1.get("final_Z_sha256") != g2.get("final_Z_sha256"):
            diffs.append(["final_Z_sha256", "", g1.get("final_Z_sha256"), g2.get("final_Z_sha256"), ""])
        n = min(len(r1), len(r2))
        for i in range(n):
            for col in ("f_outside", "a_outside", "quantization_l2"):
                v1, v2 = r1[i].get(col, ""), r2[i].get(col, "")
                if v1 != v2:
                    try:
                        dd = abs(float(v1) - float(v2))
                    except ValueError:
                        dd = ""
                    diffs.append([col, r1[i]["step"], v1, v2, dd])
        if len(r1) != len(r2):
            diffs.append(["n_saved_rows", "", len(r1), len(r2), ""])
        if diffs:
            n_mismatch += 1
            with open(DIAG / f"reproducibility_diff_{base}.csv", "w", newline="") as fh:
                w = csv.writer(fh); w.writerow(["field", "step", "exec1", "exec2", "abs_diff"]); w.writerows(diffs)

    # reproducibility summary line
    with open(DIAG / "reproducibility_overview.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["n_base_runs", "n_mismatch"]); w.writerow([len(bases), n_mismatch])
    print(f"[aggregate] summary CSV 4種 + 再現性: base={len(bases)} 不一致={n_mismatch}")


if __name__ == "__main__":
    main()
