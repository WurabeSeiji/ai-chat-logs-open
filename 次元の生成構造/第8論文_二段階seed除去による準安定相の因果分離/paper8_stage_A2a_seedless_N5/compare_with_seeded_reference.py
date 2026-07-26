#!/usr/bin/env python3
"""Stage A1b seedあり基準とA2a seedなし系列を絶対step・f水準で比較する。"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RAW = HERE / "raw"
PROCESSED = HERE / "processed"
CONFIG = json.loads((HERE / "config_locked.json").read_text(encoding="utf-8"))
EXPECTED = json.loads((HERE / "expected_hashes.json").read_text(encoding="utf-8"))
RUN = CONFIG["run_ids"][0]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def write_csv_md(stem: str, columns: list[str], rows: list[dict]) -> None:
    write_csv(PROCESSED / f"{stem}.csv", columns, rows)
    md = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for row in rows:
        md.append("| " + " | ".join(str(row.get(c, "")) for c in columns) + " |")
    (PROCESSED / f"{stem}.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def f17(value) -> str:
    if value == "" or value is None:
        return ""
    return format(float(value), ".17e")


def main() -> None:
    gate_path = PROCESSED / "exec_comparison_summary.json"
    if not gate_path.is_file():
        raise SystemExit("EXECUTION_FAILED: compare_execs.pyの結果がない")
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if gate.get("status") != "COMPARED" or not gate.get("numerical_health_passed"):
        raise SystemExit("EXECUTION_FAILED: exec比較または数値健全性gate失敗")
    for group in ("stage_a0_inputs", "stage_a1b_inputs"):
        for item in EXPECTED[group].values():
            path = REPO / item["path"]
            if not path.is_file() or sha256(path) != item["sha256"]:
                raise SystemExit(f"SOURCE_MISMATCH: {path}")

    PROCESSED.mkdir(parents=True, exist_ok=True)
    a1b_dir = HERE.parent / "paper7_N5_transition_anatomy" / "processed"
    a0_f_path = REPO / EXPECTED["stage_a0_inputs"]["fcurve_N00005_delta1e-15_seed0.csv"]["path"]
    seeded_fp = read_csv(a1b_dir / "f_first_passage_levels.csv")
    seeded_rates = read_csv(a1b_dir / "f_decade_growth_rates.csv")
    seeded_q = read_csv(a1b_dir / "first_passage_nearest_q_records.csv")
    seeded_occ = read_csv(a1b_dir / "first_passage_nearest_occupation_records.csv")
    seedless_passages = read_csv(RAW / RUN / "first_passage_measurements.csv")
    seedless_f = read_csv(RAW / RUN / "f_timeseries.csv")

    if not (len(seeded_fp) == len(seeded_q) == len(seeded_occ) == len(seedless_passages) == 31):
        raise SystemExit("EXECUTION_FAILED: f水準表の行数対応が不明")

    f_by_step = {int(r["step"]): float(r["f"]) for r in seedless_f}
    seedless_fp_rows = []
    for raw in seedless_passages:
        found = raw["status"] == "found"
        step = int(raw["first_passage_step"]) if found else None
        previous = step - 1 if found and step > 0 else None
        seedless_fp_rows.append({
            "level_index": raw["level_index"],
            "level": raw["level"],
            "level_label": raw["level_label"],
            "status": raw["status"],
            "first_passage_step": "" if step is None else step,
            "f_at_first_passage": raw["f_at_first_passage"],
            "previous_step": "" if previous is None else previous,
            "f_at_previous_step": "" if previous is None else f17(f_by_step[previous]),
            "q3_over_q1": raw["q3_over_q1"],
            "q4_over_q1": raw["q4_over_q1"],
            "direction_3_occupation": raw["direction_3_occupation"],
            "direction_4_occupation": raw["direction_4_occupation"],
            "kernel_occupation": raw["kernel_occupation"],
        })
    fp_columns = list(seedless_fp_rows[0])
    write_csv_md("seedless_first_passage_levels", fp_columns, seedless_fp_rows)

    seedless_rate_rows = []
    for i in range(len(seedless_fp_rows) - 1):
        lower, upper = seedless_fp_rows[i], seedless_fp_rows[i + 1]
        found = lower["status"] == upper["status"] == "found"
        delta_step = int(upper["first_passage_step"]) - int(lower["first_passage_step"]) if found else None
        log_delta = math.log(float(upper["level"]) / float(lower["level"]))
        rate = log_delta / delta_step if found and delta_step and delta_step > 0 else None
        seedless_rate_rows.append({
            "pair_index": i,
            "lower_level": lower["level"],
            "upper_level": upper["level"],
            "lower_level_label": lower["level_label"],
            "upper_level_label": upper["level_label"],
            "status": "found" if rate is not None else "not_found",
            "lower_first_passage_step": lower["first_passage_step"],
            "upper_first_passage_step": upper["first_passage_step"],
            "step_difference": "" if delta_step is None else delta_step,
            "threshold_log_amplitude_difference": f17(log_delta),
            "mean_exponential_rate_per_step": "" if rate is None else f17(rate),
        })
    rate_columns = list(seedless_rate_rows[0])
    write_csv_md("seedless_decade_growth_rates", rate_columns, seedless_rate_rows)

    by_level_rows: list[dict] = []
    for i, (sl, sd, sq, so) in enumerate(zip(seedless_passages, seeded_fp, seeded_q, seeded_occ)):
        level_values = [float(sl["level"]), float(sd["level"]), float(sq["level"]), float(so["level"])]
        if max(level_values) - min(level_values) > max(level_values) * 1e-14:
            raise SystemExit(f"EXECUTION_FAILED: level mapping mismatch at row {i}")
        sl_found = sl["status"] == "found"
        sd_found = sd["status"] == "found"
        q_fields = ("q3_over_q1", "q4_over_q1")
        occ_fields = ("direction_3_occupation", "direction_4_occupation", "kernel_occupation")
        row = {
            "level_index": i,
            "level": sl["level"],
            "level_label": sl["level_label"],
            "seeded_status": sd["status"],
            "seeded_first_passage_step": sd["first_passage_step"],
            "seedless_status": sl["status"],
            "seedless_first_passage_step": sl["first_passage_step"],
            "seedless_minus_seeded_first_passage_step": (
                int(sl["first_passage_step"]) - int(sd["first_passage_step"])
                if sl_found and sd_found else ""
            ),
            "seedless_measurement_step": sl["first_passage_step"],
            "seeded_q_before_step": sq["before_or_at_step"],
            "seeded_q_after_step": sq["after_or_at_step"],
            "seeded_occupation_before_step": so["before_or_at_step"],
            "seeded_occupation_after_step": so["after_or_at_step"],
        }
        for field in q_fields:
            before = float(sq[f"before_{field}"])
            after = float(sq[f"after_{field}"])
            value = float(sl[field]) if sl_found else None
            row[f"seeded_{field}_before"] = f17(before)
            row[f"seeded_{field}_after"] = f17(after)
            row[f"seedless_{field}_exact"] = "" if value is None else f17(value)
            row[f"seedless_minus_seeded_{field}_before"] = "" if value is None else f17(value - before)
            row[f"seedless_minus_seeded_{field}_after"] = "" if value is None else f17(value - after)
        for field in occ_fields:
            before = float(so[f"before_{field}"])
            after = float(so[f"after_{field}"])
            value = float(sl[field]) if sl_found else None
            row[f"seeded_{field}_before"] = f17(before)
            row[f"seeded_{field}_after"] = f17(after)
            row[f"seedless_{field}_exact"] = "" if value is None else f17(value)
            row[f"seedless_minus_seeded_{field}_before"] = "" if value is None else f17(value - before)
            row[f"seedless_minus_seeded_{field}_after"] = "" if value is None else f17(value - after)
        by_level_rows.append(row)
    by_level_columns = list(by_level_rows[0])
    write_csv_md("seeded_vs_seedless_by_f_level", by_level_columns, by_level_rows)

    growth_rows = []
    if len(seeded_rates) != len(seedless_rate_rows):
        raise SystemExit("EXECUTION_FAILED: growth-rate row mapping mismatch")
    for sl, sd in zip(seedless_rate_rows, seeded_rates):
        seedless_rate = float(sl["mean_exponential_rate_per_step"]) if sl["mean_exponential_rate_per_step"] else None
        seeded_rate = float(sd["mean_exponential_rate_per_step"]) if sd["status"] == "found" else None
        growth_rows.append({
            "pair_index": sl["pair_index"],
            "lower_level": sl["lower_level"],
            "upper_level": sl["upper_level"],
            "seeded_step_difference": sd["step_difference"],
            "seedless_step_difference": sl["step_difference"],
            "seedless_minus_seeded_step_difference": (
                int(sl["step_difference"]) - int(sd["step_difference"])
                if sl["step_difference"] != "" and sd["step_difference"] != "" else ""
            ),
            "seeded_mean_exponential_rate_per_step": "" if seeded_rate is None else f17(seeded_rate),
            "seedless_mean_exponential_rate_per_step": "" if seedless_rate is None else f17(seedless_rate),
            "seedless_minus_seeded_rate": (
                "" if seedless_rate is None or seeded_rate is None else f17(seedless_rate - seeded_rate)
            ),
        })
    growth_columns = list(growth_rows[0])
    write_csv_md("seeded_vs_seedless_growth_rate", growth_columns, growth_rows)

    seeded_f_rows = read_csv(a0_f_path)
    seeded_map = {int(r["tau"]): float(r["f"]) for r in seeded_f_rows}
    seedless_map = {int(r["step"]): float(r["f"]) for r in seedless_f}
    align_level = float(CONFIG["time_alignment_level"])
    seeded_align = next(step for step in sorted(seeded_map) if seeded_map[step] >= align_level)
    seedless_align = next(step for step in sorted(seedless_map) if seedless_map[step] >= align_level)
    relative_min = max(min(seeded_map) - seeded_align, min(seedless_map) - seedless_align)
    relative_max = min(max(seeded_map) - seeded_align, max(seedless_map) - seedless_align)
    aligned_rows = []
    log_diffs = []
    f_diffs = []
    for rel in range(relative_min, relative_max + 1):
        seeded_f_value = seeded_map[rel + seeded_align]
        seedless_f_value = seedless_map[rel + seedless_align]
        log_seeded = math.log10(seeded_f_value) if seeded_f_value > 0 else float("nan")
        log_seedless = math.log10(seedless_f_value) if seedless_f_value > 0 else float("nan")
        f_diff = seedless_f_value - seeded_f_value
        log_diff = log_seedless - log_seeded
        aligned_rows.append({
            "relative_step_from_f_ge_1e-12": rel,
            "seeded_absolute_step": rel + seeded_align,
            "seedless_absolute_step": rel + seedless_align,
            "seeded_f": f17(seeded_f_value),
            "seedless_f": f17(seedless_f_value),
            "seedless_minus_seeded_f": f17(f_diff),
            "seeded_log10_f": f17(log_seeded),
            "seedless_log10_f": f17(log_seedless),
            "seedless_minus_seeded_log10_f": f17(log_diff),
        })
        if np.isfinite(log_diff):
            log_diffs.append(log_diff)
        f_diffs.append(f_diff)
    write_csv(PROCESSED / "time_aligned_f_comparison.csv", list(aligned_rows[0]), aligned_rows)

    summary = {
        "stage": "A2a",
        "status": "COMPARISON_TABLES_COMPLETE",
        "seeded_reference": "Stage A0 reproduced f + Stage A1b actual saved q/occupation brackets",
        "interpolation_used": False,
        "seeded_q_occupation_policy": "actual before_or_at and after_or_at records are both retained",
        "alignment_level": align_level,
        "seeded_alignment_step": seeded_align,
        "seedless_alignment_step": seedless_align,
        "seedless_minus_seeded_alignment_step": seedless_align - seeded_align,
        "aligned_relative_step_range": [relative_min, relative_max],
        "aligned_max_absolute_f_difference": float(np.max(np.abs(f_diffs))),
        "aligned_max_absolute_log10_f_difference": float(np.max(np.abs(log_diffs))),
        "aligned_log10_f_rmse": float(np.sqrt(np.mean(np.square(log_diffs)))),
        "seedless_initial_f": float(seedless_f[0]["f"]),
        "seeded_initial_f": float(seeded_f_rows[0]["f"]),
        "seedless_levels_found": sum(r["status"] == "found" for r in seedless_passages),
        "seeded_levels_found": sum(r["status"] == "found" for r in seeded_fp),
    }
    (PROCESSED / "seeded_comparison_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(
        "COMPARISON_TABLES_COMPLETE: "
        f"alignment seeded={seeded_align}, seedless={seedless_align}, "
        f"offset={seedless_align - seeded_align}"
    )


if __name__ == "__main__":
    main()
