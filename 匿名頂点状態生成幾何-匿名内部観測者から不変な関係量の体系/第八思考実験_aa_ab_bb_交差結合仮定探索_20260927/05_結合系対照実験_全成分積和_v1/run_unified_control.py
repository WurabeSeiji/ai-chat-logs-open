#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""結合系の対照実験を、外から与えられた実行パターンで実行する。

実行パターン（辞書、または JSON ファイル）
  pattern_id              保存先フォルダの名前になる
  charge_multiple_n       電荷倍数 n
  sign_a, sign_b          a と b の電荷の符号
  macrosteps              実行する macrostep の数。状態の値では止めない
  echo_every_macrosteps   進み具合を表示する間隔
  rows_per_part           1 つの HDF5 ファイルに入れる行数

保存
  初期状態と、macrostep ごとの行を、1 行も省かずに HDF5 に保存する。
  1 行は step と、関係ごとの P, E, U_re, U_im, H_re, H_im, Q, N, C, D。
  microstep の途中の状態は保存しない。
"""
from __future__ import annotations
import hashlib
import json
import time
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

import unified_engine as eng
import unified_init as init

HERE = Path(__file__).resolve().parent
REQUIRED = ("pattern_id", "charge_multiple_n", "sign_a", "sign_b", "macrosteps",
            "echo_every_macrosteps", "rows_per_part")


def echo(text):
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {text}", flush=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def hms(seconds):
    seconds = int(round(seconds))
    return f"{seconds // 3600:d}:{seconds % 3600 // 60:02d}:{seconds % 60:02d}"


def write_part(out, index, blocks, columns, rows_per_part):
    block = np.concatenate(blocks)
    part = out / f"raw_macro_part{index:03d}.h5"
    with h5py.File(part, "w") as hf:
        hf.create_dataset("raw_macro_trajectory", data=block, dtype="f8",
                          chunks=(min(rows_per_part, block.shape[0]), block.shape[1]),
                          compression="lzf", shuffle=True)
        hf.attrs["columns_json"] = json.dumps(columns)
        hf.attrs["start_step"] = int(block[0, 0])
        hf.attrs["row_count"] = int(block.shape[0])
        hf.attrs["format_note"] = "Every macrostep row preserved; serialization is outside the state update."
    size = part.stat().st_size
    echo(f"saved {part.name}  rows {block.shape[0]:,}  steps {int(block[0, 0]):,}-{int(block[-1, 0]):,}  {size / 1e6:.1f} MB")
    return {"file": part.name, "start_step": int(block[0, 0]), "row_count": int(block.shape[0]),
            "bytes": int(size), "sha256": sha256(part)}


def run(pattern, out=None):
    for key in REQUIRED:
        if key not in pattern:
            raise SystemExit(f"実行パターンに {key} がない")
    macrosteps = int(pattern["macrosteps"])
    echo_every = int(pattern["echo_every_macrosteps"])
    rows_per_part = int(pattern["rows_per_part"])
    if macrosteps < 1 or echo_every < 1 or rows_per_part < 1:
        raise SystemExit("macrosteps, echo_every_macrosteps, rows_per_part は 1 以上")
    out = Path(out) if out is not None else HERE / "results" / str(pattern["pattern_id"])
    out.mkdir(parents=True, exist_ok=True)
    if any(out.glob("raw_macro_part*.h5")):
        raise SystemExit(f"保存先に raw_macro_part*.h5 が既にある: {out}")

    psi, names, metas, normalization = init.initial_state(pattern)
    kind, arg = init.action_definition(psi.shape[0])
    save_index, saved_names = init.saved_components(names)
    columns = ["step"] + saved_names
    n = psi.shape[0]
    psi_next = np.empty(n, dtype=np.float64)
    S = np.empty((n, n), dtype=np.float64)
    p_columns = [i for i, c in enumerate(columns) if c.endswith(".P")]

    rows_total = macrosteps + 1
    raw_bytes = rows_total * len(columns) * 8
    (out / "run_definition.json").write_text(json.dumps({
        "pattern": pattern,
        "normalization": normalization,
        "initializer": "paper6_init_G_q0_c1.init_case",
        "supplement_meta": metas,
        "state_dimension": n,
        "S_components": n * n,
        "S_component_kinds": init.describe(kind),
        "microsteps_per_macrostep": eng.MICROSTEPS_PER_MACRO,
        "stop": "macrostep count given from outside; no stop by state value",
        "nan_rule": "a product term that is NaN is not added",
        "saved_columns": columns,
        "rows_to_save": rows_total,
        "component_names": names,
        "initial_state": [float(v) for v in psi],
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    echo(f"pattern {pattern['pattern_id']}  n={pattern['charge_multiple_n']}  sign_a={pattern['sign_a']}  sign_b={pattern['sign_b']}")
    for m in metas:
        echo(f"  relation {m['case_id']}  C={m['C']}  D={m['D']}  N={m['N']}")
    echo(f"state components {n}  S components {n * n}  microsteps per macrostep {eng.MICROSTEPS_PER_MACRO}")
    echo(f"macrosteps {macrosteps:,}  rows to save {rows_total:,}  columns {len(columns)}  "
         f"raw size {raw_bytes / 1e6:.1f} MB (before compression)")
    echo(f"output {out}")

    echo("compiling ...")
    t_c = time.time()
    eng.generate_chunk(psi.copy(), psi_next, S, kind, arg, save_index, 0, 1, 0)   # 1 行だけ作る。状態は進めない
    echo(f"compiled in {time.time() - t_c:.1f} s.  start")

    t0 = time.time()
    parts, blocks, rows_in_part, rows_saved = [], [], 0, 0
    step0, skipped_total, first_skip = 0, 0, None
    finished = False
    while not finished:
        chunk_rows = min(echo_every, rows_per_part - rows_in_part)
        rows, count, skipped, first, finished = eng.generate_chunk(
            psi, psi_next, S, kind, arg, save_index, step0, chunk_rows, macrosteps)
        blocks.append(rows[:count].copy())
        rows_in_part += count
        rows_saved += count
        skipped_total += int(skipped)
        if first_skip is None and first >= 0:
            first_skip = int(first)
            echo(f"NaN term appeared first at macrostep {first_skip:,}; NaN terms are not added")
        step0 += count
        last = rows[count - 1]
        done = int(last[0])
        elapsed = time.time() - t0
        rate = done / elapsed if elapsed > 0 and done > 0 else 0.0
        remain = (macrosteps - done) / rate if rate > 0 else float("nan")
        p_text = " ".join(f"{columns[i]}={last[i]:.9g}" for i in p_columns)
        echo(f"macrostep {done:>10,} / {macrosteps:,} ({100.0 * done / macrosteps:6.2f}%)  "
             f"elapsed {hms(elapsed)}  {rate:,.0f} macro/s  remaining {hms(remain) if rate > 0 else '-'}  "
             f"rows {rows_saved:,}  NaN terms not added {skipped_total:,}  {p_text}")
        if rows_in_part >= rows_per_part or finished:
            parts.append(write_part(out, len(parts), blocks, columns, rows_per_part))
            blocks, rows_in_part = [], 0

    np.save(out / "final_full_state.npy", psi)
    summary = {
        "pattern_id": pattern["pattern_id"],
        "status": "completed",
        "macro_steps": macrosteps,
        "microsteps": macrosteps * eng.MICROSTEPS_PER_MACRO,
        "rows_saved": rows_saved,
        "columns": len(columns),
        "nan_terms_not_added": skipped_total,
        "first_macrostep_with_nan_term": first_skip,
        "final_state": {name: float(v) for name, v in zip(names, psi)},
        "raw_macro_parts": [p["file"] for p in parts],
        "raw_bytes_on_disk": int(sum(p["bytes"] for p in parts)),
        "elapsed_seconds": time.time() - t0,
    }
    (out / "generator_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (out / "part_manifest.json").write_text(json.dumps({"pattern_id": pattern["pattern_id"], "parts": parts},
                                                       ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    echo(f"completed  rows {rows_saved:,}  files {len(parts)}  "
         f"{summary['raw_bytes_on_disk'] / 1e6:.1f} MB on disk  elapsed {hms(summary['elapsed_seconds'])}")
    return summary
