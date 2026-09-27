#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の状態生成器を読み込み、初期状態を補遺の初期化プログラムで作って実行する。

初期状態
  補遺 paper6_init_G_q0_c1.py の init_case(case_id) が作る 41 成分をそのまま使う。
  init_case は case_relations と init_state_from_relations を呼ぶ。

第六の生成器 run_strict_charged8_5cases.py で変えるもの
  core.init_state  第六の run_case が初期状態を受け取る関数。補遺が作った初期状態を返すものに替える。
  core.CASES       実行条件（補遺の init_case が返した lambda_A, lambda_B, C, D）
  core.HERE        all_generator_summaries.json の保存先
  core.OUTROOT     ケース別フォルダの保存先

transition、macrostep、generate_chunk、保存は第六のプログラムのまま。
ケースごとに第六の run_case を呼び、第六の main() と同じ形式で要約を書く。
"""
from __future__ import annotations
import json
from wrapper_common import HERE, PROGRAMS, case_table, load_program, load_supplement


def main() -> None:
    supplement = load_supplement()
    table = case_table()
    core = load_program("run_strict_charged8_5cases.py")

    # 第六のプログラムは読み込み時に「自分の隣の cases/」を作る。使わないので、空なら消す。
    unused = PROGRAMS / "cases"
    if unused.is_dir() and not any(unused.iterdir()):
        unused.rmdir()

    core.HERE = HERE
    core.OUTROOT = HERE / "cases"
    core.OUTROOT.mkdir(exist_ok=True)
    core.CASES = tuple((case_id, la, lb, c0, d0) for (case_id, n, sa, sb, la, lb, c0, d0) in table)

    allsum = []
    records = []
    for c in core.CASES:
        case_id = c[0]
        z, meta = supplement.init_case(case_id)      # 補遺の初期化プログラムが作る初期状態
        calls = []

        def init_state(p0=50.0, n0=0.25, c0=1.09, d0=0.36, _z=z, _calls=calls):
            _calls.append({"p0": p0, "n0": n0, "c0": c0, "d0": d0})
            return _z

        core.init_state = init_state
        print("condition", c, flush=True)
        allsum.append(core.run_case(*c))
        if len(calls) != 1:
            raise RuntimeError(f"{case_id}: run_case called init_state {len(calls)} times")
        print(allsum[-1]["case_id"], allsum[-1]["macro_steps"], allsum[-1]["first_nonfinite_time_step"], flush=True)
        records.append({"case_id": case_id, "initializer": "paper6_init_G_q0_c1.init_case",
                        "init_state_calls_in_run_case": len(calls), "supplement_meta": meta})

    # 第六の main() と同じ形式
    (HERE / "all_generator_summaries.json").write_text(json.dumps(allsum, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    # 実行条件の記録（補遺のプログラムが持つ値と、補遺の init_case が返した値を書き出すだけ）
    (HERE / "CASE_TABLE.json").write_text(json.dumps({
        "normalization": {"G": supplement.G, "q0": supplement.Q0, "c": supplement.C_LIGHT,
                          "M_over_q0": str(supplement.M_OVER_Q0_EXACT)},
        "charge_of_a": "+n*q0",
        "charge_of_b": "-n*q0",
        "condition_table_given_to_supplement": [list(x) for x in supplement.CASES],
        "cases": records,
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
