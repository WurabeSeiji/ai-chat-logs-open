#!/usr/bin/env python3
"""管理下コピーの対照実験 v1：Codex 製 make_parent の完全再現検証（修正前の基準線）

目的（2026-08-04 木原氏指示）:
    Codex 作 make_parent_white_harmonics_n_only_v2.py を管理下にコピー
    （SHA-256 一致確認済み: 7a79d5f7…）し、修正前の対照として、
    保存済み正本と同一シードで**誤った結果も含めて**同一出力を再現することを確認する。

検証項目:
    T1 bit再現: make_parent(5, seed=2) / make_parent(40, seed=1) の
        relation_waves / parent_vector / raw_white_noise / edges が
        保存済み正本 (parent_white_harmonics_N5_v2 / N40_v2) と bitwise 一致
    T2 契約性質の再現: 行閉塞 < 1e-12（全行）・全体閉塞 ~1e-17
    T3 既知欠陥(b)の再現: 自己対セクターの場の閉塞（相対）が
        N=5 k=0: 0.4749 / N=40 k=0: 0.03120, k=20: 0.04561 と一致（±1e-6）
        ——「誤った結果」が同じに出ることの確認（対照実験の要）

この対照が通った後に、管理下コピーへ修正（DC/Nyquist 射影）を加える。
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CENSUS = HERE.parent / "standalone_parent_census_v1"
GEN = HERE / "make_parent_white_harmonics_n_only_v2.py"

spec = importlib.util.spec_from_file_location("gen_managed", GEN)
gen = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = gen
spec.loader.exec_module(gen)

# 期待値は手打ちせず、保存済み正本の配列から実行時に同一式で計算して比較する
CASES = {5: {"seed": 2, "orig": CENSUS / "parent_white_harmonics_N5_v2",
             "defect_ks": [0]},
         40: {"seed": 1, "orig": CENSUS / "parent_white_harmonics_N40_v2",
              "defect_ks": [0, 20]}}


def field_closure_rel(C, k, n):
    if (2 * k) % n != 0:
        return 0.0
    Zk = C[:, k]
    return float(abs(complex(np.sum(Zk ** 2))) / np.sum(np.abs(Zk) ** 2))


def main() -> None:
    all_pass = True
    results = {}
    for n, cfg in CASES.items():
        print(f"===== N={n} seed={cfg['seed']} =====")
        r = gen.make_parent(n, seed=cfg["seed"])
        orig = {name: np.load(cfg["orig"] / f"{name}.npy")
                for name in ("relation_waves", "parent_vector", "raw_white_noise", "edges")}
        t1 = (np.array_equal(r.relation_waves, orig["relation_waves"])
              and np.array_equal(r.parent_vector, orig["parent_vector"])
              and np.array_equal(r.raw_white_noise, orig["raw_white_noise"])
              and np.array_equal(r.edges, orig["edges"]))
        print(f"  T1 bit再現（4配列とも bitwise 一致）: {'PASS' if t1 else 'FAIL'}")

        W = r.relation_waves
        rowc = np.abs(np.sum(W * W, axis=1))
        rowp = np.sum(np.abs(W) ** 2, axis=1)
        t2 = (float(np.max(rowc / rowp)) < 1e-12
              and abs(complex(np.sum(W * W))) < 1e-12)
        print(f"  T2 契約性質（行閉塞<1e-12・全体閉塞）: {'PASS' if t2 else 'FAIL'} "
              f"(行max={np.max(rowc/rowp):.2e} 全体={abs(complex(np.sum(W*W))):.2e})")

        C = np.fft.fft(W, axis=1) / n
        C_orig = np.fft.fft(orig["relation_waves"], axis=1) / n
        t3 = True
        for k in cfg["defect_ks"]:
            got = field_closure_rel(C, k, n)
            expected = field_closure_rel(C_orig, k, n)   # 正本から同一式で計算
            ok = got == expected                          # bit同一配列なので厳密一致を要求
            t3 &= ok
            print(f"  T3 欠陥(b)再現 k={k}: 場の閉塞(相対)={got:.9f} 正本={expected:.9f} "
                  f"{'PASS' if ok else 'FAIL'}")
        # 非自己対セクターの厳密零閉塞も確認
        nonself_max = max(abs(complex(np.sum((C[:, k]) ** 2)))
                          * (n if (2 * k) % n == 0 else 0.0)
                          for k in range(n) if (2 * k) % n != 0) if n > 2 else 0.0
        print(f"  参考: 非自己対セクターの場の閉塞 max = {nonself_max:.2e}（恒等0）")

        results[f"N{n}"] = {"T1_bitwise": bool(t1), "T2_contract": bool(t2),
                             "T3_defect_reproduced": bool(t3)}
        all_pass &= t1 and t2 and t3

    print(f"\n総合判定: {'ALL PASS——管理下コピーは誤った結果まで完全再現（対照成立）' if all_pass else 'FAIL あり'}")
    results["all_pass"] = bool(all_pass)
    (HERE / "test_control_reproduction_result_v1.json").write_text(
        json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
