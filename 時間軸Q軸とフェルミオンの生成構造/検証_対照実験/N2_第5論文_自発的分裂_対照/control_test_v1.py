#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
対照テスト: コピーした第5論文コードが公開ベースラインを再現することの確認

規約: このシリーズの論文が別シリーズの結果に依拠する場合、そのプログラムを本フォルダに
コピーし、まず公開結果を再現できることを確認してから使う（改変前の基準点を固定）。

手順:
  1. 同梱の公開ベースライン spontaneous_splitting_result_v1/summary_v1.json を退避。
  2. run_spontaneous_splitting_preliminary_v1.py を新規実行（出力が上書きされる）。
  3. 新規出力を退避したベースラインと全数値フィールドで照合。
  4. 完全一致（または機械精度一致）を PASS とする。

実行: python3 control_test_v1.py
  （注意: これは同梱の公開結果を一度上書きするので、実行前にリポジトリの
   コミット済みベースラインが正本。差分ゼロなら上書き後も内容は同一。）
"""

import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.join(HERE, "nbody_spontaneous_splitting_reproduction_v1")
BASELINE = os.path.join(PKG, "spontaneous_splitting_result_v1", "summary_v1.json")
SCRIPT = os.path.join(PKG, "run_spontaneous_splitting_preliminary_v1.py")


def walk_compare(a, b, path="", out=None):
    if out is None:
        out = {"n": 0, "maxdiff": 0.0}
    if isinstance(a, dict):
        for k in a:
            walk_compare(a[k], b[k], f"{path}.{k}", out)
    elif isinstance(a, list):
        for i, (x, y) in enumerate(zip(a, b)):
            walk_compare(x, y, f"{path}[{i}]", out)
    elif isinstance(a, (int, float)):
        out["n"] += 1
        out["maxdiff"] = max(out["maxdiff"], abs(float(a) - float(b)))
    return out


def main():
    with open(BASELINE) as f:
        published = json.load(f)
    backup = os.path.join(HERE, "_published_baseline_backup.json")
    shutil.copy(BASELINE, backup)

    print("公開ベースラインを退避 →", backup)
    print("第5論文コードを新規実行中...")
    subprocess.run([sys.executable, SCRIPT], cwd=PKG, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    with open(BASELINE) as f:
        fresh = json.load(f)

    res = walk_compare(published["runs"], fresh["runs"])
    print(f"照合した数値フィールド数: {res['n']}")
    print(f"公開ベースラインとの最大絶対差: {res['maxdiff']:.3e}")
    verdict = "PASS（完全一致）" if res["maxdiff"] == 0.0 else (
        "PASS（機械精度一致）" if res["maxdiff"] < 1e-9 else "FAIL")
    print("対照テスト:", verdict)
    return 0 if res["maxdiff"] < 1e-9 else 1


if __name__ == "__main__":
    sys.exit(main())
