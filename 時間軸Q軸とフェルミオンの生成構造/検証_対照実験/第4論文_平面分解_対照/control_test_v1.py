#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
対照テスト: コピーした第4論文（平面分解読出し）コードが公開ベースラインを再現することの確認

N1 §5.4 が遂及説明する「あらゆる回転平面が AB 二体単位と同型」という観測を、本シリーズ内で
局所再現するための対照テスト。手順は第5論文対照テストと同形式。

実行: python3 control_test_v1.py
"""

import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RESULT = os.path.join(HERE, "nbody_plane_decomposition_readout_result_v1",
                      "plane_decomposition_readout_result_v1.json")
SCRIPT = os.path.join(HERE, "run_nbody_plane_decomposition_readout_preliminary_v1.py")


def walk(a, b, o):
    if isinstance(a, dict):
        for k in a:
            walk(a[k], b[k], o)
    elif isinstance(a, list):
        for x, y in zip(a, b):
            walk(x, y, o)
    elif isinstance(a, (int, float)):
        o[0] += 1
        o[1] = max(o[1], abs(float(a) - float(b)))
    return o


def main():
    with open(RESULT) as f:
        published = json.load(f)
    shutil.copy(RESULT, os.path.join(HERE, "_p4_baseline_backup.json"))
    print("第4論文コードを新規実行中...")
    subprocess.run([sys.executable, SCRIPT], cwd=HERE, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(RESULT) as f:
        fresh = json.load(f)
    o = walk(published, fresh, [0, 0.0])
    print(f"照合した数値フィールド数: {o[0]}")
    print(f"公開ベースラインとの最大絶対差: {o[1]:.3e}")
    verdict = "PASS（完全一致）" if o[1] == 0.0 else (
        "PASS（機械精度一致）" if o[1] < 1e-9 else "FAIL")
    print("対照テスト:", verdict)
    os.remove(os.path.join(HERE, "_p4_baseline_backup.json"))
    return 0 if o[1] < 1e-9 else 1


if __name__ == "__main__":
    sys.exit(main())
