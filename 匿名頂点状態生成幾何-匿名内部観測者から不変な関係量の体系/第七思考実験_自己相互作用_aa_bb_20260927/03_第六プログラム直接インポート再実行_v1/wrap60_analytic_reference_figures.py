#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の独立解析基準の図化プログラムを、入力と保存先だけを変えて実行する。

第六のプログラムは、保存済みの解析基準 CSV を読んで 8 枚の図を描く。軌道の再計算はしない。
ここでは、ケースごとに --input と --outdir を与えてコマンドとして呼ぶ。
"""
from __future__ import annotations
import subprocess
import sys
from wrapper_common import HERE, PROGRAMS, SELF_CASES

plotter = PROGRAMS / "plot_charged_binary_analytic_reference_v1.py"
root = HERE / "analytic_reference"


def main() -> None:
    for name, *_rest in SELF_CASES:
        cmd = [sys.executable, str(plotter),
               '--input', str(root / name / 'charged_binary_analytic_reference_v1_raw.csv'),
               '--outdir', str(root / name / 'figures')]
        p = subprocess.run(cmd, capture_output=True, text=True)
        (root / name / 'plot_stdout.txt').write_text(p.stdout, encoding='utf-8')
        (root / name / 'plot_stderr.txt').write_text(p.stderr, encoding='utf-8')
        print(name, 'success' if p.returncode == 0 else 'failed', flush=True)
        if p.returncode != 0:
            raise RuntimeError(f"{name}: plot failed\n{p.stderr}")


if __name__ == "__main__":
    main()
