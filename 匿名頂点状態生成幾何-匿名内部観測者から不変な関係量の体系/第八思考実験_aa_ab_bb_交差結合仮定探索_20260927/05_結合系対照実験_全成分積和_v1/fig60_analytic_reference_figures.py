#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の独立解析基準の図化プログラムを、入力と保存先だけを変えて実行する（ラッパー）。

  python3 fig60_analytic_reference_figures.py <pattern_id>

描かれるのは解析基準の値であり、結合系の観測 row データではない。
"""
from __future__ import annotations
import subprocess
import sys
from fig_common import PROGRAMS, figure_dir, relation_conditions

plotter = PROGRAMS / "plot_charged_binary_analytic_reference_v1.py"


def main() -> None:
    pattern_id = sys.argv[1]
    root = figure_dir(pattern_id) / "analytic_reference"
    for rel, name, *_r in relation_conditions(pattern_id):
        cmd = [sys.executable, str(plotter),
               '--input', str(root / name / 'charged_binary_analytic_reference_v1_raw.csv'),
               '--outdir', str(root / name / 'figures')]
        p = subprocess.run(cmd, capture_output=True, text=True)
        (root / name / 'plot_stdout.txt').write_text(p.stdout, encoding='utf-8')
        (root / name / 'plot_stderr.txt').write_text(p.stderr, encoding='utf-8')
        print(name, 'success' if p.returncode == 0 else 'failed', flush=True)
        if p.returncode != 0:
            raise RuntimeError(p.stderr)


if __name__ == "__main__":
    main()
