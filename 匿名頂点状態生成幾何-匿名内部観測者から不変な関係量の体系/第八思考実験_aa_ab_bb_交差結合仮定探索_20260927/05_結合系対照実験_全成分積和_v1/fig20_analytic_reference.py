#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の独立解析基準 solver を、条件と保存先だけを変えて実行する（ラッパー）。

  python3 fig20_analytic_reference.py <pattern_id>
"""
from __future__ import annotations
import json
import subprocess
import sys
from fig_common import PROGRAMS, figure_dir, relation_conditions

solver = PROGRAMS / "solve_charged_binary_analytic_reference_v1.py"


def main() -> None:
    pattern_id = sys.argv[1]
    outroot = figure_dir(pattern_id) / "analytic_reference"
    results = []
    for rel, name, la, lb, c0, d0 in relation_conditions(pattern_id):
        outdir = outroot / name
        outdir.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(solver), '--outdir', str(outdir), '--lambda-A', str(la), '--lambda-B', str(lb)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        rec = {'case': name, 'lambda_A': la, 'lambda_B': lb, 'Z': 1 - la * lb,
               'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr,
               'status': 'success' if p.returncode == 0 else 'failed_current_method'}
        (outdir / 'run_stdout.txt').write_text(p.stdout, encoding='utf-8')
        (outdir / 'run_stderr.txt').write_text(p.stderr, encoding='utf-8')
        (outdir / 'case_run_record.json').write_text(json.dumps(rec, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        results.append(rec)
        print(name, rec['status'], 'Z=', rec['Z'], flush=True)
        if p.returncode != 0:
            raise RuntimeError(p.stderr)
    (outroot / 'batch_summary.json').write_text(json.dumps(results, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


if __name__ == "__main__":
    main()
