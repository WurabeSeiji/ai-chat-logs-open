#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の独立解析基準 solver を、条件と保存先だけを変えて実行する。

第六の run_all_current_method_cases.py と同じ手順（solver をコマンドとして呼び、
標準出力・標準エラー・実行記録を保存する）で、solver のパス、保存先、条件だけが違う。
"""
from __future__ import annotations
import json
import subprocess
import sys
from wrapper_common import HERE, PROGRAMS, case_table

solver = PROGRAMS / "solve_charged_binary_analytic_reference_v1.py"
outroot = HERE / "analytic_reference"


def main() -> None:
    results = []
    for name, n, sa, sb, la, lb, c0, d0 in case_table():
        outdir = outroot / name
        outdir.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(solver), '--outdir', str(outdir), '--lambda-A', str(la), '--lambda-B', str(lb)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        rec = {
            'case': name, 'integer_multiple_n': n, 'lambda_A': la, 'lambda_B': lb,
            'Q': la * lb, 'Z': 1 - la * lb, 'abs_Fc_over_Fg': abs(la * lb),
            'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr,
            'status': 'success' if p.returncode == 0 else 'failed_current_method'
        }
        (outdir / 'run_stdout.txt').write_text(p.stdout, encoding='utf-8')
        (outdir / 'run_stderr.txt').write_text(p.stderr, encoding='utf-8')
        (outdir / 'case_run_record.json').write_text(json.dumps(rec, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        results.append(rec)
        print(name, rec['status'], 'Z=', rec['Z'], flush=True)
    (outroot / 'batch_summary.json').write_text(json.dumps(results, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


if __name__ == "__main__":
    main()
