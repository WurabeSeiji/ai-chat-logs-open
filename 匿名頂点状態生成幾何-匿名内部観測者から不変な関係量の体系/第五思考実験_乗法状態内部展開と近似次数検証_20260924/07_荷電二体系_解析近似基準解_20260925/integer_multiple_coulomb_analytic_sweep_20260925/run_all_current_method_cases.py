#!/usr/bin/env python3
import subprocess, json, sys
from pathlib import Path

solver = Path('/mnt/data/solve_charged_binary_analytic_reference_v1.py')
outroot = Path('/mnt/data/charged_binary_integer_multiple_analytic_sweep')

cases = [
    ('n1_attractive', 1, +0.30, -0.30),
    ('n1_repulsive',  1, +0.30, +0.30),
    ('n3_attractive', 3, +0.90, -0.90),
    ('n3_repulsive',  3, +0.90, +0.90),
    ('n4_attractive', 4, +1.20, -1.20),
    ('n4_repulsive',  4, +1.20, +1.20),
]

results=[]
for name,n,la,lb in cases:
    outdir=outroot/name
    outdir.mkdir(parents=True, exist_ok=True)
    cmd=[sys.executable, str(solver), '--outdir', str(outdir), '--lambda-A', str(la), '--lambda-B', str(lb)]
    p=subprocess.run(cmd, capture_output=True, text=True)
    rec={
        'case':name,'integer_multiple_n':n,'lambda_A':la,'lambda_B':lb,
        'Q':la*lb,'Z':1-la*lb,'abs_Fc_over_Fg':abs(la*lb),
        'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr,
        'status':'success' if p.returncode==0 else 'failed_current_method'
    }
    (outdir/'run_stdout.txt').write_text(p.stdout, encoding='utf-8')
    (outdir/'run_stderr.txt').write_text(p.stderr, encoding='utf-8')
    (outdir/'case_run_record.json').write_text(json.dumps(rec, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    results.append(rec)
    print(name, rec['status'], 'Z=', rec['Z'])

(outroot/'batch_summary.json').write_text(json.dumps(results, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
