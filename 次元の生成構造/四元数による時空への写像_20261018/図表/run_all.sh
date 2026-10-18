#!/bin/bash
# 図と表を作り直し、SHA256SUMS.txt を更新する
set -euo pipefail
cd "$(dirname "$0")"
python3 make_figures.py
python3 verify_tables.py
python3 - <<'PY'
import hashlib, pathlib
root = pathlib.Path('.')
files = sorted(p for p in root.rglob('*') if p.is_file() and p.name != 'SHA256SUMS.txt' and '__pycache__' not in p.parts)
lines = [f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.as_posix()}" for p in files]
(root / 'SHA256SUMS.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')
print('SHA256SUMS.txt', len(lines), 'files')
PY
