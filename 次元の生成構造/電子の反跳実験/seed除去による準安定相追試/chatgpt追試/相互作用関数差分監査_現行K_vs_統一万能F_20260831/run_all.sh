#!/bin/bash
set -e; set -o pipefail
cd "$(dirname "$0")"
python3 program/audit_checks.py 2>&1 | grep -v Warning | tee results/audit.log
find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
echo done
