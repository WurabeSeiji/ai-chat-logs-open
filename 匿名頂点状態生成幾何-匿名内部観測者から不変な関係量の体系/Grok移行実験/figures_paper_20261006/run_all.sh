#!/bin/bash
# 論文用の図の再現：../exchange_engine_20261006/run_all.sh（JSON を含む全出力、約 5 分）→ make_figures.py（約 4 分）→ SHA256SUMS。
# 同時に複数の重い Python を走らせない（Accelerate が全コアを使うので、並走させると 10 倍以上遅くなる）。
set -euo pipefail
cd "$(dirname "$0")"
if [ "${SKIP_ENGINE:-0}" != "1" ]; then
  bash ../exchange_engine_20261006/run_all.sh
fi
python3 make_figures.py > make_figures.log
shasum -a 256 make_figures.py run_all.sh README.md figures_index.md make_figures.log figures/*.png data/*.json > SHA256SUMS
echo "done: $(wc -l < SHA256SUMS) files in SHA256SUMS"
