#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""実行パターンを与えて、結合系の対照実験を実行するラッパー。

  python3 wrap_run_pattern.py patterns/<パターン>.json

実行パターンは JSON ファイルで与える。このラッパーは、パターンを読んで実行プログラムに渡すだけ。
保存先は results/<pattern_id>/ 。
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import run_unified_control


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: python3 wrap_run_pattern.py patterns/<pattern>.json")
    pattern = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    run_unified_control.run(pattern)


if __name__ == "__main__":
    main()
