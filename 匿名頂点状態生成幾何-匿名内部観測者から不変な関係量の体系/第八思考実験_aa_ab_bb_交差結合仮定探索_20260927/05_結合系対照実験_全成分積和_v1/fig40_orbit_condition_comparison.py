#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の条件比較図プログラムを、条件と保存先だけを変えて実行する（ラッパー）。

  python3 fig40_orbit_condition_comparison.py <pattern_id>

変えるのは CASES、PANEL_TITLES、STRICT_FILES（条件の表）と、保存先の引数だけ。
図の題名とファイル名の "five" は第六のプログラムの固定文字列で、変えていない。
"""
from __future__ import annotations
import shutil
import sys
from fig_common import figure_dir, load_program, relation_conditions


def main() -> None:
    pattern_id = sys.argv[1]
    base = figure_dir(pattern_id)
    out = base / "orbit_condition_comparison"
    table = relation_conditions(pattern_id)
    gen = load_program("generate_paper6_orbit_condition_comparison_v1.py")
    gen.CASES = [(cid, la, lb) for (rel, cid, la, lb, c0, d0) in table]
    gen.PANEL_TITLES = {cid: f"{cid}: lambda_A={la:+.2f}, lambda_B={lb:+.2f}" for (rel, cid, la, lb, c0, d0) in table}
    gen.STRICT_FILES = {cid: f"{cid}_figure01_orbit_overlay.svg" for (rel, cid, *_r) in table}
    strict = out / "strict_sources"
    strict.mkdir(parents=True, exist_ok=True)
    for rel, cid, *_r in table:
        shutil.copyfile(base / cid / "figure01_orbit_overlay.svg", strict / gen.STRICT_FILES[cid])
    sys.argv = [sys.argv[0], "--outdir", str(out / "figures"), "--strict-svg-dir", str(strict)]
    gen.main()


if __name__ == "__main__":
    main()
