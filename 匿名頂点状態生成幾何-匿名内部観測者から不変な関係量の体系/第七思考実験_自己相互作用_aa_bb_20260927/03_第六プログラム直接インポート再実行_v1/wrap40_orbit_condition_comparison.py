#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の条件比較図プログラムを、条件と保存先だけを変えて実行する。

変えるのは CASES、PANEL_TITLES、STRICT_FILES（どれも条件の表）と、保存先の引数だけ。
図の題名とファイル名に入っている "five" は第六のプログラムの固定文字列で、変えていない。
"""
from __future__ import annotations
import shutil
import sys
from wrapper_common import HERE, case_table, load_program

OUT = HERE / "orbit_condition_comparison"


def main() -> None:
    table = case_table()
    gen = load_program("generate_paper6_orbit_condition_comparison_v1.py")
    gen.CASES = [(cid, la, lb) for (cid, n, sa, sb, la, lb, c0, d0) in table]
    gen.PANEL_TITLES = {cid: f"{cid}: lambda_A={la:+.2f}, lambda_B={lb:+.2f}"
                        for (cid, n, sa, sb, la, lb, c0, d0) in table}
    gen.STRICT_FILES = {cid: f"{cid}_figure01_orbit_overlay.svg" for (cid, *_rest) in table}

    # 第六と同じく、ケース別の重ね描き図を strict_sources/ に集める（中身は変えない）
    strict = OUT / "strict_sources"
    strict.mkdir(parents=True, exist_ok=True)
    for cid, *_rest in table:
        shutil.copyfile(HERE / "cases" / cid / "figures" / "figure01_orbit_overlay.svg",
                        strict / gen.STRICT_FILES[cid])

    sys.argv = [sys.argv[0], "--outdir", str(OUT / "figures"), "--strict-svg-dir", str(strict)]
    gen.main()


if __name__ == "__main__":
    main()
