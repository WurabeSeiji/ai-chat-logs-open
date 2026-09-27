#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六の後処理プログラム 2 本を、条件と保存先だけを変えて実行する。

1. postprocess_strict_5case.py
     解析基準との比較（reference_comparison.json）、ケース別の図、分析 md、SHA256SUMS.json
2. postprocess_strict_charged8_5cases.py（第六の論文 11.3 節に載っている図化プログラム）
     ケース別の図 6 枚。reference_comparison.json を読むので、1. の後に実行する。

2 本は同じ figures/ に書く。名前が同じ 4 枚は 2. の出力で上書きされる。
最後に、1. と同じ方法で SHA256SUMS.json を作り直す。
"""
from __future__ import annotations
import json
from wrapper_common import HERE, SELF_CASES, load_program

IDS = [c[0] for c in SELF_CASES]


def main() -> None:
    pp5 = load_program("postprocess_strict_5case.py")
    pp5.BASE = HERE
    pp5.CASES = HERE / "cases"
    pp5.ANA = HERE / "analytic_reference"
    pp5.IDS = list(IDS)
    pp5.main()
    print("postprocess_strict_5case.py done", flush=True)

    ppc = load_program("postprocess_strict_charged8_5cases.py")
    ppc.BASE = HERE
    ppc.CASES = HERE / "cases"
    ppc.ANABASE = HERE / "analytic_reference"
    ppc.CASE_IDS = list(IDS)
    ppc.main()
    print("postprocess_strict_charged8_5cases.py done", flush=True)

    # postprocess_strict_5case.py の main() 末尾と同じ処理（図が増えたので作り直す）
    for cid in IDS:
        cdir = pp5.CASES / cid
        entries = []
        for p in sorted(cdir.rglob('*')):
            if p.is_file() and p.suffix.lower() != '.png' and p.name != 'SHA256SUMS.json':
                entries.append({'path': str(p.relative_to(cdir)), 'sha256': pp5.sha256(p), 'size': p.stat().st_size})
        (cdir / 'SHA256SUMS.json').write_text(json.dumps(entries, ensure_ascii=False, indent=2) + '\n')
    print("SHA256SUMS.json refreshed", flush=True)


if __name__ == "__main__":
    main()
