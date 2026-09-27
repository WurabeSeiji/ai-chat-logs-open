#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SHA-256 の一覧を作る。

1. ケース別の SHA256SUMS.json
   第六の postprocess_strict_5case.py の main() 末尾と同じ方法（PNG と SHA256SUMS.json 自身を除く）。
   図を作る処理がすべて終わった後に作り直す。
2. フォルダ全体の SHA256SUMS.txt
   __pycache__ と SHA256SUMS.txt 自身を除く全ファイル。
"""
from __future__ import annotations
import json
from wrapper_common import HERE, SELF_CASES, load_program


def main() -> None:
    pp5 = load_program("postprocess_strict_5case.py")
    for cid, *_rest in SELF_CASES:
        cdir = HERE / "cases" / cid
        entries = []
        for p in sorted(cdir.rglob('*')):
            if p.is_file() and p.suffix.lower() != '.png' and p.name != 'SHA256SUMS.json':
                entries.append({'path': str(p.relative_to(cdir)), 'sha256': pp5.sha256(p), 'size': p.stat().st_size})
        (cdir / 'SHA256SUMS.json').write_text(json.dumps(entries, ensure_ascii=False, indent=2) + '\n')
        print(cid, len(entries), flush=True)

    lines = []
    for p in sorted(HERE.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.name in ('SHA256SUMS.txt', '.DS_Store'):
            continue
        lines.append(f"{pp5.sha256(p)}  {p.relative_to(HERE).as_posix()}")
    (HERE / 'SHA256SUMS.txt').write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("SHA256SUMS.txt", len(lines), flush=True)


if __name__ == "__main__":
    main()
