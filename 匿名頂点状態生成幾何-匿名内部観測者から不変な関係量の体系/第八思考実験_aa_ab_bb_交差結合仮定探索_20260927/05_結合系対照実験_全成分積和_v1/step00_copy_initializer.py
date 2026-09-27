#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六思考実験 補遺の初期化プログラムを、変更せずにコピーする。

第六・第七のプログラムのうち、この実験で使うのは初期化だけである。
状態の更新は unified_engine.py に書き直したものを使う。
"""
from __future__ import annotations
import hashlib, json, shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
SERIES = HERE.parents[1]
DEST = HERE / "paper6_supplement_program"
SOURCE = ("第六思考実験_荷電8状態_重力クーロン複数相互作用_20260926/"
          "04_補講_Gq0c1_初期化同値検証_20260927/paper6_init_G_q0_c1.py")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    DEST.mkdir(exist_ok=True)
    src = SERIES / SOURCE
    dst = DEST / src.name
    shutil.copyfile(src, dst)
    a, b = sha256(src), sha256(dst)
    if a != b:
        raise RuntimeError("copy mismatch")
    (DEST / "COPY_MANIFEST.json").write_text(json.dumps([{
        "role": "初期化 G=q0=c=1（第六 補遺）", "source": SOURCE, "copy": f"paper6_supplement_program/{dst.name}",
        "bytes": dst.stat().st_size, "sha256_source": a, "sha256_copy": b, "identical": a == b,
    }], ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"{b[:16]}  {dst.stat().st_size} B  {dst.name}", flush=True)


if __name__ == "__main__":
    main()
