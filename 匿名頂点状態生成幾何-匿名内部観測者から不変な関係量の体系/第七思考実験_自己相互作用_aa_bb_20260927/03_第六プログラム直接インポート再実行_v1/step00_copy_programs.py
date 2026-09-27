#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""第六思考実験と補遺のプログラムを、変更せずに paper6_programs/ へコピーする。

コピー元とコピー先の SHA-256 を paper6_programs/COPY_MANIFEST.json に記録する。
コピーしたプログラムは書き換えない。条件と保存先はラッパー側で与える。
"""
from __future__ import annotations
import hashlib, json, shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
SERIES = HERE.parents[1]          # 匿名頂点状態生成幾何-… フォルダ
DEST = HERE / "paper6_programs"

F5 = "第五思考実験_乗法状態内部展開と近似次数検証_20260924"
F6 = "第六思考実験_荷電8状態_重力クーロン複数相互作用_20260926"

SOURCES = (
    # (役割, コピー元の相対パス)
    ("状態生成器（第六 5 条件）",
     f"{F5}/13_charged8_integer_multiple_5case_strict_numeric_20260925/00_CODE_AND_RULES/run_strict_charged8_5cases.py"),
    ("ケース別の図（第六 5 条件）",
     f"{F5}/13_charged8_integer_multiple_5case_strict_numeric_20260925/00_CODE_AND_RULES/postprocess_strict_charged8_5cases.py"),
    ("解析基準との比較・ケース別の図・分析（第六 5 条件 staging）",
     f"{F5}/13_charged8_integer_multiple_5case_strict_numeric_20260925_staging_duplicate/00_CODE/postprocess_strict_5case.py"),
    ("独立解析基準の solver（第六）",
     f"{F5}/07_荷電二体系_解析近似基準解_20260925/solve_charged_binary_analytic_reference_v1.py"),
    ("独立解析基準の図（第六）",
     f"{F5}/07_荷電二体系_解析近似基準解_20260925/plot_charged_binary_analytic_reference_v1.py"),
    ("条件比較の図（第六）",
     f"{F6}/01_軌道条件比較図/generate_paper6_orbit_condition_comparison_v1.py"),
    ("対数半径の軌道図（第六 基準ケース）",
     f"{F5}/11_荷電8状態_厳格一回写像実験_20260925/01_CODE/plot_first_10_orbits_log_radial_scale.py"),
    ("初期化 G=q0=c=1（第六 補遺）",
     f"{F6}/04_補講_Gq0c1_初期化同値検証_20260927/paper6_init_G_q0_c1.py"),
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main() -> None:
    DEST.mkdir(exist_ok=True)
    records = []
    for role, rel in SOURCES:
        src = SERIES / rel
        dst = DEST / src.name
        shutil.copyfile(src, dst)
        a, b = sha256(src), sha256(dst)
        if a != b:
            raise RuntimeError(f"copy mismatch: {rel}")
        records.append({"role": role, "source": rel, "copy": f"paper6_programs/{dst.name}",
                        "bytes": dst.stat().st_size, "sha256_source": a, "sha256_copy": b, "identical": a == b})
        print(f"{b[:16]}  {dst.stat().st_size:>6} B  {dst.name}", flush=True)
    (DEST / "COPY_MANIFEST.json").write_text(json.dumps(records, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
