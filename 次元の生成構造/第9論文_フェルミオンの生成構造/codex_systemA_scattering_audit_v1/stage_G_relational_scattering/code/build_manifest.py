"""Finalize Stage G after gates, outputs, and reference integrity pass."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STAGE_ROOT / "data"
REPORT_DIR = STAGE_ROOT / "reports"
MANIFEST = STAGE_ROOT / "manifest.json"
REQUIRED_FILES = (
    "code/relational_scattering.py",
    "code/system_A_stage_G_copy.py",
    "code/run_stage_G_unit_tests.py",
    "code/run_stage_G_repeated_systemA.py",
    "code/run_stage_G_31_series.py",
    "code/analyze_stage_G.py",
    "code/test_stage_G.py",
    "code/build_manifest.py",
    "source_copy/README.md",
    "data/reference_hashes_before.json",
    "data/reference_hashes_after.json",
    "data/reference_hash_comparison.json",
    "data/stage_G_unit_test_results.csv",
    "data/stage_G_collision_results.csv",
    "data/stage_G_run_summary.csv",
    "data/stage_G_31_series_results.csv",
    "data/stage_G_correlation_results.csv",
    "data/stage_G_numerical_residuals.csv",
    "data/stage_G_summary.json",
    "data/stage_G_final_states.npz",
    "figures/Gamma_and_R_eff_by_collision.png",
    "figures/C0_C1_relational_L_exchange.png",
    "figures/C0_C1_relational_N_eff_exchange.png",
    "figures/relation_vs_localization_difference.png",
    "figures/cycle_and_return_error_comparison.png",
    "figures/return_error_31_series.png",
    "reports/00_reversed_C1_degeneracy_proof.md",
    "reports/01_relational_candidate_definition.md",
    "reports/02_unit_test_report.md",
    "reports/03_existing_systemA_comparison.md",
    "reports/04_dynamic_relation_analysis.md",
    "reports/05_31_series_comparison.md",
    "reports/06_numerical_invariants.md",
    "reports/Stage_G_report.md",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def excluded(path: Path) -> bool:
    relative = path.relative_to(STAGE_ROOT)
    return (
        path == MANIFEST
        or "__pycache__" in relative.parts
        or ".matplotlib-cache" in relative.parts
        or path.name == ".DS_Store"
        or path.suffix == ".zip"
    )


def main() -> None:
    missing = [
        relative
        for relative in REQUIRED_FILES
        if not (STAGE_ROOT / relative).is_file()
    ]
    if missing:
        raise SystemExit(f"missing Stage G files: {missing}")
    comparison = json.loads(
        (DATA_DIR / "reference_hash_comparison.json").read_text(
            encoding="utf-8"
        )
    )
    if not comparison["all_references_unchanged"]:
        raise SystemExit("Stage G reference integrity failed")
    unit = json.loads(
        (STAGE_ROOT / "logs" / "stage_G_unit_test_run.json").read_text(
            encoding="utf-8"
        )
    )
    gate = json.loads(
        (STAGE_ROOT / "logs" / "stage_G_C0_reproduction.json").read_text(
            encoding="utf-8"
        )
    )
    if unit["status"] != "pass" or not gate["passed"]:
        raise SystemExit("Stage G gate is not pass")

    summary_path = DATA_DIR / "stage_G_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["status"] = "complete"
    summary["reference_integrity"] = {
        "all_references_unchanged": True,
        "comparison_count": comparison["comparison_count"],
    }
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    report_path = REPORT_DIR / "Stage_G_report.md"
    report = report_path.read_text(encoding="utf-8")
    if "\n## 完了状態\n" not in report:
        report += f"""

## 完了状態

参照正本{comparison["comparison_count"]}ファイルのpath・size・mtime・SHA-256は事前・事後で全件一致した。

```text
Stage G 完了。
既存System A / System B原本は変更していない。
新規実装はrelational_C1一候補のみ。
Candidate 2・3の追加実装は行っていない。
N体系へは組み込んでいない。
論文本文は変更していない。
人間の承認待ち。
```
"""
    report_path.write_text(report, encoding="utf-8")

    files = []
    for path in sorted(STAGE_ROOT.rglob("*")):
        if not path.is_file() or excluded(path):
            continue
        stat = path.stat()
        files.append(
            {
                "path": path.relative_to(STAGE_ROOT).as_posix(),
                "size": stat.st_size,
                "mtime": stat.st_mtime,
                "sha256": sha256(path),
            }
        )
    manifest = {
        "stage": "G",
        "status": "complete",
        "scope": "new files under stage_G_relational_scattering only",
        "new_candidate": "relational_C1",
        "additional_candidates": [],
        "central_decision": summary["central_question"]["decision"],
        "reference_integrity": True,
        "reference_count": comparison["comparison_count"],
        "file_count": len(files),
        "total_size": sum(row["size"] for row in files),
        "files": files,
    }
    MANIFEST.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "central_decision": manifest["central_decision"],
                "reference_integrity": manifest[
                    "reference_integrity"
                ],
                "file_count": manifest["file_count"],
                "total_size": manifest["total_size"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
