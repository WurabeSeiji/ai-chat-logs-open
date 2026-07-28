"""Finalize Stage F only after reference integrity and required outputs pass."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


STAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STAGE_ROOT / "data"
REPORT_DIR = STAGE_ROOT / "reports"
MANIFEST = STAGE_ROOT / "manifest.json"

REQUIRED_FILES = (
    "code/system_A_stage_F_copy.py",
    "code/state_dependent_scattering.py",
    "code/parity_demodulation.py",
    "code/run_stage_F_reproduction_gate.py",
    "code/run_stage_F_repeated_comparison.py",
    "code/run_stage_F_R_sweep.py",
    "code/run_stage_F_31_series_check.py",
    "code/analyze_stage_F_cycles.py",
    "code/build_stage_F_manifest.py",
    "code/test_stage_F.py",
    "source_copy/README.md",
    "source_copy/run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1_ORIGINAL_SNAPSHOT.py",
    "source_copy/run_system_A_localization_exchange_R_sweep_preliminary_v1_ORIGINAL_SNAPSHOT.py",
    "data/reference_hashes_before.json",
    "data/reference_hashes_after.json",
    "data/reference_hash_comparison.json",
    "data/stage_F_reproduction_gate.csv",
    "data/stage_F_repeated_collision_results.csv",
    "data/stage_F_run_summary.csv",
    "data/stage_F_R_sweep_results.csv",
    "data/stage_F_cycle_metrics.csv",
    "data/stage_F_31_series_results.csv",
    "data/stage_F_numerical_residuals.csv",
    "data/stage_F_summary.json",
    "data/stage_F_final_states.npz",
    "figures/C0_vs_reversed_C1_L_exchange.png",
    "figures/C0_vs_reversed_C1_N_eff_exchange.png",
    "figures/parity_and_R_eff_by_collision.png",
    "figures/R_sweep_minimum_differences.png",
    "figures/cycle_period_comparison.png",
    "figures/normalization_vs_raw_update.png",
    "figures/return_error_31_series.png",
    "reports/00_scope_and_source_hashes.md",
    "reports/01_C0_reproduction_gate.md",
    "reports/02_integration_definition.md",
    "reports/03_repeated_scattering_comparison.md",
    "reports/04_R_sweep_comparison.md",
    "reports/05_cycle_and_quasistability_analysis.md",
    "reports/06_31_series_check.md",
    "reports/07_normalization_effect.md",
    "reports/08_numerical_invariants.md",
    "reports/Stage_F_report.md",
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
        or ".matplotlib" in relative.parts
        or (
            relative.parts
            and relative.parts[0] == "source_copy"
            and any(
                part.endswith("_result_v1") for part in relative.parts
            )
        )
        or path.name == ".DS_Store"
    )


def main() -> None:
    missing = [
        relative
        for relative in REQUIRED_FILES
        if not (STAGE_ROOT / relative).is_file()
    ]
    if missing:
        raise SystemExit(f"missing required Stage F files: {missing}")

    comparison = json.loads(
        (DATA_DIR / "reference_hash_comparison.json").read_text(
            encoding="utf-8"
        )
    )
    if not comparison["all_references_unchanged"]:
        raise SystemExit("reference integrity failed")

    gate = json.loads(
        (DATA_DIR / "stage_F_reproduction_gate_summary.json").read_text(
            encoding="utf-8"
        )
    )
    if gate["status"] != "pass":
        raise SystemExit("C0 reproduction gate is not pass")

    summary_path = DATA_DIR / "stage_F_summary.json"
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

    scope_path = REPORT_DIR / "00_scope_and_source_hashes.md"
    scope_text = scope_path.read_text(encoding="utf-8")
    if "\n## 事後照合\n" not in scope_text:
        scope_text += (
            "\n## 事後照合\n\n"
            f"{comparison['comparison_count']}参照すべてについて、"
            "path・size・mtime・SHA-256が事前値と一致した。"
            " `all_references_unchanged=true`。\n"
        )
    scope_path.write_text(scope_text, encoding="utf-8")

    report_path = REPORT_DIR / "Stage_F_report.md"
    report_text = report_path.read_text(encoding="utf-8")
    if "\n## 完了状態\n" not in report_text:
        report_text += """

## 完了状態

参照正本9ファイルのpath・size・mtime・SHA-256は事前・事後で全件一致した。

```text
Stage F 完了。
既存System A / System B原本は変更していない。
反転Candidate 1は独立実験コピーにのみ統合した。
Candidate 2・3は追加していない。
N体系へは組み込んでいない。
論文本文は変更していない。
人間の承認待ち。
```
"""
    report_path.write_text(report_text, encoding="utf-8")

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
        "stage": "F",
        "scope": (
            "new files under stage_F_original_systemA_integration only"
        ),
        "status": "complete",
        "C0_reproduction_gate": "pass",
        "reference_integrity": True,
        "reference_count": comparison["comparison_count"],
        "file_count": len(files),
        "total_size": sum(item["size"] for item in files),
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
