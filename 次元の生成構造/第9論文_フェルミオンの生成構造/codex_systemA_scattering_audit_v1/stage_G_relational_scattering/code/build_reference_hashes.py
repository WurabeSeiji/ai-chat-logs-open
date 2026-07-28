"""Record immutable Stage G references before and after execution."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


STAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[5]
DATA_DIR = STAGE_ROOT / "data"
F_ROOT = (
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "codex_systemA_scattering_audit_v1/"
    "stage_F_original_systemA_integration/"
)
REFERENCE_PATHS = (
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "Codex向け_StageG_関係依存散乱の実装指示.md",
    F_ROOT + "code/parity_demodulation.py",
    F_ROOT + "code/state_dependent_scattering.py",
    F_ROOT + "code/system_A_stage_F_copy.py",
    F_ROOT + "code/run_stage_F_reproduction_gate.py",
    F_ROOT + "code/run_stage_F_31_series_check.py",
    F_ROOT + "data/stage_F_reproduction_gate_summary.json",
    F_ROOT + "data/stage_F_repeated_collision_results.csv",
    F_ROOT + "data/stage_F_31_series_results.csv",
    F_ROOT + "data/stage_F_summary.json",
    F_ROOT + "reports/Stage_F_report.md",
    F_ROOT + "manifest.json",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def snapshot() -> dict:
    rows = []
    for relative in REFERENCE_PATHS:
        path = REPO_ROOT / relative
        stat = path.stat()
        rows.append(
            {
                "path": relative,
                "size": stat.st_size,
                "mtime": stat.st_mtime,
                "sha256": sha256(path),
            }
        )
    return {"reference_count": len(rows), "references": rows}


def compare(before: dict, after: dict) -> dict:
    left = {row["path"]: row for row in before["references"]}
    right = {row["path"]: row for row in after["references"]}
    comparisons = []
    for path in sorted(set(left) | set(right)):
        unchanged = left.get(path) == right.get(path)
        comparisons.append(
            {
                "path": path,
                "unchanged": unchanged,
                "before": left.get(path),
                "after": right.get(path),
            }
        )
    return {
        "all_references_unchanged": all(
            row["unchanged"] for row in comparisons
        ),
        "comparison_count": len(comparisons),
        "comparisons": comparisons,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("before", "after"))
    args = parser.parse_args()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    current = snapshot()
    (DATA_DIR / f"reference_hashes_{args.phase}.json").write_text(
        json.dumps(current, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    if args.phase == "after":
        before = json.loads(
            (DATA_DIR / "reference_hashes_before.json").read_text(
                encoding="utf-8"
            )
        )
        result = compare(before, current)
        (DATA_DIR / "reference_hash_comparison.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        if not result["all_references_unchanged"]:
            raise SystemExit("Stage G reference mutation detected")


if __name__ == "__main__":
    main()
