"""Record and compare immutable references used by Stage F."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


STAGE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[5]
DATA_DIR = STAGE_ROOT / "data"

REFERENCE_PATHS = (
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "Codex向け_StageF_反転Candidate1の既存SystemA限定統合指示.md",
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "対照実験_波束収縮_実行環境_v1/20260713/"
    "run_exchange_scattering_matrix_fermionic_localization_transfer_preliminary_v1.py",
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "対照実験_波束収縮_実行環境_v1/20260715/"
    "run_system_A_localization_exchange_R_sweep_preliminary_v1.py",
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "対照実験_波束収縮_実行環境_v1/parity_suite_v1/"
    "01_odd_kernel_N_1_2_3_5_15_63/output/"
    "odd_kernel_N_1_2_3_5_15_63_rows_v1.csv",
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "対照実験_波束収縮_実行環境_v1/parity_suite_v1/"
    "01_odd_kernel_N_1_2_3_5_15_63/output/"
    "odd_kernel_N_1_2_3_5_15_63_result_v1.json",
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "対照実験_波束収縮_実行環境_v1/parity_suite_v1/"
    "01_odd_kernel_N_1_2_3_5_15_63/output/"
    "odd_kernel_N_1_2_3_5_15_63_report_v1.md",
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "対照実験_波束収縮_実行環境_v1/20260715/instrument_check_v1/"
    "base_run/system_A_custom_packet_A1_B1-2_R0-0p697177927_C256_rows_v1.csv",
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "対照実験_波束収縮_実行環境_v1/20260715/instrument_check_v1/"
    "base_run/system_A_custom_packet_A1_B1-2_R0-0p697177927_C256_result_v1.json",
    "次元の生成構造/第9論文_フェルミオンの生成構造/"
    "対照実験_波束収縮_実行環境_v1/20260715/instrument_check_v1/"
    "base_run/system_A_custom_packet_A1_B1-2_R0-0p697177927_C256_report_v1.md",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def snapshot() -> dict:
    records = []
    for relative in REFERENCE_PATHS:
        path = REPO_ROOT / relative
        stat = path.stat()
        records.append(
            {
                "path": relative,
                "size": stat.st_size,
                "mtime": stat.st_mtime,
                "sha256": sha256(path),
            }
        )
    return {"reference_count": len(records), "references": records}


def compare(before: dict, after: dict) -> dict:
    left = {item["path"]: item for item in before["references"]}
    right = {item["path"]: item for item in after["references"]}
    paths = sorted(set(left) | set(right))
    rows = []
    for path in paths:
        before_item = left.get(path)
        after_item = right.get(path)
        unchanged = before_item == after_item
        rows.append(
            {
                "path": path,
                "unchanged": unchanged,
                "before": before_item,
                "after": after_item,
            }
        )
    return {
        "all_references_unchanged": all(row["unchanged"] for row in rows),
        "comparison_count": len(rows),
        "comparisons": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("before", "after"))
    args = parser.parse_args()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    current = snapshot()
    output = DATA_DIR / f"reference_hashes_{args.phase}.json"
    output.write_text(
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
            raise SystemExit("reference mutation detected")


if __name__ == "__main__":
    main()
