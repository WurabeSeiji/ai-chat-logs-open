#!/usr/bin/env python3
"""論文7の現行二段階seed条件を変えずに図4を隔離再生成し、既存成果物と照合する。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ENGINE = (
    REPO
    / "時間軸Q軸とフェルミオンの生成構造"
    / "検証_対照実験"
    / "第5論文原本_自発的分裂予備実験_v1"
)
V2 = ENGINE / "exact_lowN_eigenspectrum_v2"
P7 = V2 / "paper7_longtime"
P7_CODE = P7 / "code"
OUTPUT = HERE / "outputs"
COMPARISON = HERE / "comparison"
NS = (5, 40, 300)

# 既存の論文7成果物に記録された実行経路をそのまま固定する。
RUN_CONFIG = {
    5: {
        "runner": "transverse",
        "eps": [1e-8, 1e-10, 1e-12, 1e-14],
        "rec": 50,
    },
    40: {
        "runner": "transverse_cached",
        "eps": [1e-8, 1e-10, 1e-12, 1e-14],
        "rec": 50,
    },
    300: {
        "runner": "transverse_cached",
        "eps": [1e-8, 1e-10, 1e-12],
        "rec": 100,
    },
}

SOURCE_FILES = {
    "engine": ENGINE / "run_n_scaling_lowrank_v1.py",
    "plane_exact": ENGINE / "run_plane_flow_exact_v1.py",
    "plane_approx": ENGINE / "run_plane_flow_approx_v1.py",
    "saturation": V2 / "code" / "run_n300_dimension_saturation_v2.py",
    "five_color": P7_CODE / "run_paper7_5color_timeseries.py",
    "transverse": P7_CODE / "run_paper7_transverse.py",
    "transverse_cached": P7_CODE / "run_paper7_transverse_cached.py",
    "figures": P7_CODE / "make_paper7_figures.py",
}

EXPECTED_SOURCE_HASHES = {
    "engine": "ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d",
    "plane_exact": "9cf28ca8c0d2ad8fac2f0f6dae045248695247c5809c21ccb2069ef91a94ab76",
    "plane_approx": "a9d247a8070d849fe989e35e00320e470968354614424c9bdceda57132d9f0fa",
    "saturation": "229938a66631057426f187ed80b17de08cfcb9107cfe509c30f5bbdcca3a03e6",
    "five_color": "fe5c7cbc33437890a5f50944cbbae1594e5f647739d4955402b578f515658503",
    "transverse": "ac1073bea329971de3ff4c2fd1588d926029a8502c21e8cc01f406acb86ad60b",
    "transverse_cached": "897a20d0e4d6a2cf23e8502e02755bc1f64a727bf9fcb2e849801dc90f469dbe",
    "figures": "273fec25c2e4c30c4561a7506ef373829e374d8b0d2eeb813bb7f4e1c06a8000",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_sources() -> dict:
    checks = {}
    for name, path in SOURCE_FILES.items():
        actual = sha256(path) if path.is_file() else None
        expected = EXPECTED_SOURCE_HASHES[name]
        checks[name] = {
            "path": str(path),
            "expected_sha256": expected,
            "actual_sha256": actual,
            "matched": actual == expected,
        }
    result = {
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "all_matched": all(item["matched"] for item in checks.values()),
        "checks": checks,
    }
    HERE.mkdir(parents=True, exist_ok=True)
    (HERE / "source_verification.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    if not result["all_matched"]:
        raise SystemExit("SOURCE_MISMATCH: 論文7原本のSHA-256が固定値と一致しない")
    print("[VERIFIED] 論文7原本8ファイルのSHA-256一致", flush=True)
    return result


def load_source(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare_import_paths() -> None:
    for path in (ENGINE, V2 / "code", P7_CODE):
        text = str(path)
        if text not in sys.path:
            sys.path.insert(0, text)


def run_one(n: int) -> dict:
    if n not in NS:
        raise SystemExit(f"N must be one of {NS}")
    verify_sources()
    target_csv = OUTPUT / "raw" / f"N{n:05d}" / "transverse_stability_timeseries.csv"
    target_meta = OUTPUT / "summary" / f"N{n:05d}_transverse_meta.json"
    if target_csv.exists() or target_meta.exists():
        raise SystemExit(f"既存対照出力を上書きしない: N={n}")

    prepare_import_paths()
    config = RUN_CONFIG[n]
    source = SOURCE_FILES[config["runner"]]
    module = load_source(f"paper7_figure4_control_N{n}", source)
    module.P7 = OUTPUT

    started = time.perf_counter()
    if config["runner"] == "transverse":
        summary = module.run(n)
    else:
        summary = module.run(n, list(config["eps"]), int(config["rec"]))
    elapsed = time.perf_counter() - started

    manifest = {
        "N": n,
        "M": n * (n - 1) // 2,
        "initial_prng_seed": 40260722 + 1000 * n,
        "initial_seed_amplitude": 1e-15,
        "transverse_prng_seed": 70000 + n,
        "transverse_direction_count": 3,
        "transverse_eps": config["eps"],
        "benettin_interval": 500,
        "record_interval": config["rec"],
        "runner": config["runner"],
        "elapsed_seconds": elapsed,
        "summary": summary,
        "output_csv": str(target_csv),
        "output_meta": str(target_meta),
        "output_csv_sha256": sha256(target_csv),
        "output_meta_sha256": sha256(target_meta),
        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_dir = HERE / "logs"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    (manifest_dir / f"run_N{n:05d}.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"[DONE] N={n} runner={config['runner']} t0={summary['t0']} "
        f"lambda_max={summary['lambda_max_for_N']:.12e} elapsed={elapsed:.3f}s",
        flush=True,
    )
    return manifest


def make_figures() -> None:
    verify_sources()
    missing = [
        OUTPUT / "raw" / f"N{n:05d}" / "transverse_stability_timeseries.csv"
        for n in NS
        if not (
            OUTPUT / "raw" / f"N{n:05d}" / "transverse_stability_timeseries.csv"
        ).is_file()
    ]
    if missing:
        raise SystemExit("横安定時系列が未生成: " + ", ".join(str(path) for path in missing))

    prepare_import_paths()
    module = load_source("paper7_figures_figure4_control", SOURCE_FILES["figures"])
    module.P7 = OUTPUT
    module.FIG = OUTPUT / "figures"
    module.FIG.mkdir(parents=True, exist_ok=True)
    module.NS = list(NS)
    module.fig_transverse()
    print("[DONE] 原本fig_transverse()で図4を再生成", flush=True)


def compare_csv(reference: Path, reproduced: Path) -> dict:
    with reference.open(encoding="utf-8") as handle:
        ref_rows = list(csv.DictReader(handle))
    with reproduced.open(encoding="utf-8") as handle:
        new_rows = list(csv.DictReader(handle))
    ref_header = list(ref_rows[0].keys()) if ref_rows else []
    new_header = list(new_rows[0].keys()) if new_rows else []
    max_abs = 0.0
    numeric_equal = len(ref_rows) == len(new_rows) and ref_header == new_header
    first_difference = None
    if numeric_equal:
        for row_index, (ref_row, new_row) in enumerate(zip(ref_rows, new_rows)):
            for key in ref_header:
                try:
                    ref_value = float(ref_row[key])
                    new_value = float(new_row[key])
                    if np.isnan(ref_value) and np.isnan(new_value):
                        difference = 0.0
                    else:
                        difference = abs(ref_value - new_value)
                except ValueError:
                    difference = 0.0 if ref_row[key] == new_row[key] else float("inf")
                max_abs = max(max_abs, difference)
                if difference != 0.0 and first_difference is None:
                    first_difference = {
                        "row": row_index,
                        "column": key,
                        "reference": ref_row[key],
                        "reproduced": new_row[key],
                    }
        numeric_equal = max_abs == 0.0
    return {
        "reference": str(reference),
        "reproduced": str(reproduced),
        "reference_sha256": sha256(reference),
        "reproduced_sha256": sha256(reproduced),
        "byte_equal": reference.read_bytes() == reproduced.read_bytes(),
        "row_count_reference": len(ref_rows),
        "row_count_reproduced": len(new_rows),
        "headers_equal": ref_header == new_header,
        "numeric_equal": numeric_equal,
        "max_abs_difference": max_abs,
        "first_difference": first_difference,
    }


def compare_json(reference: Path, reproduced: Path) -> dict:
    ref_data = json.loads(reference.read_text(encoding="utf-8"))
    new_data = json.loads(reproduced.read_text(encoding="utf-8"))
    return {
        "reference": str(reference),
        "reproduced": str(reproduced),
        "reference_sha256": sha256(reference),
        "reproduced_sha256": sha256(reproduced),
        "byte_equal": reference.read_bytes() == reproduced.read_bytes(),
        "semantic_equal": ref_data == new_data,
    }


def compare_png(reference: Path, reproduced: Path) -> dict:
    import matplotlib.image as mpimg

    ref_image = mpimg.imread(reference)
    new_image = mpimg.imread(reproduced)
    shape_equal = ref_image.shape == new_image.shape
    pixel_equal = shape_equal and np.array_equal(ref_image, new_image)
    max_pixel_difference = (
        float(np.max(np.abs(ref_image.astype(float) - new_image.astype(float))))
        if shape_equal
        else None
    )
    return {
        "reference": str(reference),
        "reproduced": str(reproduced),
        "reference_sha256": sha256(reference),
        "reproduced_sha256": sha256(reproduced),
        "byte_equal": reference.read_bytes() == reproduced.read_bytes(),
        "shape_reference": list(ref_image.shape),
        "shape_reproduced": list(new_image.shape),
        "pixel_equal": bool(pixel_equal),
        "max_pixel_difference": max_pixel_difference,
    }


def normalized_svg(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"<dc:date>.*?</dc:date>", "<dc:date>NORMALIZED</dc:date>", text)
    generated_ids = []
    for match in re.finditer(r'(?:id="|xlink:href="#)([mp][0-9a-f]{8,})(?=")', text):
        value = match.group(1)
        if value not in generated_ids:
            generated_ids.append(value)
    for index, value in enumerate(generated_ids):
        text = text.replace(value, f"GENERATED_ID_{index:04d}")
    return text


def compare_svg(reference: Path, reproduced: Path) -> dict:
    ref_normalized = normalized_svg(reference)
    new_normalized = normalized_svg(reproduced)
    return {
        "reference": str(reference),
        "reproduced": str(reproduced),
        "reference_sha256": sha256(reference),
        "reproduced_sha256": sha256(reproduced),
        "byte_equal": reference.read_bytes() == reproduced.read_bytes(),
        "normalized_text_equal": ref_normalized == new_normalized,
    }


def compare_outputs() -> dict:
    verify_sources()
    results: dict[str, dict] = {"csv": {}, "json": {}, "png": {}, "svg": {}}
    for n in NS:
        key = f"N{n:05d}"
        results["csv"][key] = compare_csv(
            P7 / "raw" / key / "transverse_stability_timeseries.csv",
            OUTPUT / "raw" / key / "transverse_stability_timeseries.csv",
        )
        results["json"][key] = compare_json(
            P7 / "summary" / f"{key}_transverse_meta.json",
            OUTPUT / "summary" / f"{key}_transverse_meta.json",
        )
        results["png"][key] = compare_png(
            P7 / "figures" / f"transverse_growth_{key}.png",
            OUTPUT / "figures" / f"transverse_growth_{key}.png",
        )
        results["svg"][key] = compare_svg(
            P7 / "figures" / f"transverse_growth_{key}.svg",
            OUTPUT / "figures" / f"transverse_growth_{key}.svg",
        )

    compare_key = "compare_N5_N40_N300"
    results["png"][compare_key] = compare_png(
        P7 / "figures" / "transverse_growth_compare_N5_N40_N300.png",
        OUTPUT / "figures" / "transverse_growth_compare_N5_N40_N300.png",
    )
    results["svg"][compare_key] = compare_svg(
        P7 / "figures" / "transverse_growth_compare_N5_N40_N300.svg",
        OUTPUT / "figures" / "transverse_growth_compare_N5_N40_N300.svg",
    )

    checks = []
    checks.extend(item["numeric_equal"] for item in results["csv"].values())
    checks.extend(item["semantic_equal"] for item in results["json"].values())
    checks.extend(item["pixel_equal"] for item in results["png"].values())
    checks.extend(item["normalized_text_equal"] for item in results["svg"].values())
    report = {
        "compared_at_utc": datetime.now(timezone.utc).isoformat(),
        "N_values": list(NS),
        "conditions": {
            "initial_prng_seeds": {str(n): 40260722 + 1000 * n for n in NS},
            "initial_seed_amplitude": 1e-15,
            "transverse_prng_seeds": {str(n): 70000 + n for n in NS},
            "transverse_direction_count": 3,
            "benettin_interval": 500,
            "run_config": RUN_CONFIG,
            "xmax": 55000,
        },
        "all_passed": all(checks),
        "results": results,
    }
    COMPARISON.mkdir(parents=True, exist_ok=True)
    report_path = COMPARISON / "control_comparison.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"[COMPARE] all_passed={report['all_passed']} report={report_path}", flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("verify")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("N", type=int)
    subparsers.add_parser("figures")
    subparsers.add_parser("compare")
    args = parser.parse_args()

    if args.command == "verify":
        verify_sources()
    elif args.command == "run":
        run_one(args.N)
    elif args.command == "figures":
        make_figures()
    elif args.command == "compare":
        compare_outputs()


if __name__ == "__main__":
    main()
