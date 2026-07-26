#!/usr/bin/env python3
"""Run only the locked N=5 Paper 7 reproduction into an isolated output tree."""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import json
import os
import platform
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path


sys.dont_write_bytecode = True

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[2]
CONFIG_PATH = PACKAGE_ROOT / "config_locked.json"
EXPECTED_PATH = PACKAGE_ROOT / "expected_hashes.json"
VERIFY_PATH = PACKAGE_ROOT / "comparison" / "source_verification.json"
REPRODUCED = PACKAGE_ROOT / "reproduced"
LOG_DIR = PACKAGE_ROOT / "logs"
RUN_LOG = LOG_DIR / "run_reproduction.log"
MANIFEST_PATH = LOG_DIR / "run_manifest.json"


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_source(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"import specを作成できない: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def assert_n5(n: int) -> None:
    if type(n) is not int or n != 5:
        raise RuntimeError(f"N=5以外の実行を拒否: {n!r}")


def verify_gate() -> dict:
    if not VERIFY_PATH.is_file():
        raise RuntimeError("verify_sources.py の成功記録がない")
    verification = json.loads(VERIFY_PATH.read_text(encoding="utf-8"))
    if verification.get("success") is not True or verification.get("locked_n") != 5:
        raise RuntimeError("verify_sources.py が成功していない")
    for item in verification["checks"]:
        path = Path(item["path"])
        if not path.is_file() or sha256(path) != item["expected_sha256"]:
            raise RuntimeError(f"検証後に原本または比較対象が変化した: {path}")
    return verification


def ensure_clean_output() -> None:
    if REPRODUCED.exists() and any(path.is_file() for path in REPRODUCED.rglob("*")):
        raise RuntimeError("reproduced/ が空でないため、上書きを避けて停止")
    REPRODUCED.mkdir(parents=True, exist_ok=True)


def record_environment() -> dict:
    import numpy

    try:
        import matplotlib

        matplotlib_version = matplotlib.__version__
    except Exception as exc:
        matplotlib_version = f"IMPORT_ERROR: {type(exc).__name__}: {exc}"
    return {
        "python_executable": sys.executable,
        "python_version": sys.version,
        "python_implementation": platform.python_implementation(),
        "numpy_version": numpy.__version__,
        "matplotlib_version": matplotlib_version,
        "platform": platform.platform(),
        "os": {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
        },
    }


def required_outputs(v2_out: Path, p7_out: Path, metastable_out: Path) -> list[Path]:
    return [
        metastable_out / "fcurve_N00005_delta1e-15_seed0.csv",
        metastable_out / "summary_N00005_delta1e-15_seed0.json",
        metastable_out / "metastable_N00005_delta1e-15_seed0.json",
        v2_out / "raw" / "N00005_dimension_saturation_v2" / "q_svd_N00005.csv",
        v2_out / "diagnostics" / "N00005_saturation.json",
        p7_out / "raw" / "N00005" / "paper7_long_timeseries.csv",
        p7_out / "raw" / "N00005" / "transverse_stability_timeseries.csv",
        p7_out / "summary" / "N00005_5color_meta.json",
        p7_out / "summary" / "N00005_transverse_meta.json",
        p7_out / "figures" / "figure1_N00005.png",
        p7_out / "figures" / "figure2_N00005_5color.png",
        p7_out / "figures" / "figure3_N00005_5color.png",
        p7_out / "figures" / "transverse_growth_N00005.png",
    ]


def main() -> int:
    started = time.perf_counter()
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    expected = json.loads(EXPECTED_PATH.read_text(encoding="utf-8"))
    assert_n5(config["locked_n"])
    if config.get("allowed_n") != [5]:
        raise RuntimeError("allowed_n が [5] ではない")
    verification = verify_gate()
    ensure_clean_output()

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    mpl_dir = LOG_DIR / "matplotlib"
    cache_dir = LOG_DIR / "cache"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_dir)
    os.environ["XDG_CACHE_HOME"] = str(cache_dir)

    engine = REPO_ROOT / "時間軸Q軸とフェルミオンの生成構造/検証_対照実験/第5論文原本_自発的分裂予備実験_v1"
    v2_code = engine / "exact_lowN_eigenspectrum_v2/code"
    p7_code = engine / "exact_lowN_eigenspectrum_v2/paper7_longtime/code"
    metastable_out = REPRODUCED / "metastable_series_result_v1"
    v2_out = REPRODUCED / "exact_lowN_eigenspectrum_v2"
    p7_out = v2_out / "paper7_longtime"
    (v2_out / "diagnostics").mkdir(parents=True, exist_ok=True)

    for path in (str(p7_code), str(v2_code), str(engine)):
        if path not in sys.path:
            sys.path.insert(0, path)

    timings = []
    manifest = {
        "stage": "A0",
        "locked_n": 5,
        "success": False,
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": "PYTHONDONTWRITEBYTECODE=1 python3 " + str(Path(__file__).resolve()),
        "environment": record_environment(),
        "source_verification": str(VERIFY_PATH),
        "source_hashes": {
            name: item["sha256"] for name, item in expected["sources"].items()
        },
        "timings": timings,
        "failure": None,
    }

    def timed(label, function):
        t0 = time.perf_counter()
        print(f"[START] {label}", flush=True)
        result = function()
        elapsed = time.perf_counter() - t0
        timings.append({"step": label, "seconds": elapsed})
        print(f"[DONE] {label}: {elapsed:.6f} s", flush=True)
        return result

    with RUN_LOG.open("w", encoding="utf-8") as log_handle:
        tee_out = Tee(sys.stdout, log_handle)
        tee_err = Tee(sys.stderr, log_handle)
        try:
            with contextlib.redirect_stdout(tee_out), contextlib.redirect_stderr(tee_err):
                meta_mod = load_source(
                    "run_metastable_series_v1",
                    engine / "run_metastable_series_v1.py",
                )
                meta_mod.RESULT_DIR = metastable_out
                meta = timed("run_metastable_series_v1.py:run_one(5)", lambda: meta_mod.run_one(5))
                if meta.get("base_summary", {}).get("n") != 5:
                    raise RuntimeError("metastable系列の返値がN=5ではない")

                sat_mod = load_source(
                    "run_n300_dimension_saturation_v2",
                    v2_code / "run_n300_dimension_saturation_v2.py",
                )
                sat_mod.BASE = v2_out
                sat = timed("run_n300_dimension_saturation_v2.py:run(5)", lambda: sat_mod.run(5))
                if sat.get("N") != 5:
                    raise RuntimeError("q/rank系列の返値がN=5ではない")

                p5_mod = load_source(
                    "run_paper7_5color_timeseries",
                    p7_code / "run_paper7_5color_timeseries.py",
                )
                p5_mod.P7 = p7_out
                p5 = timed("run_paper7_5color_timeseries.py:run(5)", lambda: p5_mod.run(5))
                if p5.get("N") != 5:
                    raise RuntimeError("5色系列の返値がN=5ではない")

                trans_mod = load_source(
                    "run_paper7_transverse",
                    p7_code / "run_paper7_transverse.py",
                )
                trans_mod.P7 = p7_out
                trans = timed("run_paper7_transverse.py:run(5)", lambda: trans_mod.run(5))
                if trans.get("N") != 5:
                    raise RuntimeError("横摂動系列の返値がN=5ではない")

                fig_mod = load_source(
                    "make_paper7_figures",
                    p7_code / "make_paper7_figures.py",
                )
                fig_mod.P7 = p7_out
                fig_mod.FIG = p7_out / "figures"
                fig_mod.FIG.mkdir(parents=True, exist_ok=True)
                fig_mod.NS = [5]
                timed("make_paper7_figures.py:N=5 functions", fig_mod.fig1)
                timed("make_paper7_figures.py:N=5 fig23", fig_mod.fig23)
                timed("make_paper7_figures.py:N=5 fig_transverse", fig_mod.fig_transverse)
                timed("make_paper7_figures.py:N=5 fig_lambda_vs_N", fig_mod.fig_lambda_vs_N)

                sat_fig_mod = load_source(
                    "make_saturation_comparison",
                    v2_code / "make_saturation_comparison.py",
                )
                sat_fig_mod.BASE = v2_out
                sat_fig_mod.NS = [5]
                timed("make_saturation_comparison.py:main() with NS=[5]", sat_fig_mod.main)

                missing = [str(path) for path in required_outputs(v2_out, p7_out, metastable_out) if not path.is_file()]
                if missing:
                    raise RuntimeError("必須再生成ファイル欠落: " + ", ".join(missing))

                for item in verification["checks"]:
                    path = Path(item["path"])
                    if not path.is_file() or sha256(path) != item["expected_sha256"]:
                        raise RuntimeError(f"実行中に原本または既存成果物が変化した: {path}")

                artifacts = sorted(
                    str(path.relative_to(PACKAGE_ROOT))
                    for path in REPRODUCED.rglob("*")
                    if path.is_file()
                )
                manifest.update(
                    {
                        "success": True,
                        "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                        "duration_seconds": time.perf_counter() - started,
                        "reproduced_files": artifacts,
                        "N_values_executed": [5],
                        "runtime_NS_overrides": {
                            "make_paper7_figures.py": [5],
                            "make_saturation_comparison.py": [5],
                        },
                    }
                )
                print("[SUCCESS] Stage A0 N=5 reproduction outputs generated", flush=True)
        except Exception as exc:
            manifest.update(
                {
                    "success": False,
                    "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                    "duration_seconds": time.perf_counter() - started,
                    "failure": {
                        "type": type(exc).__name__,
                        "message": str(exc),
                        "traceback": traceback.format_exc(),
                    },
                }
            )
            print(f"STOP: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
            MANIFEST_PATH.write_text(
                json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            return 1

    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
