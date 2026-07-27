#!/usr/bin/env python3
"""無seed自然軌道を一度だけ実行し、図3個別占有列と図4用自然観測量を同時記録する。"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PAPER8 = HERE.parent
ENGINE = (
    REPO
    / "時間軸Q軸とフェルミオンの生成構造"
    / "検証_対照実験"
    / "第5論文原本_自発的分裂予備実験_v1"
)
V2 = ENGINE / "exact_lowN_eigenspectrum_v2"
P7 = V2 / "paper7_longtime"
P7_CODE = P7 / "code"
PAPER8_CODE = PAPER8 / "code"
OUTPUT = HERE / "outputs"
NS = (5, 40, 300)
XMAX = 55000
GUARD = 3000
SIG_REL = 1e-6
Q_REL_TAU = 1e-8
SAMPLE = {5: 25, 40: 25, 300: 100}
FMT = "%.10e"

SOURCE_FILES = {
    "engine": ENGINE / "run_n_scaling_lowrank_v1.py",
    "plane_exact": ENGINE / "run_plane_flow_exact_v1.py",
    "plane_approx": ENGINE / "run_plane_flow_approx_v1.py",
    "saturation": V2 / "code" / "run_n300_dimension_saturation_v2.py",
    "five_color": P7_CODE / "run_paper7_5color_timeseries.py",
    "figures": P7_CODE / "make_paper7_figures.py",
    "preliminary_seed_ablation": PAPER8_CODE / "run_preliminary_seed_ablation_v1.py",
}

EXPECTED_SOURCE_HASHES = {
    "engine": "ba0fc19b03caf06d16e97e2cd5da499ed2ed8e95288f6dbf2062bedcaf11176d",
    "plane_exact": "9cf28ca8c0d2ad8fac2f0f6dae045248695247c5809c21ccb2069ef91a94ab76",
    "plane_approx": "a9d247a8070d849fe989e35e00320e470968354614424c9bdceda57132d9f0fa",
    "saturation": "229938a66631057426f187ed80b17de08cfcb9107cfe509c30f5bbdcca3a03e6",
    "five_color": "fe5c7cbc33437890a5f50944cbbae1594e5f647739d4955402b578f515658503",
    "figures": "273fec25c2e4c30c4561a7506ef373829e374d8b0d2eeb813bb7f4e1c06a8000",
    "preliminary_seed_ablation": "75a10a5b951302bef2ba77cf363066fe8cbdce9a7f1a0af70ccd2df4e520b1d8",
}

BASE_COLUMNS = [
    "step",
    "time",
    "N",
    "condition",
    "initial_seed_enabled",
    "metastable_seed_enabled",
    "initial_seed_amplitude",
    "metastable_seed_amplitude",
    "parent_plane_occupation",
    "f_outside_parent",
    "q1",
    "q2",
    "q3",
    "q4",
    "rank_Q",
    "dominant_plane_occupation",
    "non_dominant_occupation",
    "kernel_occupation",
    "residual_occupation",
    "norm_Z",
    "dagger_norm_error",
    "zero_square_real",
    "zero_square_imag",
    "zero_square_abs",
    "projection_closure_error",
    "crossing_detected",
    "metastable_start_detected",
]

ADDED_COLUMNS = [
    "direction_1_occupation",
    "direction_2_occupation",
    "direction_3_occupation",
    "direction_4_occupation",
    "other_rotating_occupation",
    "occupation_sum",
    "plane_1_occupation",
    "plane_2_occupation",
    "splitting_fraction",
    "crossing_flag",
    "norm_error",
    "conservation_error",
]


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
        checks[name] = {
            "path": str(path),
            "expected_sha256": EXPECTED_SOURCE_HASHES[name],
            "actual_sha256": actual,
            "matched": actual == EXPECTED_SOURCE_HASHES[name],
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
        raise SystemExit("SOURCE_MISMATCH: 再利用元のSHA-256が固定値と一致しない")
    print("[VERIFIED] 再利用元7ファイルを確認", flush=True)
    return result


def prepare_import_paths() -> None:
    for path in (ENGINE, V2 / "code", P7_CODE, PAPER8_CODE):
        text = str(path)
        if text not in sys.path:
            sys.path.insert(0, text)


def load_source(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def qsv4(B0: np.ndarray, Bdom: np.ndarray) -> np.ndarray:
    Q4 = np.column_stack([B0, Bdom])
    eigenvalues = np.clip(np.linalg.eigvalsh(Q4.T @ Q4)[::-1], 0, None)
    return np.sqrt(eigenvalues)


def build_seedless(n: int):
    """条件Aと同じ乱数消費順。kernel seed gを生成せず、Z0=vとする。"""
    prepare_import_paths()
    from run_n_scaling_lowrank_v1 import LowRankSystem, make_parent
    from run_n300_dimension_saturation_v2 import dominant_plane, gram_reduce
    from run_plane_flow_approx_v1 import parent_plane_split_approx
    from run_plane_flow_exact_v1 import parent_plane_split_exact

    sys_lr = LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000 * n)
    v, residual, sig = make_parent(sys_lr, rng, iters=1200, tol=1e-12)
    if n <= 40:
        _, B_p1, B_rot, _ = parent_plane_split_exact(sys_lr, v)
        method = "exact"
    else:
        _, B_p1, B_rot, _, _ = parent_plane_split_approx(sys_lr, v, SIG_REL)
        method = "approx"
    _, B0, _, _, _ = dominant_plane(sys_lr, gram_reduce(sys_lr, v))
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    Z0 = v.copy()
    wp = rng.normal(size=sys_lr.m)
    return sys_lr, v, B_p1, B_rot, B0, p, q, Z0, wp, residual, sig, method


def run_one(n: int) -> dict:
    if n not in NS:
        raise SystemExit(f"N must be one of {NS}")
    verify_sources()
    target_csv = OUTPUT / "raw" / f"N{n:05d}" / "paper7_long_timeseries.csv"
    target_meta = OUTPUT / "summary" / f"N{n:05d}_5color_meta.json"
    if target_csv.exists() or target_meta.exists():
        raise SystemExit(f"既存の無seed出力を上書きしない: N={n}")

    prepare_import_paths()
    from run_n300_dimension_saturation_v2 import dominant_plane, gram_reduce
    from run_paper7_5color_timeseries import align_2d, occ, s4_new_dirs

    (
        sys_lr,
        v,
        B_p1,
        B_rot,
        B0,
        p,
        q,
        Z,
        wp,
        parent_residual,
        parent_sigma,
        method,
    ) = build_seedless(n)
    started = time.perf_counter()
    sample_every = SAMPLE[n]
    crossing = None
    metastable_start = None
    f_prev = None
    max_occupation_closure_error = 0.0
    max_norm_error = 0.0
    max_zero_square_abs = 0.0
    target_csv.parent.mkdir(parents=True, exist_ok=True)

    with target_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BASE_COLUMNS + ADDED_COLUMNS)
        writer.writeheader()
        t = 0
        while True:
            Z_perp = Z - p * (p @ Z) - q * (q @ Z)
            f_step = float(np.real(np.conj(Z_perp) @ Z_perp)) / float(
                np.real(np.conj(Z) @ Z)
            )
            if crossing is None and f_step > 0.05:
                crossing = t
                metastable_start = crossing + GUARD

            if t % sample_every == 0 or t == XMAX:
                totZ = float(np.real(np.conj(Z) @ Z))
                E_P1 = occ(B_p1, Z)
                E_other = occ(B_rot, Z)
                E_ker = totZ - E_P1 - E_other
                f = 1.0 - E_P1 / totZ

                gr = gram_reduce(sys_lr, Z)
                _, Bdom, _, _, _ = dominant_plane(sys_lr, gr)
                E_dom = occ(Bdom, Z)
                qs = qsv4(B0, Bdom)
                rankQ = int(np.sum(qs > Q_REL_TAU * qs[0]))
                ztz = complex(Z @ Z)
                closure = abs(totZ - E_P1 - E_other - E_ker)

                e34 = s4_new_dirs(B0, Bdom)
                projected = B_rot @ (B_rot.T @ e34)
                fq, _ = np.linalg.qr(projected)
                f34 = fq[:, :2] if fq.shape[1] >= 2 else fq
                f34 = align_2d(f_prev, f34)
                f_prev = f34
                E_d3 = occ(f34[:, [0]], Z)
                E_d4 = occ(f34[:, [1]], Z) if f34.shape[1] > 1 else 0.0
                E_remaining = max(0.0, E_other - E_d3 - E_d4)
                E_a1 = occ(B_p1[:, [0]], Z)
                E_a2 = occ(B_p1[:, [1]], Z)
                occupation_sum = (
                    E_P1 + E_d3 + E_d4 + E_remaining + E_ker
                ) / totZ
                occupation_closure_error = abs(occupation_sum - 1.0)
                norm_error = abs(totZ - 1.0)

                row = {
                    "step": t,
                    "time": t,
                    "N": n,
                    "condition": "A",
                    "initial_seed_enabled": 0,
                    "metastable_seed_enabled": 0,
                    "initial_seed_amplitude": FMT % 0.0,
                    "metastable_seed_amplitude": FMT % 0.0,
                    "parent_plane_occupation": FMT % (E_P1 / totZ),
                    "f_outside_parent": FMT % f,
                    "q1": FMT % qs[0],
                    "q2": FMT % qs[1],
                    "q3": FMT % qs[2],
                    "q4": FMT % qs[3],
                    "rank_Q": rankQ,
                    "dominant_plane_occupation": FMT % (E_dom / totZ),
                    "non_dominant_occupation": FMT % (E_other / totZ),
                    "kernel_occupation": FMT % (E_ker / totZ),
                    "residual_occupation": FMT % (closure / totZ),
                    "norm_Z": FMT % np.sqrt(totZ),
                    "dagger_norm_error": FMT % norm_error,
                    "zero_square_real": FMT % ztz.real,
                    "zero_square_imag": FMT % ztz.imag,
                    "zero_square_abs": FMT % abs(ztz),
                    "projection_closure_error": FMT % (closure / totZ),
                    "crossing_detected": int(crossing is not None),
                    "metastable_start_detected": int(
                        metastable_start is not None and t >= metastable_start
                    ),
                    "direction_1_occupation": FMT % (E_a1 / totZ),
                    "direction_2_occupation": FMT % (E_a2 / totZ),
                    "direction_3_occupation": FMT % (E_d3 / totZ),
                    "direction_4_occupation": FMT % (E_d4 / totZ),
                    "other_rotating_occupation": FMT % (E_remaining / totZ),
                    "occupation_sum": FMT % occupation_sum,
                    "plane_1_occupation": FMT % (E_P1 / totZ),
                    "plane_2_occupation": FMT % ((E_d3 + E_d4) / totZ),
                    "splitting_fraction": FMT % f,
                    "crossing_flag": int(crossing is not None),
                    "norm_error": FMT % norm_error,
                    "conservation_error": FMT % norm_error,
                }
                writer.writerow(row)
                max_occupation_closure_error = max(
                    max_occupation_closure_error, occupation_closure_error
                )
                max_norm_error = max(max_norm_error, norm_error)
                max_zero_square_abs = max(max_zero_square_abs, abs(ztz))

            if t >= XMAX:
                break
            sys_lr.set_theta(np.angle(Z))
            sigma_estimate, wp = sys_lr.sigma_max_power(wp)
            Z = sys_lr.cayley_step(Z, sigma_estimate)
            t += 1

    elapsed = time.perf_counter() - started
    summary = {
        "N": n,
        "M": int(sys_lr.m),
        "condition": "A_seedless_natural",
        "initial_seed_enabled": False,
        "metastable_seed_enabled": False,
        "initial_seed_amplitude": 0.0,
        "metastable_seed_amplitude": 0.0,
        "benettin_enabled": False,
        "state_feedback_from_observation": False,
        "initial_state_rule": "Z0 = v.copy(); kernel seed is not generated",
        "parent_prng_seed": 40260722 + 1000 * n,
        "crossing": crossing,
        "metastable_start": metastable_start,
        "xmax": XMAX,
        "sample_every": sample_every,
        "method_parent_basis": method,
        "parent_residual": float(parent_residual),
        "parent_sigma": [float(value) for value in parent_sigma],
        "max_occupation_closure_error": max_occupation_closure_error,
        "max_norm_error": max_norm_error,
        "max_zero_square_abs": max_zero_square_abs,
        "elapsed_seconds": elapsed,
        "csv_sha256": sha256(target_csv),
    }
    target_meta.parent.mkdir(parents=True, exist_ok=True)
    target_meta.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    logs = HERE / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    (logs / f"run_N{n:05d}.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"[DONE] N={n} seedless crossing={crossing} "
        f"metastable_start={metastable_start} elapsed={elapsed:.3f}s",
        flush=True,
    )
    return summary


def compare_existing_condition_a(n: int) -> dict:
    reproduced = OUTPUT / "raw" / f"N{n:05d}" / "paper7_long_timeseries.csv"
    reference = PAPER8 / "raw" / f"N{n:05d}" / "condition_A_no_seed.csv"
    if not reproduced.is_file() or not reference.is_file():
        raise SystemExit("比較対象CSVが不足")
    with reference.open(encoding="utf-8") as handle:
        ref_rows = list(csv.DictReader(handle))
    with reproduced.open(encoding="utf-8") as handle:
        new_rows = list(csv.DictReader(handle))
    first_difference = None
    equal = len(ref_rows) == len(new_rows)
    if equal:
        for row_index, (ref_row, new_row) in enumerate(zip(ref_rows, new_rows)):
            for column in BASE_COLUMNS:
                if ref_row[column] != new_row[column]:
                    first_difference = {
                        "row": row_index,
                        "column": column,
                        "reference": ref_row[column],
                        "reproduced": new_row[column],
                    }
                    equal = False
                    break
            if not equal:
                break

    max_added_closure_error = 0.0
    for row in new_rows:
        components = (
            float(row["direction_1_occupation"])
            + float(row["direction_2_occupation"])
            + float(row["direction_3_occupation"])
            + float(row["direction_4_occupation"])
            + float(row["other_rotating_occupation"])
            + float(row["kernel_occupation"])
        )
        max_added_closure_error = max(max_added_closure_error, abs(components - 1.0))
    result = {
        "N": n,
        "reference": str(reference),
        "reproduced": str(reproduced),
        "base_column_count": len(BASE_COLUMNS),
        "added_column_count": len(ADDED_COLUMNS),
        "row_count_reference": len(ref_rows),
        "row_count_reproduced": len(new_rows),
        "base_columns_string_equal": equal,
        "first_difference": first_difference,
        "max_added_component_sum_error_from_csv": max_added_closure_error,
    }
    comparison_dir = HERE / "comparison"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    path = comparison_dir / f"compare_existing_A_N{n:05d}.json"
    path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"[COMPARE] N={n} base_columns_equal={equal} "
        f"max_added_sum_error={max_added_closure_error:.3e}",
        flush=True,
    )
    return result


def make_figures() -> None:
    verify_sources()
    available = [
        n
        for n in NS
        if (
            OUTPUT / "raw" / f"N{n:05d}" / "paper7_long_timeseries.csv"
        ).is_file()
        and (OUTPUT / "summary" / f"N{n:05d}_5color_meta.json").is_file()
    ]
    if not available:
        raise SystemExit("描画可能な無seed時系列がない")

    prepare_import_paths()
    original_figures = load_source("paper7_seedless_original_figures", SOURCE_FILES["figures"])
    original_figures.P7 = OUTPUT
    original_figures.FIG = OUTPUT / "figures"
    original_figures.FIG.mkdir(parents=True, exist_ok=True)
    original_figures.NS = list(available)
    original_figures.fig23()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def load_rows(n: int):
        path = OUTPUT / "raw" / f"N{n:05d}" / "paper7_long_timeseries.csv"
        with path.open(encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    fig, axes = plt.subplots(
        len(available),
        1,
        figsize=(11, 3.4 * len(available)),
        sharex=True,
        squeeze=False,
    )
    for ax, n in zip(axes[:, 0], available):
        rows = load_rows(n)
        t = np.array([float(row["time"]) for row in rows])
        f = np.array([float(row["f_outside_parent"]) for row in rows])
        q3 = np.array([float(row["q3"]) for row in rows])
        q4 = np.array([float(row["q4"]) for row in rows])
        meta = json.loads(
            (OUTPUT / "summary" / f"N{n:05d}_5color_meta.json").read_text(
                encoding="utf-8"
            )
        )
        ax.plot(t, f, color="black", lw=0.9, label="f")
        ax.plot(t, q3, color="#F58518", lw=0.8, label="q3")
        ax.plot(t, q4, color="#E45756", lw=0.8, label="q4")
        ax.axvline(meta["crossing"], color="black", ls=":", lw=0.7)
        ax.axvline(meta["metastable_start"], color="#54A24B", ls="--", lw=0.7)
        ax.set_xlim(0, XMAX)
        ax.set_xticks(np.arange(0, XMAX + 1, 5000))
        ax.set_ylabel(f"N={n}")
    axes[0, 0].legend(fontsize=7, loc="upper right")
    axes[-1, 0].set_xlabel("step (absolute)")
    fig.suptitle(
        "Seedless natural metastable observables "
        "(no perturbation, no Benettin feedback)"
    )
    fig.tight_layout()
    figure_dir = OUTPUT / "figures"
    fig.savefig(figure_dir / "figure4_seedless_natural_f_q3_q4_compare.png", dpi=130)
    fig.savefig(figure_dir / "figure4_seedless_natural_f_q3_q4_compare.svg")
    plt.close(fig)

    # 振幅が小さい波形の形だけを調べるmin-max表示。
    # N=300未再実行時は、共通27列が同じ既存条件A CSVを読み取り専用で使う。
    passive_available = [
        n
        for n in NS
        if (
            OUTPUT / "raw" / f"N{n:05d}" / "paper7_long_timeseries.csv"
        ).is_file()
        or (PAPER8 / "raw" / f"N{n:05d}" / "condition_A_no_seed.csv").is_file()
    ]

    def load_passive_rows(n: int):
        current = OUTPUT / "raw" / f"N{n:05d}" / "paper7_long_timeseries.csv"
        fallback = PAPER8 / "raw" / f"N{n:05d}" / "condition_A_no_seed.csv"
        source = current if current.is_file() else fallback
        with source.open(encoding="utf-8") as handle:
            return list(csv.DictReader(handle)), source

    def passive_times(n: int):
        current = OUTPUT / "summary" / f"N{n:05d}_5color_meta.json"
        if current.is_file():
            meta = json.loads(current.read_text(encoding="utf-8"))
            return int(meta["crossing"]), int(meta["metastable_start"])
        fallback = PAPER8 / "diagnostics" / f"N{n:05d}_condition_A.json"
        meta = json.loads(fallback.read_text(encoding="utf-8"))
        return int(meta["crossing_step"]), int(meta["metastable_start_step"])

    def minmax(values: np.ndarray):
        low = float(np.min(values))
        high = float(np.max(values))
        amplitude = high - low
        if amplitude == 0.0:
            return np.zeros_like(values), low, high, amplitude
        return (values - low) / amplitude, low, high, amplitude

    range_report = {
        "normalization": "(x - min_window) / (max_window - min_window)",
        "warning": "min-max normalization shows shape; monotonic relaxation is not oscillation",
        "windows": {},
    }
    for window_name, fixed_start in (
        ("metastable", None),
        ("late_50000_55000", 50000),
    ):
        fig, axes = plt.subplots(
            len(passive_available),
            1,
            figsize=(11, 3.4 * len(passive_available)),
            sharex=True,
            squeeze=False,
        )
        range_report["windows"][window_name] = {}
        for ax, n in zip(axes[:, 0], passive_available):
            rows, source = load_passive_rows(n)
            crossing, metastable_start = passive_times(n)
            t_all = np.array([float(row["time"]) for row in rows])
            start = metastable_start if fixed_start is None else fixed_start
            mask = t_all >= start
            t = t_all[mask]
            report_n = {
                "source_csv": str(source),
                "window_start": int(start),
                "window_end": XMAX,
                "crossing": crossing,
                "metastable_start": metastable_start,
                "observables": {},
            }
            for key, label, color in (
                ("f_outside_parent", "f", "black"),
                ("q3", "q3", "#F58518"),
                ("q4", "q4", "#E45756"),
            ):
                values = np.array([float(row[key]) for row in rows])[mask]
                scaled, low, high, amplitude = minmax(values)
                differences = np.diff(values)
                nonzero = differences[differences != 0.0]
                signs = np.sign(nonzero)
                turning_count = (
                    int(np.sum(signs[1:] != signs[:-1])) if len(signs) > 1 else 0
                )
                ax.plot(t, scaled, color=color, lw=0.8, label=label)
                report_n["observables"][key] = {
                    "min": low,
                    "max": high,
                    "peak_to_peak": amplitude,
                    "turning_count": turning_count,
                }
            amplitudes = report_n["observables"]
            ax.text(
                0.01,
                0.04,
                (
                    f"Δf={amplitudes['f_outside_parent']['peak_to_peak']:.3e}, "
                    f"Δq3={amplitudes['q3']['peak_to_peak']:.3e}, "
                    f"Δq4={amplitudes['q4']['peak_to_peak']:.3e}"
                ),
                transform=ax.transAxes,
                fontsize=7,
                va="bottom",
            )
            ax.set_ylim(-0.05, 1.05)
            ax.set_ylabel(f"N={n}\nmin-max")
            range_report["windows"][window_name][str(n)] = report_n
        axes[0, 0].legend(fontsize=7, loc="upper right")
        axes[-1, 0].set_xlabel("step (absolute)")
        axes[-1, 0].set_xlim(
            min(item["window_start"] for item in range_report["windows"][window_name].values()),
            XMAX,
        )
        fig.suptitle(
            "Seedless natural observables — per-series min-max normalization "
            f"({window_name})"
        )
        fig.tight_layout()
        tag = f"figure4_seedless_natural_minmax_{window_name}_compare"
        fig.savefig(figure_dir / f"{tag}.png", dpi=130)
        fig.savefig(figure_dir / f"{tag}.svg")
        plt.close(fig)

    summary_dir = OUTPUT / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    (summary_dir / "figure4_seedless_natural_minmax_ranges.json").write_text(
        json.dumps(range_report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    # N=40,300の漸近曲線を中心に、残差だけを20倍した表示。
    from scipy.optimize import curve_fit

    asymptotic_start = 20000
    residual_amplification = 20.0
    asymptotic_ns = [n for n in (40, 300) if n in passive_available]
    fit_report = {
        "fit_window": [asymptotic_start, XMAX],
        "model": "c + a * exp(-(t - fit_start) / tau)",
        "display_coordinate": (
            "amplification * (data - asymptotic_curve); "
            "zero line represents the fitted asymptotic curve"
        ),
        "amplification": residual_amplification,
        "N": {},
    }

    def exponential_asymptote(u, c, a, tau):
        return c + a * np.exp(-u / tau)

    fig, axes = plt.subplots(
        len(asymptotic_ns),
        3,
        figsize=(15, 4.2 * len(asymptotic_ns)),
        sharex=True,
        squeeze=False,
    )
    observables = (
        ("f_outside_parent", "f", "black"),
        ("q3", "q3", "#F58518"),
        ("q4", "q4", "#E45756"),
    )
    for row_index, n in enumerate(asymptotic_ns):
        rows, source = load_passive_rows(n)
        t_all = np.array([float(row["time"]) for row in rows])
        mask = t_all >= asymptotic_start
        t = t_all[mask]
        u = t - asymptotic_start
        fit_report["N"][str(n)] = {"source_csv": str(source), "observables": {}}
        for column_index, (key, label, color) in enumerate(observables):
            ax = axes[row_index, column_index]
            values = np.array([float(row[key]) for row in rows])[mask]
            initial = (values[-1], values[0] - values[-1], 5000.0)
            parameters, _ = curve_fit(
                exponential_asymptote,
                u,
                values,
                p0=initial,
                bounds=([-np.inf, -np.inf, 100.0], [np.inf, np.inf, 1e6]),
                maxfev=100000,
            )
            trend = exponential_asymptote(u, *parameters)
            residual = values - trend
            displayed_residual = residual_amplification * residual
            residual_ss = float(np.sum(residual**2))
            total_ss = float(np.sum((values - np.mean(values)) ** 2))
            r_squared = 1.0 - residual_ss / total_ss if total_ss else 1.0
            differences = np.diff(residual)
            nonzero = differences[differences != 0.0]
            signs = np.sign(nonzero)
            turning_count = (
                int(np.sum(signs[1:] != signs[:-1])) if len(signs) > 1 else 0
            )
            residual_peak_to_peak = float(np.ptp(residual))
            fit_report["N"][str(n)]["observables"][key] = {
                "c": float(parameters[0]),
                "a": float(parameters[1]),
                "tau": float(parameters[2]),
                "r_squared": r_squared,
                "residual_std": float(np.std(residual)),
                "residual_peak_to_peak": residual_peak_to_peak,
                "residual_turning_count": turning_count,
                "displayed_peak_to_peak_about_trend": (
                    residual_amplification * residual_peak_to_peak
                ),
            }
            ax.axhline(
                0.0,
                color="#2F2F2F",
                ls="--",
                lw=1.0,
                label="fitted asymptote (=0)",
            )
            ax.plot(
                t,
                residual,
                color="#B8B8B8",
                lw=0.7,
                label="actual residual (×1)",
            )
            ax.plot(
                t,
                displayed_residual,
                color=color,
                lw=0.9,
                label=f"residual ×{residual_amplification:.0f}",
            )
            ax.set_title(f"N={n}  {label}")
            ax.text(
                0.02,
                0.04,
                (
                    f"τ={parameters[2]:.1f}, R²={r_squared:.7f}\n"
                    f"residual Δ={residual_peak_to_peak:.3e}, "
                    f"turns={turning_count}"
                ),
                transform=ax.transAxes,
                fontsize=8,
                va="bottom",
            )
            ax.grid(alpha=0.15)
            if column_index == 0:
                ax.set_ylabel(
                    f"{residual_amplification:.0f} × "
                    "(data − fitted asymptote)"
                )
    axes[0, 0].legend(fontsize=7, loc="best")
    for ax in axes[-1, :]:
        ax.set_xlabel("step (absolute)")
        ax.set_xlim(asymptotic_start, XMAX)
    fig.suptitle(
        "Seedless natural trajectory: residual around fitted asymptote "
        f"(display ×{residual_amplification:.0f}; dashed line = asymptote)"
    )
    fig.tight_layout()
    asymptotic_tag = "figure4_seedless_asymptotic_centered_x20_N40_N300"
    fig.savefig(figure_dir / f"{asymptotic_tag}.png", dpi=150)
    fig.savefig(figure_dir / f"{asymptotic_tag}.svg")
    plt.close(fig)
    (summary_dir / "figure4_seedless_asymptotic_centered_x20_fit.json").write_text(
        json.dumps(fit_report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"[FIGURES] 図3 N={available}; 無介入自然図4 N={passive_available} を生成",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("verify")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("N", type=int)
    compare_parser = subparsers.add_parser("compare-existing")
    compare_parser.add_argument("N", type=int)
    subparsers.add_parser("figures")
    args = parser.parse_args()

    if args.command == "verify":
        verify_sources()
    elif args.command == "run":
        run_one(args.N)
    elif args.command == "compare-existing":
        compare_existing_condition_a(args.N)
    elif args.command == "figures":
        make_figures()


if __name__ == "__main__":
    main()
