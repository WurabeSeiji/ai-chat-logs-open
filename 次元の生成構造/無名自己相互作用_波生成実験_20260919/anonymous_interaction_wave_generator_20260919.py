#!/usr/bin/env python3
"""
Anonymous complex-wave generation experiment (2026-09-19)

Physical/model input: one initial complex wave z1 = a + i b.

Generation rule (applied literally; no algebraic simplification in the update):
    T_ij = (conj(z_i) * z_j) * z_j
for every currently existing ordered pair (i, j), INCLUDING i == j.
The next wave is the negative of the full left-hand sum:
    z_new = - sum_i sum_j T_ij
Old waves are kept; exactly one new wave is appended per generation.

No Delta-tau, no adjacency matrix, no particle labels, no exclusion of self-interaction.
The max generation count and overflow guard are numerical/computational limits only,
not physical parameters of the rule.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Tuple

import matplotlib.pyplot as plt


@dataclass
class StepRecord:
    generation: int
    wave_count: int
    new_re: float
    new_im: float
    new_abs: float
    new_phase: float
    left_sum_re: float
    left_sum_im: float
    ordered_pair_terms: int
    finite: bool


def literal_left_sum(waves: List[complex]) -> complex:
    """Literal O(N^2) ordered-pair sum. Deliberately NOT algebraically simplified."""
    total = 0.0 + 0.0j
    n = len(waves)
    for i in range(n):
        zi_conj = waves[i].conjugate()
        for j in range(n):
            zj = waves[j]
            term = (zi_conj * zj) * zj
            total += term
    return total


def run_generation(a: float, b: float, max_generations: int = 200,
                   stop_abs: float = 1e120) -> Tuple[List[complex], List[StepRecord], str]:
    if max_generations < 1:
        raise ValueError("max_generations must be >= 1")

    waves: List[complex] = [complex(a, b)]
    records: List[StepRecord] = [
        StepRecord(
            generation=1,
            wave_count=1,
            new_re=a,
            new_im=b,
            new_abs=abs(waves[0]),
            new_phase=math.atan2(b, a),
            left_sum_re=float('nan'),
            left_sum_im=float('nan'),
            ordered_pair_terms=0,
            finite=math.isfinite(a) and math.isfinite(b),
        )
    ]

    stop_reason = "max_generations"
    while len(waves) < max_generations:
        n = len(waves)
        left_sum = literal_left_sum(waves)
        z_new = -left_sum
        finite = math.isfinite(z_new.real) and math.isfinite(z_new.imag)
        mag = abs(z_new) if finite else float('inf')
        phase = math.atan2(z_new.imag, z_new.real) if finite else float('nan')

        records.append(
            StepRecord(
                generation=n + 1,
                wave_count=n + 1,
                new_re=z_new.real,
                new_im=z_new.imag,
                new_abs=mag,
                new_phase=phase,
                left_sum_re=left_sum.real,
                left_sum_im=left_sum.imag,
                ordered_pair_terms=n * n,
                finite=finite,
            )
        )

        if not finite:
            stop_reason = "non_finite"
            break

        waves.append(z_new)

        if mag > stop_abs:
            stop_reason = f"abs>{stop_abs:g}"
            break

    return waves, records, stop_reason


def write_csv(path: str, records: List[StepRecord]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fields = list(asdict(records[0]).keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in records:
            w.writerow(asdict(r))


def safe_case_name(a: float, b: float) -> str:
    def f(x: float) -> str:
        s = f"{x:.12g}"
        return s.replace("-", "m").replace("+", "p").replace(".", "d")
    return f"a_{f(a)}__b_{f(b)}"


def plot_single(outdir: str, a: float, b: float, records: List[StepRecord]) -> None:
    gens = [r.generation for r in records]
    mags = [r.new_abs for r in records]
    plot_cap = 1e40
    res = [r.new_re for r in records]
    ims = [r.new_im for r in records]
    phases = [r.new_phase for r in records]

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    lin_pairs = [(g,m) for g,m in zip(gens,mags) if math.isfinite(m) and m <= plot_cap]
    if lin_pairs:
        ax.plot([p[0] for p in lin_pairs], [p[1] for p in lin_pairs], marker="o", markersize=2, linewidth=1)
    ax.set_xlabel("generation / wave index")
    ax.set_ylabel("|z_new|")
    ax.set_title(f"Generated-wave amplitude: z1={a:g}+i({b:g})")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "amplitude_linear.png"), dpi=180)
    plt.close(fig)

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    nonzero = [(g, m) for g, m in zip(gens, mags) if m > 0 and math.isfinite(m) and m <= plot_cap]
    if nonzero:
        ax.plot([x[0] for x in nonzero], [x[1] for x in nonzero], marker="o", markersize=2, linewidth=1)
        ax.set_yscale("log")
    ax.set_xlabel("generation / wave index")
    ax.set_ylabel("|z_new| (log scale)")
    ax.set_title(f"Generated-wave amplitude (log): z1={a:g}+i({b:g})")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "amplitude_log.png"), dpi=180)
    plt.close(fig)

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.plot(res, ims, marker="o", markersize=3, linewidth=1)
    ax.axhline(0, linewidth=0.7)
    ax.axvline(0, linewidth=0.7)
    ax.set_xlabel("Re(z_new)")
    ax.set_ylabel("Im(z_new)")
    ax.set_title("Generated waves in complex plane")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "complex_plane.png"), dpi=180)
    plt.close(fig)

    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)
    ax.plot(gens, phases, marker="o", markersize=2, linewidth=1)
    ax.set_xlabel("generation / wave index")
    ax.set_ylabel("arg(z_new) [rad]")
    ax.set_title("Generated-wave phase")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "phase.png"), dpi=180)
    plt.close(fig)


def run_case(base_outdir: str, a: float, b: float, max_generations: int, stop_abs: float) -> Dict[str, Any]:
    case = safe_case_name(a, b)
    outdir = os.path.join(base_outdir, case)
    os.makedirs(outdir, exist_ok=True)
    waves, records, stop_reason = run_generation(a, b, max_generations=max_generations, stop_abs=stop_abs)
    write_csv(os.path.join(outdir, "generated_waves.csv"), records)
    plot_single(outdir, a, b, records)

    nz = [r.new_abs for r in records if r.new_abs > 0 and math.isfinite(r.new_abs)]
    summary = {
        "a": a,
        "b": b,
        "initial_abs": abs(complex(a, b)),
        "initial_phase": math.atan2(b, a),
        "record_count": len(records),
        "last_wave_abs": records[-1].new_abs,
        "max_finite_abs": max(nz) if nz else 0.0,
        "stop_reason": stop_reason,
        "rule": "z_new = -sum_i sum_j ((conj(z_i)*z_j)*z_j), including i=j; append one wave; keep all prior waves",
        "literal_double_loop": True,
    }
    with open(os.path.join(outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    return {"case": case, **summary, "records": records}


def run_suite(base_outdir: str, max_generations: int, stop_abs: float) -> None:
    s2 = math.sqrt(2.0)
    cases = [
        (1.0, 0.0, "unit_real"),
        (0.0, 1.0, "unit_imag"),
        (1.0 / s2, 0.0, "subcritical_real_1_over_sqrt2"),
        (0.0, 1.0 / s2, "subcritical_imag_i_over_sqrt2"),
        (s2, 0.0, "supercritical_real_sqrt2"),
        (0.0, s2, "supercritical_imag_i_sqrt2"),
        (1.001, 0.0, "nearcritical_real_1p001"),
        (0.999, 0.0, "nearcritical_real_0p999"),
        (1.0 / s2, 1.0 / s2, "unit_phase45"),
        (1.001 / s2, 1.001 / s2, "nearcritical_phase45_1p001"),
        (0.999 / s2, 0.999 / s2, "nearcritical_phase45_0p999"),
    ]

    suite_rows = []
    all_results = []
    for a, b, label in cases:
        result = run_case(base_outdir, a, b, max_generations, stop_abs)
        result["label"] = label
        all_results.append(result)
        suite_rows.append({k: v for k, v in result.items() if k not in ("records",)})

    with open(os.path.join(base_outdir, "suite_summary.csv"), "w", newline="", encoding="utf-8") as f:
        fields = ["label", "case", "a", "b", "initial_abs", "initial_phase", "record_count", "last_wave_abs", "max_finite_abs", "stop_reason", "literal_double_loop", "rule"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in suite_rows:
            w.writerow({k: row[k] for k in fields})

    # Comparison plot: all cases, log amplitude.
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    for result in all_results:
        rr = result["records"]
        xs = [r.generation for r in rr if r.new_abs > 0 and math.isfinite(r.new_abs) and r.new_abs <= 1e40]
        ys = [r.new_abs for r in rr if r.new_abs > 0 and math.isfinite(r.new_abs) and r.new_abs <= 1e40]
        if xs:
            ax.plot(xs, ys, linewidth=1.2, label=result["label"])
    ax.set_yscale("log")
    ax.set_xlabel("generation / wave index")
    ax.set_ylabel("|z_new| (log scale)")
    ax.set_title("Literal anonymous self+cross interaction: parameter comparison")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(base_outdir, "suite_amplitude_log.png"), dpi=200)
    plt.close(fig)

    # Near-critical linear-scale comparison.
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    wanted = {"nearcritical_real_1p001", "nearcritical_real_0p999", "nearcritical_phase45_1p001", "nearcritical_phase45_0p999", "unit_real", "unit_phase45"}
    for result in all_results:
        if result["label"] not in wanted:
            continue
        rr = result["records"]
        xs = [r.generation for r in rr if math.isfinite(r.new_abs)]
        ys = [r.new_abs for r in rr if math.isfinite(r.new_abs)]
        ax.plot(xs, ys, linewidth=1.2, label=result["label"])
    ax.set_xlabel("generation / wave index")
    ax.set_ylabel("|z_new|")
    ax.set_title("Near |z1|=1: generated-wave amplitude")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(base_outdir, "suite_nearcritical_linear.png"), dpi=200)
    plt.close(fig)

    # Near-critical zoom: omit the first two O(1) waves so the small generated branch is visible.
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    wanted_zoom = {"nearcritical_real_1p001", "nearcritical_real_0p999", "nearcritical_phase45_1p001", "nearcritical_phase45_0p999"}
    for result in all_results:
        if result["label"] not in wanted_zoom:
            continue
        rr = [r for r in result["records"] if 3 <= r.generation <= 80 and math.isfinite(r.new_abs)]
        ax.plot([r.generation for r in rr], [r.new_abs for r in rr], linewidth=1.4, label=result["label"])
    ax.set_xlabel("generation / wave index")
    ax.set_ylabel("|z_new|")
    ax.set_title("Near |z1|=1: small generated branch (generation 3-80)")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(base_outdir, "suite_nearcritical_zoom.png"), dpi=200)
    plt.close(fig)

    # Complex-plane rays for mixed/pure phase cases.
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)
    wanted2 = {"nearcritical_real_1p001", "nearcritical_phase45_1p001", "unit_imag", "subcritical_imag_i_over_sqrt2"}
    for result in all_results:
        if result["label"] not in wanted2:
            continue
        rr = result["records"]
        xs = [r.new_re for r in rr if math.isfinite(r.new_re) and abs(r.new_re) < 1e6 and abs(r.new_im) < 1e6]
        ys = [r.new_im for r in rr if math.isfinite(r.new_im) and abs(r.new_re) < 1e6 and abs(r.new_im) < 1e6]
        ax.plot(xs, ys, marker="o", markersize=2, linewidth=1, label=result["label"])
    ax.axhline(0, linewidth=0.7)
    ax.axvline(0, linewidth=0.7)
    ax.set_xlabel("Re(z_new)")
    ax.set_ylabel("Im(z_new)")
    ax.set_title("Generated waves in the complex plane")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(base_outdir, "suite_complex_plane.png"), dpi=200)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description="Literal anonymous complex-wave generation experiment")
    p.add_argument("a", nargs="?", type=float, help="initial real part a")
    p.add_argument("b", nargs="?", type=float, help="initial imaginary coefficient b (z1=a+i b)")
    p.add_argument("--max-generations", type=int, default=200, help="numerical run limit; not a physical parameter")
    p.add_argument("--stop-abs", type=float, default=1e120, help="overflow guard; not a physical parameter")
    p.add_argument("--outdir", default="results_single", help="output directory")
    p.add_argument("--suite", action="store_true", help="run the built-in parameter test suite")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    if args.suite:
        run_suite(args.outdir, args.max_generations, args.stop_abs)
        return
    if args.a is None or args.b is None:
        p.error("provide a and b, or use --suite")
    result = run_case(args.outdir, args.a, args.b, args.max_generations, args.stop_abs)
    print(json.dumps({k: v for k, v in result.items() if k != "records"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
