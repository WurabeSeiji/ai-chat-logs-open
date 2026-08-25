#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot n=64 round12 brute-force results.

Reads the CSV produced by:
    n64_round12_simulation.py

Outputs:
    n64_round12_results.png
    n64_round12_results.svg

Usage:
    python3 n64_round12_plot.py

Optional:
    python3 n64_round12_plot.py \
        --input n64_round12_results.csv \
        --png n64_round12_results.png \
        --svg n64_round12_results.svg
"""

import argparse
import csv
from pathlib import Path
import math

import matplotlib.pyplot as plt


def read_results(path: Path):
    rows = []

    with path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            rows.append({
                "step": int(row["step"]),
                "state_count": int(row["state_count"]),
                "max_harmonic": int(row["max_harmonic"]),
                "zero_closure": float(row["zero_closure"]),
            })

    return rows


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        default="n64_round12_results.csv",
    )
    parser.add_argument(
        "--png",
        default="n64_round12_results.png",
    )
    parser.add_argument(
        "--svg",
        default="n64_round12_results.svg",
    )
    parser.add_argument(
        "--title",
        default=(
            "Original n=64 model, round(real,12)/round(imag,12) "
            "WRITE deduplication"
        ),
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    rows = read_results(input_path)

    if not rows:
        raise RuntimeError(f"No data rows found in {input_path}")

    steps = [r["step"] for r in rows]
    counts = [r["state_count"] for r in rows]
    closures = [r["zero_closure"] for r in rows]

    fig, ax_left = plt.subplots(figsize=(10.5, 6.0))

    line_states = ax_left.plot(
        steps,
        counts,
        marker="o",
        markersize=3,
        color="#1f77b4",
        label="Distinct states by rounded complex value",
    )

    ax_left.set_xlabel("Processing step")
    ax_left.set_ylabel("Distinct state count")
    ax_left.set_xlim(min(steps), max(steps))
    ax_left.grid(True, alpha=0.30)

    ax_right = ax_left.twinx()

    line_closure = ax_right.plot(
        steps,
        closures,
        marker="s",
        markersize=2.5,
        color="#d62728",
        label=r"Zero-closure observation  $|\sum z^2|/\sum |z|^2$",
    )

    ax_right.set_ylabel(
        r"Zero-closure observation  $|\sum z^2|/\sum |z|^2$"
    )

    lines = line_states + line_closure
    labels = [line.get_label() for line in lines]
    ax_left.legend(lines, labels, loc="best")

    ax_left.set_title(args.title)

    fig.tight_layout()

    png_path = Path(args.png)
    svg_path = Path(args.svg)

    fig.savefig(
        png_path,
        dpi=180,
        bbox_inches="tight",
    )
    fig.savefig(
        svg_path,
        bbox_inches="tight",
    )

    print(f"Read: {input_path}")
    print(f"Rows: {len(rows)}")
    print(f"Last completed step: {rows[-1]['step']}")
    print(f"PNG: {png_path}")
    print(f"SVG: {svg_path}")

    plt.show()


if __name__ == "__main__":
    main()
