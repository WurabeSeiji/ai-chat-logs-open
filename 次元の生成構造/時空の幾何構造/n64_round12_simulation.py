#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
n=64 simple brute-force experiment
==================================

Current experiment definition (kept exactly as discussed):

- Start from:
      Wave(harmonic=1, value=1+0j)

- Base phase:
      omega = exp(2*pi*i/n), n=64

- Rotation:
      Wave(h, z) -> Wave(h, omega**h * z)

- Pair read:
      Wave(ha, a), Wave(hb, b)
          -> Wave(ha+hb, a*b)

- Synchronous step:
    1. Read only the current snapshot.
    2. Write rotate(w) for every current wave.
    3. For every unordered pair including self-pairs, write
       rotate(pair_read(a,b)).

- WRITE deduplication:
    A candidate is NOT written if another candidate in the same next-state
    already has the same rounded complex value key:

        (round(value.real, 12), round(value.imag, 12))

    The first Wave encountered for each rounded complex value is retained.

- Sweep:
      step = 0 ... n*2 = 128

This is intentionally the simple brute-force formulation.
No Fourier aggregation, no modulo-state compression, no distinct-harmonic
optimization, and no physical reinterpretation are added here.

The program writes one CSV row per completed step, so plotting can be done
separately even if the run is interrupted.

Usage:
    python3 n64_round12_simulation.py

Optional:
    python3 n64_round12_simulation.py --max-step 128
    python3 n64_round12_simulation.py --output results.csv

WARNING:
The all-pairs operation is O(M^2) in the number of retained states.
The run can become extremely expensive.
"""

from dataclasses import dataclass
import argparse
import cmath
import csv
import math
import time
from pathlib import Path


N = 1
DEFAULT_N = 64


@dataclass(frozen=True)
class Wave:
    harmonic: int
    value: complex


def make_omega(n: int) -> complex:
    return cmath.exp(2j * math.pi / n)


def rotate(w: Wave, omega: complex) -> Wave:
    return Wave(
        w.harmonic,
        (omega ** w.harmonic) * w.value
    )


def pair_read(a: Wave, b: Wave) -> Wave:
    return Wave(
        a.harmonic + b.harmonic,
        a.value * b.value
    )


def value_key(z: complex):
    # Current agreed experimental key.
    return (
        round(z.real, 12),
        round(z.imag, 12),
    )


def synchronous_step(state, omega):
    """
    Current-state snapshot -> all writes for next state.

    Deduplication is ONLY by rounded complex value.
    """
    snapshot = tuple(state)
    writes = {}

    def write_if_new(w: Wave):
        key = value_key(w.value)
        if key not in writes:
            writes[key] = w

    # Single-state writes.
    for w in snapshot:
        write_if_new(rotate(w, omega))

    # Unordered pair writes, including self-pairs.
    for i, a in enumerate(snapshot):
        for j in range(i, len(snapshot)):
            candidate = rotate(pair_read(a, snapshot[j]), omega)
            write_if_new(candidate)

    return tuple(writes.values())


def zero_closure_observation(state):
    """
    Observation only:
        |sum z^2| / sum |z|^2
    """
    q = sum((w.value * w.value for w in state), 0j)
    den = sum(abs(w.value) ** 2 for w in state)
    return abs(q) / den if den else float("nan")


def write_csv_header(path: Path):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "step",
            "state_count",
            "max_harmonic",
            "max_abs_value",
            "zero_closure",
            "step_seconds",
            "elapsed_seconds",
        ])


def append_csv_row(path: Path, row):
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(row)
        f.flush()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=DEFAULT_N)
    parser.add_argument(
        "--max-step",
        type=int,
        default=None,
        help="Default: n*2",
    )
    parser.add_argument(
        "--output",
        default="n64_round12_results.csv",
        help="CSV output path",
    )
    args = parser.parse_args()

    n = args.n
    max_step = args.max_step if args.max_step is not None else n * 2
    omega = make_omega(n)

    output_path = Path(args.output)
    write_csv_header(output_path)

    state = (Wave(1, 1.0 + 0.0j),)

    run_start = time.perf_counter()
    previous_step_seconds = 0.0

    print(f"n = {n}")
    print(f"max_step = {max_step}")
    print(f"output = {output_path}")
    print()

    for step in range(max_step + 1):
        closure = zero_closure_observation(state)
        state_count = len(state)
        max_harmonic = max(w.harmonic for w in state)
        max_abs_value = max(abs(w.value) for w in state)
        elapsed = time.perf_counter() - run_start

        append_csv_row(
            output_path,
            [
                step,
                state_count,
                max_harmonic,
                repr(max_abs_value),
                repr(closure),
                repr(previous_step_seconds),
                repr(elapsed),
            ],
        )

        print(
            f"step={step:4d}  "
            f"states={state_count:10d}  "
            f"max_harmonic={max_harmonic:12d}  "
            f"max|z|={max_abs_value:.12g}  "
            f"closure={closure:.12g}  "
            f"last_step={previous_step_seconds:.3f}s  "
            f"elapsed={elapsed:.3f}s",
            flush=True,
        )

        if step >= max_step:
            break

        t0 = time.perf_counter()
        state = synchronous_step(state, omega)
        previous_step_seconds = time.perf_counter() - t0

    total = time.perf_counter() - run_start
    print()
    print(f"Completed through step {max_step}.")
    print(f"Total elapsed: {total:.3f} s")


if __name__ == "__main__":
    main()
