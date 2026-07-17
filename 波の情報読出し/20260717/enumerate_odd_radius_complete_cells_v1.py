#!/usr/bin/env python3
"""Enumerate integer-centred unit 4-cubes fully contained in 4-balls.

For a cell centred at c in Z^4, complete containment in the radius-r ball is
equivalent to

    sum_j (|c_j| + 1/2)^2 <= r^2.

The implementation evaluates the equivalent integer inequality

    sum_j (2|c_j| + 1)^2 <= 4 r^2

so the counts and boundary-contact tests contain no floating-point rounding.
"""

from __future__ import annotations

import argparse
import csv
from itertools import product
from pathlib import Path


DEFAULT_RADII = (1, 3, 5, 7)
EXPECTED_COUNTS = {1: 1, 3: 137, 5: 1545, 7: 7281}


def containment_numerator(cell: tuple[int, int, int, int]) -> int:
    """Return four times the maximum squared radius of the unit 4-cube."""
    return sum((2 * abs(coordinate) + 1) ** 2 for coordinate in cell)


def enumerate_complete_cells(radius: int) -> tuple[int, int]:
    """Return total contained cells and cells touching the radius-r shell."""
    if radius <= 0:
        raise ValueError("radius must be a positive integer")

    radius_numerator = 4 * radius * radius
    total = 0
    boundary_contacts = 0
    coordinate_range = range(-radius, radius + 1)

    for cell in product(coordinate_range, repeat=4):
        numerator = containment_numerator(cell)
        if numerator <= radius_numerator:
            total += 1
            if numerator == radius_numerator:
                boundary_contacts += 1

    return total, boundary_contacts


def parse_radii(value: str) -> tuple[int, ...]:
    radii = tuple(int(item.strip()) for item in value.split(",") if item.strip())
    if not radii:
        raise argparse.ArgumentTypeError("at least one radius is required")
    if any(radius <= 0 for radius in radii):
        raise argparse.ArgumentTypeError("all radii must be positive integers")
    return radii


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--radii",
        type=parse_radii,
        default=DEFAULT_RADII,
        help="comma-separated positive integer radii (default: 1,3,5,7)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent
        / "odd_radius_cell_enumeration_v1_results",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, int | str]] = []
    for radius in args.radii:
        total, boundary_contacts = enumerate_complete_cells(radius)
        expected = EXPECTED_COUNTS.get(radius)
        verified = expected is None or total == expected
        if not verified:
            raise AssertionError(
                f"radius {radius}: expected {expected} cells, obtained {total}"
            )

        rows.append(
            {
                "k_for_r_eq_2k_plus_1": (radius - 1) // 2
                if radius % 2 == 1
                else "",
                "radius": radius,
                "complete_cell_count": total,
                "boundary_contact_cell_count": boundary_contacts,
                "expected_count": expected if expected is not None else "",
                "verified": "true" if verified else "false",
            }
        )

    output_path = args.output_dir / "odd_radius_cell_counts_v1.csv"
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"r={row['radius']}: N_cell={row['complete_cell_count']}, "
            f"boundary_contacts={row['boundary_contact_cell_count']}"
        )
    print(output_path)


if __name__ == "__main__":
    main()
