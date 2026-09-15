#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate raw complex-state data for the U^4 = I, base + one harmonic experiment.

Definitions
-----------
U = i
z_base     = U^k
z_harmonic = U^(m*k)
z_sum      = z_base + z_harmonic

m = 2, 3

Forward order:
0, 1, 2, 3

Reverse order:
0, 3, 2, 1

Outputs
-------
K4_harmonic_complex_raw_states_m2_m3_forward_reverse_20260915.csv
K4_harmonic_complex_raw_states_m2_20260915.csv
K4_harmonic_complex_raw_states_m3_20260915.csv
K4_harmonic_complex_raw_states_README_20260915.md
"""

from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
U = 1j

def clean(x: float) -> float:
    x = float(np.round(x, 15))
    if abs(x) < 1e-14:
        return 0.0
    return x

def build_rows():
    rows = []

    for m in (2, 3):
        for direction, seq in (
            ("forward", [0, 1, 2, 3]),
            ("reverse", [0, 3, 2, 1]),
        ):
            for step, k in enumerate(seq):
                z_base = U ** k
                z_harm = U ** (m * k)
                z_sum = z_base + z_harm

                rows.append({
                    "m": m,
                    "direction": direction,
                    "step": step,
                    "k": k,
                    "theta_rad": k * np.pi / 2,
                    "theta_over_pi": k / 2,
                    "base_re": clean(np.real(z_base)),
                    "base_im": clean(np.imag(z_base)),
                    "harmonic_re": clean(np.real(z_harm)),
                    "harmonic_im": clean(np.imag(z_harm)),
                    "sum_re": clean(np.real(z_sum)),
                    "sum_im": clean(np.imag(z_sum)),
                })

    return pd.DataFrame(rows)

def write_outputs(df: pd.DataFrame):
    all_path = OUT / "K4_harmonic_complex_raw_states_m2_m3_forward_reverse_20260915.csv"
    m2_path = OUT / "K4_harmonic_complex_raw_states_m2_20260915.csv"
    m3_path = OUT / "K4_harmonic_complex_raw_states_m3_20260915.csv"
    readme_path = OUT / "K4_harmonic_complex_raw_states_README_20260915.md"

    df.to_csv(all_path, index=False)
    df[df["m"] == 2].to_csv(m2_path, index=False)
    df[df["m"] == 3].to_csv(m3_path, index=False)

    readme_path.write_text(
        """# K4 harmonic complex raw states

Raw complex-state data for the U^4=I, base + one harmonic experiment.

## Definitions

U = i

z_base = U^k

z_harmonic = U^(m k)

z_sum = z_base + z_harmonic

m = 2, 3

Forward state order:
0, 1, 2, 3

Reverse state order:
0, 3, 2, 1

## Columns

- m
- direction
- step
- k
- theta_rad
- theta_over_pi
- base_re
- base_im
- harmonic_re
- harmonic_im
- sum_re
- sum_im

`step` is only the row/order index in the selected sequence.
`k` is the discrete phase-state index.

No interpolation, inference, randomization, or fitted values are used.
All complex states are generated directly from U=i.
""",
        encoding="utf-8",
    )

    return [all_path, m2_path, m3_path, readme_path]

def main():
    df = build_rows()
    outputs = write_outputs(df)

    print(df.to_string(index=False))
    print("\nGenerated:")
    for p in outputs:
        print(p)

if __name__ == "__main__":
    main()
