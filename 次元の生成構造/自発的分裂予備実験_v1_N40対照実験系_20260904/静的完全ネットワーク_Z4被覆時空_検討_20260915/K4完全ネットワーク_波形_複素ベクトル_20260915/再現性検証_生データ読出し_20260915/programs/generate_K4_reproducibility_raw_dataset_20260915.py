#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate the complete raw dataset needed to reproduce the K4 harmonic study.

This script intentionally preserves two numerically distinct but mathematically
equivalent state conventions used by the original programs:

1) static figures:
   U = exp(2*pi*i/4)

2) HTML / MP4 animation:
   U = 1j

It also stores the exact sampled wave arrays used by:
- static PNG/SVG figures: 2000 points, theta in [-2*pi, 2*pi]
- HTML animation: 800 points, theta in [0, 2*pi]
- MP4 animation: 1200 points, theta in [0, 2*pi]

All CSVs are written with 17 significant digits.
"""

from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "raw_data"
OUT.mkdir(parents=True, exist_ok=True)

def write_states(U, filename):
    rows = []
    for m in (2, 3):
        for direction, seq in (
            ("forward", [0,1,2,3]),
            ("reverse", [0,3,2,1]),
        ):
            for step, k in enumerate(seq):
                z1 = U**k
                zm = U**(m*k)
                zs = z1 + zm
                rows.append({
                    "m": m,
                    "direction": direction,
                    "step": step,
                    "k": k,
                    "theta_rad": k*np.pi/2,
                    "theta_over_pi": k/2,
                    "base_re": float(np.real(z1)),
                    "base_im": float(np.imag(z1)),
                    "harmonic_re": float(np.real(zm)),
                    "harmonic_im": float(np.imag(zm)),
                    "sum_re": float(np.real(zs)),
                    "sum_im": float(np.imag(zs)),
                })
    pd.DataFrame(rows).to_csv(OUT/filename,index=False,float_format="%.17g")

write_states(np.exp(2j*np.pi/4), "complex_states_static_exp.csv")
write_states(1j, "complex_states_animation_exact_i.csv")

# canonical alias = animation state data
pd.read_csv(
    OUT/"complex_states_animation_exact_i.csv",
    float_precision="round_trip"
).to_csv(
    OUT/"complex_states.csv",
    index=False,
    float_format="%.17g"
)

for m in (2,3):
    # static figures
    th = np.linspace(-2*np.pi, 2*np.pi, 2000)
    pd.DataFrame({
        "sample_index": np.arange(len(th)),
        "theta_rad": th,
        "base": np.cos(th),
        "harmonic": np.cos(m*th),
        "composite": np.cos(th)+np.cos(m*th),
    }).to_csv(
        OUT/f"wave_samples_static_2000_m{m}.csv",
        index=False,
        float_format="%.17g"
    )

    # HTML
    th = np.linspace(0, 2*np.pi, 800)
    pd.DataFrame({
        "sample_index": np.arange(len(th)),
        "theta_rad": th,
        "base": np.cos(th),
        "harmonic": np.cos(m*th),
        "composite": np.cos(th)+np.cos(m*th),
    }).to_csv(
        OUT/f"wave_samples_html_800_m{m}.csv",
        index=False,
        float_format="%.17g"
    )

    # MP4
    th = np.linspace(0, 2*np.pi, 1200)
    pd.DataFrame({
        "sample_index": np.arange(len(th)),
        "theta_rad": th,
        "base": np.cos(th),
        "harmonic": np.cos(m*th),
        "composite": np.cos(th)+np.cos(m*th),
    }).to_csv(
        OUT/f"wave_samples_mp4_1200_m{m}.csv",
        index=False,
        float_format="%.17g"
    )

rep = []
for m in (2,3):
    th = np.pi/4
    z1 = np.exp(1j*th)
    zm = np.exp(1j*m*th)
    rep.append({
        "m": m,
        "theta_rad": th,
        "base_re": float(z1.real),
        "base_im": float(z1.imag),
        "harmonic_re": float(zm.real),
        "harmonic_im": float(zm.imag),
    })
pd.DataFrame(rep).to_csv(
    OUT/"representative_vectors_theta_pi4.csv",
    index=False,
    float_format="%.17g"
)

manifest = pd.DataFrame([
    ["complex_states_static_exp.csv","Static K4/vector state raw data using exp(2*pi*i/4)"],
    ["complex_states_animation_exact_i.csv","HTML/MP4 discrete state raw data using U=1j"],
    ["complex_states.csv","Canonical alias of animation exact-i state data"],
    ["wave_samples_static_2000_m2.csv","Static waveform samples m=2"],
    ["wave_samples_static_2000_m3.csv","Static waveform samples m=3"],
    ["wave_samples_html_800_m2.csv","HTML waveform samples m=2"],
    ["wave_samples_html_800_m3.csv","HTML waveform samples m=3"],
    ["wave_samples_mp4_1200_m2.csv","MP4 waveform samples m=2"],
    ["wave_samples_mp4_1200_m3.csv","MP4 waveform samples m=3"],
    ["representative_vectors_theta_pi4.csv","Representative theta=pi/4 complex vectors"],
], columns=["file","purpose"])
manifest.to_csv(OUT/"RAW_DATA_MANIFEST.csv",index=False)

print("Generated raw dataset in", OUT)
