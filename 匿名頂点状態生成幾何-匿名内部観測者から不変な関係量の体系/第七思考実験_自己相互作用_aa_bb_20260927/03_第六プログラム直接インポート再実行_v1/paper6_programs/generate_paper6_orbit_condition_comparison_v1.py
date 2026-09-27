#!/usr/bin/env python3
"""Paper 6 orbit-condition comparison figures.

This script produces comparison figures for the five charged-binary cases already
used in the strict eight-state experiment.  It does NOT rerun or modify the strict
state generator.

Outputs
-------
1) figure01_five_case_reference_orbits_same_scale.svg/png
   Independent analytic LO reference orbits, five cases, identical x/y scale.
2) figure02_five_case_strict_generator_overlay_panels.svg/png
   Montage of the existing per-case strict-generator vs analytic-reference overlay
   SVGs.  This is a presentation-only composition; source SVGs are unchanged.
3) figure03_five_case_radius_vs_cycles.svg/png
   Radius as a function of accumulated orbital cycles for all five analytic cases.
4) figure04_five_case_first_10_orbits_same_scale.svg/png
   First ten analytic orbits for each case, same x/y scale.
5) paper6_orbit_condition_case_summary.csv

The analytic formulas are the same closed-form LO benchmark used in
07_荷電二体系_解析近似基準解_20260925:
  Z = 1 - lambda_A lambda_B
  D = (lambda_A-lambda_B)^2
  dr/dt = -A/r^2 - B/r^3
  A = (4/3) mu D Z
  B = (64/5) mu Z^2
  omega^2 = Z/r^3
with G=M=c=1, mA=mB=0.5, mu=0.25, r0=50, rf=20.

The strict montage expects the five existing figure01_orbit_overlay.svg files.
Provide them via --strict-svg-dir with filenames listed in STRICT_FILES below, or
use --strict-svg explicitly five times in CASES order.
"""
from __future__ import annotations

import argparse
import base64
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

TWO_PI = 2.0 * math.pi
MU = 0.25
R0 = 50.0
RF = 20.0

CASES = [
    ("n1_attractive", +0.30, -0.30),
    ("n1_repulsive", +0.30, +0.30),
    ("n3_attractive", +0.90, -0.90),
    ("n3_repulsive", +0.90, +0.90),
    ("n4_attractive", +1.20, -1.20),
]

PANEL_TITLES = {
    "n1_attractive": r"n1 attractive: $\lambda_A=+0.30,\ \lambda_B=-0.30$",
    "n1_repulsive": r"n1 repulsive: $\lambda_A=+0.30,\ \lambda_B=+0.30$",
    "n3_attractive": r"n3 attractive: $\lambda_A=+0.90,\ \lambda_B=-0.90$",
    "n3_repulsive": r"n3 repulsive: $\lambda_A=+0.90,\ \lambda_B=+0.90$",
    "n4_attractive": r"n4 attractive: $\lambda_A=+1.20,\ \lambda_B=-1.20$",
}

# When --strict-svg-dir points to a source directory containing these files.
STRICT_FILES = {
    "n1_attractive": "n1_attractive_figure01_orbit_overlay.svg",
    "n1_repulsive": "n1_repulsive_figure01_orbit_overlay.svg",
    "n3_attractive": "n3_attractive_figure01_orbit_overlay.svg",
    "n3_repulsive": "n3_repulsive_figure01_orbit_overlay.svg",
    "n4_attractive": "n4_attractive_figure01_orbit_overlay.svg",
}


def primitive_time(r: np.ndarray, A: float, B: float) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    if abs(A) < 1e-30:
        return r**4 / (4.0 * B)
    return (
        r**3 / (3.0 * A)
        - B * r**2 / (2.0 * A**2)
        + B**2 * r / A**3
        - B**3 * np.log(A * r + B) / A**4
    )


def primitive_phase(r: np.ndarray, A: float, B: float, Z: float) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    if abs(A) < 1e-30:
        return math.sqrt(Z) * (2.0 / (5.0 * B)) * r**2.5
    sr = np.sqrt(r)
    term = (
        2.0 * r**1.5 / (3.0 * A)
        - 2.0 * B * sr / A**2
        + 2.0 * B**1.5 / A**2.5 * np.arctan(np.sqrt(A * r / B))
    )
    return math.sqrt(Z) * term


def make_reference(lambda_A: float, lambda_B: float, base_rows: int = 200001):
    Z = 1.0 - lambda_A * lambda_B
    D = (lambda_A - lambda_B) ** 2
    if Z <= 0:
        raise ValueError(f"Z={Z} <= 0; outside the circular reference domain")
    A = (4.0 / 3.0) * MU * D * Z
    B = (64.0 / 5.0) * MU * Z**2

    # Dense monotone r grid only for closed-form evaluation/interpolation.
    r_dense = np.linspace(R0, RF, base_rows)
    Ft0 = float(primitive_time(np.array([R0]), A, B)[0])
    Fp0 = float(primitive_phase(np.array([R0]), A, B, Z)[0])
    t_dense = Ft0 - primitive_time(r_dense, A, B)
    phi_dense = Fp0 - primitive_phase(r_dense, A, B, Z)
    cycles_final = float(phi_dense[-1] / TWO_PI)

    # Plot uniformly in phase so high-cycle cases are not aliased by a uniform-r grid.
    # 16 points/orbit, bounded for SVG size while keeping the spiral readable.
    n_phase = int(max(8000, min(600000, math.ceil(cycles_final * 128.0))))
    phi = np.linspace(0.0, float(phi_dense[-1]), n_phase)
    # np.interp requires increasing xp. phi_dense is increasing and r_dense decreasing.
    r = np.interp(phi, phi_dense, r_dense)
    x = r * np.cos(phi)
    y = r * np.sin(phi)

    return {
        "Z": Z,
        "D": D,
        "A": A,
        "B": B,
        "cycles_final": cycles_final,
        "t_final": float(t_dense[-1]),
        "phi_final": float(phi_dense[-1]),
        "phi": phi,
        "r": r,
        "x": x,
        "y": y,
        "phi_dense": phi_dense,
        "r_dense": r_dense,
    }


def save_both(fig, outdir: Path, stem: str):
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_reference_panels(refs: dict, outdir: Path):
    fig, axes = plt.subplots(3, 2, figsize=(11, 15), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, (case, la, lb) in zip(axes, CASES):
        d = refs[case]
        ax.plot(d["x"], d["y"], lw=0.55, rasterized=(d["cycles_final"] > 1000))
        ax.scatter([d["x"][0], d["x"][-1]], [d["y"][0], d["y"][-1]], s=18)
        ax.set_title(f"{case}\nC={d['Z']:.2f}, D={d['D']:.2f}, cycles={d['cycles_final']:.1f}")
        ax.set_aspect("equal", "box")
        ax.set_xlim(-52, 52)
        ax.set_ylim(-52, 52)
        ax.grid(True, alpha=0.20)
        ax.set_xlabel("relative x / M")
        ax.set_ylabel("relative y / M")
    axes[-1].axis("off")
    fig.suptitle("Five charged-binary conditions: analytic reference orbits (same spatial scale)", y=0.995)
    fig.tight_layout()
    save_both(fig, outdir, "figure01_five_case_reference_orbits_same_scale")


def plot_radius_vs_cycles(refs: dict, outdir: Path):
    fig, ax = plt.subplots(figsize=(9, 6))
    for case, _, _ in CASES:
        d = refs[case]
        cycles = d["phi_dense"] / TWO_PI
        # Thin by index for compact vector output while preserving monotonic curve.
        step = max(1, len(cycles) // 12000)
        ax.plot(cycles[::step], d["r_dense"][::step], label=f"{case} (C={d['Z']:.2f}, D={d['D']:.2f})")
    ax.set_xlabel("accumulated orbital cycles")
    ax.set_ylabel("r / M")
    ax.set_title("Inspiral rate by orbital cycle for the five charged-binary conditions")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_both(fig, outdir, "figure03_five_case_radius_vs_cycles")


def plot_first_10(refs: dict, outdir: Path):
    fig, axes = plt.subplots(3, 2, figsize=(11, 15), sharex=True, sharey=True)
    axes = axes.ravel()
    phi_end = 10.0 * TWO_PI
    for ax, (case, _, _) in zip(axes, CASES):
        d = refs[case]
        stop = min(phi_end, d["phi_dense"][-1])
        phi = np.linspace(0.0, stop, 4000)
        r = np.interp(phi, d["phi_dense"], d["r_dense"])
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        dr10 = R0 - float(r[-1])
        ax.plot(x, y, lw=0.8)
        ax.set_title(f"{case}\nfirst 10 cycles: $\\Delta r$={dr10:.4g}")
        ax.set_aspect("equal", "box")
        ax.set_xlim(-52, 52)
        ax.set_ylim(-52, 52)
        ax.grid(True, alpha=0.20)
        ax.set_xlabel("relative x / M")
        ax.set_ylabel("relative y / M")
    axes[-1].axis("off")
    fig.suptitle("First ten analytic orbits for the five charged-binary conditions (same spatial scale)", y=0.995)
    fig.tight_layout()
    save_both(fig, outdir, "figure04_five_case_first_10_orbits_same_scale")


def compose_svg_montage(svg_paths: list[Path], out_svg: Path):
    """Embed five existing SVGs in a 3x2 SVG montage without modifying them."""
    W, H = 1200, 1650
    margin = 20
    cell_w = (W - 3 * margin) / 2
    cell_h = (H - 4 * margin) / 3
    positions = [
        (margin, margin),
        (2 * margin + cell_w, margin),
        (margin, 2 * margin + cell_h),
        (2 * margin + cell_w, 2 * margin + cell_h),
        (margin, 3 * margin + 2 * cell_h),
    ]
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
        '<text x="600" y="18" text-anchor="middle" font-family="sans-serif" font-size="16">Five strict 8-state orbit overlays: analytic reference vs generator</text>',
    ]
    for (case, _, _), src, (x, y) in zip(CASES, svg_paths, positions):
        payload = base64.b64encode(src.read_bytes()).decode("ascii")
        label_y = y + 18
        image_y = y + 24
        image_h = cell_h - 24
        lines.append(f'<text x="{x + cell_w/2:.1f}" y="{label_y:.1f}" text-anchor="middle" font-family="sans-serif" font-size="13">{case}</text>')
        lines.append(
            f'<image x="{x:.1f}" y="{image_y:.1f}" width="{cell_w:.1f}" height="{image_h:.1f}" '
            f'preserveAspectRatio="xMidYMid meet" href="data:image/svg+xml;base64,{payload}"/>'
        )
    lines.append('</svg>')
    out_svg.write_text("\n".join(lines), encoding="utf-8")


def rasterize_svg(svg_path: Path, png_path: Path):
    try:
        import cairosvg
        cairosvg.svg2png(url=str(svg_path), write_to=str(png_path), output_width=1600)
    except Exception as e:
        print(f"warning: PNG preview not generated for {svg_path.name}: {e}")


def write_summary(refs: dict, outdir: Path):
    p = outdir / "paper6_orbit_condition_case_summary.csv"
    with p.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["case", "lambda_A", "lambda_B", "C_Z", "D_delta_lambda_sq", "A_em", "B_gw", "cycles_r50_to_r20", "t_final"])
        for case, la, lb in CASES:
            d = refs[case]
            w.writerow([case, la, lb, f"{d['Z']:.17g}", f"{d['D']:.17g}", f"{d['A']:.17g}", f"{d['B']:.17g}", f"{d['cycles_final']:.17g}", f"{d['t_final']:.17g}"])


def resolve_strict_paths(args) -> list[Path] | None:
    if args.strict_svg:
        if len(args.strict_svg) != 5:
            raise SystemExit("--strict-svg must be supplied exactly five times, in CASES order")
        paths = [Path(x) for x in args.strict_svg]
    elif args.strict_svg_dir:
        root = Path(args.strict_svg_dir)
        paths = [root / STRICT_FILES[c] for c, _, _ in CASES]
    else:
        return None
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise SystemExit("missing strict SVG source(s): " + ", ".join(map(str, missing)))
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent / "figures")
    ap.add_argument("--strict-svg-dir", type=Path, default=None)
    ap.add_argument("--strict-svg", action="append", default=[])
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    refs = {case: make_reference(la, lb) for case, la, lb in CASES}
    plot_reference_panels(refs, args.outdir)
    plot_radius_vs_cycles(refs, args.outdir)
    plot_first_10(refs, args.outdir)
    write_summary(refs, args.outdir)

    strict_paths = resolve_strict_paths(args)
    if strict_paths:
        out_svg = args.outdir / "figure02_five_case_strict_generator_overlay_panels.svg"
        compose_svg_montage(strict_paths, out_svg)
        rasterize_svg(out_svg, args.outdir / "figure02_five_case_strict_generator_overlay_panels.png")
    else:
        print("strict montage skipped: supply --strict-svg-dir or five --strict-svg arguments")

    print("created:")
    for p in sorted(args.outdir.iterdir()):
        print(" ", p)


if __name__ == "__main__":
    main()
