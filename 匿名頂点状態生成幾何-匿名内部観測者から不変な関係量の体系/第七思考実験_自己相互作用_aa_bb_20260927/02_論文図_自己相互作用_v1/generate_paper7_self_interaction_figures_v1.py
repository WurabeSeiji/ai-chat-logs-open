#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate Paper 7 self-interaction figures (SVG + PNG).

This script is intentionally narrow.  It does not change the Paper 6 transition
map.  It visualizes the independent aa/bb self-interaction cases obtained from
G=q0=c=1, M/q0=20/3, with D_self=0 and C_self=1-lambda_self^2.

Inputs expected next to this script or under --data-root:
  results/self_n1/sampled_trajectory.csv
  results/self_n3_fresh/sampled_trajectory_prefix100k.csv (optional check only)

The complete n=1 and n=3 orbit curves are evaluated from the same leading-order
closed-form benchmark used by Paper 6.  For n=1, the fresh strict-generator sample
is overlaid.  For n=3, the fresh 100k-step prefix is overlaid as an independent
check; the complete strict n=3 trajectory is already identical to Paper 6
n3_repulsive because the initial full state and transition(z) are identical.

Outputs are written in both SVG and PNG for reproducibility/presentation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

TWO_PI = 2.0 * math.pi
MU = 0.25
R0 = 50.0
RF = 20.0
STEPS_PER_ORBIT = 4000
HSTEP = TWO_PI / STEPS_PER_ORBIT

CASES = [
    ("self_n1", 1, 0.30, 0.91, 0.0),
    ("self_n3", 3, 0.90, 0.19, 0.0),
]


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


def primitive_phase(r: np.ndarray, A: float, B: float, C: float) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    if abs(A) < 1e-30:
        return math.sqrt(C) * (2.0 / (5.0 * B)) * r**2.5
    sr = np.sqrt(r)
    term = (
        2.0 * r**1.5 / (3.0 * A)
        - 2.0 * B * sr / A**2
        + 2.0 * B**1.5 / A**2.5 * np.arctan(np.sqrt(A * r / B))
    )
    return math.sqrt(C) * term


def make_reference(C: float, D: float, base_rows: int = 200001):
    if C <= 0:
        raise ValueError(f"C={C} <= 0 outside current real quasicircular domain")
    A = (4.0 / 3.0) * MU * D * C
    B = (64.0 / 5.0) * MU * C**2
    r_dense = np.linspace(R0, RF, base_rows)
    Ft0 = float(primitive_time(np.array([R0]), A, B)[0])
    Fp0 = float(primitive_phase(np.array([R0]), A, B, C)[0])
    t_dense = Ft0 - primitive_time(r_dense, A, B)
    phi_dense = Fp0 - primitive_phase(r_dense, A, B, C)
    cycles_final = float(phi_dense[-1] / TWO_PI)
    n_phase = int(max(12000, min(600000, math.ceil(cycles_final * 128.0))))
    phi = np.linspace(0.0, float(phi_dense[-1]), n_phase)
    r = np.interp(phi, phi_dense, r_dense)
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    pgw_dense = (32.0 / 5.0) * MU**2 * C**3 / (r_dense**5)
    pem_dense = (2.0 / 3.0) * MU**2 * D * C**2 / (r_dense**4)
    return dict(C=C, D=D, A=A, B=B, r_dense=r_dense, t_dense=t_dense,
                phi_dense=phi_dense, cycles_final=cycles_final, phi=phi, r=r,
                x=x, y=y, pgw_dense=pgw_dense, pem_dense=pem_dense)


def save_both(fig, outdir: Path, stem: str):
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def load_csv(path: Path):
    if not path.exists():
        return None
    return np.genfromtxt(path, delimiter=",", names=True)


def plot_orbits(refs, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5), sharex=True, sharey=True)
    for ax, (name, n, lam, C, D) in zip(axes, CASES):
        d = refs[name]
        ax.plot(d["x"], d["y"], lw=0.55, rasterized=(d["cycles_final"] > 1000))
        ax.scatter([d["x"][0], d["x"][-1]], [d["y"][0], d["y"][-1]], s=18)
        ax.set_title(f"aa = bb, n={n}: C={C:.2f}, D={D:.2f}\ncycles={d['cycles_final']:.3f}")
        ax.set_aspect("equal", "box")
        ax.set_xlim(-52, 52); ax.set_ylim(-52, 52)
        ax.grid(True, alpha=0.20)
        ax.set_xlabel("x readout")
        ax.set_ylabel("y readout")
    fig.suptitle("Independent self-interaction orbits: aa and bb are identical")
    fig.tight_layout()
    save_both(fig, outdir, "figure01_self_orbits_n1_n3_same_scale")


def plot_radius_cycles(refs, outdir):
    fig, ax = plt.subplots(figsize=(9, 6))
    for name, n, lam, C, D in CASES:
        d = refs[name]
        cycles = d["phi_dense"] / TWO_PI
        step = max(1, len(cycles)//12000)
        ax.plot(cycles[::step], d["r_dense"][::step], label=f"n={n}, C={C:.2f}, D=0")
    ax.set_xlabel("accumulated orbital cycles")
    ax.set_ylabel("P (radius readout)")
    ax.set_title("Self-interaction inspiral: radius versus accumulated cycles")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    save_both(fig, outdir, "figure02_self_radius_vs_cycles")


def plot_powers(refs, outdir):
    fig, axes = plt.subplots(2, 1, figsize=(9, 9), sharex=True)
    ax=axes[0]
    for name, n, lam, C, D in CASES:
        d=refs[name]
        ax.plot(d["r_dense"], d["pgw_dense"], label=f"GW n={n}, C={C:.2f}")
    ax.set_yscale("log")
    ax.set_ylabel("leading GW quadrupole power")
    ax.set_title("Self-interaction radiation in the Paper 6 approximation")
    ax.grid(True, alpha=0.25)
    ax.legend()
    ax=axes[1]
    for name, n, lam, C, D in CASES:
        d=refs[name]
        ax.plot(d["r_dense"], d["pem_dense"], label=f"EM dipole n={n}")
    ax.set_xlabel("P (radius readout)")
    ax.set_ylabel("leading EM dipole power")
    ax.set_ylim(-1e-18, 1e-18)
    ax.grid(True, alpha=0.25)
    ax.legend()
    axes[0].invert_xaxis(); axes[1].invert_xaxis()
    fig.tight_layout()
    save_both(fig, outdir, "figure03_self_gw_and_em_power")


def plot_generator_checks(refs, data_root, outdir):
    n1 = load_csv(data_root / "results/self_n1/sampled_trajectory.csv")
    n3 = load_csv(data_root / "results/self_n3_fresh/sampled_trajectory_prefix100k.csv")
    fig, axes = plt.subplots(2, 1, figsize=(9, 9))
    for ax, arr, name, n in [(axes[0],n1,"self_n1",1),(axes[1],n3,"self_n3",3)]:
        d=refs[name]
        if arr is not None:
            cyc = arr["step"] / STEPS_PER_ORBIT
            ax.scatter(cyc, arr["P"], s=5, label="strict generator sampled rows")
            phi = cyc * TWO_PI
            rref = np.interp(phi, d["phi_dense"], d["r_dense"])
            ax.plot(cyc, rref, lw=0.8, label="Paper 6 closed-form reference")
            resid = np.max(np.abs(arr["P"]-rref))
            ax.set_title(f"n={n} generator check; max sampled |ΔP|={resid:.3e}")
        else:
            ax.text(0.5,0.5,"sample file not found",transform=ax.transAxes,ha="center")
        ax.set_xlabel("accumulated cycles")
        ax.set_ylabel("P")
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.tight_layout()
    save_both(fig, outdir, "figure04_self_generator_reference_checks")


def write_summary(refs, outdir):
    with (outdir/"paper7_self_figure_summary.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(["case","n","lambda_self_abs","C","D","cycles_r50_to_r20","t_final","P_EM_zero"])
        for name,n,lam,C,D in CASES:
            d=refs[name]
            w.writerow([name,n,lam,C,D,f"{d['cycles_final']:.17g}",f"{d['t_dense'][-1]:.17g}",True])


def sha256(path: Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()


def write_manifest(outdir: Path):
    files=[]
    for p in sorted(outdir.iterdir()):
        if p.is_file(): files.append({"file":p.name,"sha256":sha256(p),"bytes":p.stat().st_size})
    (outdir/"FIGURE_SHA256SUMS.json").write_text(json.dumps(files,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--data-root",default=str(Path(__file__).resolve().parent))
    ap.add_argument("--outdir",default="figures")
    args=ap.parse_args()
    data_root=Path(args.data_root)
    outdir=Path(args.outdir)
    if not outdir.is_absolute(): outdir=data_root/outdir
    outdir.mkdir(parents=True,exist_ok=True)
    refs={name:make_reference(C,D) for name,n,lam,C,D in CASES}
    plot_orbits(refs,outdir)
    plot_radius_cycles(refs,outdir)
    plot_powers(refs,outdir)
    plot_generator_checks(refs,data_root,outdir)
    write_summary(refs,outdir)
    write_manifest(outdir)
    print(outdir)

if __name__ == "__main__":
    main()
