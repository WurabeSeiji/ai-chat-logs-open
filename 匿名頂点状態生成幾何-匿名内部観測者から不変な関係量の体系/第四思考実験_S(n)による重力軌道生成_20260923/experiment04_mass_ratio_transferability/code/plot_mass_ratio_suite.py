#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot Experiment 04 from saved raw CSV rows only; no dynamics recomputed."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def read_metrics(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def save(fig, figs, stem):
    fig.tight_layout()
    fig.savefig(figs / f"{stem}.png", dpi=180)
    fig.savefig(figs / f"{stem}.svg")
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--base", default=None); args=ap.parse_args()
    base=Path(args.base) if args.base else Path(__file__).resolve().parent.parent
    raw, summary, figs=base/"raw", base/"summary", base/"figures"
    figs.mkdir(parents=True, exist_ok=True)
    cfg=json.load(open(summary/"cases.json", encoding="utf-8"))
    cases=cfg["cases"]
    rows=read_metrics(summary/"case_metrics.csv")

    fig, axs=plt.subplots(2,2, figsize=(11.5,10.0))
    for ax,c in zip(axs.ravel(), cases):
        cid=c["case_id"]
        pn=read_csv(raw/f"pn_plot_{cid}.csv"); sn=read_csv(raw/f"sn_{cid}.csv")
        ax.plot(pn["x"],pn["y"],lw=1.15,label="Paper 3 PN")
        ax.plot(sn["x"],sn["y"],lw=1.0,ls="--",label="same S(n)")
        ax.plot([0],[0],marker="+",ms=6)
        ax.set_aspect("equal",adjustable="datalim")
        ax.set_title(f"q={c['q_mass_ratio']:g}, nu={c['nu']:.5f}")
        ax.set_xlabel("x [M]"); ax.set_ylabel("y [M]"); ax.grid(True,alpha=.2)
    axs[0,0].legend(fontsize=8)
    fig.suptitle("Experiment 04: same S(n) rule across mass ratios (p=60, e=0.45)", fontsize=13)
    fig.tight_layout(rect=(0,0,1,.97))
    fig.savefig(figs/"experiment04_mass_ratio_orbits.png",dpi=180)
    fig.savefig(figs/"experiment04_mass_ratio_orbits.svg")
    plt.close(fig)

    q=np.array([float(r["q_mass_ratio"]) for r in rows])
    nu=np.array([float(r["nu"]) for r in rows])
    rms=np.array([float(r["radial_rms_over_a"]) for r in rows])
    dprec=np.array([float(r["precession_error_rad"]) for r in rows])
    sh=np.array([float(r["sn_peri_shrink_fraction"])-float(r["pn_peri_shrink_fraction"]) for r in rows])

    fig,ax=plt.subplots(figsize=(7.2,5.2)); ax.plot(q,100*rms,marker="o")
    ax.set_xscale("log"); ax.set_xlabel("mass ratio q=m_A/m_B"); ax.set_ylabel("radial RMS / a [%]")
    ax.set_title("Transferability error vs mass ratio"); ax.grid(True,alpha=.25); save(fig,figs,"experiment04_rms_vs_q")

    fig,ax=plt.subplots(figsize=(7.2,5.2)); ax.plot(nu,dprec,marker="o")
    ax.set_xlabel("symmetric mass ratio nu"); ax.set_ylabel("S(n) - PN precession [rad/orbit]")
    ax.set_title("Precession residual vs symmetric mass ratio"); ax.grid(True,alpha=.25); save(fig,figs,"experiment04_precession_vs_nu")

    fig,ax=plt.subplots(figsize=(7.2,5.2)); ax.plot(nu,sh,marker="o")
    ax.set_xlabel("symmetric mass ratio nu"); ax.set_ylabel("S(n) - PN periapsis shrink fraction")
    ax.set_title("Dissipative shrink residual vs symmetric mass ratio"); ax.grid(True,alpha=.25); save(fig,figs,"experiment04_shrink_vs_nu")
    print("wrote 4 figure pairs")

if __name__=="__main__": main()
