#!/usr/bin/env python3
"""Plot the saved raw CSV only; does not recompute the orbital solution."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def load_csv(path: Path):
    with path.open(newline='', encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    out={}
    for k in rows[0]:
        if k=='index':
            out[k]=np.array([int(r[k]) for r in rows])
        else:
            out[k]=np.array([float(r[k]) for r in rows], dtype=float)
    return out


def save_both(fig, outdir: Path, stem: str):
    fig.savefig(outdir/(stem+'.svg'), bbox_inches='tight')
    fig.savefig(outdir/(stem+'.png'), dpi=220, bbox_inches='tight')
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input', type=Path, default=Path(__file__).resolve().parent/'charged_binary_analytic_reference_v1_raw.csv')
    ap.add_argument('--outdir', type=Path, default=Path(__file__).resolve().parent/'figures')
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    d=load_csv(args.input)

    # 1 Relative orbit spiral.
    fig,ax=plt.subplots(figsize=(7,7))
    ax.plot(d['x_rel'],d['y_rel'],lw=0.75)
    ax.scatter([d['x_rel'][0],d['x_rel'][-1]],[d['y_rel'][0],d['y_rel'][-1]],s=22)
    ax.set_aspect('equal','box'); ax.set_xlabel('relative x / M'); ax.set_ylabel('relative y / M')
    ax.set_title('Charged binary analytic reference: relative inspiral')
    ax.grid(True,alpha=0.25)
    save_both(fig,args.outdir,'figure01_relative_orbit_spiral')

    # 2 Two body COM trajectories.
    fig,ax=plt.subplots(figsize=(7,7))
    ax.plot(d['x_A'],d['y_A'],lw=0.75,label='body A')
    ax.plot(d['x_B'],d['y_B'],lw=0.75,label='body B')
    ax.scatter([0],[0],s=18,label='center of mass')
    ax.set_aspect('equal','box'); ax.set_xlabel('x / M'); ax.set_ylabel('y / M')
    ax.set_title('Center-of-mass trajectories')
    ax.legend(); ax.grid(True,alpha=0.25)
    save_both(fig,args.outdir,'figure02_two_body_com_orbits')

    # 3 Radius vs time.
    fig,ax=plt.subplots(figsize=(8,5))
    ax.plot(d['t'],d['r'])
    ax.set_xlabel('t / M'); ax.set_ylabel('r / M'); ax.set_title('Adiabatic separation decay')
    ax.grid(True,alpha=0.25)
    save_both(fig,args.outdir,'figure03_radius_vs_time')

    # 4 Radiation powers.
    fig,ax=plt.subplots(figsize=(8,5))
    ax.plot(d['t'],d['P_gw'],label='GW quadrupole power')
    ax.plot(d['t'],d['P_em'],label='EM dipole power')
    ax.plot(d['t'],d['P_total'],label='total',lw=1.2)
    ax.set_yscale('log'); ax.set_xlabel('t / M'); ax.set_ylabel('power (geometric units)')
    ax.set_title('Radiation power'); ax.legend(); ax.grid(True,alpha=0.25)
    save_both(fig,args.outdir,'figure04_radiation_power')

    # 5 Cumulative emitted energies.
    fig,ax=plt.subplots(figsize=(8,5))
    ax.plot(d['t'],d['E_gw_cumulative'],label='GW emitted energy')
    ax.plot(d['t'],d['E_em_cumulative'],label='EM emitted energy')
    ax.plot(d['t'],d['E_total_cumulative'],label='total emitted energy',lw=1.2)
    ax.plot(d['t'],d['E_balance_expected'],'--',label='orbital-energy loss check')
    ax.set_xlabel('t / M'); ax.set_ylabel('energy / M'); ax.set_title('Cumulative radiated energy')
    ax.legend(); ax.grid(True,alpha=0.25)
    save_both(fig,args.outdir,'figure05_cumulative_radiated_energy')

    # 6 EM/GW power ratio against radius; crossing is visible at ratio=1.
    fig,ax=plt.subplots(figsize=(8,5))
    ax.plot(d['r'],d['P_em_over_P_gw'])
    ax.axhline(1.0,ls='--',lw=1.0)
    ax.invert_xaxis(); ax.set_xlabel('r / M (inspiral direction ->)'); ax.set_ylabel('P_EM / P_GW')
    ax.set_title('Competition between EM dipole and GW quadrupole radiation')
    ax.grid(True,alpha=0.25)
    save_both(fig,args.outdir,'figure06_em_to_gw_power_ratio')

    # 7 1PN conservative diagnostic, kept separate from LO trajectory generation.
    fig,ax=plt.subplots(figsize=(8,5))
    ax.plot(d['r'],100*d['rel_r_1pn_correction'],label='radius correction at fixed omega')
    ax.plot(d['r'],100*d['rel_omega_1pn_correction'],label='omega correction at fixed r')
    ax.invert_xaxis(); ax.set_xlabel('LO r / M (inspiral direction ->)'); ax.set_ylabel('relative diagnostic correction (%)')
    ax.set_title('Charged 1PN Kepler-law diagnostic (not fed back)')
    ax.legend(); ax.grid(True,alpha=0.25)
    save_both(fig,args.outdir,'figure07_1pn_kepler_diagnostic')

    # 8 Radiation frequencies.
    fig,ax=plt.subplots(figsize=(8,5))
    ax.plot(d['t'],d['f_em_dipole'],label='EM dipole frequency = orbital frequency')
    ax.plot(d['t'],d['f_gw_quadrupole'],label='GW quadrupole frequency = 2 x orbital')
    ax.set_xlabel('t / M'); ax.set_ylabel('frequency (1/M)'); ax.set_title('Dominant radiation frequencies')
    ax.legend(); ax.grid(True,alpha=0.25)
    save_both(fig,args.outdir,'figure08_radiation_frequencies')

    print(f'created figures in {args.outdir}')

if __name__=='__main__':
    main()
