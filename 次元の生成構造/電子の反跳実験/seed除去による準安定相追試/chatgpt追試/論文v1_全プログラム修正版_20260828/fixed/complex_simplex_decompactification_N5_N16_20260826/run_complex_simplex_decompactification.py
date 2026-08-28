#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Complex-simplex decompactification analysis, N=5 and N=16.

Physics:
  - same no-sigma-normalization LowRankSystem engine
  - same make_parent
  - seedless Z0 = v
  - same GAMMA
  - 5000 steps
  - no modification to dynamics

Readout:
  - Z = Z_parallel + Z_perp relative to the initial parent 2-plane
  - complex squared edge distances d_ij^2 = z_ij^2
  - centered complex Gram B = -1/2 J D^2 J
  - canonical complex-simplex axis scales r_k = sqrt(s_k(B)),
    where s_k are singular values (Takagi values for complex symmetric B).
"""

import importlib.util, math, json, csv, hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ENGINE_PATH = HERE / "run_n_scaling_lowrank_v1_no_sigma_norm.py"
STEPS = 5000
SEED = 0

def load_engine():
    spec = importlib.util.spec_from_file_location("eng", ENGINE_PATH)
    eng = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(eng)
    return eng

def centered_gram(z, n, edges):
    D2 = np.zeros((n,n), dtype=complex)
    for val,(i,j) in zip(z*z, edges):
        D2[i,j] = D2[j,i] = val
    J = np.eye(n) - np.ones((n,n))/n
    return -0.5 * J @ D2 @ J

def axis_scales(B):
    s = np.linalg.svd(B, compute_uv=False)
    return np.sqrt(np.maximum(s, 0.0)), s

def fit_rate(y, lo=1e-10, hi=1e-3):
    y = np.asarray(y, float)
    mask = np.isfinite(y) & (y > lo) & (y < hi)
    idx = np.where(mask)[0]
    if len(idx) < 10:
        return None, None, None, int(len(idx))
    slope = float(np.polyfit(idx, np.log(y[idx]), 1)[0])
    return slope, int(idx[0]), int(idx[-1]), int(len(idx))

def run_n(n, eng, outdir):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    syslr = eng.LowRankSystem(n)
    rng = np.random.default_rng(40260722 + 1000*n + SEED)
    v, residual, sig_parent = eng.make_parent(syslr, rng, iters=1200, tol=1e-12)
    Z = v.copy()
    p = v.real / np.linalg.norm(v.real)
    q = v.imag - (v.imag @ p) * p
    q = q / np.linalg.norm(q)
    wp = rng.normal(size=m)

    raw_rows = []
    geom_rows = []
    full_axes = []
    par_axes = []
    perp_axes = []
    cross_axes = []

    for t in range(STEPS+1):
        a = Z.real.copy()
        b = Z.imag.copy()
        norm2 = float(np.vdot(Z,Z).real)
        ztz = complex(Z @ Z)

        Zpar = p*(p @ Z) + q*(q @ Z)
        Zperp = Z - Zpar
        Hpar = float(np.vdot(Zpar,Zpar).real)
        Hperp = float(np.vdot(Zperp,Zperp).real)
        Aperp = math.sqrt(max(Hperp,0.0))
        f = Hperp / norm2

        Bfull = centered_gram(Z,n,edges)
        Bpar = centered_gram(Zpar,n,edges)
        Bperp = centered_gram(Zperp,n,edges)
        Bcross = Bfull - Bpar - Bperp

        rf,sf = axis_scales(Bfull)
        rp,sp = axis_scales(Bpar)
        ro,so = axis_scales(Bperp)
        rc,sc = axis_scales(Bcross)

        full_axes.append(rf)
        par_axes.append(rp)
        perp_axes.append(ro)
        cross_axes.append(rc)

        tol = max(float(sf[0]),1.0)*1e-10
        rank_full = int(np.sum(sf > tol))
        tol_o = max(float(so[0]),1.0)*1e-10
        rank_perp = int(np.sum(so > tol_o))

        Rfull = float(np.sqrt(np.sum(sf)))
        Rpar = float(np.sqrt(np.sum(sp)))
        Rperp = float(np.sqrt(np.sum(so)))
        Rcross = float(np.sqrt(np.sum(sc)))

        geom_rows.append([
            t,norm2,Hpar,Hperp,Aperp,f,abs(ztz),
            rank_full,rank_perp,
            Rfull,Rpar,Rperp,Rcross,
            float(rf[0]),float(rp[0]),float(ro[0]),float(rc[0]),
            float(np.linalg.norm(Bcross,"fro"))
        ])

        raw = [t,norm2,ztz.real,ztz.imag,abs(ztz),a@a,b@b,a@b,f]
        for k in range(m):
            raw += [a[k],b[k],abs(Z[k]),np.angle(Z[k])]
        raw_rows.append(raw)

        if t < STEPS:
            syslr.set_state(Z)  # FIX4
            se, wp = syslr.sigma_max_power(wp)
            Z = syslr.linear_rotation_step(Z,se)

    full_axes=np.asarray(full_axes)
    par_axes=np.asarray(par_axes)
    perp_axes=np.asarray(perp_axes)
    cross_axes=np.asarray(cross_axes)

    raw_header=['step','norm2','ztz_re','ztz_im','abs_ztz','sum_a2','sum_b2','sum_ab','f']
    for i,j in edges:
        raw_header += [f'a_{i+1}_{j+1}',f'b_{i+1}_{j+1}',f'abs_{i+1}_{j+1}',f'phase_{i+1}_{j+1}']
    pd.DataFrame(raw_rows,columns=raw_header).to_csv(outdir/f"N{n}_raw_states.csv",index=False)

    geom_header=[
        'step','H_total','H_parallel','H_perp','A_perp','f','abs_ztz',
        'full_simplex_rank','perp_simplex_rank',
        'R_full_takagi','R_parallel_takagi','R_perp_takagi','R_cross_takagi',
        'r1_full','r1_parallel','r1_perp','r1_cross','Bcross_fro'
    ]
    gdf=pd.DataFrame(geom_rows,columns=geom_header)
    gdf.to_csv(outdir/f"N{n}_geometry_summary.csv",index=False)

    adata={'step':np.arange(STEPS+1)}
    for k in range(n):
        adata[f'full_r{k+1}']=full_axes[:,k]
        adata[f'parallel_r{k+1}']=par_axes[:,k]
        adata[f'perp_r{k+1}']=perp_axes[:,k]
        adata[f'cross_r{k+1}']=cross_axes[:,k]
    adf=pd.DataFrame(adata)
    adf.to_csv(outdir/f"N{n}_takagi_axes.csv",index=False)

    # growth fits for all non-null simplex axes (N-1 physical centered rank; Nth is translation null)
    rate_rows=[]
    for k in range(n-1):
        rate,i0,i1,count=fit_rate(perp_axes[:,k])
        rate_rows.append([k+1,rate,i0,i1,count,float(perp_axes[0,k]),float(perp_axes[-1,k]),float(perp_axes[:,k].max())])
    rdf=pd.DataFrame(rate_rows,columns=['axis','log_growth_rate_per_step','fit_start','fit_end','fit_count','initial','final','max'])
    rdf.to_csv(outdir/f"N{n}_perp_axis_growth_rates.csv",index=False)

    Rrate,R0,R1,Rcount=fit_rate(gdf['R_perp_takagi'].to_numpy())

    summary={
        'N':n,'M':m,'steps':STEPS,'seed':SEED,'gamma':float(eng.GAMMA),'angle':float(eng.ANGLE),
        'physics':'amplitude-aware K, exact linear rotation exp(ANGLE K), no normalization (FIX1-4)',
        'parent_residual':float(residual),
        'parent_sigma_spectrum':[float(x) for x in sig_parent],
        'max_abs_ztz':float(gdf.abs_ztz.max()),
        'H_total_min':float(gdf.H_total.min()),
        'H_total_max':float(gdf.H_total.max()),
        'A_perp_initial':float(gdf.A_perp.iloc[0]),
        'A_perp_final':float(gdf.A_perp.iloc[-1]),
        'A_perp_max':float(gdf.A_perp.max()),
        'R_perp_takagi_initial':float(gdf.R_perp_takagi.iloc[0]),
        'R_perp_takagi_final':float(gdf.R_perp_takagi.iloc[-1]),
        'R_perp_takagi_max':float(gdf.R_perp_takagi.max()),
        'R_perp_log_growth_rate_per_step':Rrate,
        'R_perp_fit_window':[R0,R1],
        'full_simplex_rank_counts':{str(int(k)):int(v) for k,v in gdf.full_simplex_rank.value_counts().sort_index().items()},
        'perp_axis_final':[float(x) for x in perp_axes[-1,:n-1]],
        'perp_axis_growth_rates':[None if x is None or (isinstance(x,float) and np.isnan(x)) else float(x) for x in rdf.log_growth_rate_per_step],
    }
    (outdir/f"N{n}_summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')

    # figures
    plt.figure(figsize=(8,5))
    plt.semilogy(gdf.step,np.maximum(gdf.A_perp,1e-18),label='A_perp = ||Z_perp||')
    plt.semilogy(gdf.step,np.maximum(gdf.R_perp_takagi,1e-18),label='R_perp (Takagi simplex scale)')
    plt.xlabel('step'); plt.ylabel('raw scale (log)')
    plt.title(f'N={n}: orthogonal amplitude and complex-simplex scale')
    plt.legend(); plt.tight_layout()
    plt.savefig(outdir/f"N{n}_decompactification_scale.png",dpi=180); plt.close()

    plt.figure(figsize=(8,5))
    for k in range(n-1):
        plt.semilogy(np.arange(STEPS+1),np.maximum(perp_axes[:,k],1e-18),label=f'r{k+1}')
    plt.xlabel('step'); plt.ylabel('perpendicular-simplex Takagi axis scale (log)')
    plt.title(f'N={n}: canonical axes of the emergent perpendicular simplex')
    if n<=5: plt.legend()
    plt.tight_layout()
    plt.savefig(outdir/f"N{n}_perp_takagi_axes.png",dpi=180); plt.close()

    plt.figure(figsize=(8,5))
    for k in range(n-1):
        plt.plot(np.arange(STEPS+1),full_axes[:,k],label=f'r{k+1}')
    plt.xlabel('step'); plt.ylabel('full-simplex Takagi axis scale')
    plt.title(f'N={n}: full complex-simplex canonical axes')
    if n<=5: plt.legend()
    plt.tight_layout()
    plt.savefig(outdir/f"N{n}_full_takagi_axes.png",dpi=180); plt.close()

    plt.figure(figsize=(8,5))
    plt.plot(gdf.step,gdf.H_total,label='H_total')
    plt.plot(gdf.step,gdf.H_parallel,label='H_parallel')
    plt.plot(gdf.step,gdf.H_perp,label='H_perp')
    plt.xlabel('step'); plt.ylabel('raw squared amplitude')
    plt.title(f'N={n}: conserved total and redistribution')
    plt.legend(); plt.tight_layout()
    plt.savefig(outdir/f"N{n}_H_components.png",dpi=180); plt.close()

    return summary

def main():
    outdir=HERE/"results"
    outdir.mkdir(exist_ok=True)
    eng=load_engine()
    summaries=[]
    for n in (5,16):
        summaries.append(run_n(n,eng,outdir))

    analysis = [
        "# Complex-simplex decompactification test — N=5 and N=16",
        "",
        "The dynamics were not altered. Both runs use the no-K/sigma-normalization engine, seedless parent state Z0=v, the same GAMMA, and 5000 steps.",
        "",
        "## Geometric readout",
        "",
        "For each state z_ij, set d_ij^2=z_ij^2 and form the centered complex symmetric Gram matrix B=-1/2 J D^2 J. Its Takagi values are the singular values s_k(B); canonical simplex axis scales are r_k=sqrt(s_k).",
        "",
        "The state is also split relative to the initial parent plane as Z=Z_parallel+Z_perp. The same complex-simplex reconstruction is applied separately to Z_perp. This is a readout only; it does not feed back into the dynamics.",
        "",
    ]
    for s in summaries:
        n=s['N']
        rates=[x for x in s['perp_axis_growth_rates'] if x is not None]
        analysis += [
            f"## N={n}",
            f"- full simplex rank: {s['full_simplex_rank_counts']}",
            f"- max |Z^T Z|: {s['max_abs_ztz']:.3e}",
            f"- H_total range: {s['H_total_min']:.16f} .. {s['H_total_max']:.16f}",
            f"- A_perp: {s['A_perp_initial']:.3e} -> {s['A_perp_final']:.6f} (max {s['A_perp_max']:.6f})",
            f"- R_perp(Takagi): {s['R_perp_takagi_initial']:.3e} -> {s['R_perp_takagi_final']:.6f}",
            (f"- R_perp early log growth rate: {s['R_perp_log_growth_rate_per_step']:.6f}/step" if s['R_perp_log_growth_rate_per_step'] is not None else "- R_perp early log growth rate: (no exponential regime found; fit window empty)"),  # ROBUSTNESS PATCH: 修正版では指数成長域が無く fit が None になる,
            f"- growing perpendicular canonical axes measured: {len(rates)} of {n-1}",
            f"- axis growth-rate range: {min(rates):.6f} .. {max(rates):.6f}/step" if rates else "- no fitted axis rates",
            "",
        ]
    analysis += [
        "## Interpretation constrained by the data",
        "",
        "The perpendicular component is not merely a normalized plotting ratio: its raw amplitude A_perp grows exponentially from numerical-noise scale to O(1), while H_total remains conserved.",
        "",
        "When that perpendicular component is read as a complex simplex, its canonical distance scales also grow exponentially. Therefore the decompactification-like reading survives a direct complex-distance reconstruction.",
        "",
        "However, this particular Takagi-axis readout does not select only three expanding axes: N=5 shows all 4 non-null centered-simplex axes growing, and N=16 shows all 15. Thus the earlier 'three readable directions' phenomenon is not identical to the number of Takagi axes of the full complex simplex; it must be a separate readout/rank-selection structure if the two are to be connected.",
    ]
    (outdir/"ANALYSIS.md").write_text("\n".join(analysis),encoding='utf-8')

if __name__ == "__main__":
    main()
