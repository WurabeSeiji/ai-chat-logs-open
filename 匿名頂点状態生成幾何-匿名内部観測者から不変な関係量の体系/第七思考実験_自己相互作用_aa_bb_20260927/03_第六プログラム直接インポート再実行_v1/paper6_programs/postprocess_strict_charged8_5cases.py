
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
import json, math
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent
CASES=BASE/"cases"
ANABASE=Path("/mnt/data/charged_binary_integer_multiple_analytic_sweep")
CASE_IDS=["n1_attractive","n1_repulsive","n3_attractive","n3_repulsive","n4_attractive"]
HSTEP=2*math.pi/4000

def parts(cdir):
    m=json.loads((cdir/"part_manifest.json").read_text())
    return [cdir/p["file"] for p in m["parts"]]

def sample_raw(paths,total_rows,limit=8000):
    stride=max(1,total_rows//limit)
    out=[]
    for p in paths:
        with h5py.File(p,"r") as hf:
            ds=hf["raw_macro_trajectory"]
            for i in range(0,ds.shape[0],250_000):
                b=ds[i:min(i+250_000,ds.shape[0])]
                st=b[:,0].astype(np.int64)
                m=(st%stride)==0
                if np.any(m): out.append(b[m])
    with h5py.File(paths[-1],"r") as hf:
        last=hf["raw_macro_trajectory"][-1:]
    if not out or int(out[-1][-1,0])!=int(last[0,0]): out.append(last)
    return np.concatenate(out)

def save(fig,figdir,name):
    fig.tight_layout()
    fig.savefig(figdir/(name+".svg"))
    fig.savefig(figdir/(name+".png"),dpi=180)
    plt.close(fig)

def plot_case(cid):
    cdir=CASES/cid
    figdir=cdir/"figures"; figdir.mkdir(exist_ok=True)
    gen=json.loads((cdir/"generator_summary.json").read_text())
    cmp=json.loads((cdir/"reference_comparison.json").read_text())
    ref=pd.read_csv(ANABASE/cid/"charged_binary_analytic_reference_v1_raw.csv",
                    usecols=["r","t","phi","x_rel","y_rel"])
    rr=ref["r"].to_numpy(); tr=ref["t"].to_numpy(); pr=ref["phi"].to_numpy()
    xr=ref["x_rel"].to_numpy(); yr=ref["y_rel"].to_numpy()
    ir=np.linspace(0,len(ref)-1,min(8000,len(ref)),dtype=int)
    g=sample_raw(parts(cdir),gen["macro_steps"]+1)
    st=g[:,0]; r=g[:,1]; t=g[:,11]; x=g[:,12]; y=g[:,13]

    fig,ax=plt.subplots(figsize=(7,7))
    ax.plot(xr[ir],yr[ir],label="Analytic reference",linewidth=1.1)
    ax.plot(x,y,"--",label="Strict state generator",linewidth=.9)
    ax.set_xlabel("x / M"); ax.set_ylabel("y / M"); ax.set_aspect("equal",adjustable="box")
    ax.grid(True,alpha=.25); ax.legend(); ax.set_title(cid+": orbit overlay")
    save(fig,figdir,"figure01_orbit_overlay")

    fig,ax=plt.subplots(figsize=(8,5))
    ax.plot(tr[ir],rr[ir],label="Analytic reference")
    mt=np.isfinite(t)
    ax.plot(t[mt],r[mt],"--",label="Generator (finite t readout)")
    ax.set_xlabel("t / M"); ax.set_ylabel("r / M"); ax.grid(True,alpha=.25); ax.legend()
    ax.set_title(cid+": radius vs time")
    save(fig,figdir,"figure02_radius_vs_time")

    mt=np.isfinite(t)&(t<=tr[-1])
    fig,ax=plt.subplots(figsize=(8,5))
    if np.any(mt):
        ri=np.interp(t[mt],tr,rr); pi=np.interp(t[mt],tr,pr)
        er=np.abs(r[mt]-ri); exy=np.hypot(x[mt]-ri*np.cos(pi),y[mt]-ri*np.sin(pi))
        ax.semilogy(t[mt],np.maximum(er,1e-18),label="|Delta r| at same t")
        ax.semilogy(t[mt],np.maximum(exy,1e-18),label="||Delta(x,y)|| at same t")
    ax.set_xlabel("t / M"); ax.set_ylabel("absolute residual"); ax.grid(True,alpha=.25); ax.legend()
    ax.set_title(cid+": finite-time residuals")
    save(fig,figdir,"figure03_reference_residuals")

    fig,ax=plt.subplots(figsize=(8,5))
    ax.plot(st,g[:,8],label="N"); ax.plot(st,g[:,9],label="C"); ax.plot(st,g[:,10],label="D")
    ax.set_xlabel("macro step"); ax.set_ylabel("state value"); ax.grid(True,alpha=.25); ax.legend()
    ax.set_title(cid+": identity-propagated states")
    save(fig,figdir,"figure04_identity_states")

    m=pd.read_csv(cdir/"raw_first_macro_microsteps.csv")
    fig,ax=plt.subplots(figsize=(8,5))
    for nm,label in [("k1p","k1_P"),("k2p","k2_P"),("k3p","k3_P"),("k4p","k4_P"),("dP","dP")]:
        ax.plot(m["micro"],m[nm],marker="o",label=label)
    ax.set_xlabel("microstep"); ax.set_ylabel("work-state value"); ax.grid(True,alpha=.25); ax.legend()
    ax.set_title(cid+": first macrostep work states")
    save(fig,figdir,"figure05_microphase_work_states")

    ph=st*HSTEP; mp=ph<=pr[-1]+1e-12
    ri=np.interp(ph[mp],pr,rr)
    er=np.abs(r[mp]-ri); exy=np.hypot(x[mp]-ri*np.cos(ph[mp]),y[mp]-ri*np.sin(ph[mp]))
    fig,ax=plt.subplots(figsize=(8,5))
    ax.semilogy(st[mp]/4000,np.maximum(er,1e-18),label="|Delta r| at same phase")
    ax.semilogy(st[mp]/4000,np.maximum(exy,1e-18),label="||Delta(x,y)|| at same phase")
    ax.set_xlabel("orbital cycles"); ax.set_ylabel("absolute residual"); ax.grid(True,alpha=.25); ax.legend()
    ax.set_title(cid+": full-run phase-domain residuals")
    save(fig,figdir,"figure06_full_phase_residuals")

    (cdir/"FIGURE_INDEX_ja.md").write_text("""# 図一覧

1. `figure01_orbit_overlay` — 解析基準と厳格生成器の軌道重ね合わせ
2. `figure02_radius_vs_time` — 半径対時間
3. `figure03_reference_residuals` — 有限な時刻読み出し範囲の同一時刻残差
4. `figure04_identity_states` — N,C,D の同一写像による保持
5. `figure05_microphase_work_states` — 最初の macrostep の明示 work state
6. `figure06_full_phase_residuals` — 全走行の同一位相残差

各図について PNG と SVG の双方を生成する。Google Drive 保存対象はユーザー指示により SVG のみ。
""",encoding="utf-8")

def main():
    for cid in CASE_IDS:
        plot_case(cid)
        print(cid)

if __name__=="__main__":
    main()
