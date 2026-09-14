#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
base=Path(__file__).resolve().parent
s=pd.read_csv(base/"highN_N2_N40_summary.csv")
fig,ax=plt.subplots(figsize=(8,5))
ax.plot(s["N"],s["sample_mean_kappa"],marker="o",markersize=3,label="sample mean")
ax.plot(s["N"],s["exact_mean_kappa"],label="exact mean")
ax.set_xlabel("N"); ax.set_ylabel("mean kappa"); ax.set_title("Generic complex-Gaussian null: mean kappa, N=2..40")
ax.legend(); ax.grid(True,alpha=0.25); fig.tight_layout(); fig.savefig(base/"mean_kappa_vs_N.png",dpi=160); fig.savefig(base/"mean_kappa_vs_N.svg"); plt.close(fig)
fig,ax=plt.subplots(figsize=(8,5))
ax.semilogy(s["N"],s["max_reconstruction_error"].clip(lower=1e-18),marker="o",markersize=3)
ax.set_xlabel("N"); ax.set_ylabel("max reconstruction error (log scale)"); ax.set_title("Numerical conditioning of Newton/root reconstruction")
ax.grid(True,alpha=0.25); fig.tight_layout(); fig.savefig(base/"reconstruction_conditioning_vs_N.png",dpi=160); fig.savefig(base/"reconstruction_conditioning_vs_N.svg"); plt.close(fig)
