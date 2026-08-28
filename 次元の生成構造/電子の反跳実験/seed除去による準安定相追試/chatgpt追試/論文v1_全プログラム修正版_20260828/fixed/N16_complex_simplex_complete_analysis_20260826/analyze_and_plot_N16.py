#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
geom=pd.read_csv(HERE/"decompact_N16_geometry_summary.csv")
glob=pd.read_csv(HERE/"N16_global_summary.csv")
final=pd.read_csv(HERE/"N16_step5000_final_edges.csv")
axes=pd.read_csv(HERE/"N16_takagi_axes.csv")
clus=pd.read_csv(HERE/"N16_triplet_cluster_counts_selected_steps.csv")

glob["r2_mean"]=(glob.r2_min+glob.r2_max)/2
glob["r2_rel_spread"]=(glob.r2_max-glob.r2_min)/glob.r2_mean

plt.figure(figsize=(9,5))
plt.semilogy(geom.step,np.maximum(geom.H_perp,1e-30),label="H_perp")
plt.semilogy(glob.step,np.maximum(glob.r2_rel_spread,1e-16),label="edge |z|^2 relative spread")
plt.xlabel("step"); plt.ylabel("log scale")
plt.title("N=16: rapid decompactification and edge-amplitude equalization")
plt.legend(); plt.tight_layout()
plt.savefig(HERE/"N16_inflation_and_equalization.png",dpi=180); plt.close()

plt.figure(figsize=(9,5))
for k in range(15):
    plt.plot(axes.step,axes[f"takagi_r{k+1}"])
plt.xlabel("step"); plt.ylabel("Takagi axis scale")
plt.title("N=16: all 15 non-null complex-simplex axes")
plt.tight_layout()
plt.savefig(HERE/"N16_all_takagi_axes.png",dpi=180); plt.close()

theta=np.mod(final.theta.to_numpy(),np.pi)/np.pi
plt.figure(figsize=(8,5))
plt.hist(theta,bins=30)
plt.xlabel("theta / pi (mod 1)"); plt.ylabel("edge count")
plt.title("N=16 step 5000: edge phase distribution")
plt.tight_layout()
plt.savefig(HERE/"N16_final_phase_histogram.png",dpi=180); plt.close()

r2=final.r2.to_numpy()
plt.figure(figsize=(6,6))
plt.scatter(final.z2_re/r2,final.z2_im/r2,s=18)
plt.xlabel("Re(z^2)/|z|^2"); plt.ylabel("Im(z^2)/|z|^2")
plt.title("N=16 step 5000: normalized complex squared distances")
plt.gca().set_aspect("equal",adjustable="box")
plt.tight_layout()
plt.savefig(HERE/"N16_final_z2_phase_circle.png",dpi=180); plt.close()

plt.figure(figsize=(9,5))
for tol in [1e-2,1e-3,1e-4,1e-6]:
    d=clus[clus.tol==tol]
    plt.plot(d.step,d.triplet_cluster_count,label=f"tol={tol:g}")
plt.xlabel("step"); plt.ylabel("number of classes")
plt.title("N=16: (a2,b2,ab) class count")
plt.legend(); plt.tight_layout()
plt.savefig(HERE/"N16_class_count_evolution.png",dpi=180); plt.close()
