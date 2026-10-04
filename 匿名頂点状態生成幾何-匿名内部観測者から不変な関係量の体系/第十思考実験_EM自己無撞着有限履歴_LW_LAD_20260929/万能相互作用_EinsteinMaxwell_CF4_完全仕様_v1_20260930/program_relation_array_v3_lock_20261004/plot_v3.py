#!/usr/bin/env python3
"""results/lock_*.csv を図にする（ξ、2θ/2π、Δb、s_dir の節点変化）。SVG と PNG。"""
import csv, glob, os, sys
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as e:
    print("matplotlib not available:", e)
    sys.exit(0)

files = sorted(glob.glob(os.path.join("results", "lock_*.csv")))
if not files:
    sys.exit(0)
fig, axes = plt.subplots(4, 1, figsize=(9, 12), sharex=True)
for p in files:
    with open(p) as fh:
        rows = list(csv.DictReader(fh))
    n = [float(r["node"]) for r in rows]
    lab = os.path.basename(p).replace("lock_", "").replace(".csv", "")
    axes[0].plot(n, [float(r["xi"]) for r in rows], lw=0.8, label=lab)
    axes[1].plot(n, [float(r["two_theta_over_2pi"]) for r in rows], lw=0.8, label=lab)
    axes[2].plot(n, [float(r["Delta_b"]) for r in rows], lw=0.6, label=lab)
    axes[3].plot(n, [max(float(r["s_dir"]), 1e-20) for r in rows], lw=0.6, label=lab)
axes[0].set_ylabel("ξ")
axes[1].set_ylabel("2θ/2π (turn per node)")
axes[2].set_ylabel("Δb = 2θ − x3b (mod 2π)")
axes[3].set_ylabel("s_dir (source factor)")
axes[3].set_yscale("log")
axes[3].set_xlabel("node")
for ax in axes:
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=7)
fig.suptitle("relation array v3: three-phase integer relation, dynamic node clock, source stops at lock")
fig.tight_layout()
fig.savefig("results/lock_v3.svg")
fig.savefig("results/lock_v3.png", dpi=130)
print("saved results/lock_v3.svg, results/lock_v3.png")
