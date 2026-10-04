#!/usr/bin/env python3
"""results/lock_*.csv を図にする（2θ/2π、|H|²、m、x3a の節点変化）。SVG と PNG。"""
import csv, glob, os, sys

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as e:  # 図は補助。無ければ事実だけ残す
    print("matplotlib not available:", e)
    sys.exit(0)

def load(p):
    with open(p) as fh:
        rows = list(csv.DictReader(fh))
    def col(k):
        return [float(r[k]) for r in rows]
    return col("node"), col("two_theta_over_2pi"), col("H2"), col("m"), col("x3a")

files = sorted(glob.glob(os.path.join("results", "lock_*.csv")))
if not files:
    sys.exit(0)
fig, axes = plt.subplots(4, 1, figsize=(9, 12), sharex=True)
for p in files:
    n, frac, H2, m, x3a = load(p)
    lab = os.path.basename(p).replace("lock_", "").replace(".csv", "")
    axes[0].plot(n, frac, lw=0.8, label=lab)
    axes[1].plot(n, H2, lw=0.6, label=lab)
    axes[2].plot(n, m, lw=0.8, label=lab)
    axes[3].plot(n, x3a, lw=0.8, label=lab)
axes[0].set_ylabel("2θ/2π (turn per node)")
axes[1].set_ylabel("|H|² (finite-history direction)")
axes[1].set_yscale("log")
axes[2].set_ylabel("m")
axes[3].set_ylabel("x3a (node clock)")
axes[3].set_xlabel("node")
for ax in axes:
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=7, ncol=2)
fig.suptitle("relation array v2: lock-in patterns A (clock) / B (orbit), windows N")
fig.tight_layout()
os.makedirs("results", exist_ok=True)
fig.savefig("results/lockin_v2.svg")
fig.savefig("results/lockin_v2.png", dpi=130)
print("saved results/lockin_v2.svg, results/lockin_v2.png")
