#!/usr/bin/env python3
"""図：ξ、八状態（a0,b0 / a1,b1 / a2 / a3,b3）、a と e、δτ と ΔT、s_dir。SVG と PNG。"""
import csv, glob, gzip, os, sys
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as ex:
    print("matplotlib not available:", ex); sys.exit(0)

files = sorted(glob.glob(os.path.join("results", "run_*.csv"))) or sorted(glob.glob(os.path.join("results", "run_*.csv.gz")))
for p in files:
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as fh:
        rows = list(csv.DictReader(fh))
    g = lambda k: [float(r[k]) for r in rows]
    n = g("node")
    lab = os.path.basename(p).replace(".csv.gz", "").replace(".csv", "")
    fig, ax = plt.subplots(7, 1, figsize=(9, 18), sharex=True)
    ax[0].plot(n, g("xi"), lw=0.6); ax[0].set_ylabel("xi")
    ax[1].plot(n, g("a0"), lw=0.6, label="a0 (ln amp A)"); ax[1].plot(n, g("b0"), lw=0.6, label="b0 (ln amp B)"); ax[1].legend(fontsize=7); ax[1].set_ylabel("length state")
    ax[2].plot(n, g("a1"), lw=0.4, label="a1 (clock face A)"); ax[2].plot(n, g("b1"), lw=0.4, label="b1 (clock face B)"); ax[2].legend(fontsize=7); ax[2].set_ylabel("clock faces (state)")
    ax[3].plot(n, g("a2"), lw=0.6, label="a2 = ln|g| (length increment)"); ax[3].legend(fontsize=7); ax[3].set_ylabel("increment state")
    ax[4].plot(n, g("a3"), lw=0.6, label="a3 (clock increment A)"); ax4b = ax[4].twinx(); ax4b.plot(n, g("b3"), lw=0.6, color="C1", label="b3"); ax[4].set_ylabel("a3"); ax4b.set_ylabel("b3")
    ax[5].plot(n, g("a_orb"), lw=0.6); ax5b = ax[5].twinx(); ax5b.plot(n, g("ecc"), lw=0.6, color="C3"); ax[5].set_ylabel("a (semi-major)"); ax5b.set_ylabel("e")
    ax[6].plot(n, [max(v, 1e-20) for v in g("s_dir")], lw=0.4); ax[6].set_yscale("log"); ax[6].set_ylabel("s_dir"); ax[6].set_xlabel("node")
    for a_ in ax: a_.grid(alpha=0.3)
    fig.suptitle(f"eight-state v2 (clock faces in state), theta0 = 23pi/124: {lab}")
    fig.tight_layout()
    fig.savefig(f"results/fig_{lab}.svg"); fig.savefig(f"results/fig_{lab}.png", dpi=120); plt.close(fig)
    print("saved", f"results/fig_{lab}.svg/.png")
