#!/usr/bin/env python3
"""図：各走行の ξ、八状態（a0..a3、b0..b3）、|G| と arg G、a と e、δτ、読み (a) の s_dir。SVG と PNG。"""
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
    fig, ax = plt.subplots(6, 1, figsize=(9, 16), sharex=True)
    ax[0].plot(n, g("xi"), lw=0.6); ax[0].set_ylabel("xi (distance)")
    for k in ("a0", "b0"): ax[1].plot(n, g(k), lw=0.6, label=k + " (ln amplitude)")
    ax[1].legend(fontsize=7); ax[1].set_ylabel("state: log amplitude")
    for k in ("a2", "a3"): ax[2].plot(n, g(k), lw=0.6, label=k + (" (ln|G|)" if k == "a2" else " (arg G)"))
    ax[2].legend(fontsize=7); ax[2].set_ylabel("state: increment")
    ax[3].plot(n, g("a_orb"), lw=0.6, label="a"); ax3b = ax[3].twinx(); ax3b.plot(n, g("ecc"), lw=0.6, color="C3"); ax[3].set_ylabel("a (semi-major)"); ax3b.set_ylabel("e")
    ax[4].plot(n, g("dtau"), lw=0.6, label="dtau (node pitch, aux state)"); ax4b = ax[4].twinx(); ax4b.plot(n, g("dT_node"), lw=0.6, color="C2"); ax[4].set_ylabel("dtau"); ax4b.set_ylabel("dT per node")
    ax[5].plot(n, [max(v, 1e-20) for v in g("s_dir")], lw=0.5); ax[5].set_yscale("log"); ax[5].set_ylabel("s_dir (reading a)"); ax[5].set_xlabel("node")
    for a_ in ax: a_.grid(alpha=0.3)
    fig.suptitle(f"eight-state exchange, theta0 = 23pi/124: {lab}")
    fig.tight_layout()
    fig.savefig(f"results/fig_{lab}.svg"); fig.savefig(f"results/fig_{lab}.png", dpi=120); plt.close(fig)
    print("saved", f"results/fig_{lab}.svg/.png")
