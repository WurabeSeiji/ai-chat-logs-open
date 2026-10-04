#!/usr/bin/env python3
"""図：各走行について ξ(t)、配分 N_p・N_q と ⟨p|q⟩、a と e、読み (a) の s_dir と Δ。SVG と PNG。"""
import csv, glob, gzip, os, sys
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception as ex:
    print("matplotlib not available:", ex); sys.exit(0)

files = sorted(glob.glob(os.path.join("results", "run_*.csv"))) or sorted(glob.glob(os.path.join("results", "run_*.csv.gz")))
for p in files:
    opener = gzip.open if p.endswith(".gz") else open
    with opener(p, "rt") as fh:
        rows = list(csv.DictReader(fh))
    g = lambda k: [float(r[k]) for r in rows]
    n = g("node"); t = g("t")
    lab = os.path.basename(p).replace(".csv.gz", "").replace(".csv", "")
    fig, ax = plt.subplots(5, 1, figsize=(9, 14), sharex=True)
    ax[0].plot(n, g("xi"), lw=0.7); ax[0].set_ylabel("xi (distance)")
    ax[1].plot(n, g("N_p"), lw=0.7, label="N_p (current)"); ax[1].plot(n, g("N_q"), lw=0.7, label="N_q (rate)")
    ax[1].plot(n, [a + b for a, b in zip(g("N_p"), g("N_q"))], lw=0.7, label="N_p+N_q (= a, invariant)"); ax[1].legend(fontsize=7); ax[1].set_ylabel("distribution (norms)")
    ax[2].plot(n, g("Re_pq"), lw=0.7, label="Re<p|q> = L/(4w) (invariant)"); ax[2].plot(n, g("Im_pq"), lw=0.7, label="Im<p|q> = -(dxi/ds)/(4w) (rotating)"); ax[2].legend(fontsize=7); ax[2].set_ylabel("<p|q>")
    ax[3].plot(n, g("a_orb"), lw=0.7, label="a"); ax3b = ax[3].twinx(); ax3b.plot(n, g("ecc"), lw=0.7, color="C3", label="e"); ax[3].set_ylabel("a (semi-major)"); ax3b.set_ylabel("e")
    ax[4].plot(n, [max(v, 1e-20) for v in g("s_dir")], lw=0.5, label="s_dir (reading a)"); ax[4].set_yscale("log"); ax[4].set_ylabel("s_dir"); ax[4].set_xlabel("node")
    for a_ in ax: a_.grid(alpha=0.3)
    fig.suptitle(f"half-angle exchange: {lab}")
    fig.tight_layout()
    fig.savefig(f"results/fig_{lab}.svg"); fig.savefig(f"results/fig_{lab}.png", dpi=120); plt.close(fig)
    print("saved", f"results/fig_{lab}.svg/.png")
