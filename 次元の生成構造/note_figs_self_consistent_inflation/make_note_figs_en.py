# -*- coding: utf-8 -*-
"""English-labelled versions of the note figures (same data as make_note_figs.py)."""
import os, sys, csv, zipfile, io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1]; ZIPDIR = sys.argv[2]
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.titlesize": 15, "axes.titleweight": "bold",
                     "axes.labelsize": 13, "mathtext.fontset": "dejavusans", "axes.spines.top": False, "axes.spines.right": False})
RED, BLUE, GREEN, GRAY, NAVY, ORANGE = "#d7263d", "#1f5fd8", "#2e8b57", "#777777", "#12245e", "#e8871e"

def read(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for k in rows[0]:
        vals = [r[k] for r in rows]
        try: out[k] = np.array([float(v) for v in vals])
        except ValueError: out[k] = vals
    return out

# ---- Fig 1: inflation, log scale, N=5 and N=16
z = zipfile.ZipFile(os.path.join(ZIPDIR, "complex_simplex_decompactification_N5_N16_20260826.zip"))
def load(n):
    rows = list(csv.DictReader(io.StringIO(z.read(n).decode())))
    return {k: np.array([float(r[k]) for r in rows]) for k in ("step", "H_perp", "H_parallel", "H_total")}
d5, d16 = load("results/N5_geometry_summary.csv"), load("results/N16_geometry_summary.csv")
fig, ax = plt.subplots(figsize=(10, 6.2), dpi=170)
ax.semilogy(d16["step"], d16["H_perp"], color=RED, lw=2.4, label=r"N = 16: amplitude$^2$ in the new directions, $H_\perp$")
ax.semilogy(d5["step"], d5["H_perp"], color=ORANGE, lw=2.0, label=r"N = 5: amplitude$^2$ in the new directions, $H_\perp$")
ax.semilogy(d16["step"], d16["H_parallel"], color=BLUE, lw=1.8, alpha=0.9, label=r"N = 16: amplitude$^2$ in the original plane, $H_\parallel$")
ax.set_ylim(1e-34, 5); ax.set_xlim(0, 800)
h0, hf = d16["H_perp"][0], d16["H_perp"][-1]
ax.axhline(h0, color=GRAY, lw=0.8, ls=":"); ax.axhline(hf, color=GRAY, lw=0.8, ls=":")
ax.annotate("", xy=(520, hf), xytext=(520, h0), arrowprops=dict(arrowstyle="<->", color=NAVY, lw=2.2))
ax.text(540, 1e-15, "about 31 orders\nof magnitude\n($10^{31}$ times)", color=NAVY, fontsize=15, weight="bold", va="center")
ax.text(200, hf * 0.12, f"N = 16 final {hf:.2f}  (N = 5: 0.4)", fontsize=11, color=RED, va="top")
ax.text(795, h0 * 0.25, "start $10^{-32}$ (N = 16) / $10^{-31}$ (N = 5)", ha="right", fontsize=11, color=RED)
ax.text(60, 1e-8, "straight line on a log scale\n= exponential growth at a constant rate", fontsize=10, color=RED, rotation=66, ha="center", va="center")
ax.set_xlabel("step"); ax.set_ylabel("squared amplitude (total = 1, log scale)")
ax.set_title("Inflation: the amplitude in the new directions climbs 31 orders of magnitude")
ax.legend(loc="lower right", frameon=False, bbox_to_anchor=(0.99, 0.08), fontsize=10); fig.tight_layout()
fig.savefig(os.path.join(OUT, "note_fig0_inflation_log_en.png")); plt.close(fig)

# ---- Fig 2: transfer (linear)
d = read(os.path.join(HERE, "pump_depletion_timeseries.csv"))
fig, ax = plt.subplots(figsize=(10, 5.6), dpi=170)
ax.plot(d["step"], d["H_parallel"], color=BLUE, lw=2.4, label=r"component in the original plane, $H_\parallel$")
ax.plot(d["step"], d["H_perp"], color=RED, lw=2.4, label=r"component moved to the new directions, $H_\perp$")
ax.plot(d["step"], d["H_total"], color=GREEN, lw=2.0, ls="--", label=r"total $H_\parallel + H_\perp$ (exactly constant)")
ax.axvspan(0, 449, color="#fff1f1", zorder=0); ax.text(225, 0.5, "rapid expansion\n(to step 449)", ha="center", color=RED, fontsize=12)
ax.text(2600, 0.72, r"$H_\perp$ = 0.678", color=RED, fontsize=13); ax.text(2600, 0.27, r"$H_\parallel$ = 0.322", color=BLUE, fontsize=13)
ax.set_xlim(0, 5000); ax.set_ylim(-0.02, 1.08); ax.set_xlabel("step"); ax.set_ylabel("squared amplitude (total = 1)")
ax.set_title("What the expansion really is: a transfer — the total never grows")
ax.legend(loc="center right", frameon=False); fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig1_transfer_en.png")); plt.close(fig)

# ---- Fig 3: onset vs residual
eps = np.array([3.8728081613e-7, 1.8154355031e-9, 5.0849521984e-11, 2.3845913288e-13]); on = np.array([72, 134, 176, 238])
x = -np.log(eps); a, b = 11.616225, -99.563139
fig, ax = plt.subplots(figsize=(10, 5.6), dpi=170)
xx = np.linspace(13, 30, 10); ax.plot(xx, a * xx + b, color=GRAY, lw=1.5, ls="--", label=f"linear fit: slope {a:.3f}  (R² = 0.99999)")
ax.plot(x, on, "o", color=RED, ms=11, label="measured (4 runs)", zorder=3)
for xi, yi, e in zip(x, on, eps):
    ax.annotate(f"initial offset {e:.1e}\n→ onset at step {int(yi)}", (xi, yi), textcoords="offset points", xytext=(12, -28), fontsize=10, color="#333")
ax.text(13.5, 225, "prediction: slope = 1 / ln μ₁ = 11.593\n(from the eigenvalue in the next figure)", fontsize=12, color=NAVY,
        bbox=dict(boxstyle="round,pad=0.5", fc="#eef2ff", ec=NAVY))
ax.set_xlabel("smallness of the initial offset, −ln ε  (further right = smaller offset)"); ax.set_ylabel("step at which the expansion starts")
ax.set_title("The smaller the offset, the later the onset — precisely logarithmic"); ax.set_xlim(13, 30); ax.set_ylim(50, 260)
ax.legend(loc="lower right", frameon=False); fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig2_onset_residual_en.png")); plt.close(fig)

# ---- Fig 4: Floquet unit circle
fl = read(os.path.join(HERE, "floquet_spectrum.csv")); m = fl["fd_eps"] == 3e-6
re_, im_, mod = fl["eig_re"][m], fl["eig_im"][m], fl["modulus"][m]
fig, ax = plt.subplots(figsize=(8.6, 7.6), dpi=170)
th = np.linspace(0, 2 * np.pi, 400); ax.plot(np.cos(th), np.sin(th), color=GRAY, lw=1.2, ls="--", label="unit circle (on it: unchanged)")
ax.axhline(0, color="#cccccc", lw=0.8); ax.axvline(0, color="#cccccc", lw=0.8)
stable = mod < 1 - 1e-6; neutral = np.abs(mod - 1) <= 1e-6; unst = mod > 1 + 1e-6
ax.plot(re_[stable], im_[stable], "o", color=BLUE, ms=9, label="shrinking directions (|μ| < 1)")
ax.plot(re_[neutral], im_[neutral], "o", color=GRAY, ms=9, mfc="white", label="neutral directions (|μ| = 1)")
ax.plot(re_[unst], im_[unst], "o", color=RED, ms=13, label="growing directions (|μ| > 1) = the fall")
ax.annotate("μ₁ = 1.0901 (double)\nthe fastest-falling 2 dimensions", (1.0901, 0), xytext=(0.45, 0.55), fontsize=12, color=RED, arrowprops=dict(arrowstyle="->", color=RED))
ax.annotate("μ₂ = 1.0526 (double)", (1.0526, 0), xytext=(0.35, -0.55), fontsize=11, color=RED, arrowprops=dict(arrowstyle="->", color=RED))
ax.set_aspect("equal"); ax.set_xlim(-1.25, 1.35); ax.set_ylim(-1.2, 1.2)
ax.set_xlabel("real part"); ax.set_ylabel("imaginary part"); ax.set_title("Computing how the tilting wave falls, from the wave itself")
ax.legend(loc="lower left", fontsize=10, frameon=False); fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig3_floquet_circle_en.png")); plt.close(fig)

# ---- Fig 5: triple consistency
fig, ax = plt.subplots(figsize=(11, 6.2), dpi=170); ax.axis("off"); ax.set_xlim(0, 11); ax.set_ylim(0, 6.2)
def box(x, y, w, h, title, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.15,rounding_size=0.25", fc="white", ec=color, lw=2.2))
    ax.text(x + w / 2, y + h - 0.45, title, ha="center", va="center", fontsize=12.5, weight="bold", color=color)
    ax.text(x + w / 2, y + h / 2 - 0.35, body, ha="center", va="center", fontsize=12, color="#222", linespacing=1.6)
box(0.3, 3.7, 3.2, 2.2, "1. Measure the time evolution", "growth rate of the expansion\n0.172513 / step", BLUE)
box(0.3, 0.4, 3.2, 2.2, "2. Vary the initial offset", "slope of onset time\n11.616", GREEN)
box(7.3, 2.0, 3.4, 2.4, "3. Compute how fast it falls", "eigenvalue of the linearization\nμ₁ = 1.090086569", RED)
ax.add_patch(FancyBboxPatch((4.2, 1.9), 2.6, 2.6, boxstyle="round,pad=0.15,rounding_size=0.3", fc="#fff8e1", ec=NAVY, lw=2.5))
ax.text(5.5, 3.85, "one number", ha="center", fontsize=14, weight="bold", color=NAVY)
ax.text(5.5, 3.05, "2 ln μ₁ = 0.172514", ha="center", fontsize=13, color=BLUE)
ax.text(5.5, 2.45, "1 / ln μ₁ = 11.593", ha="center", fontsize=13, color=GREEN)
for (x0, y0, x1, y1, c) in [(3.5, 4.8, 4.25, 3.6, BLUE), (3.5, 1.5, 4.25, 2.7, GREEN), (7.3, 3.2, 6.8, 3.2, RED)]:
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=22, color=c, lw=2.2))
ax.text(5.5, 0.55, "Three independent measurements explained by one eigenvalue (agreement 0.0006 % and 0.2 %)", ha="center", fontsize=12, color="#333")
ax.set_title("Triple consistency — the onset is the linear instability of a tilting pure wave", fontsize=14.5, weight="bold", color=NAVY)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig4_triple_en.png")); plt.close(fig)

# ---- Fig 6: spectral entropy
e = read(os.path.join(HERE, "N5_spectral_entropy_timeseries.csv"))
fig, ax = plt.subplots(figsize=(10, 5.4), dpi=170)
ax.plot(e["step"], e["entropy_over_lnM"], color=NAVY, lw=2.4)
ax.axhline(1.0, color=GRAY, lw=1, ls="--"); ax.text(4950, 1.0008, "perfect equipartition = 1", ha="right", fontsize=11, color=GRAY)
ax.annotate("dips once during\nthe expansion (step 375)", (375, 0.97322), xytext=(1200, 0.978), fontsize=11, arrowprops=dict(arrowstyle="->", color="#333"))
ax.annotate("1.000000 at step 5000", (5000, 1.0), xytext=(3300, 0.9905), fontsize=11, arrowprops=dict(arrowstyle="->", color="#333"))
ax.set_xlim(0, 5000); ax.set_ylim(0.972, 1.003); ax.set_xlabel("step"); ax.set_ylabel("evenness of the amplitudes,  S / ln M")
ax.set_title("After it stops: the amplitudes of the 10 waves become perfectly equal")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig5_entropy_en.png")); plt.close(fig)

# ---- Fig 7: K5 and pyramid
CLASS = {"A+": ["12", "13", "45"], "A-": ["14", "15", "23"], "B+": ["24", "35"], "B-": ["25", "34"]}
EDGE = {e_: c for c, es in CLASS.items() for e_ in es}
STYLE = {"A+": (RED, "-"), "A-": (RED, "--"), "B+": (BLUE, "-"), "B-": (BLUE, "--")}
K5 = {"1": (0.0, 1.0), "2": (-0.95, 0.31), "3": (-0.59, -0.81), "4": (0.59, -0.81), "5": (0.95, 0.31)}
PYR = {"1": (0, 0, 1.25), "2": (-1, -1, 0), "4": (1, -1, 0), "3": (1, 1, 0), "5": (-1, 1, 0)}
def proj(p): x, y, z_ = p; return (x * 0.9 + 0.55 * y, z_ * 0.85 + 0.35 * y - 0.075 * x)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5.8), dpi=170)
for e_, c in EDGE.items():
    col, ls = STYLE[c]; p, q = K5[e_[0]], K5[e_[1]]
    a1.plot([p[0], q[0]], [p[1], q[1]], color=col, ls=ls, lw=2.6)
    a1.text((p[0] + q[0]) / 2 * 1.1, (p[1] + q[1]) / 2 * 1.1, e_, color=col, fontsize=10, ha="center", va="center", bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))
    P, Q = proj(PYR[e_[0]]), proj(PYR[e_[1]]); diag = e_ in ("23", "45")
    a2.plot([P[0], Q[0]], [P[1], Q[1]], color=col, ls=ls, lw=2.6 if not diag else 1.8, alpha=1 if not diag else 0.7)
for v, (x, y) in K5.items():
    a1.add_patch(plt.Circle((x, y), 0.12, fc="white", ec="black", lw=1.4, zorder=4)); a1.text(x, y, v, ha="center", va="center", fontsize=12, zorder=5)
for v, pt in PYR.items():
    x, y = proj(pt); a2.add_patch(plt.Circle((x, y), 0.12, fc=RED if v == "1" else BLUE, ec="black", lw=0.8, zorder=4)); a2.text(x, y, v, ha="center", va="center", fontsize=11, color="white", weight="bold", zorder=5)
h = [Line2D([0], [0], color=STYLE[c][0], ls=STYLE[c][1], lw=2.6, label=lbl) for c, lbl in [("A+", "A+ (3 edges)"), ("A-", "A− (3 edges) = sign flip of A+"), ("B+", "B+ (2 edges)"), ("B-", "B− (2 edges) = sign flip of B+")]]
a1.legend(handles=h, loc="lower center", fontsize=10, frameon=False, ncol=2, bbox_to_anchor=(0.5, -0.02))
for a in (a1, a2): a.set_aspect("equal"); a.axis("off")
a1.set_xlim(-1.3, 1.3); a1.set_ylim(-1.5, 1.25); a2.set_xlim(-1.9, 1.9); a2.set_ylim(-1.35, 1.6)
a1.set_title("The 10 relations of 5 bodies split into 4 classes", fontsize=13); a2.set_title("The same 10 read as a square pyramid\n(8 outer edges + 2 base diagonals)", fontsize=13)
fig.suptitle("The settled shape for N = 5:  3 + 3 + 2 + 2", fontsize=15, weight="bold", color=NAVY)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig7_k5_pyramid_en.png")); plt.close(fig)

# ---- Fig 8: seed sweep
sd = read(os.path.join(HERE, "N5_moduli_seed_sweep.csv"))
fig, ax = plt.subplots(figsize=(10, 5.2), dpi=170)
ax.axhline(0, color=GRAY, lw=1, ls="--"); ax.bar(sd["seed"], sd["relative_phase_mod_pi_rad"], color=BLUE, width=0.55)
for s_, v in zip(sd["seed"], sd["relative_phase_mod_pi_rad"]):
    ax.text(s_, v + (0.012 if v >= 0 else -0.03), f"{v:+.3f}", ha="center", fontsize=10, color="#333")
ax.text(0.98, 0.95, "All 8 runs: the same 4 classes (3+3+2+2), all amplitudes 0.1.\nOnly the relative phase between the two distance families differs.",
        transform=ax.transAxes, va="top", ha="right", fontsize=11.5, bbox=dict(boxstyle="round,pad=0.5", fc="#f5f6fc", ec=NAVY))
ax.set_xlabel("random seed (used to build the parent state)"); ax.set_ylabel("relative phase of the two distance families [rad]"); ax.set_ylim(-0.16, 0.42)
ax.set_title("What is fixed, and what is not"); fig.tight_layout(); fig.savefig(os.path.join(OUT, "note_fig6_seed_phase_en.png")); plt.close(fig)
print("done:", sorted(os.listdir(OUT)))
