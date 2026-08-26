# -*- coding: utf-8 -*-
"""English re-creation of the three Japanese N=5 infographics.
All numbers are taken from the analysis package CSVs (step 5000)."""
import os, sys, math
import csv
import numpy as np
import matplotlib

def read_csv(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for k in rows[0].keys():
        vals = [r[k] for r in rows]
        try:
            out[k] = np.array([float(v) for v in vals])
        except ValueError:
            out[k] = vals
    return out
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10.5,
                     "axes.titleweight": "bold", "mathtext.fontset": "dejavusans"})

# ------------------------------------------------------------------ data
cls = read_csv(os.path.join(HERE, "N5_step5000_final_complex_distance_classes.csv"))
_g  = read_csv(os.path.join(HERE, "N5_step5000_four_group_summary.csv"))
ts  = read_csv(os.path.join(HERE, "N5_inflation_vs_ordering_timeseries.csv"))
class _Grp:  # minimal .loc[i, col] accessor
    class _Loc:
        def __getitem__(self, key): return _g[key[1]][key[0]]
    loc = _Loc()
grp = _Grp()

CLASS = {"A+": ["12", "13", "45"], "A-": ["14", "15", "23"], "B+": ["24", "35"], "B-": ["25", "34"]}
EDGE_CLASS = {e: c for c, es in CLASS.items() for e in es}
PH = {"A+": 0.321107, "A-": 0.821107, "B+": 0.348702, "B-": 0.848702}   # theta_ij / pi (mod 1)
DELTA = 0.027595
RED, BLUE = "#e8003d", "#1f4fd8"
STYLE = {"A+": (RED, "-"), "A-": (RED, "--"), "B+": (BLUE, "-"), "B-": (BLUE, "--")}
NAVY = "#12245e"

def lab(c):  # class label with subscript
    return {"A+": r"$A_+$", "A-": r"$A_-$", "B+": r"$B_+$", "B-": r"$B_-$"}[c]

# ------------------------------------------------------------------ helpers
def panel(ax, title=None, color="#555555", lw=0.8):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_edgecolor(color); s.set_linewidth(lw)
    if title:
        ax.set_title(title, loc="left", color=NAVY, pad=4)

def textbox(ax, lines, x=0.03, y=0.95, fs=8.6, color="black", box=None, va="top", lsp=1.45):
    t = ax.text(x, y, "\n".join(lines), transform=ax.transAxes, fontsize=fs, va=va, ha="left",
                color=color, linespacing=lsp)
    if box:
        t.set_bbox(dict(boxstyle="round,pad=0.5", fc=box, ec="#8890c0", lw=0.8))
    return t

def table(ax, header, rows, colw, bbox=(0.02, 0.02, 0.96, 0.96), fs=8.2, rowcolors=None, hcolor="#eef0f8"):
    tb = ax.table(cellText=rows, colLabels=header, colWidths=colw, bbox=bbox, cellLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(fs)
    for (r, c), cell in tb.get_celld().items():
        cell.set_edgecolor("#777777"); cell.set_linewidth(0.6)
        if r == 0:
            cell.set_facecolor(hcolor); cell.set_text_props(weight="bold")
        elif rowcolors and rowcolors[r - 1]:
            cell.set_text_props(color=rowcolors[r - 1])
    return tb

K5POS = {"1": (0.0, 1.0), "2": (-0.95, 0.31), "3": (-0.59, -0.81), "4": (0.59, -0.81), "5": (0.95, 0.31)}

def draw_k5(ax, labels=True, node_r=0.11):
    for e, c in EDGE_CLASS.items():
        col, ls = STYLE[c]
        p, q = K5POS[e[0]], K5POS[e[1]]
        ax.plot([p[0], q[0]], [p[1], q[1]], color=col, ls=ls, lw=2.0, zorder=1)
        if labels:
            mx, my = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
            # push label slightly away from centre to reduce overlap
            ax.text(mx * 1.1, my * 1.1, e, color=col, fontsize=8.5, ha="center", va="center", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
    for v, (x, y) in K5POS.items():
        ax.add_patch(plt.Circle((x, y), node_r, fc="white", ec="black", lw=1.2, zorder=4))
        ax.text(x, y, v, ha="center", va="center", fontsize=9.5, zorder=5)
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-2.05, 1.25); ax.set_aspect("equal", adjustable="datalim")

def k5_legend(ax, loc="upper right", fs=8):
    h = [Line2D([0], [0], color=STYLE[c][0], ls=STYLE[c][1], lw=2, label=f"{lab(c)}  ({len(CLASS[c])} edges)") for c in CLASS]
    ax.legend(handles=h, loc=loc, fontsize=fs, title="edge colour (class)", title_fontsize=8, framealpha=0.95)

# square pyramid: apex 1, base 4-cycle 2-4-3-5
PYR = {"1": (0, 0, 1.25), "2": (-1, -1, 0), "4": (1, -1, 0), "3": (1, 1, 0), "5": (-1, 1, 0)}

def project(p, view):
    x, y, z = p
    if view == "top":      return (x, y)
    if view == "front":    return (x + 0.18 * y, z + 0.12 * y)
    if view == "oblique":  return (x * 0.9 + 0.55 * y, z * 0.85 + 0.35 * y - 0.25 * x * 0.3)
    return (x, z)

def draw_pyramid(ax, view, labels=True, lw=1.9, title=None):
    for e, c in EDGE_CLASS.items():
        col, ls = STYLE[c]
        p, q = project(PYR[e[0]], view), project(PYR[e[1]], view)
        base_diag = e in ("23", "45")
        ax.plot([p[0], q[0]], [p[1], q[1]], color=col, ls=ls, lw=lw * (0.8 if base_diag else 1), alpha=0.75 if base_diag else 1, zorder=1)
        if labels:
            ax.text((p[0] + q[0]) / 2, (p[1] + q[1]) / 2, e, color=col, fontsize=7.5, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8), zorder=3)
    for v, pt in PYR.items():
        x, y = project(pt, view)
        ax.add_patch(plt.Circle((x, y), 0.11, fc="#c0392b" if v == "1" else "#2c5fd8", ec="black", lw=0.8, zorder=4))
        ax.text(x, y, v, ha="center", va="center", fontsize=8, color="white", zorder=5, weight="bold")
    pts = np.array([project(p, view) for p in PYR.values()])
    cx, cy = pts.mean(axis=0); span = max(np.ptp(pts[:, 0]), np.ptp(pts[:, 1])) * 0.62
    ax.set_xlim(cx - span, cx + span); ax.set_ylim(cy - span, cy + span)
    ax.set_aspect("equal", adjustable="datalim")
    if title: ax.set_title(title, fontsize=8.5, color=NAVY, pad=2)

def zplane(ax, title=None, fs=8):
    r = math.sqrt(0.1)
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(r * np.cos(th), r * np.sin(th), color="#999999", lw=0.8)
    ax.axhline(0, color="#bbbbbb", lw=0.6); ax.axvline(0, color="#bbbbbb", lw=0.6)
    for k in range(len(cls["edge"])):
        c = EDGE_CLASS[cls["edge"][k].replace("-", "")]
        col = STYLE[c][0]
        ax.plot(cls["a"][k], cls["b"][k], "o", color=col, ms=7, mfc=col if "+" in c else "white", mew=1.6, zorder=3)
    for i, (c, ang) in enumerate([("A+", 0.321107), ("A-", 0.821107), ("B+", 0.348702), ("B-", 0.848702)]):
        ax.plot([-0.70], [0.44 - i * 0.075], "o", color=STYLE[c][0], ms=6, mfc=STYLE[c][0] if "+" in c else "white", mew=1.5, clip_on=False)
        ax.text(-0.66, 0.44 - i * 0.075, f"{lab(c)} ({len(CLASS[c])})  θ/π ≈ {ang:.4f}", color=STYLE[c][0], fontsize=fs, va="center")
    ax.text(r * 1.04, 0.012, "Re z", fontsize=7, color="#666"); ax.text(0.012, r * 1.06, "Im z", fontsize=7, color="#666")
    ax.text(-0.7, -0.46, r"radius $\sqrt{0.1}\simeq0.316228$;  $z$ and $-z$ give the same $z^2$",
            fontsize=6.6, color="#666")
    ax.set_xlim(-0.74, 0.46); ax.set_ylim(-0.5, 0.5); ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([]); ax.set_yticks([])
    if title: ax.set_title(title, loc="left", color=NAVY)

# common tables --------------------------------------------------------
GROUP_ROWS = [[lab(c), f"{PH[c]:.6f}" + ("" if "+" in c else f"  (= {lab(c[0]+'+')} + 0.5)"), str(len(CLASS[c])),
               ", ".join(CLASS[c]), ("+" if "+" in c else "−") + c[0]] for c in CLASS]
GROUP_HDR = ["class", "phase θᵢⱼ/π (mod 1)", "count", "edges", "sign of z²"]
GROUP_W = [0.11, 0.37, 0.11, 0.24, 0.17]
GROUP_COLORS = [RED, RED, BLUE, BLUE]
PHASE_ROWS = [["0", "8", "within the same class\n(3C2)+(3C2)+(2C2)+(2C2) = 8"],
              [f"{DELTA:.6f}\n(= δ)", "12", "A₊–B₊ : 6 pairs\nA₋–B₋ : 6 pairs"],
              [f"{0.5-DELTA:.6f}\n(= 0.5 − δ)", "12", "A₊–B₋ : 6 pairs\nA₋–B₊ : 6 pairs"],
              ["0.5", "13", "A₊–A₋ : 9 pairs\nB₊–B₋ : 4 pairs"],
              ["total", "45", "(10C2) = 45"]]
PHASE_HDR = ["|Δθ|/π (mod 1)", "count", "pairs involved"]

def family_boxes(ax, y=0.22):
    ax.text(0.03, y, "Family A (6 edges): {12, 13, 14, 15, 23, 45}\n= 4 edges from vertex 1 + (23, 45)\n→ contains the 2 base diagonals",
            transform=ax.transAxes, fontsize=7.0, color=RED, va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#fff0f4", ec=RED, lw=0.7))
    ax.text(0.545, y, "Family B (4 edges): {24, 25, 34, 35}\n= closed 4-cycle 2-4-3-5-2\n→ outer boundary of the base",
            transform=ax.transAxes, fontsize=7.0, color=BLUE, va="top",
            bbox=dict(boxstyle="round,pad=0.4", fc="#eef2ff", ec=BLUE, lw=0.7))

# ================================================================== FIGURE A
FACES = [("1-2-3", "{12, 13, 23}", "A (same-phase family)", RED), ("1-2-4", "{12, 14, 24}", "mixed (A–B)", BLUE),
         ("1-2-5", "{12, 15, 25}", "mixed (A–B)", BLUE), ("1-3-4", "{13, 14, 34}", "mixed (A–B)", BLUE),
         ("1-3-5", "{13, 15, 35}", "mixed (A–B)", BLUE), ("1-4-5", "{14, 15, 45}", "A (same-phase family)", RED),
         ("2-3-4", "{23, 24, 34}", "mixed (A–B)", BLUE), ("2-3-5", "{23, 25, 35}", "mixed (A–B)", BLUE),
         ("2-4-5", "{24, 25, 45}", "mixed (A–B)", BLUE), ("3-4-5", "{34, 35, 45}", "mixed (A–B)", BLUE)]

def figure_A():
    fig = plt.figure(figsize=(15.36, 10.4), dpi=200); fig.patch.set_facecolor("white")
    fig.text(0.012, 0.982, "N = 5: structure of the complex relation distances (step 5000)", fontsize=15, weight="bold", color=NAVY, va="top")
    fig.text(0.012, 0.952, "Phase structure engraved on the regular 4-simplex formed by the 10 relations", fontsize=10.5, color="#333", va="top")
    gs = fig.add_gridspec(nrows=3, ncols=3, left=0.012, right=0.988, top=0.905, bottom=0.015,
                          width_ratios=[1.05, 1.0, 1.0], height_ratios=[1.0, 0.78, 0.5], hspace=0.22, wspace=0.07)
    # 1 + 2
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "1. Basic data (step 5000)")
    textbox(ax, ["• number of vertices N = 5      number of relations M = 10",
                 r"• $|z_{ij}|^2 = a^2 + b^2 \simeq 0.1$  (equal for all 10 edges)",
                 "• Gram matrix rank = 4  →  the 5 points close in 4 dimensions"], y=0.95, fs=8.6)
    textbox(ax, ["All distance magnitudes are equal → outer shell of a regular 4-simplex (regular 5-cell)"], y=0.66, fs=8.0, color=NAVY, box="#eef0fb")
    ax2 = ax.inset_axes([0.0, 0.0, 1.0, 0.5]); ax2.axis("off")
    ax2.set_title("2. Four-class split of the 10 relations (phase θᵢⱼ and count)", loc="left", color=NAVY, fontsize=9.3, weight="bold")
    table(ax2, GROUP_HDR, GROUP_ROWS, GROUP_W, bbox=(0.02, 0.05, 0.96, 0.92), rowcolors=GROUP_COLORS, fs=7.8)
    # 3 K5
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "3. Edge colouring on K₅ (class layout)")
    draw_k5(ax); k5_legend(ax, loc="upper right", fs=7.3); family_boxes(ax)
    # 4 + 5
    ax = fig.add_subplot(gs[0, 2]); panel(ax, "4. Phase differences |Δθ|/π (45 pairs)")
    table(ax, PHASE_HDR, PHASE_ROWS, [0.3, 0.14, 0.56], bbox=(0.03, 0.42, 0.94, 0.55), fs=7.8)
    textbox(ax, ["5. Symmetry",
                 "• magnitudes |zᵢⱼ| alone: full S₅ symmetry (120 permutations)",
                 "• vertex permutations preserving the complex zᵢⱼ exactly: 2",
                 "   (the identity and the global reflection)",
                 "• allowing the global sign flip z² → −z² as one overall phase: 4"], y=0.35, fs=7.9)
    # notes
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "Notes on the four classes")
    textbox(ax, [r"$A_\pm$ :  $z^2 = \pm A$    ($|z^2| = 0.1000000001$)",
                 r"$B_\pm$ :  $z^2 = \pm B$    ($|z^2| = 0.1000000001$)",
                 f"phase offset between families A and B:  δ = {DELTA:.6f} π",
                 r"$A_-$ is $A_+$ shifted by π/2;  $B_-$ is $B_+$ shifted by π/2",
                 "",
                 "Per-class means at step 5000:",
                 f"  A₊: a² = {grp.loc[0,'mean_a2']:.6f},  b² = {grp.loc[0,'mean_b2']:.6f},  ab = {grp.loc[0,'mean_ab']:+.6f}",
                 f"  A₋: a² = {grp.loc[1,'mean_a2']:.6f},  b² = {grp.loc[1,'mean_b2']:.6f},  ab = {grp.loc[1,'mean_ab']:+.6f}",
                 f"  B₊: a² = {grp.loc[2,'mean_a2']:.6f},  b² = {grp.loc[2,'mean_b2']:.6f},  ab = {grp.loc[2,'mean_ab']:+.6f}",
                 f"  B₋: a² = {grp.loc[3,'mean_a2']:.6f},  b² = {grp.loc[3,'mean_b2']:.6f},  ab = {grp.loc[3,'mean_ab']:+.6f}",
                 "",
                 r"(A₊, A₋: $a^2$ and $b^2$ exchanged, sign of $ab$ flipped — likewise B₊, B₋)"], y=0.94, fs=8.0, lsp=1.5)
    # 6 projections
    sub = gs[1, 1].subgridspec(1, 3, wspace=0.05)
    for i, (v, t, lb) in enumerate([("oblique", "(1) oblique view", False), ("front", "(2) front view", False), ("top", "(3) top view (1 up)", True)]):
        ax = fig.add_subplot(sub[0, i]); panel(ax, None); draw_pyramid(ax, v, labels=lb, lw=1.7)
        ax.set_title(t, fontsize=8.2, color=NAVY, pad=2)
        if i == 0:
            ax.text(0.0, 1.13, "6. Symmetric 3D projection (three views)", transform=ax.transAxes, fontsize=10.5, weight="bold", color=NAVY)
            ax.text(0.03, 0.03, "all outer faces (triangles)\nare equilateral", transform=ax.transAxes, fontsize=6.9, color="#333", va="bottom")
    ax = fig.add_subplot(gs[1, 2]); panel(ax, "Triangular faces classified by edge colour")
    for i, (f, es, k, col) in enumerate(FACES):
        ax.text(0.04, 0.93 - i * 0.093, f"▲ {f} : {es}   →   {k}", transform=ax.transAxes, fontsize=8.2, color=col, va="center")
    # bottom
    ax = fig.add_subplot(gs[2, 0]); panel(ax, "Summary")
    textbox(ax, ["• The 5 points form a 4-dimensional simplex (rank = 4) and all 10 edge",
                 "   lengths are equal (regular outer shell).",
                 "• On top of that, the complex phases split into the four classes 3+3+2+2:",
                 "   A₊ / A₋ give 6 edges and B₊ / B₋ give 4 edges — an internal structure.",
                 "• The 4 edges of family B form a closed 4-cycle on the vertices {2, 3, 4, 5}."], y=0.92, fs=7.9, lsp=1.4)
    ax = fig.add_subplot(gs[2, 1]); panel(ax, "Essential structure")
    textbox(ax, ["• absolute distance |z| : fully symmetric (regular 4-simplex)",
                 "• complex phase : breaks the symmetry spontaneously and is organised",
                 "   into two basic phases ±A and ±B in an orthogonal (π/2) relation"], y=0.9, fs=8.0, color="#333", box="#fffbe6", lsp=1.4)
    ax = fig.add_subplot(gs[2, 2]); panel(ax, "Legend (edge colours)")
    for i, (c, t) in enumerate([("A+", "(phase θ)"), ("A-", "(θ + π/2)"), ("B+", "(θ + δ)"), ("B-", "(θ + δ + π/2)")]):
        col, ls = STYLE[c]
        ax.plot([0.06, 0.22], [0.82 - i * 0.2] * 2, color=col, ls=ls, lw=2.2, transform=ax.transAxes)
        ax.text(0.27, 0.82 - i * 0.2, f"{lab(c)}  {t}", transform=ax.transAxes, fontsize=8.5, va="center", color=col)
    fig.savefig(os.path.join(OUT, "N5_complex_distance_structure_overview_en.png"), dpi=200, facecolor="white"); plt.close(fig)

# ================================================================== FIGURE B
def figure_B():
    fig = plt.figure(figsize=(15.36, 10.4), dpi=200); fig.patch.set_facecolor("white")
    fig.add_artist(Rectangle((0.0, 0.935), 1.0, 0.065, transform=fig.transFigure, fc=NAVY, ec="none"))
    fig.text(0.012, 0.988, "N = 5 complex simplex: complete analysis set (final state at step 5000)", fontsize=13.5, weight="bold", color="white", va="top")
    fig.text(0.012, 0.958, "— structure, interpretation and time evolution —          date: 2026-08-26", fontsize=9, color="white", va="top")
    gs = fig.add_gridspec(nrows=3, ncols=3, left=0.035, right=0.988, top=0.905, bottom=0.05,
                          width_ratios=[1.05, 1.1, 1.0], height_ratios=[1.0, 0.85, 0.75], hspace=0.28, wspace=0.09)
    ax = fig.add_subplot(gs[0, 0]); panel(ax, "1. Summary of the final state (step 5000)")
    textbox(ax, ["• number of vertices N = 5,  number of relations M = 10",
                 r"• complex distance $z_{ij} = a_{ij} + i\,b_{ij}$",
                 r"• $|z_{ij}|^2 = a^2 + b^2 \simeq 0.1$  (equal for all 10 edges)",
                 "• Gram rank = 4  →  5 points close in 4 dimensions (4-simplex)",
                 "• the 10 edges split by phase into 4 classes: 3 + 3 + 2 + 2",
                 "• |Δθ|/π takes 4 discrete values: 0, δ, 0.5 − δ, 0.5",
                 f"     δ = {DELTA:.6f} π"], y=0.95, fs=8.3)
    ax2 = ax.inset_axes([0.0, 0.0, 1.0, 0.42]); ax2.axis("off")
    ax2.set_title("2. Four-class split of the 10 edges (phase θ/π of the distance)", loc="left", color=NAVY, fontsize=9.3, weight="bold")
    table(ax2, GROUP_HDR, GROUP_ROWS, GROUP_W, bbox=(0.02, 0.05, 0.96, 0.9), rowcolors=GROUP_COLORS, fs=7.6)
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "3. Edge colouring on K₅ (class layout)")
    draw_k5(ax); k5_legend(ax, loc="upper right", fs=7.3); family_boxes(ax)
    ax = fig.add_subplot(gs[0, 2]); panel(ax, "4. Phase-difference statistics |Δθ|/π (45 pairs)")
    table(ax, PHASE_HDR, PHASE_ROWS, [0.3, 0.14, 0.56], bbox=(0.03, 0.03, 0.94, 0.94), fs=7.9)
    ax = fig.add_subplot(gs[1, 0]); panel(ax, "5. Reading the 3D structure (why it looks 3D)")
    textbox(ax, ["• all absolute distances are equal",
                 "   → outer shell = regular 4-simplex (4-dimensional)",
                 "• but the phase structure splits the edges naturally into 8 + 2",
                 "• with vertex 1 as apex and {2, 4, 3, 5} as base this is a square pyramid:",
                 "   – outer edges (8): 12, 13, 14, 15, 24, 43, 35, 52",
                 "   – base diagonals (2): 23, 45",
                 "• hence: a 3D structure is engraved on the outer shell of the 4-simplex"], y=0.94, fs=8.1, lsp=1.55)
    sub = gs[1, 1].subgridspec(1, 3, wspace=0.05)
    for i, (v, t) in enumerate([("top", "(1) top view (1 up)"), ("front", "(2) front view"), ("oblique", "(3) oblique view")]):
        ax = fig.add_subplot(sub[0, i]); panel(ax, None); draw_pyramid(ax, v, labels=(v == "top"), lw=1.6)
        ax.set_title(t, fontsize=8, color=NAVY, pad=2)
        if i == 0:
            ax.text(0.0, 1.2, "6. 3D visualisation by projection (square-pyramid reading)", transform=ax.transAxes, fontsize=10.5, weight="bold", color=NAVY)
    ax = fig.add_subplot(gs[1, 2]); panel(ax, "8. The four classes in the complex plane"); zplane(ax)
    # 7 time evolution
    fig.text(0.035, 0.305, "7. Rapid expansion and metastable evolution (time-evolution view)", fontsize=10.5, weight="bold", color=NAVY, va="bottom")
    sub = gs[2, :].subgridspec(1, 3, wspace=0.2, width_ratios=[1, 1, 1.15])
    ax = fig.add_subplot(sub[0, 0])
    ax.loglog(np.maximum(ts["step"], 1), np.maximum(ts["H_perp"], 1e-30), color="#c0392b", lw=1.6)
    ax.axvspan(1, 449, color="#fde8e8", alpha=0.6); ax.axvspan(449, 6000, color="#e8f5e9", alpha=0.6)
    ax.axvline(430, color="#888", ls=":", lw=0.8); ax.axvline(449, color="#888", ls=":", lw=0.8)
    ax.text(2.5, 3e-2, "rapid expansion", color="#c0392b", fontsize=8); ax.text(650, 3e-2, "metastable stage", color="#2e7d32", fontsize=8)
    ax.text(1.6, 3e-7, "H⊥ reaches 95 % of its final value at step 430\nand 99 % at step 449 (dotted lines)", fontsize=6.8, color="#444")
    ax.set_ylim(1e-10, 3); ax.set_xlim(1, 6000); ax.set_xlabel("step (log axis)", fontsize=8); ax.tick_params(labelsize=7)
    ax.set_title("(1) transverse component H⊥ (outside the parent plane)", loc="left", color=NAVY, fontsize=8.8)
    ax = fig.add_subplot(sub[0, 1])
    ax.loglog(np.maximum(ts["step"], 1), np.maximum(ts["four_group_error"], 1e-12), color="#c0392b", lw=1.4, label="deviation from the four-class pattern")
    ax.axvline(449, color="#2e7d32", ls="--", lw=0.8); ax.text(500, 3e-1, "ordering continues", color="#2e7d32", fontsize=8)
    ax.set_ylim(1e-10, 3); ax.set_xlim(1, 6000); ax.set_xlabel("step (log axis)", fontsize=8); ax.tick_params(labelsize=7)
    ax.legend(fontsize=7, loc="lower left"); ax.set_title("(2) four-class agreement error", loc="left", color=NAVY, fontsize=8.8)
    ax = fig.add_subplot(sub[0, 2]); panel(ax, "(3) representative steps and indicators")
    textbox(ax, ["spatial scale (rapid expansion), H⊥ relative to its final value:",
                 "    50 %: step 357    90 %: step 395    95 %: step 430    99 %: step 449",
                 "phase ordering (four-class error stays below the threshold from this step on):",
                 "    < 10⁻¹ : 814      < 10⁻² : 1531     < 10⁻³ : 2099",
                 "    < 10⁻⁴ : 2627     < 10⁻⁶ : 3791     < 10⁻⁸ : 4923"], y=0.94, fs=7.6, lsp=1.5)
    textbox(ax, ["Stage 1 (rapid expansion): the transverse scale H⊥ grows exponentially and",
                 "   the overall size is fixed (essentially complete by step ≈ 400).",
                 "Stage 2 (ordering in the metastable stage): with the scale fixed, the internal",
                 "   phase/distance structure slowly self-organises into the stable 3+3+2+2",
                 "   pattern (completed over thousands of steps)."], y=0.4, fs=7.4, color="#333", box="#f1f3fb", lsp=1.4)
    fig.savefig(os.path.join(OUT, "N5_complete_simplex_analysis_overview_en.png"), dpi=200, facecolor="white"); plt.close(fig)

# ================================================================== FIGURE C
def figure_C():
    fig = plt.figure(figsize=(15.36, 10.4), dpi=200); fig.patch.set_facecolor("white")
    fig.text(0.012, 0.982, "N = 5   complex distance structure of the relation waves (final state at step 5000)", fontsize=14, weight="bold", color=NAVY, va="top")
    gs = fig.add_gridspec(nrows=3, ncols=3, left=0.012, right=0.988, top=0.935, bottom=0.015,
                          width_ratios=[1.25, 1.0, 1.0], height_ratios=[1.0, 0.85, 0.55], hspace=0.2, wspace=0.07)
    ax = fig.add_subplot(gs[0, 0]); panel(ax, None)
    textbox(ax, ["Basic results",
                 "• number of vertices N = 5,  number of relations M = 10",
                 r"• $|z_{ij}|^2 = a^2 + b^2 \simeq 0.1$  (equal for all 10 edges)",
                 "• Gram rank = 4  →  5 points in 4 dimensions (complex)",
                 "• the 10 edges split by phase into 4 classes: 3+3+2+2",
                 "• |Δθ|/π takes 4 discrete values: 0, δ, 0.5 − δ, 0.5",
                 f"   (δ = {DELTA:.6f}…)"], y=0.96, fs=8.0, color=NAVY, box="#f5f6fc", lsp=1.4)
    ax2 = ax.inset_axes([0.0, 0.0, 1.0, 0.52]); ax2.axis("off")
    ax2.set_title("1. Class assignment of the 10 edges (z = a + ib)", loc="left", color=NAVY, fontsize=9.3, weight="bold")
    rows = [[lab(c), f"{PH[c]:.6f}" + ("" if "+" in c else f" (= {c[0]}₊ + 0.5)"), "0.1",
             f"{grp.loc[i,'mean_a2']:.6f}", f"{grp.loc[i,'mean_b2']:.6f}", f"{grp.loc[i,'mean_ab']:+.6f}", str(len(CLASS[c])), "{" + ", ".join(CLASS[c]) + "}"]
            for i, c in enumerate(CLASS)]
    table(ax2, ["class", "phase θ/π (mod 1)", "|z|²", "a²", "b²", "ab", "count", "edges"], rows,
          [0.08, 0.27, 0.07, 0.12, 0.12, 0.13, 0.08, 0.18], bbox=(0.0, 0.3, 1.0, 0.66), rowcolors=GROUP_COLORS, fs=7.2)
    ax2.text(0.0, 0.24, "A₊, A₋ differ in phase by π/2  (a² and b² are exchanged, the sign of ab flips);  B₊, B₋ also differ by π/2\n"
                        f"families A (6 edges) and B (4 edges) differ by δ = {DELTA:.6f} π  or  π/2 − δ", transform=ax2.transAxes, fontsize=7.3, va="top", color="#333")
    ax = fig.add_subplot(gs[0, 1]); panel(ax, "2. Edge colouring on K₅ (class layout)")
    draw_k5(ax); k5_legend(ax, loc="upper right", fs=7.3); family_boxes(ax)
    ax = fig.add_subplot(gs[0, 2]); panel(ax, "3. Phase differences |Δθ|/π (45 pairs)")
    table(ax, PHASE_HDR, PHASE_ROWS, [0.3, 0.14, 0.56], bbox=(0.03, 0.03, 0.94, 0.94), fs=7.9)
    # row 1: pyramid reading (5 panels)
    sub = gs[1, :].subgridspec(1, 6, wspace=0.06, width_ratios=[1.15, 0.85, 1.15, 0.75, 0.75, 0.75])
    ax = fig.add_subplot(sub[0, 0]); panel(ax, "4. Square pyramid — (1) vertices"); draw_pyramid(ax, "oblique", labels=True)
    ax.legend(handles=[Line2D([0], [0], color="k", lw=1.5, label="outer edges (8)"), Line2D([0], [0], color="k", lw=1.2, ls="--", alpha=0.6, label="base diagonals (2)")],
              loc="lower right", fontsize=6.8)
    ax = fig.add_subplot(sub[0, 1]); panel(ax, "(2) edge classes")
    L = [("outer edges (8)", "#333", 0.93), ("12  (A₊)", RED, 0.86), ("13  (A₊)", RED, 0.80), ("14  (A₋)", RED, 0.74), ("15  (A₋)", RED, 0.68),
         ("24  (B₊)", BLUE, 0.59), ("43  (B₋)", BLUE, 0.53), ("35  (B₊)", BLUE, 0.47), ("52  (B₋)", BLUE, 0.41),
         ("base diagonals (2)", "#333", 0.30), ("23  (A₋)", RED, 0.22), ("45  (A₊)", RED, 0.16)]
    for t, col, y in L:
        ax.text(0.06, y, t, transform=ax.transAxes, fontsize=8, color=col, va="center", weight="bold" if col == "#333" else "normal")
    ax.text(0.46, 0.77, "} apex 1 → base\n   (4 edges)", transform=ax.transAxes, fontsize=7.4, color="#333", va="center")
    ax.text(0.46, 0.50, "} base 4-cycle\n   (outer boundary)", transform=ax.transAxes, fontsize=7.4, color="#333", va="center")
    ax.text(0.46, 0.19, "} diagonals of the\n   base square", transform=ax.transAxes, fontsize=7.4, color="#333", va="center")
    ax = fig.add_subplot(sub[0, 2]); panel(ax, "(3) colour-coded pyramid"); draw_pyramid(ax, "oblique", labels=False)
    for i, (v, t) in enumerate([("top", "5. Projections — (1) top view"), ("front", "(2) front view"), ("oblique", "(3) oblique view")]):
        ax = fig.add_subplot(sub[0, 3 + i]); panel(ax, None); draw_pyramid(ax, v, labels=False, lw=1.5); ax.set_title(t, fontsize=7.8, color=NAVY, pad=2, loc="left")
    ax = fig.add_subplot(gs[2, 0]); panel(ax, "6. Summary of symmetries")
    textbox(ax, ["• looking only at the magnitudes, all 10 edges are equal",
                 "   → high symmetry of a regular 4-simplex (S₅)",
                 "• the phase structure reduces this symmetry spontaneously to the 4 classes 3+3+2+2",
                 "• in 3D this reads naturally as 'outer shell of a square pyramid (8) + base diagonals (2)'",
                 "• the base 4-cycle (2-4-3-5) and its diagonals (23, 45) are distinguished by phase"], y=0.92, fs=7.8, lsp=1.45)
    ax = fig.add_subplot(gs[2, 1]); panel(ax, "7. Complex distances in the z-plane"); zplane(ax, fs=7)
    ax = fig.add_subplot(gs[2, 2]); panel(ax, "8. Numerical verification (step 5000)")
    r2 = cls["r2"]
    textbox(ax, [f"• range of |zᵢⱼ|² : {r2.min():.12f} – {r2.max():.12f}",
                 f"   maximum deviation from 0.1 : ≈ ±{max(abs(r2.max()-0.1), abs(r2.min()-0.1)):.1e}",
                 "• Gram matrix rank : 4  (numerical error ~ 10⁻¹⁴)",
                 f"• representative phase difference δ/π : {DELTA:.6f}  (mode)"], y=0.9, fs=8.2, lsp=1.6)
    fig.savefig(os.path.join(OUT, "N5_square_pyramid_interpretation_en.png"), dpi=200, facecolor="white"); plt.close(fig)

if __name__ == "__main__":
    figure_A(); figure_B(); figure_C()
    print("done ->", OUT)
