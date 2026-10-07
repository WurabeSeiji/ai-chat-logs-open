#!/usr/bin/env python3
"""English version of the note figure (general audience, schematic, no data).

Left: the picture in space. Two bodies (proton and electron) at the centre, a shell at a great distance that
      absorbs all radiation. A sphere in the rest frame of the absorber; an ellipse r(θ) = t₀/(1 − β cos θ)
      when cut by the simultaneity surface of a frame in which the pair moves.
Right: the picture in spacetime. Each emission event on the worldline of the pair travels along the surface
      of its own future light cone (the null surface), lands on future infinity (the absorber) and becomes a
      record (time, direction, frequency, polarization). The background radiation at 2.7255 K comes in from
      past infinity.
Output: note_fig_lightcone_absorber_en.png
Same layout as make_note_fig_lightcone_ja.py (the Japanese version); the paper's figure 16 redrawn for a
general audience.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "figure.dpi": 150,
                     "savefig.dpi": 200, "axes.unicode_minus": False})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 6.4))

# ------------------------------------------------------------------ left: space
ax1.set_aspect("equal")
ax1.axis("off")
ax1.set_xlim(-1.5, 2.05)
ax1.set_ylim(-1.74, 1.42)
ax1.set_anchor("N")
ax1.set_title("In space:\nthe two bodies and the complete absorber at infinity", fontsize=12, pad=8)

ax1.add_patch(Circle((0, 0), 1.0, fill=False, color="C0", lw=3.0))
beta, t0 = 0.30, 1.0
a = t0 / (1 - beta**2)
b = a * np.sqrt(1 - beta**2)
ax1.add_patch(Ellipse((a * beta, 0), 2 * a, 2 * b, fill=False, color="C3", lw=2.6, ls="--"))

ax1.add_patch(Circle((0, 0), 0.14, fill=False, color="0.4", lw=0.9, ls=":"))
ax1.scatter([0], [0], s=160, color="0.15", zorder=6)
ax1.scatter([0.14 * np.cos(0.9)], [0.14 * np.sin(0.9)], s=40, color="C2", zorder=6)
ax1.text(0.02, -0.21, "proton", fontsize=10, ha="center", va="top")
ax1.text(0.17, 0.17, "electron", fontsize=10, ha="left", va="bottom", color="C2")

for ang in (0.35, 1.25, 2.3, 3.4, 4.3, 5.3):
    ax1.add_patch(FancyArrowPatch((0.25 * np.cos(ang), 0.25 * np.sin(ang)),
                                  (0.90 * np.cos(ang), 0.90 * np.sin(ang)),
                                  arrowstyle="-|>", mutation_scale=14, lw=1.3, color="0.45"))
ax1.text(-0.97, 0.62, "radiation (light and\ngravitational waves)\ngoes out and is absorbed\nby the absorber",
         fontsize=9.5, color="0.3", ha="right", va="center")

for ang in (0.0, 1.9, 4.6):
    s = np.linspace(0.98, 0.42, 60)
    wig = 0.025 * np.sin(2 * np.pi * (s - 0.42) / 0.14)
    x = s * np.cos(ang) - wig * np.sin(ang)
    y = s * np.sin(ang) + wig * np.cos(ang)
    ax1.plot(x, y, color="C1", lw=1.3)
    ax1.add_patch(FancyArrowPatch((x[-2], y[-2]), (x[-1], y[-1]), arrowstyle="-|>",
                                  mutation_scale=12, lw=1.0, color="C1"))
ax1.text(1.27, -0.62, "background radiation\n(2.7255 K) comes in", fontsize=9.5, color="C1", ha="left", va="center")

ax1.text(0.0, 1.07, "complete absorber at infinity\n(a sphere in the absorber's rest frame)", fontsize=10,
         color="C0", ha="center", va="bottom")
ax1.text(a * beta + 0.0, -b - 0.06,
         "an ellipse in a frame where the pair moves at speed v\n(the geometry of the Doppler shift and aberration)",
         fontsize=10, color="C3", ha="center", va="top")
ax1.text(0.2, -1.71,
         "The motion of the pair is read out only from the record of emission and absorption.\n"
         "No coordinate grid and no metric are held in the state.",
         fontsize=9.2, color="0.3", ha="center", va="bottom")

# ------------------------------------------------------------------ right: spacetime
ax2.set_aspect("equal")
ax2.axis("off")
ax2.set_xlim(-0.80, 1.70)
ax2.set_ylim(-1.46, 1.30)
ax2.set_anchor("N")
ax2.set_title("In spacetime:\neach emission event lands on the absorber on its own light cone", fontsize=12, pad=8)

ax2.plot([0, 0], [-1.0, 1.0], color="k", lw=2.6)
ax2.text(-0.03, 1.03, "worldline of the pair (e.g. a hydrogen atom)", fontsize=9.5, ha="right", va="bottom")

ax2.plot([1.0, 0.0], [0.0, 1.0], color="C0", lw=3.2)
ax2.annotate("future infinity = the absorber\n(surface of the future light cone, the null surface)",
             xy=(0.24, 0.76), xytext=(0.50, 0.98), fontsize=9.8, color="C0", ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color="C0", lw=0.8))
ax2.plot([1.0, 0.0], [0.0, -1.0], color="C1", lw=3.2)
ax2.annotate("past infinity = incoming background radiation\n(last scattering, 2.7255 K)",
             xy=(0.24, -0.76), xytext=(0.50, -0.98), fontsize=9.8, color="C1", ha="left", va="center",
             arrowprops=dict(arrowstyle="-", color="C1", lw=0.8))

events = ((-0.78, "emission 1", "u₁"), (-0.36, "emission 2", "u₂"), (0.08, "emission 3", "u₃"))
for te, lab, u in events:
    s = (1 - te) / 2
    ax2.plot([0, s], [te, te + s], color="C0", lw=1.1)
    ax2.scatter([0], [te], color="k", s=22, zorder=5)
    ax2.scatter([s], [te + s], color="C0", s=30, zorder=5)
    ax2.text(-0.04, te, lab, fontsize=9.5, ha="right", va="center")
    ax2.text(s + 0.04, te + s - 0.01, "record (time %s, direction,\nfrequency, polarization)" % u,
             fontsize=9.0, va="center", ha="left")

for ta in (-0.9, -0.56, -0.14):
    s = (1 + ta) / 2
    ax2.plot([s, 0], [ta - s, ta], color="C1", lw=0.9, ls="--")

ax2.annotate("", xy=(-0.12, 0.64), xytext=(-0.12, 0.26), arrowprops=dict(arrowstyle="-|>", lw=1.2, color="0.3"))
ax2.text(-0.17, 0.45, "between emissions:\nno readable change", fontsize=9.5, color="0.3", ha="right", va="center")

ax2.annotate("", xy=(-0.62, -0.78), xytext=(-0.62, -1.08), arrowprops=dict(arrowstyle="-|>", lw=1.0, color="0.4"))
ax2.annotate("", xy=(-0.32, -1.08), xytext=(-0.62, -1.08), arrowprops=dict(arrowstyle="-|>", lw=1.0, color="0.4"))
ax2.text(-0.60, -0.74, "time", fontsize=9.5, color="0.4", ha="center", va="bottom")
ax2.text(-0.28, -1.08, "space", fontsize=9.5, color="0.4", ha="left", va="center")

ax2.text(-0.78, -1.44,
         "Light travels along the 45-degree lines (the surface of the light cone).\n"
         "A record is carried unchanged along that line; the Doppler shift appears\n"
         "as the mismatch of the clocks at its two ends.",
         fontsize=9.2, color="0.3", ha="left", va="bottom")

fig.tight_layout()
out = HERE / "note_fig_lightcone_absorber_en.png"
fig.savefig(out)
print("saved", out)
