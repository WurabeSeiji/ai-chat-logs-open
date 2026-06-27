#!/usr/bin/env python3
"""
Normalized Intensity Plot for Odd Modes (n=1,3,5,7,9)
- Complex amplitude Z_n normalized |Z_n|=1
- Central phase = 0, ± considered but I = cos²(n φ)
- Figure x-axis spans ±360° (for future multi-path convenience)
- Curves DRAWN ONLY in central ±90° range
- Outside ±90°: deliberately left BLANK (no curves)
- All text/captions in English
- Saves PNG (300dpi), SVG (vector), and this script itself
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Use a clean sans-serif font
rcParams['font.family'] = 'DejaVu Sans'
rcParams['mathtext.fontset'] = 'cm'  # for nice math if needed

# Parameters
phi = np.linspace(-360.0, 360.0, 14401)  # 0.05° resolution for smooth high-n curves
odd_ns = [1, 3, 5, 7, 9]
mode_labels = [f"n={n}  (λ = 1/{n})" for n in odd_ns]

# Nice distinguishable colors (tab10 palette, slightly desaturated for clarity)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # blue, orange, green, red, purple

fig, ax = plt.subplots(figsize=(18, 9), facecolor='white')  # wide figure for ±360° span

for idx, n in enumerate(odd_ns):
    # Compute intensity I(φ) = cos²(n φ)   [φ in degrees]
    phi_rad = np.deg2rad(n * phi)
    I = np.cos(phi_rad) ** 2
    
    # CRITICAL: Only keep values inside central ±90° ; set NaN outside → no curve drawn
    I_plot = np.where(np.abs(phi) <= 90.0, I, np.nan)
    
    # Plot the line only where defined (central region)
    # Pure line plot for maximum clarity of oscillations (especially higher n)
    ax.plot(phi, I_plot, 
            label=mode_labels[idx], 
            color=colors[idx], 
            linewidth=2.4, 
            solid_capstyle='round',
            zorder=3)

# Mark the closed-system boundary (±90°)
ax.axvline(x=-90, color='#555555', linestyle='--', linewidth=1.8, alpha=0.85, zorder=4)
ax.axvline(x=90, color='#555555', linestyle='--', linewidth=1.8, alpha=0.85, zorder=4)

# Very light background highlight ONLY for the plotted central region
ax.axvspan(-90, 90, alpha=0.06, color='#4a90d9', zorder=0)

# Axis limits and labels
ax.set_xlim(-360, 360)
ax.set_ylim(-0.02, 1.08)
ax.set_xlabel("Phase φ (degrees)   [Reference scale: λ=1  →  central span of ±90°  (180° closed system)]", 
              fontsize=13, labelpad=10)
ax.set_ylabel("Normalized Intensity   I(φ) = |ψ(φ)|²   [maximum height canonically set to 1.0]", 
              fontsize=13, labelpad=10)

# Title (English)
ax.set_title(
    "Normalized Intensity Distributions for Odd Harmonic Modes\n"
    "in Closed Boundaryless System (Formal Mathematical Model)\n"
    "Zₙ = cos(nπ) + i sin(nπ) = −1 (odd n)   |   Iₙ(φ) = cos²(n φ)   |   Central phase = 0°",
    fontsize=15, pad=18, linespacing=1.3
)

# Legend
legend = ax.legend(loc='upper right', fontsize=11, framealpha=0.95, 
                   edgecolor='gray', fancybox=True, shadow=False)
legend.get_frame().set_linewidth(0.8)

# Grid
ax.grid(True, which='major', alpha=0.25, linestyle='-', linewidth=0.7)
ax.grid(True, which='minor', alpha=0.12, linestyle=':', linewidth=0.5)
ax.minorticks_on()

# Informative text box (English) explaining the plotting rule
textstr = (
    "PLOTTING RULE (as specified):\n"
    "• Figure x-axis width: full ±360° (prepared for future double-slit / multi-path phase shifts)\n"
    "• Curves drawn EXCLUSIVELY inside central ±90° range\n"
    "• Outside ±90° → deliberately BLANK (no curves, no repetition)\n"
    "• 180° closed system assumed; ±90° marks the boundary\n"
    "• All modes share the same reference phase scale (λ=1 defines ±90°)"
)
props = dict(boxstyle='round,pad=0.6', facecolor='#fff8dc', edgecolor='#8b7355', alpha=0.92)
ax.text(0.015, 0.985, textstr, transform=ax.transAxes, fontsize=9.5,
        verticalalignment='top', bbox=props, family='sans-serif',
        linespacing=1.35)

# Footer note
fig.text(0.5, 0.01, 
         "Formal model only — physical interpretations ignored.  |  "
         "Internal data computed to ±360° but display strictly clipped to central ±90° as requested.",
         ha='center', fontsize=9, style='italic', color='#444444')

plt.tight_layout(rect=[0, 0.03, 1, 1])  # leave space for footer

# === Save outputs ===
output_dir = "/home/workdir/artifacts/"
base_name = "odd_modes_intensity_central90_only"

png_path = f"{output_dir}{base_name}.png"
svg_path = f"{output_dir}{base_name}.svg"

plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(svg_path, bbox_inches='tight', facecolor='white', edgecolor='none')

print(f"✓ PNG saved:  {png_path}")
print(f"✓ SVG saved:  {svg_path}")
print("✓ Plot generated with curves ONLY in central ±90°; outer regions left blank.")
print("Ready for download.")