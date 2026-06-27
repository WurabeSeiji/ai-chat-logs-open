#!/usr/bin/env python3
"""
Correct Coherent Interference Plot for Two Sources at ±15° Phase Offset
Odd discrete modes n=1,3,5,7,9 (λ=1, 1/3, 1/5, 1/7, 1/9)
Full horizontal axis ±360° with central closed-system region |φ|≤90° highlighted.
All captions in English as requested.
Uses correct formulation: ±15° phase offset is OUTSIDE n multiplier.
ψ(φ) = Σ_n [cos(nφ - α) + cos(nφ + α)] with α=15°
Then I(φ) = [ψ(φ)]² normalized to max=1.0
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os

# Parameters
modes = [1, 3, 5, 7, 9]
n_points = 8001  # High resolution for fine structure from high n
phi_deg = np.linspace(-360.0, 360.0, n_points)
phi = np.deg2rad(phi_deg)
alpha_deg = 15.0
alpha = np.deg2rad(alpha_deg)

# Correct coherent amplitude sum (amplitudes first, then intensity)
psi = np.zeros_like(phi)
for n in modes:
    # ±15° is common source position phase offset, NOT multiplied by n
    psi += np.cos(n * phi - alpha) + np.cos(n * phi + alpha)

I = psi ** 2
I_norm = I / np.max(I)

# Plot setup
fig, ax = plt.subplots(figsize=(16, 9), dpi=150)

# Main curve - full range ±360°
ax.plot(phi_deg, I_norm, color='#1f77b4', linewidth=1.2, label='Normalized Intensity I(φ)')

# Highlight central closed-system region |φ| ≤ 90° (the "180° closed system" core)
ax.axvspan(-90, 90, alpha=0.15, color='green', label='Central closed-system region (|φ| ≤ 90°)')
ax.axvline(x=-90, color='red', linestyle='--', linewidth=1.5, alpha=0.8)
ax.axvline(x=90, color='red', linestyle='--', linewidth=1.5, alpha=0.8)

# Annotations for central region
ax.annotate('', xy=(90, 0.95), xytext=(-90, 0.95),
            arrowprops=dict(arrowstyle='<->', color='darkgreen', lw=2))
ax.text(0, 0.97, 'Central Core Region\n(±90° reference for closed system)', 
        ha='center', va='bottom', fontsize=10, color='darkgreen', fontweight='bold')

# Title and labels (ALL IN ENGLISH)
ax.set_title(
    'Coherent Two-Source Interference\n'
    'Sources at ±15° Phase Offset around Center (Phase Difference = 30°)\n'
    'Odd Modes n=1,3,5,7,9  (λ = 1, 1/3, 1/5, 1/7, 1/9)  |  Full Range ±360°',
    fontsize=14, fontweight='bold', pad=20
)
ax.set_xlabel('Phase φ (degrees)', fontsize=12)
ax.set_ylabel('Normalized Intensity  I(φ) / max[I]', fontsize=12)
ax.set_xlim(-360, 360)
ax.set_ylim(0, 1.05)
ax.grid(True, alpha=0.3, linestyle=':')

# Setup description text box (English)
textstr = (
    'Setup (Formal, Phase-Only Model):\n'
    '• 5 odd modes simultaneously from each source (coherent)\n'
    '• Two sources placed symmetrically at α = ±15° around φ=0\n'
    '• Correct formula:  ψ(φ) = Σₙ [cos(nφ − α) + cos(nφ + α)]   (α=15°)\n'
    '• Equivalent to:     ψ(φ) = 2·cos(15°) · Σₙ cos(nφ)\n'
    '• Intensity: I(φ) = [ψ(φ)]²   → normalized max=1.0\n'
    '• ±15° is COMMON position offset for ALL modes (n-independent)\n'
    '• Full ±360° shown; central ±90° is the reference closed-system core\n'
    '• Coherent superposition of 10 waves (5 modes × 2 sources) at amplitude level'
)
props = dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.9, edgecolor='gray')
ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
        verticalalignment='top', bbox=props, family='monospace')

# Legend
ax.legend(loc='upper right', fontsize=9, framealpha=0.95)

# Add note about symmetry
ax.text(0.98, 0.02, 'Strict left-right symmetry about φ=0\n(Central peak is single maximum)',
        transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
        bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.8))

plt.tight_layout()

# Save outputs
output_dir = '/home/workdir/artifacts'
png_path = os.path.join(output_dir, 'two_source_coherent_interference_corrected_v4.png')
svg_path = os.path.join(output_dir, 'two_source_coherent_interference_corrected_v4.svg')
hires_png_path = os.path.join(output_dir, 'two_source_coherent_interference_corrected_v4_hires.png')

plt.savefig(png_path, dpi=150, bbox_inches='tight')
plt.savefig(svg_path, format='svg', bbox_inches='tight')
plt.savefig(hires_png_path, dpi=400, bbox_inches='tight')

print(f"Files saved successfully:")
print(f"  PNG (standard): {png_path}")
print(f"  SVG (vector):   {svg_path}")
print(f"  PNG (hires 400dpi): {hires_png_path}")
print(f"  Python script:  /home/workdir/artifacts/two_source_coherent_interference_corrected_v4.py")

# Also print some verification values
print("\nVerification at key points (after norm):")
for p in [0, 15, 30, 45, 90, 180, 270, 360]:
    idx = np.argmin(np.abs(phi_deg - p))
    print(f"  φ = {p:6.1f}° : I_norm = {I_norm[idx]:.6f}")

plt.close()