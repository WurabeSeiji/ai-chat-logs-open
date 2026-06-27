#!/usr/bin/env python3
"""
Normalized Intensity Interference Pattern
Two coherent sources with n=5 (λ=1/5), placed at ±15° phase centers,
with 30° relative phase difference.
Complex amplitude representation: Z = cosθ + i sinθ with |Z|=1 canonical.
Intensity I = |Ψ|^2 normalized so max wave amplitude =1, thus I max=1.
Full x-axis ±360° to display all interference fringes (no clipping).
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters
n = 5                    # mode index, λ = 1/n
delta_deg = 15.0         # source placement offset from center (degrees)
alpha_deg = 30.0         # relative phase difference between the two sources (degrees)
phi_min, phi_max = -360.0, 360.0
num_points = 4000        # high resolution for dense fringes

# Observation phase variable (degrees)
phi_deg = np.linspace(phi_min, phi_max, num_points)
phi_rad = np.deg2rad(phi_deg)

# Convert offsets to radians
delta_rad = np.deg2rad(delta_deg)
alpha_rad = np.deg2rad(alpha_deg)

# Wave from source 1 (centered at -delta, i.e. offset +delta in argument)
# Using real part of Z: ψ = Re( exp(i * n * (phi - phi_c)) ) = cos( n * (phi - phi_c) )
psi1 = np.cos( n * (phi_rad + delta_rad) )

# Wave from source 2 (centered at +delta), with additional relative phase alpha
psi2 = np.cos( n * (phi_rad - delta_rad) + alpha_rad )

# Total coherent superposition (complex amplitudes add, here real representation)
psi_total = psi1 + psi2

# Canonical normalization: max wave height (amplitude) to 1.0
max_amp = np.max(np.abs(psi_total))
if max_amp > 0:
    psi_norm = psi_total / max_amp
else:
    psi_norm = psi_total

# Intensity I ∝ |ψ|^2 , normalized so max I = 1.0
I = psi_norm ** 2

# Create figure
fig, ax = plt.subplots(figsize=(14, 6), dpi=150)

# Plot the interference pattern
ax.plot(phi_deg, I, color='#1f77b4', linewidth=1.2, label='Normalized Intensity I(φ)')

# Highlight central ±90° region with light background (for reference, but pattern shown full range)
ax.axvspan(-90, 90, alpha=0.08, color='yellow', label='Central ±90° region (for reference)')

# Add vertical lines at ±90°
ax.axvline(x=-90, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
ax.axvline(x=90, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)

# Labels and title (ALL IN ENGLISH as requested)
ax.set_xlabel('Phase φ (degrees)', fontsize=12)
ax.set_ylabel('Normalized Intensity I(φ)  [max = 1.0]', fontsize=12)
ax.set_title(
    'Interference Pattern from Two Coherent Sources (Complex Amplitude Z = cosθ + i sinθ, |Z|=1)\n'
    'n=5 (λ=1/5), Sources placed at ±15° around center, Relative phase difference = 30°\n'
    'Full ±360° range — Interference fringes visible across entire domain (no clipping outside ±90°)',
    fontsize=11, pad=10
)

# Set axis limits
ax.set_xlim(phi_min, phi_max)
ax.set_ylim(0.0, 1.05)

# Grid and legend
ax.grid(True, which='major', linestyle=':', alpha=0.6)
ax.legend(loc='upper right', fontsize=9)

# Add annotation box explaining key points
textstr = '\n'.join([
    'Key settings (formal, no physical mapping):',
    '• Single mode n=5 only',
    '• Each source: ψ = cos(n × (φ ∓ 15°))',
    '• Source 2 has extra +30° phase shift',
    '• Ψ = ψ₁ + ψ₂ , normalized max|Ψ|=1',
    '• I(φ) = [Ψ_norm]²   (from |Z| canonical)',
    '• X-axis full ±360° to show all fringes',
    '• Central ±90° indicated by dashed lines (yellow tint)',
    '',
    'This demonstrates coherent self-interference',
    'of the two-source system in the abstract phase space.'
])
props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.85)
ax.text(0.015, 0.98, textstr, transform=ax.transAxes, fontsize=8,
        verticalalignment='top', bbox=props, family='monospace')

# Tight layout
plt.tight_layout()

# Save outputs
output_dir = '/home/workdir/artifacts/'
png_path = output_dir + 'double_source_interference_n5.png'
svg_path = output_dir + 'double_source_interference_n5.svg'

plt.savefig(png_path, dpi=150, bbox_inches='tight', facecolor='white')
plt.savefig(svg_path, bbox_inches='tight', facecolor='white')

print(f"PNG saved to: {png_path}")
print(f"SVG saved to: {svg_path}")
print("Plot generated successfully with full ±360° range and interference fringes.")
