import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Coherent Self-Interference Plot for Odd Modes
# Z_n = cos(n φ) + i sin(n φ), |Z_n| = 1 canonical normalization
# Closed boundaryless system, central phase = 0
# ============================================================

# High-resolution phase data: internal up to ±360°
x_deg = np.linspace(-360.0, 360.0, 8000)
x_rad = np.deg2rad(x_deg)

# Odd modes (λ = 1, 1/3, 1/5, 1/7, 1/9)
odd_ns = [1, 3, 5, 7, 9]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # distinct colors

# Individual mode intensities: I_n(φ) = cos²(n φ)
Is_masked = {}
for n in odd_ns:
    theta = n * x_rad
    I = np.cos(theta) ** 2
    # Draw ONLY in central ±90°; outer regions blank (NaN)
    I_masked = np.where(np.abs(x_deg) <= 90.0, I, np.nan)
    Is_masked[n] = I_masked

# Total coherent self-interference
# Ψ(φ) = Σ Z_n(φ) = Σ exp(i n φ)   [coherent sum, phase-aligned at center]
total_psi = np.zeros_like(x_rad, dtype=complex)
for n in odd_ns:
    total_psi += np.exp(1j * n * x_rad)

I_total = np.abs(total_psi) ** 2

# Normalize using max in the physical central region only
central_mask = np.abs(x_deg) <= 90.0
max_I_total = np.nanmax(I_total[central_mask])
I_total_norm = I_total / max_I_total

# Mask total to central ±90° only
I_total_masked = np.where(central_mask, I_total_norm, np.nan)

# ============================================================
# Plotting
# ============================================================
fig, ax = plt.subplots(figsize=(16, 9), facecolor='white')

# Light shade for the physical observation region (central ±90°)
ax.axvspan(-90, 90, alpha=0.06, color='steelblue', zorder=0)

# Dashed vertical lines marking closed system boundary
ax.axvline(x=-90, color='#555555', linestyle='--', linewidth=1.5, alpha=0.85, zorder=1)
ax.axvline(x=90, color='#555555', linestyle='--', linewidth=1.5, alpha=0.85, zorder=1)

# Plot individual mode curves (thinner)
for idx, n in enumerate(odd_ns):
    ax.plot(x_deg, Is_masked[n], color=colors[idx], linewidth=1.6, 
            label=f'n={n}  (λ=1/{n})', zorder=2)

# Plot TOTAL coherent self-interference (prominent, thick)
ax.plot(x_deg, I_total_masked, color='#000000', linewidth=3.2, 
        linestyle='-', label='TOTAL Coherent Self-Interference  |Σ Z_n|² (normalized)', zorder=3)

# Styling
ax.set_xlim(-360, 360)
ax.set_ylim(-0.02, 1.08)
ax.set_xlabel('Phase φ (degrees)   [Reference scale based on λ=1]\n'
              'Curves drawn ONLY in central ±90° region  •  Outer regions intentionally left blank (closed 180° system rule)',
              fontsize=10.5)
ax.set_ylabel('Normalized Intensity   I(φ)   (maximum = 1.0)', fontsize=11)

ax.set_title('Coherent Self-Interference Wave Expression\n'
             'Odd Harmonic Modes with Canonical Complex Amplitude  Z = cosθ + i sinθ   (|Z| = 1)\n'
             'Closed Boundaryless System  •  Central Phase = 0  •  ± Signs Considered',
             fontsize=13, pad=12, fontweight='medium')

ax.legend(loc='upper right', fontsize=9, framealpha=0.92, edgecolor='gray')
ax.grid(True, which='major', linestyle='-', alpha=0.25, zorder=0)
ax.grid(True, which='minor', linestyle=':', alpha=0.15, zorder=0)
ax.minorticks_on()

# Annotations
ax.text(0.0, 0.97, 'Physical observation range (closed system):  ±90°', 
        ha='center', va='top', fontsize=9, color='#333333',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white', alpha=0.85, edgecolor='#888888', linewidth=0.8))

ax.text(0.0, 0.90, 'Outer regions (±90° to ±360°): blank — no curve repetition shown', 
        ha='center', va='top', fontsize=8, color='#555555', style='italic')

ax.text(270, 0.55, 'Internal phase data\nprepared up to ±360°\n(for future multi-turn\nor extended analysis)', 
        ha='center', va='center', fontsize=8, color='#666666', alpha=0.75,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#fffde7', alpha=0.6, edgecolor='none'))

# Bottom note
fig.text(0.5, 0.015, 
         'All captions in English.  Formal treatment only (no physical mapping).  '
         'Coherent superposition of odd-n modes with phase alignment at center.  '
         'I_n(φ) = cos²(nφ)   |   Total I(φ) = |Σ exp(i n φ)|² / max',
         ha='center', va='bottom', fontsize=7.5, style='italic', color='#444444')

plt.tight_layout(rect=[0, 0.04, 1, 0.98])

# Save outputs
png_path = '/home/workdir/artifacts/coherent_self_interference_odd_modes.png'
svg_path = '/home/workdir/artifacts/coherent_self_interference_odd_modes.svg'
plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(svg_path, bbox_inches='tight', facecolor='white', edgecolor='none')
print(f"PNG saved: {png_path}")
print(f"SVG saved: {svg_path}")
print("Plot generation complete.")
