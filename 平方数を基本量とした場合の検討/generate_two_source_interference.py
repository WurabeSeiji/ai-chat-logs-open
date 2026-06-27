#!/usr/bin/env python3
"""
Generate 5 separate plots for coherent two-source interference patterns.
Each plot corresponds to one odd mode n=1,3,5,7,9 (λ=1,1/3,1/5,1/7,1/9).
Two coherent sources (same λ) placed symmetrically at ±15° phase positions around center 0°.
Phase difference between sources: 30°.
Interference calculated strictly from amplitude sum, then intensity normalized max=1.
X-axis: full ±360° with interference fringes displayed across entire range.
All captions and labels in English.
Saves PNG (high dpi) and SVG (vector) for each n, plus this script.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import os

output_dir = "/home/workdir/artifacts"
os.makedirs(output_dir, exist_ok=True)

# Odd modes n and corresponding wavelengths
ns = [1, 3, 5, 7, 9]
wavelengths = {1: "1", 3: "1/3", 5: "1/5", 7: "1/7", 9: "1/9"}

for n in ns:
    # High resolution for smooth high-n oscillations
    phi = np.linspace(-360.0, 360.0, 12000)
    
    alpha1 = -15.0  # degrees
    alpha2 = 15.0   # degrees
    
    # Phase arguments (radians)
    arg1 = n * (phi - alpha1) * np.pi / 180.0
    arg2 = n * (phi - alpha2) * np.pi / 180.0
    
    # Real-valued wave amplitudes (consistent with canonical |Z|=1 -> cos representation)
    psi1 = np.cos(arg1)
    psi2 = np.cos(arg2)
    
    # Total wave from coherent superposition
    psi_total = psi1 + psi2
    
    # Intensity
    I = psi_total ** 2
    
    # Normalize to maximum wave height = 1 (canonical normalization)
    max_I = np.max(I)
    I_norm = I / max_I if max_I > 0 else I
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 5.5))
    
    # Main interference curve (full ±360°)
    ax.plot(phi, I_norm, color='#1f77b4', linewidth=0.9, 
            label=f'Interference intensity I(φ)  [n={n}]')
    
    # Highlight central ±90° region (reference to closed 180° system)
    ax.axvspan(-90, 90, alpha=0.08, color='green', zorder=0)
    ax.axvline(x=90, color='red', linestyle='--', linewidth=1.2, alpha=0.75, 
               label='Closed-system reference boundary (±90°)')
    ax.axvline(x=-90, color='red', linestyle='--', linewidth=1.2, alpha=0.75)
    
    # Axes and limits
    ax.set_xlim(-360, 360)
    ax.set_ylim(0, 1.08)
    
    # Ticks
    ax.set_xticks(np.arange(-360, 361, 90))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(30))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    
    # Labels and title (all English)
    ax.set_xlabel('Phase φ (degrees)', fontsize=12)
    ax.set_ylabel('Normalized Intensity  I(φ)   [maximum = 1.0]', fontsize=12)
    
    ax.set_title(
        f'Coherent Two-Source Interference Fringes (Self-Interfered Canonical Wave)\n'
        f'Odd Mode n={n}   |   λ = {wavelengths[n]}   |   Sources placed at ±15° (phase diff = 30°)\n'
        f'Full phase range ±360° shown   |   Symmetric around center 0°   |   Max wave height canonically normalized to 1',
        fontsize=11, pad=12, linespacing=1.4
    )
    
    # Grid
    ax.grid(True, which='major', alpha=0.35, linestyle='-')
    ax.grid(True, which='minor', alpha=0.15, linestyle=':')
    
    # Explanatory text box (English)
    textstr = (
        'Setup (formal, no physical mapping):\n'
        '• Two coherent sources, identical λ (same n), placed symmetrically at phase positions α = ±15° around 0°.\n'
        '• Each source contributes canonical wave: ψₖ(φ) = cos[ n · (φ − αₖ) ]   (|Z| = 1 normalized)\n'
        '• Total amplitude: ψ(φ) = ψ₁(φ) + ψ₂(φ)\n'
        '• Observed intensity: I(φ) = [ψ(φ)]²   →   normalized so peak = 1.0\n'
        '• Due to symmetric placement, pattern remains centered; fringes visible across entire ±360° range.'
    )
    props = dict(boxstyle='round,pad=0.4', facecolor='lightyellow', alpha=0.92, edgecolor='gray')
    ax.text(0.015, 0.97, textstr, transform=ax.transAxes, fontsize=7.5,
            verticalalignment='top', horizontalalignment='left', bbox=props,
            family='monospace', linespacing=1.25)
    
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    
    # Tight layout
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Save PNG and SVG
    base_name = f'two_source_interference_n{n}'
    png_path = os.path.join(output_dir, f'{base_name}.png')
    svg_path = os.path.join(output_dir, f'{base_name}.svg')
    
    fig.savefig(png_path, dpi=220, bbox_inches='tight', facecolor='white', edgecolor='none')
    fig.savefig(svg_path, bbox_inches='tight', facecolor='white', edgecolor='none')
    
    plt.close(fig)
    print(f"Generated: {png_path}")
    print(f"Generated: {svg_path}")

print("\n=== All 5 interference plots (PNG + SVG) generated successfully ===")
print(f"Output directory: {output_dir}")
print("Files ready for download.")