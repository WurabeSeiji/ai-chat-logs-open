import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 日本語フォント設定
jp_fonts = [f.name for f in fm.fontManager.ttflist if any(k in f.name for k in ['Noto Sans CJK', 'IPAexGothic', 'IPAGothic', 'Takao', 'VL Gothic', 'Meiryo', 'Hiragino', 'Yu Gothic'])]
if jp_fonts:
    plt.rcParams['font.family'] = jp_fonts[0]
    use_jp = True
else:
    use_jp = False

phi_deg = np.linspace(-180, 180, 5000)
phi = np.deg2rad(phi_deg)

def S10(x):
    return np.sum([np.cos(n * x) for n in range(1, 11)], axis=0)

f_fast = S10(10 * phi)
shifts = np.arange(-180, 181, 36)
fig, axes = plt.subplots(len(shifts), 1, figsize=(11, 15), sharex=True, sharey=True)
amp_factor = np.sqrt(0.5)

for idx, w in enumerate(shifts):
    ax = axes[idx]
    w_rad = np.deg2rad(w)
    
    window_sum = S10(phi - w_rad) + amp_factor * S10(phi + w_rad)
    Y_squared = (f_fast * window_sum) ** 2
    
    ax.plot(phi_deg, Y_squared, color='#d62728', linewidth=1.5)
    ax.plot(phi_deg, f_fast ** 2, color='gray', linestyle=':', alpha=0.3, linewidth=1.0)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_xlim(-180, 180)
    ax.set_ylim(0, 44000)
    ax.set_yticks([0, 10000, 40000])
    ax.set_yticklabels(['0', '10k', '40k'], fontsize=8)
    
    if w == 0:
        status = 'ω = 0° (Central Crossing: ~29.1k)'
    elif abs(w) == 180:
        status = f'ω = {w}° (Boundary: ~29.1k)'
    elif w < 0:
        status = f'ω = {w}° (Before: Left 10k, Right 5k)'
    else:
        status = f'ω = +{w}° (After: Left 5k, Right 10k)'
        
    ax.text(-176, 32000, status, fontsize=9.5, fontweight='bold', bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=1))
    
    if idx == len(shifts) - 1:
        ax.set_xticks(np.arange(-180, 181, 36))
        ax.set_xticklabels([f'{x}°' for x in np.arange(-180, 181, 36)], fontsize=9)
        ax.set_xlabel('Phase [Degrees] (φ)', fontsize=11)

fig.suptitle('Asymmetric Transmission Simulation (Right Wave Amp 0.707)', fontsize=12.5, fontweight='bold', y=0.985)
fig.text(0.015, 0.5, 'Y²', va='center', rotation='vertical', fontsize=11)
plt.subplots_adjust(hspace=0.25, top=0.94, bottom=0.045, left=0.08, right=0.96)
plt.savefig('harmonics_asymmetric_crossing_sim.png', dpi=150)
plt.savefig('harmonics_asymmetric_crossing_sim.svg')
plt.close()
