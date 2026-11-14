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

phi_deg = np.linspace(-180, 180, 6000)
phi = np.deg2rad(phi_deg)

def S10(x):
    return np.sum([np.cos(n * x) for n in range(1, 11)], axis=0)

f_fast = S10(10 * phi)
shifts = np.arange(-180, 181, 36)
fig, axes = plt.subplots(len(shifts), 1, figsize=(11, 14.5), sharex=True, sharey=True)

for idx, w in enumerate(shifts):
    ax = axes[idx]
    w_rad = np.deg2rad(w)
    window_sum = S10(phi - w_rad) + S10(phi + w_rad)
    Y_squared = (f_fast * window_sum) ** 2
    
    ax.plot(phi_deg, Y_squared, color='#d62728', linewidth=1.6)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_xlim(-180, 180)
    ax.set_ylim(0, 44000)
    ax.set_yticks([0, 10000, 40000])
    ax.set_yticklabels(['0', '10k', '40k'], fontsize=8)
    
    if w == 0:
        status = 'ω = 0° (Central Crossing: 40k)'
    elif abs(w) == 180:
        status = f'ω = {w}° (Boundary Crossing: 40k)'
    elif w < 0:
        status = f'ω = {w}° (Before: 10k x 2)'
    else:
        status = f'ω = +{w}° (After: 10k x 2)'
        
    ax.text(-176, 31000, status, fontsize=9.5, fontweight='bold', bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=1.5))
    
    if idx == len(shifts) - 1:
        ax.set_xticks(np.arange(-180, 181, 36))
        ax.set_xticklabels([f'{x}°' for x in np.arange(-180, 181, 36)], fontsize=9)
        ax.set_xlabel('Phase [Degrees] (φ)', fontsize=11)

fig.suptitle('Crossing Simulation of Counter-Propagating Peaks (Squared Amplitude)', fontsize=12.5, fontweight='bold', y=0.98)
fig.text(0.02, 0.5, 'Y²', va='center', rotation='vertical', fontsize=11)
plt.subplots_adjust(hspace=0.25, top=0.93, bottom=0.05, left=0.08, right=0.96)
plt.savefig('harmonics_crossing_11panels.png', dpi=150)
plt.savefig('harmonics_crossing_11panels.svg')
plt.close()
