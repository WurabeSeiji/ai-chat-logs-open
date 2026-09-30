from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch

OUT = Path(__file__).resolve().parent

# Japanese font if available
candidates = [
    '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',
    '/usr/share/fonts/opentype/noto/NotoSansCJKJP-Regular.otf',
]
for p in candidates:
    if Path(p).exists():
        font_manager.fontManager.addfont(p)
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=p).get_name()
        break
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['svg.fonttype'] = 'none'

DPI = 220

def save(fig, stem):
    fig.savefig(OUT / f'{stem}.png', dpi=DPI, bbox_inches='tight')
    fig.savefig(OUT / f'{stem}.svg', bbox_inches='tight')
    plt.close(fig)

# ----------------------------------------------------------------------
# Figure 1: monochromatic vs harmonics / intra-period resolution
# ----------------------------------------------------------------------
phi = np.linspace(0, 2*np.pi, 2000)
mono = np.cos(phi)
M = 9
# Equal-amplitude cosine harmonics normalized by M. This is illustrative, not an optimized pulse.
comb = np.sum([np.cos(m*phi) for m in range(1, M+1)], axis=0) / M

fig = plt.figure(figsize=(11.5, 6.8))
gs = fig.add_gridspec(2, 1, hspace=0.42)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0])

ax1.plot(phi/(2*np.pi), mono, lw=2)
ax1.set_title('単色波：1周期内部の構造は位相1本で決まる', fontsize=15)
ax1.set_xlim(0,1); ax1.set_ylim(-1.25,1.25)
ax1.set_ylabel('振幅')
ax1.set_xticks([0,0.25,0.5,0.75,1.0])
ax1.set_xticklabels(['0','T/4','T/2','3T/4','T'])
ax1.grid(alpha=0.25)
ax1.annotate('ベース周期 T', xy=(0.5,-1.1), ha='center', fontsize=12)

ax2.plot(phi/(2*np.pi), comb, lw=2)
ax2.set_title(f'倍音を含む波：m=1…{M} の干渉で1周期内部に細かな読出し構造が生じる', fontsize=15)
ax2.set_xlim(0,1)
ax2.set_ylabel('正規化振幅')
ax2.set_xlabel('ベース周期内部の位置')
ax2.set_xticks(np.linspace(0,1,10))
ax2.set_xticklabels([f'{i}/9 T' if i not in (0,9) else ('0' if i==0 else 'T') for i in range(10)], rotation=0)
ax2.grid(alpha=0.25)
for x in np.linspace(0,1,10):
    ax2.axvline(x, lw=0.6, alpha=0.18)
ax2.text(0.5, ax2.get_ylim()[0]*0.78,
         '倍音帯域 ↑  →  一周期内部の区別可能な位相位置 ↑',
         ha='center', va='center', fontsize=13)
fig.suptitle('図1　単色波と倍音波：一周期内部の時間分解能', fontsize=18, y=0.99)
save(fig, 'figure01_monochromatic_vs_harmonic_time_resolution')

# ----------------------------------------------------------------------
# Figure 2: cyclic FIFO / phase delay memory
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10.5, 8.0))
ax.set_aspect('equal'); ax.axis('off')
R = 2.8
circle = Circle((0,0), R, fill=False, lw=2)
ax.add_patch(circle)
N = 8
angles = np.linspace(np.pi/2, np.pi/2-2*np.pi, N, endpoint=False)
labels = [r'$X(\tau)$'] + [rf'$X(\tau-{j}\delta\tau)$' for j in range(1,N)]
for i,(ang,lab) in enumerate(zip(angles,labels)):
    x,y = R*np.cos(ang), R*np.sin(ang)
    ax.add_patch(Circle((x,y), 0.20, fill=True, alpha=0.15))
    ax.text(x*1.18,y*1.18,lab,ha='center',va='center',fontsize=12)
    ang2 = angles[(i+1)%N]
    p1 = ((R-0.18)*np.cos(ang-0.08),(R-0.18)*np.sin(ang-0.08))
    p2 = ((R-0.18)*np.cos(ang2+0.08),(R-0.18)*np.sin(ang2+0.08))
    ax.add_patch(FancyArrowPatch(p1,p2,arrowstyle='-|>',mutation_scale=13,lw=1.2,connectionstyle='arc3,rad=-0.17'))
ax.text(0,0.45,'一周期 T を',ha='center',va='center',fontsize=16)
ax.text(0,0.0,'循環履歴レジスタとして読む',ha='center',va='center',fontsize=17,weight='bold')
ax.text(0,-0.5,r'$T=N\,\delta\tau$',ha='center',va='center',fontsize=15)

# conventional FIFO below
start_y=-4.25
ax.text(-4.35,start_y+0.95,'通常のFIFO',fontsize=14,weight='bold',ha='left')
for j in range(5):
    x=-3.1+j*1.55
    box=FancyBboxPatch((x,start_y),1.15,0.72,boxstyle='round,pad=0.05',fill=False,lw=1.5)
    ax.add_patch(box)
    ax.text(x+0.575,start_y+0.36,f'cell {j}',ha='center',va='center',fontsize=10)
    if j<4:
        ax.add_patch(FancyArrowPatch((x+1.15,start_y+0.36),(x+1.5,start_y+0.36),arrowstyle='-|>',mutation_scale=12,lw=1.2))
ax.text(-4.35,start_y-0.55,'本モデル：固定セルではなく、干渉で識別された位相位置を「phase cell」として読む',fontsize=12,ha='left')
ax.set_xlim(-5.2,5.2); ax.set_ylim(-5.3,4.2)
ax.set_title('図2　周期波を循環FIFO（phase-delay memory）として読む模式図', fontsize=18, pad=16)
save(fig, 'figure02_wave_period_fifo_phase_delay_memory')

# ----------------------------------------------------------------------
# Figure 3: frequency comb <-> time-domain localization
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(12, 7.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1,1.25], hspace=0.46, wspace=0.28)
Ms = [1,3,5,15]
for col, Mv in enumerate([1,15]):
    ax = fig.add_subplot(gs[0,col])
    freqs = np.arange(1,Mv+1)
    amps = np.ones_like(freqs)
    markerline, stemlines, baseline = ax.stem(freqs, amps, basefmt=' ')
    plt.setp(stemlines, linewidth=1.2)
    plt.setp(markerline, markersize=4)
    ax.set_ylim(0,1.25)
    ax.set_xlim(0.5, max(2,Mv+0.5))
    ax.set_title(f'周波数領域：M={Mv}')
    ax.set_xlabel(r'周波数次数  $m\omega_0$')
    ax.set_ylabel('相対振幅')
    ax.grid(alpha=0.2)

ax = fig.add_subplot(gs[1,:])
for Mv in Ms:
    y=np.sum([np.cos(m*phi) for m in range(1,Mv+1)],axis=0)/Mv
    ax.plot(phi/(2*np.pi),y,lw=1.6,label=f'M={Mv}')
ax.set_xlim(0,1)
ax.set_xlabel('ベース周期内部 t/T')
ax.set_ylabel('正規化合成振幅')
ax.set_title('時間領域：帯域（倍音数）を増やすと局所構造が細くなる')
ax.legend(ncol=4,loc='upper right')
ax.grid(alpha=0.25)
fig.suptitle('図3　周波数コムと時間領域構造の対応', fontsize=18, y=0.98)
save(fig, 'figure03_frequency_comb_to_time_domain_resolution')

# ----------------------------------------------------------------------
# Figure 4: odd/even decomposition and readout gauge
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11.5,7.2))
ax.axis('off')
ax.set_xlim(0,10); ax.set_ylim(0,10)
ax.text(5,9.5,'図4　奇数系列・偶数系列の分解：4π carrier と 2π 読出しゲージ',ha='center',fontsize=18,weight='bold')
ax.text(5,8.55,r'$V=e^{i\phi/2},\qquad U=V^2=e^{i\phi}$',ha='center',fontsize=18)

# Odd block
odd = FancyBboxPatch((0.7,4.2),4.0,3.2,boxstyle='round,pad=0.18',fill=False,lw=1.7)
ax.add_patch(odd)
ax.text(2.7,6.95,'奇数半整数系列',ha='center',fontsize=15,weight='bold')
ax.text(2.7,6.2,r'$V,\;V^3,\;V^5,\;\ldots$',ha='center',fontsize=17)
ax.text(2.7,5.45,r'$V^{2k+1}=V\,U^k$',ha='center',fontsize=17)
ax.text(2.7,4.75,'共通因子 V が 4π 二重被覆を保持',ha='center',fontsize=12)

# Even block
even = FancyBboxPatch((5.3,4.2),4.0,3.2,boxstyle='round,pad=0.18',fill=False,lw=1.7)
ax.add_patch(even)
ax.text(7.3,6.95,'偶数系列',ha='center',fontsize=15,weight='bold')
ax.text(7.3,6.2,r'$V^2,\;V^4,\;V^6,\;\ldots$',ha='center',fontsize=17)
ax.text(7.3,5.45,r'$V^{2k}=U^k$',ha='center',fontsize=17)
ax.text(7.3,4.75,'2π 内部読出しゲージだけが残る',ha='center',fontsize=12)

ax.add_patch(FancyArrowPatch((2.7,4.0),(5.0,2.6),arrowstyle='-|>',mutation_scale=16,lw=1.4))
ax.add_patch(FancyArrowPatch((7.3,4.0),(5.0,2.6),arrowstyle='-|>',mutation_scale=16,lw=1.4))
ax.text(5,2.35,r'ベース波との干渉で $U^k$ を読み出す',ha='center',fontsize=15,weight='bold')
ax.text(5,1.7,'→ 一周期内部の位相位置（履歴セル）の識別へ',ha='center',fontsize=14)
ax.text(5,0.85,'解釈：奇数系列 = 4π carrier × intra-period readout gauge',ha='center',fontsize=13)
save(fig, 'figure04_odd_even_harmonic_roles_and_readout_gauge')

# README
readme = f'''# 倍音干渉による有限履歴読出し：説明図 v1\n\n生成スクリプト: `generate_harmonic_history_figures_v1.py`\n\n## 生成図\n\n1. `figure01_monochromatic_vs_harmonic_time_resolution`\n   - 単色波と倍音波の一周期内部時間分解能の比較。\n2. `figure02_wave_period_fifo_phase_delay_memory`\n   - 周期波を循環FIFO / phase-delay memory として読む模式図。\n3. `figure03_frequency_comb_to_time_domain_resolution`\n   - 周波数コムの帯域増加と時間領域局所構造の対応。\n4. `figure04_odd_even_harmonic_roles_and_readout_gauge`\n   - 奇数半整数系列と偶数系列の役割分解。\n\n## 再現\n\n```bash\npython generate_harmonic_history_figures_v1.py\n```\n\nPython依存: numpy, matplotlib。\nすべての図は同じスクリプトから PNG / SVG を同時生成する。\n\n## 注意\n\n図1・図3は解析式から生成した可視化であり、実測データではない。\n図2・図4は概念構造を説明する模式図である。\n'''
(OUT/'README_ja.md').write_text(readme,encoding='utf-8')
print('generated:')
for p in sorted(OUT.iterdir()):
    print(p.name)
