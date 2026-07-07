# -*- coding: utf-8 -*-
"""図4: 一様に回る点を真横から見ると端に溜まる（逆正弦則）"""
import numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = 'Hiragino Sans'
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

# ---------- (a) 円周上の等間隔点を u 軸へ射影 ----------
ax = axes[0]
t = np.linspace(0, 2 * np.pi, 300)
ax.plot(np.cos(t), np.sin(t), color='lightgray', lw=1.2)
th = np.deg2rad(np.arange(5, 360, 10))  # 等間隔 36 点
ax.plot(np.cos(th), np.sin(th), 'o', color='C0', ms=5, label='等間隔に回る点')
u_axis_y = -1.55
for a in th[np.sin(th) >= 0]:  # 上半分だけ射影線（見やすさ）
    ax.plot([np.cos(a), np.cos(a)], [np.sin(a), u_axis_y + 0.05],
            color='C0', lw=0.5, ls=':', alpha=0.55)
ax.axhline(u_axis_y, color='gray', lw=1)
ax.plot(np.cos(th), np.full_like(th, u_axis_y), '|', color='C3', ms=14, mew=1.6)
ax.text(0, u_axis_y - 0.22, '真横から見た位置 $u=\\cos\\theta$（端ほど密）',
        ha='center', fontsize=10, color='C3')
ax.text(1.13, u_axis_y + 0.13, '$+1$', ha='center', fontsize=9)
ax.text(-1.13, u_axis_y + 0.13, '$-1$', ha='center', fontsize=9)
ax.set_aspect('equal')
ax.set_xlim(-1.7, 1.7)
ax.set_ylim(-2.1, 1.4)
ax.axis('off')
ax.set_title('(a) 円周上を一様に回る点を、\n真横（$u$ 軸）へ射影する', fontsize=11)

# ---------- (b) 溜まり方の形 = 逆正弦則 ----------
ax = axes[1]
th_fine = np.linspace(0, 2 * np.pi, 36001)[:-1]
u_samples = np.cos(th_fine)
ax.hist(u_samples, bins=40, density=True, color='C0', alpha=0.45,
        label='ヒストグラム（一様な $\\theta$ から）')
u = np.linspace(-0.999, 0.999, 500)
ax.plot(u, 1 / (np.pi * np.sqrt(1 - u ** 2)), color='C3', lw=2.2,
        label='$p(u)=\\dfrac{1}{\\pi\\sqrt{1-u^2}}$')
ax.set_xlabel('$u = \\cos\\theta$', fontsize=11)
ax.set_ylabel('頻度（密度）', fontsize=10)
ax.set_ylim(0, 2.4)
ax.legend(fontsize=10, loc='upper center')
ax.set_title('(b) 端 $u=\\pm1$ に溜まる——逆正弦則', fontsize=11)

plt.tight_layout()
plt.savefig('fig_text_ch1_arcsine.png', dpi=160, bbox_inches='tight')
plt.savefig('fig_text_ch1_arcsine.svg', bbox_inches='tight')
print('saved')
