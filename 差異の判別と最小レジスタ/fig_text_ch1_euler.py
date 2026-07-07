# -*- coding: utf-8 -*-
"""図: e^{iθ} 記法の幾何 — (a) 積=角の足し算 (b) 絶対値は掛け算・偏角は足し算 (c) 共役=実軸折り返し"""
import numpy as np
import matplotlib
matplotlib.rcParams['font.family'] = 'Hiragino Sans'
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.8))

def setup(ax, lim):
    ax.axhline(0, color='gray', lw=0.8, zorder=0)
    ax.axvline(0, color='gray', lw=0.8, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel('実部', fontsize=10)
    ax.set_ylabel('虚部', fontsize=10)
    ax.tick_params(labelsize=8)

t = np.linspace(0, 2 * np.pi, 300)
th1, th2 = np.deg2rad(40), np.deg2rad(70)

# ---------- (a) 単位円上: 積 = 角の足し算 ----------
ax = axes[0]
setup(ax, 1.45)
ax.plot(np.cos(t), np.sin(t), color='lightgray', lw=1.2, zorder=0)
items = [(th1, 'C0', '$e^{i40°}$'),
         (th2, 'C2', '$e^{i70°}$'),
         (th1 + th2, 'C3', '積 $e^{i110°}$')]
for th, c, lab in items:
    ax.annotate('', xy=(np.cos(th), np.sin(th)), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=c, lw=2.2))
    ax.text(1.22 * np.cos(th), 1.22 * np.sin(th), lab, color=c,
            ha='center', va='center', fontsize=10)
ax.add_patch(Arc((0, 0), 0.75, 0.75, theta1=0, theta2=40, color='C0', lw=2))
ax.add_patch(Arc((0, 0), 0.95, 0.95, theta1=40, theta2=110, color='C2', lw=2))
ax.text(0.45, 0.12, '40°', color='C0', fontsize=9)
ax.text(0.16, 0.52, '+70°', color='C2', fontsize=9)
ax.set_title('(a) 掛け算＝回転の合成\n40° + 70° = 110°', fontsize=11)

# ---------- (b) 一般: 絶対値は掛け算・偏角は足し算 ----------
ax = axes[1]
setup(ax, 2.15)
z1 = 1.5 * np.exp(1j * th1)
z2 = 1.2 * np.exp(1j * th2)
z3 = z1 * z2
for r in (1.2, 1.5, 1.8):
    ax.plot(r * np.cos(t), r * np.sin(t), color='lightgray', lw=0.7, ls=':', zorder=0)
for z, c, lw in [(z1, 'C0', 2.2), (z2, 'C2', 2.2), (z3, 'C3', 2.6)]:
    ax.annotate('', xy=(z.real, z.imag), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=c, lw=lw))
ax.text(z1.real + 0.42, z1.imag - 0.05, '$Z_1 = 1.5\\,e^{i40°}$', color='C0', fontsize=9, ha='center')
ax.text(z2.real + 0.55, z2.imag + 0.18, '$Z_2 = 1.2\\,e^{i70°}$', color='C2', fontsize=9, ha='center')
ax.text(z3.real - 0.15, z3.imag + 0.22, '$Z_1 Z_2 = 1.8\\,e^{i110°}$', color='C3', fontsize=9, ha='center')
ax.set_title('(b) 絶対値は掛け算 1.5×1.2=1.8\n偏角は足し算 40°+70°=110°', fontsize=11)

# ---------- (c) 共役 = 実軸折り返し ----------
ax = axes[2]
setup(ax, 1.85)
thc = np.deg2rad(35)
z = 1.4 * np.exp(1j * thc)
ax.plot(1.4 * np.cos(t), 1.4 * np.sin(t), color='lightgray', lw=1.0, zorder=0)
ax.annotate('', xy=(z.real, z.imag), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color='C0', lw=2.4))
ax.annotate('', xy=(z.real, -z.imag), xytext=(0, 0),
            arrowprops=dict(arrowstyle='->', color='C4', lw=2.4))
ax.text(z.real + 0.12, z.imag + 0.16, '$Z = \\rho\\,e^{i\\theta}$', color='C0', fontsize=11)
ax.text(z.real + 0.12, -z.imag - 0.20, '$\\bar{Z} = \\rho\\,e^{-i\\theta}$', color='C4', fontsize=11)
ax.add_patch(Arc((0, 0), 0.8, 0.8, theta1=0, theta2=35, color='C0', lw=1.6))
ax.add_patch(Arc((0, 0), 0.8, 0.8, theta1=-35, theta2=0, color='C4', lw=1.6))
ax.text(0.48, 0.14, '$+\\theta$', color='C0', fontsize=10)
ax.text(0.48, -0.24, '$-\\theta$', color='C4', fontsize=10)
ax.annotate('', xy=(z.real, -z.imag + 0.06), xytext=(z.real, z.imag - 0.06),
            arrowprops=dict(arrowstyle='<->', color='gray', lw=1.2, ls='--'))
ax.text(z.real + 0.08, 0.0, '実軸で\n折り返し', color='gray', fontsize=9, va='center')
ax.set_title('(c) 共役 $\\bar{Z}$＝実軸で折り返し\n（逆向きに回る相棒）', fontsize=11)

plt.tight_layout()
plt.savefig('fig_text_ch1_euler.png', dpi=160, bbox_inches='tight')
plt.savefig('fig_text_ch1_euler.svg', bbox_inches='tight')
print('saved')
