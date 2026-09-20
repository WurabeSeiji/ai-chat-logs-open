import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib import font_manager as fm
import numpy as np

BASE = Path('/mnt/data/第一思考実験_等価原理_20260920')
FIG = BASE / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

FONT = '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
fp = fm.FontProperties(fname=FONT)


def rotation_matrix(n: int) -> np.ndarray:
    theta = 2 * math.pi / n
    return np.array([
        [math.cos(theta), -math.sin(theta)],
        [math.sin(theta),  math.cos(theta)],
    ])


def orbit_points(n: int, R: float = 1.0):
    """Return S_n and points X_0,...,X_n for X_{k+1}=S_n X_k."""
    S = rotation_matrix(n)
    X = np.array([R, 0.0])
    pts = [X.copy()]
    for _ in range(n):
        X = S @ X
        X = np.where(np.abs(X) < 1e-12, 0.0, X)
        pts.append(X.copy())
    return S, np.array(pts)


def add_label(ax, x, y, s, size=12, **kwargs):
    ax.text(x, y, s, fontproperties=fp, fontsize=size, **kwargs)


def matrix_block_text(n: int) -> str:
    if n == 2:
        return (
            'S2 = [[-1, 0], [0, -1]]\n'
            'θ = 2π/2 = π\n'
            'X_k = S2^k X0,   X0 = (1,0)^T\n'
            'z_k = a_k + i b_k = exp(i k pi)\n'
            'z1 = -1,   z2 = 1'
        )
    elif n == 6:
        return (
            'S6 = [[1/2, -sqrt(3)/2], [sqrt(3)/2, 1/2]]\n'
            'θ = 2π/6 = π/3\n'
            'X_k = S6^k X0,   X0 = (1,0)^T\n'
            'z_k = a_k + i b_k = exp(i k pi/3)\n'
            'z1,...,z6 は単位円上の正六角形をなす'
        )
    return ''


def plot_orbit(n: int, R: float = 1.0):
    S, pts = orbit_points(n, R)

    fig, ax = plt.subplots(figsize=(8.8, 8.8), dpi=220)
    ax.set_aspect('equal')
    lim = 1.48 * R
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, linewidth=1.0)
    ax.axvline(0, linewidth=1.0)

    circle = Circle((0, 0), R, fill=False, linewidth=2.0, linestyle='--')
    ax.add_patch(circle)

    x0, y0 = pts[0]
    ax.scatter([x0], [y0], s=70, marker='s')
    add_label(ax, x0 + 0.05, y0 - 0.14, '初期点 $z_0=(1,0)$', size=11)

    for k in range(1, n + 1):
        x_prev, y_prev = pts[k - 1]
        x, y = pts[k]
        arr = FancyArrowPatch((x_prev, y_prev), (x, y), arrowstyle='->', mutation_scale=14, linewidth=1.6)
        ax.add_patch(arr)
        ax.scatter([x], [y], s=55)
        offset_x = 0.04 if x <= 1.1 else -0.1
        offset_y = 0.04 if y <= 1.1 else -0.1
        add_label(ax, x + offset_x, y + offset_y, f'$z_{k}$', size=12)

    ax.set_xlabel('Re$(z)=a$', fontproperties=fp, fontsize=12)
    ax.set_ylabel('Im$(z)=b$', fontproperties=fp, fontsize=12)
    ax.set_title(f'離散位相相互作用の軌道（n={n}, R={R:g}）', fontproperties=fp, fontsize=15)

    ax.text(-1.43 * R, -1.23 * R, matrix_block_text(n), fontsize=11,
            fontproperties=fp,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.88))

    caption = (
        f'図：$S_{n}$ を配列（行列）として定義し、$X_0=(1,0)^T$ に順次作用させた軌道。\n'
        f'すべての点で $a_k^2+b_k^2=R^2=1$ が保たれる。'
    )
    add_label(ax, -1.43 * R, -1.50 * R, caption, size=10)

    stem = f'discrete_phase_orbit_n{n}_ja'
    png = FIG / f'{stem}.png'
    svg = FIG / f'{stem}.svg'
    fig.savefig(png, bbox_inches='tight')
    fig.savefig(svg, bbox_inches='tight')
    plt.close(fig)
    return png, svg, S, pts


def main():
    outputs = []
    for n in (2, 6):
        png, svg, S, pts = plot_orbit(n)
        outputs.extend([str(png), str(svg)])
        print(f'n={n}')
        print('S=')
        print(np.array2string(S, precision=6, suppress_small=True))
        print('points=')
        print(np.array2string(pts, precision=6, suppress_small=True))
    print('\nGenerated files:')
    print('\n'.join(outputs))


if __name__ == '__main__':
    main()
