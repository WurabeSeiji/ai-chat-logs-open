import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Arc
from matplotlib import font_manager as fm
from pathlib import Path

BASE = Path('/mnt/data/第一思考実験_等価原理_20260920')
FIG = BASE / 'figures'
FIG.mkdir(parents=True, exist_ok=True)

JP_FONT = '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
EN_FONT = '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'

jp_prop = fm.FontProperties(fname=JP_FONT)
en_prop = fm.FontProperties(fname=EN_FONT)


def add_text(ax, x, y, s, fontprops, size=12, **kwargs):
    ax.text(x, y, s, fontproperties=fontprops, fontsize=size, **kwargs)


def draw_diagram(lang='ja'):
    fig, ax = plt.subplots(figsize=(10, 8), dpi=200)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_xlim(-2.2, 7.2)
    ax.set_ylim(-3.8, 3.8)

    if lang == 'ja':
        fp = jp_prop
        title = '等価原理の思考実験：重力加速度と遠心加速度の釣り合い'
        center_label = '中心質量 A'
        point_label = '質点 B'
        radius_label = '曲率半径 R'
        omega_label = '角速度 ω'
        grav_label = '重力加速度  g(R)=GM/R²'
        cent_label = '遠心加速度  a_c=ω²R'
        balance = '釣り合い条件：  g(R)=a_c'
        balance2 = 'すなわち  GM/R² = ω²R'
        note = '回転系では、外向きの遠心加速度と内向きの重力加速度がつり合う。'
        note2 = '本論文では、このとき保たれる R を曲率半径として読む。'
        caption = '図1　中心 A と円周上の B、半径 R、角速度 ω、遠心加速度と重力加速度の釣り合い'
    else:
        fp = en_prop
        title = 'Thought Experiment for the Equivalence Principle'
        center_label = 'Central mass A'
        point_label = 'Point mass B'
        radius_label = 'Curvature radius R'
        omega_label = 'Angular velocity ω'
        grav_label = 'Gravitational acceleration  g(R)=GM/R²'
        cent_label = 'Centrifugal acceleration  a_c=ω²R'
        balance = 'Balance condition:  g(R)=a_c'
        balance2 = 'that is,  GM/R² = ω²R'
        note = 'In the rotating frame, the outward centrifugal acceleration'
        note2 = 'balances the inward gravitational acceleration.'
        caption = 'Figure 1. A at the center, B on the orbit, radius R, angular velocity ω, and the balance of gravity and centrifugal acceleration.'

    # Geometry
    center = (0, 0)
    R = 4.0
    theta = math.radians(35)
    B = (R * math.cos(theta), R * math.sin(theta))

    # Orbit and bodies
    orbit = Circle(center, R, fill=False, linewidth=2.0)
    bodyA = Circle(center, 0.18, fill=True)
    bodyB = Circle(B, 0.14, fill=True)
    ax.add_patch(orbit)
    ax.add_patch(bodyA)
    ax.add_patch(bodyB)

    # radius line
    ax.plot([center[0], B[0]], [center[1], B[1]], linewidth=1.8)
    mid = ((center[0]+B[0])/2, (center[1]+B[1])/2)
    add_text(ax, mid[0]-0.2, mid[1]+0.2, radius_label, fp, size=12)

    # labels
    add_text(ax, center[0]-0.55, center[1]-0.55, center_label, fp, size=12)
    add_text(ax, B[0]+0.15, B[1]+0.15, point_label, fp, size=12)

    # gravitational acceleration arrow (inward)
    grav_end = (B[0] - 1.6 * math.cos(theta), B[1] - 1.6 * math.sin(theta))
    grav = FancyArrowPatch(B, grav_end, arrowstyle='->', mutation_scale=16, linewidth=2.2)
    ax.add_patch(grav)
    add_text(ax, grav_end[0]-2.2, grav_end[1]-0.2, grav_label, fp, size=12)

    # centrifugal arrow (outward)
    cent_end = (B[0] + 1.6 * math.cos(theta), B[1] + 1.6 * math.sin(theta))
    cent = FancyArrowPatch(B, cent_end, arrowstyle='->', mutation_scale=16, linewidth=2.2)
    ax.add_patch(cent)
    add_text(ax, cent_end[0]-0.1, cent_end[1]-0.1, cent_label, fp, size=12)

    # tangent direction and omega label
    tangent_theta = theta + math.pi/2
    tang_end = (B[0] + 1.4 * math.cos(tangent_theta), B[1] + 1.4 * math.sin(tangent_theta))
    tang = FancyArrowPatch(B, tang_end, arrowstyle='->', mutation_scale=16, linewidth=1.8)
    ax.add_patch(tang)
    add_text(ax, tang_end[0]+0.05, tang_end[1], omega_label, fp, size=12)

    # rotation arc near center
    arc = Arc(center, 2.0, 2.0, angle=0, theta1=10, theta2=80, linewidth=1.6)
    ax.add_patch(arc)
    arrow_head_start = (1.0 * math.cos(math.radians(75)), 1.0 * math.sin(math.radians(75)))
    arrow_head_end = (1.0 * math.cos(math.radians(80)), 1.0 * math.sin(math.radians(80)))
    arc_arrow = FancyArrowPatch(arrow_head_start, arrow_head_end, arrowstyle='->', mutation_scale=14, linewidth=1.6)
    ax.add_patch(arc_arrow)

    # formulas / note
    add_text(ax, -2.0, -2.2, balance, fp, size=15)
    add_text(ax, -2.0, -2.75, balance2, fp, size=15)
    add_text(ax, -2.0, -3.25, note, fp, size=11)
    add_text(ax, -2.0, -3.6, note2, fp, size=11)
    add_text(ax, -2.0, 3.35, title, fp, size=16)
    add_text(ax, -2.0, -3.95, caption, fp, size=10)

    suffix = 'ja' if lang == 'ja' else 'en'
    png_path = FIG / f'equivalence_principle_centrifugal_gravity_{suffix}.png'
    svg_path = FIG / f'equivalence_principle_centrifugal_gravity_{suffix}.svg'
    fig.savefig(png_path, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close(fig)
    return png_path, svg_path


def main():
    outputs = []
    outputs.extend(draw_diagram('ja'))
    outputs.extend(draw_diagram('en'))
    print('\n'.join(str(p) for p in outputs))


if __name__ == '__main__':
    main()
