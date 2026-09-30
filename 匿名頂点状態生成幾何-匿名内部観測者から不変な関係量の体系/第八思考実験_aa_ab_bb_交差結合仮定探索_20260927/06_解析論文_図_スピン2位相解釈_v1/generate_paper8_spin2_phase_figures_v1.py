from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Arc
from matplotlib import font_manager

OUT = Path(__file__).resolve().parent
FONT = 'Noto Sans CJK JP'
plt.rcParams['font.family'] = FONT
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['axes.unicode_minus'] = False


def save(fig, stem):
    fig.savefig(OUT / f'{stem}.svg', bbox_inches='tight')
    fig.savefig(OUT / f'{stem}.png', dpi=240, bbox_inches='tight')
    plt.close(fig)


def add_box(ax, xy, wh, title, body='', lw=1.6, fc='white'):
    x, y = xy; w, h = wh
    p = FancyBboxPatch((x,y), w,h, boxstyle='round,pad=0.02,rounding_size=0.025',
                       linewidth=lw, edgecolor='black', facecolor=fc)
    ax.add_patch(p)
    ax.text(x+w/2, y+h*0.68, title, ha='center', va='center', fontsize=15, weight='bold')
    if body:
        ax.text(x+w/2, y+h*0.34, body, ha='center', va='center', fontsize=11, linespacing=1.35)
    return p


def arrow(ax, p1, p2, text=None, rad=0.0, style='-|>'):
    a = FancyArrowPatch(p1,p2,arrowstyle=style,mutation_scale=14,linewidth=1.5,
                        connectionstyle=f'arc3,rad={rad}',color='black')
    ax.add_patch(a)
    if text:
        xm=(p1[0]+p2[0])/2; ym=(p1[1]+p2[1])/2
        ax.text(xm,ym+0.025,text,ha='center',va='bottom',fontsize=10)


def figure1():
    fig, ax = plt.subplots(figsize=(13,7.5))
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    ax.text(0.5,0.95,'図1　第7論文から第8論文へ：交差項と放射安定性の未解決問題',
            ha='center',va='center',fontsize=19,weight='bold')

    # left: state vector
    add_box(ax,(0.05,0.58),(0.22,0.25),r'$\Psi=(aa,\,ab,\,bb)^T$',
            '3つの関係状態を\n1本の長い状態ベクトルへ')

    # center: block diagonal S
    add_box(ax,(0.36,0.56),(0.28,0.29),r'$S_{\rm diag}$',
            'S = diag(S_aa, S_ab, S_bb)\n\n交差ブロック = 0')
    arrow(ax,(0.27,0.70),(0.36,0.70),'作用')

    # right: same as independent
    add_box(ax,(0.73,0.58),(0.22,0.25),'個別相互作用と同一',
            '$aa,ab,bb$ を\n別々に走らせた結果と一致')
    arrow(ax,(0.64,0.70),(0.73,0.70),'第7論文で確認')

    # unresolved cross terms
    ax.text(0.50,0.47,'しかし',ha='center',va='center',fontsize=14,weight='bold')
    add_box(ax,(0.28,0.26),(0.20,0.14),'未決定 1','非対角ブロック\n（交差相互作用）')
    add_box(ax,(0.52,0.26),(0.20,0.14),'未決定 2','放射損失を持つため\n安定状態が未説明')

    # radiation arrows
    for x, label in [(0.17,'aa'),(0.50,'ab'),(0.83,'bb')]:
        y0=0.14
        ax.add_patch(Circle((x,y0),0.035,fill=False,lw=1.4))
        ax.text(x,y0,label,ha='center',va='center',fontsize=12,weight='bold')
        arrow(ax,(x,y0+0.035),(x,0.22))
    ax.text(0.17,0.085,'GW',ha='center',fontsize=11)
    ax.text(0.50,0.085,'GW + EM',ha='center',fontsize=11)
    ax.text(0.83,0.085,'GW',ha='center',fontsize=11)

    ax.text(0.50,0.03,
            '第8論文の出発点：交差項を適切に設定すれば、放射を内部交換として閉じられるのではないか？',
            ha='center',va='bottom',fontsize=13,weight='bold')
    save(fig,'figure01_paper7_to_paper8_cross_terms_and_radiation_problem')


def figure2():
    fig, ax = plt.subplots(figsize=(13,7.5))
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off')
    ax.text(0.5,0.95,'図2　自己関係 $aa$ の再解釈：空間半径ではなく発展方向の位相勾配',
            ha='center',va='center',fontsize=19,weight='bold')

    # left timeline
    ax.text(0.17,0.85,'1点でも発展方向がある',ha='center',fontsize=14,weight='bold')
    xs=np.linspace(0.06,0.30,5); y=0.66
    for i,x in enumerate(xs):
        ax.add_patch(Circle((x,y),0.028,fill=False,lw=1.6))
        ax.text(x,y,f'$a_{i}$',ha='center',va='center',fontsize=11)
        if i < len(xs)-1:
            arrow(ax,(x+0.028,y),(xs[i+1]-0.028,y))
    ax.text(0.18,0.57,r'$a_{k+1}=e^{i\Delta\phi}\,a_k$',ha='center',fontsize=15)
    ax.text(0.18,0.50,r'$\Delta\phi=\arg(a_{k+1}/a_k)$',ha='center',fontsize=15)
    ax.text(0.18,0.43,'空間的な広がりがなくても\n処理ステップ方向の位相傾斜は定義できる',
            ha='center',va='center',fontsize=11)

    # middle complex plane illustration
    ax.text(0.50,0.85,'自己位相の回転',ha='center',fontsize=14,weight='bold')
    cx,cy,R=0.50,0.64,0.13
    ax.add_patch(Circle((cx,cy),R,fill=False,lw=1.4))
    for th in [20,55,90]:
        t=np.deg2rad(th)
        arrow(ax,(cx,cy),(cx+R*np.cos(t),cy+R*np.sin(t)))
    ax.add_patch(Arc((cx,cy),0.17,0.17,theta1=20,theta2=55,lw=1.2))
    ax.text(cx+0.10,cy+0.07,r'$\Delta\phi$',fontsize=12)
    ax.text(cx,0.43,'絶対角ではなく\n連続する自己状態間の位相差',ha='center',fontsize=11)

    # right mapping to spin2
    ax.text(0.82,0.85,'二次自己関係',ha='center',fontsize=14,weight='bold')
    add_box(ax,(0.71,0.61),(0.22,0.15),r'$a\sim e^{i\phi}$','一次位相')
    arrow(ax,(0.82,0.61),(0.82,0.52),'自己関係 $aa=a^2$')
    add_box(ax,(0.71,0.36),(0.22,0.15),r'$aa=a^2\sim e^{i2\phi}$','位相が2倍で変換')
    ax.text(0.82,0.26,r'$aa\rightarrow e^{i2\theta}aa$',ha='center',fontsize=16,weight='bold')
    ax.text(0.82,0.19,'weight 2 / spin-2 型の変換則',ha='center',fontsize=12)

    # bottom conclusion
    ax.plot([0.05,0.95],[0.12,0.12],lw=1.0,color='black')
    ax.text(0.50,0.065,
            '再解釈：$aa$ の「軌道」は自己距離ではなく、発展方向の位相回転として読める',
            ha='center',va='center',fontsize=14,weight='bold')
    save(fig,'figure02_self_relation_phase_gradient_spin2_interpretation')


if __name__ == '__main__':
    figure1()
    figure2()
    print('generated:')
    for p in sorted(OUT.glob('figure0*')):
        print(p.name)
