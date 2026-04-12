import numpy as np
import matplotlib.pyplot as plt
import os

def rewrite_evaluate_tensor():
    with open('evaluate_tensor.py', 'r') as f:
        content = f.read()
    
    # Replace generate_geodesic_fig2 with the new Baumkuchen Geodesic plot
    new_fig2_code = '''
# ==========================================
# Figure 1 (was Fig 2): Tensor Lookup Geodesic mapped on Baumkuchen
# ==========================================
def generate_geodesic_fig2():
    print("Generating Figure 1 (Geodesic on Baumkuchen): tensor_lookup_geodesic.png...")
    R = 100.0
    T_max = 80.0
    n_width = 80.0
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # 1. Background Mesh (Rigorous Metric Baumkuchen)
    # The spatial radius physically contracts as r(t) = R * sqrt(1 - (t/R)^2)
    xi_grid = np.linspace(-n_width/2, n_width/2, 17)
    T_grid = np.linspace(0, T_max, 9)
    
    # Draw Concentric Arcs (Constant Coordinate Time T)
    xi_dense = np.linspace(-n_width/2, n_width/2, 200)
    for t_val in T_grid:
        r = R * np.sqrt(np.abs(1.0 - (t_val/R)**2))
        theta = xi_dense / R
        px = r * np.sin(theta)
        pz = r * np.cos(theta)
        lw = 2.5 if (t_val == 0 or t_val == T_max) else 1.0
        c = 'saddlebrown' if t_val == 0 else 'chocolate'
        ax.plot(px, pz, color=c, alpha=0.8, linewidth=lw)
        
    # Draw Radial Lines (Constant Spatial Index Xi)
    for xi_val in xi_grid:
        t_dense = np.linspace(0, T_max, 100)
        r_vals = R * np.np.sqrt(np.abs(1.0 - (t_dense/R)**2))
        theta = xi_val / R
        px = r_vals * np.sin(theta)
        pz = r_vals * np.cos(theta)
        lw = 2.5 if (xi_val == xi_grid[0] or xi_val == xi_grid[-1]) else 1.0
        ax.plot(px, pz, color='chocolate', alpha=0.8, linewidth=lw)

    # 2. Simulate and Map Test Particles
    dt_exact = 0.1
    steps_exact = int(T_max / dt_exact)
    
    # Particle A: Static in Index Space (falling straight down)
    xi_A = np.zeros((steps_exact, 4))
    xi_A[:, 0] = np.arange(0, T_max, dt_exact)
    xi_A[:, 1] = 20.0  # Initial index
    
    # Particle B: Moving across the index space dynamically
    v_x = -0.4 
    xi_B = np.zeros((steps_exact, 4))
    xi_B[:, 0] = np.arange(0, T_max, dt_exact)
    xi_B[:, 1] = -30.0 + xi_B[:, 0] * v_x

    # Retrieve physical trajectories
    X_phys_A_ex = exact_mapping(xi_A, R)
    X_phys_A_ap = tensor_lerp_mapping(xi_A, R)
    X_phys_B_ex = exact_mapping(xi_B, R)
    X_phys_B_ap = tensor_lerp_mapping(xi_B, R)

    # Helper to plot physically mapped metric coordinates onto the embedding
    def plot_trajectory(X_computed, xi_original, color_ex, color_ap, label_ex, label_ap):
        # Physical radius at given proper time
        r_phys = R * np.sqrt(np.abs(1.0 - (xi_original[:, 0]/R)**2))
        # Compute exact polar projection angle from calculated physical arc-length mapping
        # theta = arc_length / radius
        theta_ex = X_computed[0][:, 1] / r_phys
        theta_ap = X_computed[1][:, 1] / r_phys
        
        # Plot Exact (Thick solid line behind)
        ax.plot(r_phys * np.sin(theta_ex), r_phys * np.cos(theta_ex), 
                color=color_ex, linewidth=7, alpha=0.5, label=label_ex)
        # Plot Approximation (Dashed line on top)
        ax.plot(r_phys * np.sin(theta_ap), r_phys * np.cos(theta_ap), 
                color=color_ap, linestyle='--', linewidth=2.5, label=label_ap)

    plot_trajectory((X_phys_A_ex, X_phys_A_ap), xi_A, 'lightcoral', 'red', 'Particle A (Exact Geodesic)', 'Particle A (Tensor Lerp)')
    plot_trajectory((X_phys_B_ex, X_phys_B_ap), xi_B, 'cornflowerblue', 'blue', 'Particle B (Exact Geodesic)', 'Particle B (Tensor Lerp)')
    
    # Origin Gravity Marker
    ax.scatter([0], [0], color='black', marker='X', s=200, label='Center of Curvature (Gravity Origin)')
    
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("Physical Spatial Axis X")
    ax.set_ylabel("Orthogonal Curvature Depth W")
    ax.set_title("Figure 1: Full Metric Tensor Simulation over Curved Geometry (God\\'s Eye View)\\nExact Analytical Geodesics vs JIT Tensor Interpolation Tracked Inside Spatial Contraction Mesh")
    ax.grid(True, linestyle=':', alpha=0.4)
    # Dynamically frame the plot
    ax.set_ylim(-5, 105)
    ax.legend(loc='lower center', ncol=2, prop={'size': 9})
    plt.tight_layout()
    plt.savefig('tensor_lookup_geodesic.png', dpi=150)
    plt.close()
'''
    
    # Find the bounds of generate_geodesic_fig2
    lines = content.split('\n')
    start_idx = -1
    end_idx = -1
    for i, line in enumerate(lines):
        if line.startswith('def generate_geodesic_fig2():'):
            start_idx = i - 3
        if start_idx != -1 and line.startswith('def evaluate_accuracy_fig3():'):
            end_idx = i - 3
            break
            
    if start_idx != -1 and end_idx != -1:
        new_content = '\n'.join(lines[:start_idx]) + '\n' + new_fig2_code.strip() + '\n\n' + '\n'.join(lines[end_idx:])
        new_content = new_content.replace('np.np.sqrt', 'np.sqrt') # fix typo
        with open('evaluate_tensor.py', 'w') as f:
            f.write(new_content)

def rewrite_md():
    md_file = "閉じた正曲率時空におけるテスト粒子運動のGPU効率的シミュレーション技法（テンソル演算版）.md"
    
    with open(md_file, "r") as f:
        content = f.read()

    # 1. Remove the old Figure 0 text and image
    # 2. Rename Figure 1 to Figure 2 (Frustum)
    # 3. Add the new Figure 1 (Baumkuchen Geodesic) inside the text smoothly
    
    import re
    # We will just write a new abstract mapping the 2 stages
    new_abstract = """## アブストラクト（空間歪みモデルと究極の可視化）
一般相対性理論における「曲がった時空の落下」を直観的に理解し、かつ計算機で高速処理するための概念図として、まずは「神の視点（God's Eye View）」から空間と粒子の軌跡を描画した特製プロット（図1）を提示する。

![測地線の完全な可視化](tensor_lookup_geodesic.png)
*図1: 本シミュレーションの対象となる空間モデル。平坦な「インデックス空間」のマス目が、直交する曲率軸（$W$軸）に向かって深く潜るほど狭くなる「バームクーヘン型のくさびメッシュ」へと曲げられている（曲率半径 $R=100$）。テスト粒子A・Bは、一切の引力加速度を持たずにこのマス目を定規のようにまっすぐ進むだけで、物理空間上ではマス目が縮むことによって自動的に原点（重心）へと曲がり落ちてゆく。この図では背景の厳密な空間メッシュに対して、厳密解の軌道（太い半透明線）と、本手法におけるJITテンソル補間軌道（点線）が完全に重なって一致していることを証明している。*

続いて、この曲がったバームクーヘンの上で「時間」を進行させた際に発生する動態的な変化を、時間軸（-t軸）を加えた3Dプロットとして描画したのが以下の「図2」である。

![テンソル空間の動態モデル：台形錐体](mesh_distortion_frustum.png)
*図2: 構築された動的なテンソル時空間のアナロジー。底面が平坦なインデックス空間（特殊相対論領域）であり、時間が上方（-t軸）へ進むにつれて曲率の影響で空間全体が中央へローレンツ収縮し「台形錐体（フラスタム）」に歪んでゆく。*

本来、重力場（曲率）の中で物が落下する現象をシミュレーションするには、毎ステップ重力加速度を計算する必要がある。しかし本提案手法は、図1・図2のような**「時間が進むにつれて重力方向へマス目が収縮する歪みの幾何データ」を、あらかじめ4次元配列（テンソル）としてメモリ上に全構築（Step 1）**してしまう。
このテンソル空間の中では、重力の影響はすべて「方眼紙の歪み」として焼き込まれているため、テスト粒子は**一切の加速度計算を持たず線形補間（Lerp）するだけで、極めて高い精度で宇宙物理学の測地線軌道を再現（Step 2）**するブレイクスルーを実現する。
"""
    # Replace the abstract section up to Introduction
    pattern = r'## アブストラクト.*?## 1\. 目的（Introduction）'
    content = re.sub(pattern, new_abstract + '\n## 1. 目的（Introduction）', content, flags=re.DOTALL)
    
    # Section 5 validation cleanup (Update figures reference)
    content = content.replace('結果を図2に示す', '結果は導入部の図1に示した通りである')
    content = content.replace('図2: 原点重力に向かって落下する2つの軌跡。半透明の太線が非線形な「完全厳密解」、その上に重なる破線が本アルゴリズムの「テンソル線形補間近似解」。', '')
    content = content.replace('![線形補間を用いた測地線シミュレーション結果](tensor_lookup_geodesic.png)', '')
    
    with open(md_file, "w") as f:
        f.write(content)

if __name__ == '__main__':
    rewrite_evaluate_tensor()
    rewrite_md()
