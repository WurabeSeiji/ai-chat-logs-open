import numpy as np
import matplotlib.pyplot as plt

def build_tensor_space(nt, nx, ny, R, c=1.0):
    T_mesh = np.zeros((nt, nx, ny, 3))
    for it in range(nt):
        v = (c * it) / R
        v_clip = min(v, 0.999 * c)
        gamma = 1.0 / np.sqrt(1.0 - (v_clip / c)**2)
        for ix in range(nx):
            for iy in range(ny):
                X_phys = ix * (1.0 / gamma) 
                Y_phys = iy * (1.0 / gamma)
                T_phys = it * gamma
                T_mesh[it, ix, iy, 0] = X_phys
                T_mesh[it, ix, iy, 1] = Y_phys
                T_mesh[it, ix, iy, 2] = T_phys
    return T_mesh

# ====== テスト条件 ======
R = 100.0   
nt = 80     
nx, ny = 20, 20
c = 1.0

print("Step 1: 物理空間テンソル T_mesh を構築中...")
T_mesh = build_tensor_space(nt, nx, ny, R, c)
print(f"-> テンソル構築完了。配列サイズ: {T_mesh.shape}")

def simulate_particle(T_mesh, start_idx, v_local):
    path_phys_X = []
    path_phys_Y = []
    
    curr_x = start_idx[0]
    curr_y = start_idx[1]
    
    for it in range(nt):
        # 1. 局所慣性系での等速直線移動
        curr_x += v_local[0]
        curr_y += v_local[1]
        
        if curr_x >= nx-1 or curr_y >= ny-1 or curr_x < 0 or curr_y < 0:
            break
            
        # 2. テンソル空間からの『線形補間（Bilinear Interpolation）』によるルックアップ
        # ※以前の単なる int() 切り捨ては、不連続なテレポート（ギザギザのノコギリ波）を生むミスでした
        ix = int(np.floor(curr_x))
        iy = int(np.floor(curr_y))
        dx = curr_x - ix
        dy = curr_y - iy
        
        # テンソルの周囲4点を参照
        P00 = T_mesh[it, ix, iy, :2]
        P10 = T_mesh[it, ix+1, iy, :2]
        P01 = T_mesh[it, ix, iy+1, :2]
        P11 = T_mesh[it, ix+1, iy+1, :2]
        
        # 平面（ポリゴン）内での位置を線形計算（完全な連続体としての折れ線軌道を生成）
        phys_pos = (P00 * (1 - dx) * (1 - dy) +
                    P10 * dx * (1 - dy) +
                    P01 * (1 - dx) * dy +
                    P11 * dx * dy)
        
        path_phys_X.append(phys_pos[0])
        path_phys_Y.append(phys_pos[1])
        
    return path_phys_X, path_phys_Y

print("Step 2: テスト粒子を配置し、線形補間ルックアップを開始...")

# 粒子A：静止（v=0）
pathA_x, pathA_y = simulate_particle(T_mesh, start_idx=(18.0, 18.0), v_local=(0.0, 0.0))
# 粒子B：等速横移動（v_x=0.2）
pathB_x, pathB_y = simulate_particle(T_mesh, start_idx=(2.0, 18.0), v_local=(0.2, 0.0))

plt.figure(figsize=(9, 8))
plt.grid(True, linestyle='--', color='gray', alpha=0.3)

plt.plot(pathA_x, pathA_y, marker='s', markersize=3, label='Particle A (Static / v=0)', color='red', alpha=0.7)
plt.scatter([pathA_x[0]], [pathA_y[0]], color='darkred', s=100, zorder=5, label='Start A')

plt.plot(pathB_x, pathB_y, marker='o', markersize=3, label='Particle B (Moving SR / v=0.2)', color='blue', alpha=0.7)
plt.scatter([pathB_x[0]], [pathB_y[0]], color='darkblue', s=100, zorder=5, label='Start B')

plt.title('Step 2 Simulation: Tensor-mapped Geodesics (Corrected Continuous Interpolation)', fontsize=14)
plt.xlabel('Physical X (Mapped via Tensor)', fontsize=12)
plt.ylabel('Physical Y (Mapped via Tensor)', fontsize=12)
plt.legend(fontsize=11)
plt.xlim(0, 20)
plt.ylim(0, 20)

plt.tight_layout()
plt.savefig('tensor_lookup_geodesic.png', dpi=150)
print("Saved tensor_lookup_geodesic.png")
