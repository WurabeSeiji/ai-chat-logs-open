import numpy as np
import matplotlib.pyplot as plt

def map_coordinates(x, y, z, t, R=100.0, C=1.0):
    # 初期速度ベクトル=0、曲率Rの方向に向けた加速度落下
    # 時間tが経過するごとに、曲率による速度vを獲得する
    v = (C * t) / R 
    
    # ローレンツ計算 (Lorentz Transformation constraints)
    # vがCを超えないようクリップし、ガンマ因子を計算
    v_clipped = np.clip(v, 0.0, 0.999 * C)
    gamma = 1.0 / np.sqrt(1.0 - (v_clipped / C)**2)
    
    # 空間の歪み（ローレンツ収縮）
    # 静止系（インデックス空間）の等間隔・正立方体メッシュは、
    # 落下速度が上がるにつれて進行方向に対して空間収縮を起こす
    X = x * np.sqrt(1.0 - (v_clipped / C)**2)
    Y = y * np.sqrt(1.0 - (v_clipped / C)**2)
    Z = z * np.sqrt(1.0 - (v_clipped / C)**2)
    
    # 時間の歪み（ローレンツ時間遅れ）
    # インデックスとしての時間tが、物理的な時間Tへ変換される
    T = t * gamma
    
    return X, Y, Z, T

R = 100.0
n = 10
C = 1.0

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 1. 時間方向へのスライス（x,y 平面）をワイヤーフレームで描画
for t in range(n):
    X_slice = np.zeros((n, n))
    Y_slice = np.zeros((n, n))
    T_slice = np.zeros((n, n))
    for x in range(n):
        for y in range(n):
            X, Y, Z_val, mapped_T = map_coordinates(x, y, 0, t, R, C)
            X_slice[x, y] = X
            Y_slice[x, y] = Y
            # ご指示の通り「-t」の進行方向としてセット
            T_slice[x, y] = -mapped_T
            
    # メッシュの横の線を描画
    ax.plot_wireframe(X_slice, Y_slice, T_slice, color='blue', alpha=0.4, linewidth=1)

# 2. 時間軸方向（縦軸）の赤い線を描画し、台形錐体（フラスタム）の傾きを可視化
corners = [(0, 0), (n-1, 0), (0, n-1), (n-1, n-1)]
for cx, cy in corners:
    X_line, Y_line, T_line = [], [], []
    for t in range(n):
        X, Y, Z_val, mapped_T = map_coordinates(cx, cy, 0, t, R, C)
        X_line.append(X)
        Y_line.append(Y)
        T_line.append(-mapped_T)
    # 四隅の柱を描画（これが上端に向かって内側に傾くことで台形錐体を形成する）
    ax.plot(X_line, Y_line, T_line, color='red', linewidth=3)
    
    # 頂点の座標をコンソールに出力して証明
    if t == 0 or t == n-1:
        print(f"Corner(x={cx}, y={cy}, t={t}) -> Mapped(X={X_line[-1]:.4f}, Y={Y_line[-1]:.4f}, -T={T_line[-1]:.4f})")

ax.set_xlabel('Mapped X Space')
ax.set_ylabel('Mapped Y Space')
ax.set_zlabel('Mapped -T (Time)')
ax.set_title('z=0 Spacetime Mesh Visualization (R=100, n=10)\nDeformation into Trapezoidal Frustum (台形錐体) via Lorentz Contraction')

# 描画範囲を少し広げて、内側に傾く様子をわかりやすくスケーリング
ax.set_xlim(-1, n)
ax.set_ylim(-1, n)

plt.tight_layout()
plt.savefig('mesh_distortion_frustum.png', dpi=150)
print("Saved mesh_distortion_frustum.png")
