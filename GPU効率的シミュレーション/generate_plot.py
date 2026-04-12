import torch
import matplotlib.pyplot as plt
import os

device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))

n = 50                          
R = 5.0                         
n_steps = 1000
dt = 0.01

N = n * n * n

X = torch.randn(N, 4, device=device)
X = X / torch.norm(X, dim=1, keepdim=True) * R

# 初期速度の振幅を少し上げて軌跡を見やすくする
V = torch.randn(N, 4, device=device) * 5.0
# [修正1] 半径Rの球に対する正しい接空間射影 (R**2で割る必要がある)
V = V - (X * V).sum(dim=1, keepdim=True) * X / (R**2)

history = torch.zeros(n_steps, N, 9, device=device)
history[0, :, 0:4] = X
history[0, :, 4:8] = V
history[0, :, 8] = 0.0

for step in range(1, n_steps):
    t = step * dt
    # [修正2] スケール因子による「毎ステップごとの正しい減衰比率」
    a_prev = 1.0 + 0.001 * (t - dt)
    a_curr = 1.0 + 0.001 * t
    
    speed = torch.norm(V, dim=1, keepdim=True)
    theta = (speed / R) * dt
    
    cos_th = torch.cos(theta)
    sin_th = torch.sin(theta)
    unit_V = V / (speed + 1e-8)
    
    # [修正3] 球面上での正しい進み方。unit_V の係数に R を掛ける必要がある
    X_new = X * cos_th + unit_V * sin_th * R
    X_new = X_new / torch.norm(X_new, dim=1, keepdim=True) * R
    
    # [修正4] 速度ベクトル自体の平行移動 (位置が回転した分だけ、接ベクトルも方向を変える)
    V_transport = - (X / R) * speed * sin_th + unit_V * speed * cos_th
    # スケール因子膨張による赤方偏移（V ∝ 1/a の減衰）
    V_new = V_transport * (a_prev / a_curr)
    
    v_pec = speed.squeeze(1) / (a_curr * R)
    dtau = dt * torch.sqrt(torch.clamp(1.0 - v_pec**2, min=1e-8))
    
    X = X_new
    V = V_new
    tau = history[step-1, :, 8] + dtau
    
    history[step, :, 0:4] = X
    history[step, :, 4:8] = V
    history[step, :, 8] = tau

print("シミュレーション完了。描画を開始します...")

num_samples = 5
history_cpu = history[:, :num_samples, 0:3].cpu().numpy()

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

colors = ['b', 'g', 'r', 'c', 'm']
for i in range(num_samples):
    ax.plot(history_cpu[:, i, 0], history_cpu[:, i, 1], history_cpu[:, i, 2], color=colors[i%len(colors)], alpha=0.7, linewidth=1.5, label=f'Particle {i}')
    ax.scatter(history_cpu[0, i, 0], history_cpu[0, i, 1], history_cpu[0, i, 2], color=colors[i%len(colors)], s=50, marker='o') # 始点
    ax.scatter(history_cpu[-1, i, 0], history_cpu[-1, i, 1], history_cpu[-1, i, 2], color=colors[i%len(colors)], s=50, marker='x') # 終点

u = torch.linspace(0, 2 * torch.pi, 30)
v = torch.linspace(0, torch.pi, 30)
x_sphere = R * torch.outer(torch.cos(u), torch.sin(v)).numpy()
y_sphere = R * torch.outer(torch.sin(u), torch.sin(v)).numpy()
z_sphere = R * torch.outer(torch.ones(30), torch.cos(v)).numpy()
ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color='k', alpha=0.1)

ax.set_title("Test Particle Trajectories in S^3 (Orthographic 3D Projection)")
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim([-R, R])
ax.set_ylim([-R, R])
ax.set_zlim([-R, R])
ax.legend()

plt.savefig('particle_tracks_3d.png', dpi=150, bbox_inches='tight')
print("3Dグラフを particle_tracks_3d.png として保存しました。")
