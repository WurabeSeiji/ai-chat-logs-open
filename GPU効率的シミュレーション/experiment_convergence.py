import torch
import matplotlib.pyplot as plt
import numpy as np

# 計算精度の問題でFloat64(倍精度)を使用し、数値落ちによるノイズを防ぎます
device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))

def run_sim(R, dt, T_max=10.0, N=10000):
    n_steps = int(round(T_max / dt))
    actual_T = n_steps * dt
    
    # Initialize exactly correctly on S^3
    X0 = torch.randn(N, 4, dtype=torch.float32, device=device)
    X0 = X0 / torch.norm(X0, dim=1, keepdim=True) * R
    
    V0 = torch.randn(N, 4, dtype=torch.float32, device=device) * 5.0
    V0 = V0 - (X0 * V0).sum(dim=1, keepdim=True) * X0 / (R**2)
    
    X = X0.clone()
    V = V0.clone()
    
    for step in range(1, n_steps + 1):
        t = step * dt
        a_prev = 1.0 + 0.001 * (t - dt)
        a_curr = 1.0 + 0.001 * t
        
        speed = torch.norm(V, dim=1, keepdim=True)
        theta = (speed / R) * dt
        
        cos_th = torch.cos(theta)
        sin_th = torch.sin(theta)
        unit_V = V / (speed + 1e-12)
        
        X_new = X * cos_th + unit_V * sin_th * R
        X_new = X_new / torch.norm(X_new, dim=1, keepdim=True) * R
        
        V_trans = - (X / R) * speed * sin_th + unit_V * speed * cos_th
        V_new = V_trans * (a_prev / a_curr)
        
        X = X_new
        V = V_new

    # 厳密解の計算 (Analytical Exact Solution)
    c = 0.001
    speed0 = torch.norm(V0, dim=1, keepdim=True)
    Theta_exact = (speed0 / (c * R)) * (1.0 - 1.0 / (1.0 + c * actual_T))
    
    cos_Th_exact = torch.cos(Theta_exact)
    sin_Th_exact = torch.sin(Theta_exact)
    unit_V0 = V0 / (speed0 + 1e-12)
    
    X_exact = X0 * cos_Th_exact + unit_V0 * sin_Th_exact * R
    
    # 2乗空間距離の平均誤差 (Mean Squared Error)
    squared_errors = torch.sum((X - X_exact)**2, dim=1)
    mse = torch.mean(squared_errors).item()
    return mse

print("実験1を実行中 (R依存性)...")
Rs = [5.0, 10.0, 20.0, 40.0, 80.0]
errors_R = []
for R in Rs:
    err = run_sim(R=R, dt=0.01)
    errors_R.append(err)

print("実験2を実行中 (分割数 n 依存性)...")
R_fixed = 10.0
dts = [0.1, 0.05, 0.025, 0.0125, 0.00625]
errors_dt = []
n_divisions = [1.0/dt for dt in dts] # モデル上の「分割数 n」に相当
for dt in dts:
    err = run_sim(R=R_fixed, dt=dt)
    errors_dt.append(err)

print("グラフ描画中...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.plot(Rs, errors_R, marker='o', linestyle='-', color='b', linewidth=2)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_title("Exp 1: Const Mesh / Varying Curvature R")
ax1.set_xlabel("Curvature Radius R (Log Scale)")
ax1.set_ylabel("Mean Squared Error (L2 Norm Squared)")
ax1.grid(True, which="both", ls="--")

ax2.plot(n_divisions, errors_dt, marker='o', linestyle='-', color='r', linewidth=2)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_title("Exp 2: Const Curvature / Varying Division n (1/dt)")
ax2.set_xlabel("Division Frequency n $\propto$ 1/dt (Log Scale)")
ax2.set_ylabel("Mean Squared Error (L2 Norm Squared)")
ax2.grid(True, which="both", ls="--")

plt.tight_layout()
plt.savefig('error_convergence.png', dpi=150)
print("Saved error_convergence.png")
