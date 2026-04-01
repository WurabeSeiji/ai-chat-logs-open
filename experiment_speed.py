import time
import torch
import matplotlib.pyplot as plt

device = torch.device("mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu"))

# MPS（Mac GPU）およびCUDA用の同期関数
def sync_device():
    if device.type == "mps":
        torch.mps.synchronize()
    elif device.type == "cuda":
        torch.cuda.synchronize()

def run_speed_benchmark(N, n_steps=1000, R=10.0, dt=0.01):
    sync_device()
    t0 = time.time()
    
    # Step 1: 初期化（行列生成、メモリ確保）
    X0 = torch.randn(N, 4, dtype=torch.float32, device=device)
    X0 = X0 / torch.norm(X0, dim=1, keepdim=True) * R
    
    V0 = torch.randn(N, 4, dtype=torch.float32, device=device) * 5.0
    V0 = V0 - (X0 * V0).sum(dim=1, keepdim=True) * X0 / (R**2)
    
    X = X0.clone()
    V = V0.clone()
    
    sync_device()
    t_step1 = time.time() - t0
    
    # Step 2: シミュレーション本番（1000ループの行列演算）
    t1 = time.time()
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
        
    sync_device()
    t_step2 = time.time() - t1
    return t_step1, t_step2

print(f"Device being used: {device}")
# 粒子数を対数的に変化させる: 1万から100万
N_list = [10000, 50000, 100000, 500000, 1000000]

t1_results = []
t2_results = []

print("Running warmup...")
run_speed_benchmark(100, n_steps=10)

for N in N_list:
    print(f"Benchmarking N={N}...")
    t1, t2 = run_speed_benchmark(N)
    t1_results.append(t1)
    t2_results.append(t2)

# 理論的な O(N) 線形スケーリング（もしCPUで直列に処理した場合の増加予測）
# 基準として N=10000 のときの時間を利用
linear_reference = [t2_results[0] * (n / N_list[0]) for n in N_list]

fig, ax = plt.subplots(figsize=(9, 6))

ax.plot(N_list, t1_results, marker='o', label='Step 1: Initialization', color='blue')
ax.plot(N_list, t2_results, marker='s', label='Step 2: Simulation Loop (GPU Actual)', color='red', linewidth=2)
ax.plot(N_list, linear_reference, marker='x', linestyle='--', label='CPU-like Linear Scaling O(N) Reference', color='gray', alpha=0.7)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Number of Test Particles N (Log Scale)')
ax.set_ylabel('Execution Time [sec] (Log Scale)')
ax.set_title('GPU Processing Time Scaling Benchmark\n(1000 Steps)')
ax.grid(True, which="both", ls="--")
ax.legend()

plt.tight_layout()
plt.savefig('speed_benchmark.png', dpi=150)
print("Saved speed_benchmark.png")
