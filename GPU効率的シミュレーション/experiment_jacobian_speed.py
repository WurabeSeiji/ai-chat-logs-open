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

def bench_methods(N, n_steps=1000, R=10.0, dt=0.01):
    sync_device()
    # 初期位置（S^3球面上）
    X = torch.randn(N, 4, dtype=torch.float32, device=device)
    X = X / torch.norm(X, dim=1, keepdim=True) * R
    
    # === 手法A（従来手法）のセットアップ ===
    Va = torch.randn(N, 4, dtype=torch.float32, device=device) * 5.0
    Va = Va - (X * Va).sum(dim=1, keepdim=True) * X / (R**2)
    Xa = X.clone()
    
    # === 手法B（提案手法：Jacobianメッシュ線形近似）のセットアップ ===
    Xb = X.clone()
    # Step 1で計算済みのJacobian行列 J を読み込んだと仮定 (N個の 4x3 行列)
    J = torch.randn(N, 4, 3, dtype=torch.float32, device=device)
    # 局所メッシュ内の3次元速度
    V_local = torch.randn(N, 3, 1, dtype=torch.float32, device=device) * 5.0
    
    # --- Benchmark 手法A: 毎回の非線形・厳密計算（sin, cos）重いループ ---
    sync_device()
    t0_a = time.time()
    for step in range(1, n_steps + 1):
        t = step * dt
        a_curr = 1.0 + 0.001 * t
        
        # 三角関数や平方根を多用する重い非線形演算
        speed = torch.norm(Va, dim=1, keepdim=True)
        theta = (speed / R) * dt
        cos_th = torch.cos(theta)
        sin_th = torch.sin(theta)
        unit_V = Va / (speed + 1e-12)
        
        Xa_new = Xa * cos_th + unit_V * sin_th * R
        Xa = Xa_new
        
        # 速度ベクトルの正確な球面平行移動と赤方偏移（減速）
        V_trans = - (Xa / R) * speed * sin_th + unit_V * speed * cos_th
        Va = V_trans * (1.0 / a_curr)
        
    sync_device()
    time_a = time.time() - t0_a
    
    # --- Benchmark 手法B: 提案のJacobian行列を用いた直交線形ループ ---
    sync_device()
    t0_b = time.time()
    for step in range(1, n_steps + 1):
        t = step * dt
        a_curr = 1.0 + 0.001 * t
        
        # ローカル速度に赤方偏移（減速）のみ適用
        V_local_t = V_local * (1.0 / a_curr)
        
        # 事前Jacobian行列との線形積（batch_matmul）のみで空間変位を高速計算
        # 非線形関数(sin, cos)を完全に排除
        dX = torch.bmm(J, V_local_t).squeeze(-1) * dt
        Xb = Xb + dX
        # （球面上に戻す規格化のみ実施）
        Xb = Xb / torch.norm(Xb, dim=1, keepdim=True) * R
        
    sync_device()
    time_b = time.time() - t0_b
    
    return time_a, time_b

print(f"Device being used: {device}")
N_list = [10000, 50000, 100000, 500000, 1000000]

timea_results = []
timeb_results = []

print("Running warmup...")
bench_methods(1000, 10)

for N in N_list:
    print(f"Benchmarking N={N}...")
    ta, tb = bench_methods(N)
    timea_results.append(ta)
    timeb_results.append(tb)
    # それぞれのNにおける速度比（何倍速くなったか）を表示
    print(f"  Method A (Exact): {ta:.4f} sec | Method B (Jacobian): {tb:.4f} sec -> Speedup: {ta/tb:.2f}x")

fig, ax = plt.subplots(figsize=(9, 6))

ax.plot(N_list, timea_results, marker='s', label='Method A: Continuous / Non-linear (sin, cos)', color='red', linewidth=2)
ax.plot(N_list, timeb_results, marker='o', label='Method B: Proposed Linear Jacobian Mesh', color='blue', linewidth=2)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Number of Test Particles N (Log Scale)')
ax.set_ylabel('Execution Time [sec] (Log Scale)')
ax.set_title('Computational Time Comparison of Simulation Loops\n(1000 Steps on GPU)')
ax.grid(True, which="both", ls="--")
ax.legend()

plt.tight_layout()
plt.savefig('speed_benchmark_jacobian.png', dpi=150)
print("Saved speed_benchmark_jacobian.png")
