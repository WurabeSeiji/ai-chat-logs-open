import numpy as np
import matplotlib.pyplot as plt
import itertools

# Physical constants
c = 1.0

def gamma(v_over_c):
    v_over_c = np.clip(v_over_c, -0.9999, 0.9999)
    return 1.0 / np.sqrt(1.0 - v_over_c**2)

def mapping(xi, R=100.0):
    t = xi[:, 0]
    v = c * t / R
    g = gamma(v / c)
    X = np.zeros_like(xi)
    X[:, 0] = t * g
    X[:, 1] = xi[:, 1] / g
    X[:, 2] = xi[:, 2] / g
    X[:, 3] = xi[:, 3] / g
    return X

def tensor_lerp_mapping(xi, R=100.0):
    X_interp = np.zeros_like(xi, dtype=np.float64)
    corners = np.array(list(itertools.product([0, 1], repeat=4)))
    idx_floor = np.floor(xi)
    w = xi - idx_floor
    
    for i in range(len(xi)):
        p_floor = idx_floor[i]
        p_w = w[i]
        val = np.zeros(4)
        for c_idx in corners:
            c_pos = p_floor + c_idx
            c_w = np.prod(np.where(c_idx == 1, p_w, 1.0 - p_w))
            c_val = mapping(c_pos.reshape(1, 4), R)[0]
            val += c_w * c_val
        X_interp[i] = val
    return X_interp

# Simulation settings
R = 100.0
dt_exact = 0.05
steps_exact = int(80 / dt_exact)

# --- Particle A (Stationary in index space, moving in time) ---
xi_A_exact = np.zeros((steps_exact, 4))
xi_A_exact[:, 0] = np.arange(0, 80, dt_exact)
xi_A_exact[:, 1] = 18.0
xi_A_exact[:, 2] = 18.0

X_phys_A_exact = mapping(xi_A_exact, R)
X_phys_A_approx = tensor_lerp_mapping(xi_A_exact, R) # Now using Lerp

# --- Particle B (Moving horizontally in index space) ---
v_x = 0.35 
xi_B_exact = np.zeros((steps_exact, 4))
xi_B_exact[:, 0] = np.arange(0, 80, dt_exact)
xi_B_exact[:, 1] = -10.0 + xi_B_exact[:, 0] * v_x
xi_B_exact[:, 2] = 18.0

X_phys_B_exact = mapping(xi_B_exact, R)
X_phys_B_approx = tensor_lerp_mapping(xi_B_exact, R) # Now using Lerp

# Plotting XY Plane
plt.figure(figsize=(10, 8))

# Draw Exact Paths First (Thick translucent lines)
plt.plot(X_phys_A_exact[:, 1], X_phys_A_exact[:, 2], 
         color='lightcoral', linestyle='-', linewidth=6, alpha=0.5, label='Particle A (Exact Geodesic)')
plt.plot(X_phys_B_exact[:, 1], X_phys_B_exact[:, 2], 
         color='cornflowerblue', linestyle='-', linewidth=6, alpha=0.5, label='Particle B (Exact Geodesic)')

# Draw Approx Paths (Dashed lines to show precise alignment)
plt.plot(X_phys_A_approx[:, 1], X_phys_A_approx[:, 2], 
         color='red', linestyle='--', linewidth=2, label='Particle A (Tensor Lerp Approximation)')
plt.plot(X_phys_B_approx[:, 1], X_phys_B_approx[:, 2], 
         color='blue', linestyle='--', linewidth=2, label='Particle B (Tensor Lerp Approximation)')

plt.scatter([0], [0], color='black', marker='X', s=200, label='Center of Curvature (Gravity)')
plt.xlabel('Physical X-coordinate')
plt.ylabel('Physical Y-coordinate')
plt.title('Figure 2: Physical Geodesic Trajectories\n(Exact Solution vs Tensor GPU Linear Interpolation)')

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='lower left')
plt.axis('equal') 
plt.tight_layout()

plt.savefig('tensor_lookup_geodesic.png', dpi=150)
print("Generated tensor_lookup_geodesic.png with Lerp.")
