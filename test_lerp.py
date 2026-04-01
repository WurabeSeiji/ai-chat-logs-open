import numpy as np
import itertools
c = 1.0

def gamma(v_over_c):
    return 1.0 / np.sqrt(1.0 - v_over_c**2)

def exact_mapping(xi, R):
    xi0 = xi[:, 0]
    v = c * xi0 / R
    g = gamma(v/c)
    
    X_phys = np.zeros_like(xi, dtype=np.float64)
    X_phys[:, 0] = xi[:, 0] * g
    X_phys[:, 1] = xi[:, 1] / g
    X_phys[:, 2] = xi[:, 2] / g
    X_phys[:, 3] = xi[:, 3] / g
    return X_phys

def tensor_lerp_mapping(xi, R):
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
            c_val = exact_mapping(c_pos.reshape(1, 4), R)[0]
            val += c_w * c_val
        X_interp[i] = val
    return X_interp

Rs = np.logspace(2, 4, 10) 
steps = 50
dt = 1.05  # so final t is 52.5

for R in Rs:
    xi_path = []
    xi = np.array([0.0, 10.0, 10.0, 10.0])
    v_local = np.array([1.0, 0.2, 0.1, 0.0])
    for s in range(steps):
        xi_path.append(xi.copy())
        xi = xi + v_local * dt
    xi_path = np.array(xi_path)
    
    X_exact = exact_mapping(xi_path, R)
    X_approx_floor = exact_mapping(np.floor(xi_path), R)
    X_approx_lerp = tensor_lerp_mapping(xi_path, R)
    
    dist_floor = np.linalg.norm(X_exact[-1] - X_approx_floor[-1])
    dist_lerp = np.linalg.norm(X_exact[-1] - X_approx_lerp[-1])
    
    print(f"R={R:.1f}, Error(Floor)={dist_floor:.4f}, Error(Lerp)={dist_lerp:.8f}")
