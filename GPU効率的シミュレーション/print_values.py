import numpy as np
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

Rs = np.logspace(np.log10(55), 4, 10) 
steps = 500
dt = 0.1 
v_local = np.array([1.0, 0.2, 0.1, 0.0])

for R in Rs:
    xi_path = []
    xi = np.array([0.0, 10.0, 10.0, 10.0])
    for s in range(steps):
        xi_path.append(xi.copy())
        xi = xi + v_local * dt
    xi_path = np.array(xi_path)
    
    X_exact = exact_mapping(xi_path, R)
    diffs_ex = np.diff(X_exact, axis=0)
    L_exact = np.sum(np.linalg.norm(diffs_ex, axis=1))
    
    X_approx = exact_mapping(np.floor(xi_path), R)
    diffs_ap = np.diff(X_approx, axis=0)
    L_approx = np.sum(np.linalg.norm(diffs_ap, axis=1))
    
    print(f"R={R:.1f}, L_exact={L_exact:.3f}, L_approx={L_approx:.3f}, Diff={abs(L_exact-L_approx):.3f}")
